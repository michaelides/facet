#' Recalculate VAR with Alternative Weights
#'
#' @keywords internal
recalculate_var_with_weights <- function(dstudy_obj, weights, dims, ci, probs, n = NULL) {

  gstudy_obj <- dstudy_obj$gstudy
  dimensions <- dims
  universe_spec <- dstudy_obj$universe
  error_spec <- dstudy_obj$error
  object_spec <- dstudy_obj$object

  is_mom <- inherits(gstudy_obj$model, "momfit")

  if (is.null(n)) {
    if (isTRUE(dstudy_obj$is_sweep)) {
      n <- as.list(gstudy_obj$facet_n)
    } else {
      n <- dstudy_obj$n
    }
  }

  if (is_mom) {
    mom_draws <- generate_mom_variance_and_covariance_draws(gstudy_obj)
    vc_draws <- mom_draws$vc_draws
    cov_draws <- mom_draws$cov_draws
    n_draws <- mom_draws$n_draws
  } else {
    vc_draws <- extract_variance_draws_from_gstudy(gstudy_obj)
    cov_draws <- extract_covariance_draws(gstudy_obj$model, dimensions)
    n_draws <- length(vc_draws[[1]][[1]])
  }

  components <- unique(unlist(lapply(vc_draws, names)))

  error_info <- identify_error_components_for_draws(
    components = components,
    universe_spec = universe_spec,
    error_spec = error_spec
  )

  n_dims <- length(dimensions)

  scale_factors <- compute_scale_factors_for_viable(components, n, universe_spec, object_spec)

  # Construct covariance arrays for Haberman formula
  # uni_cov_draws: universe-score covariance (G-study, UNSCALED)
  # total_rel_cov_draws: observed covariance (D-study, SCALED)
  # total_abs_cov_draws: observed covariance (D-study, SCALED)
  uni_cov_draws <- array(0, dim = c(n_draws, n_dims, n_dims))
  total_rel_cov_draws <- array(0, dim = c(n_draws, n_dims, n_dims))
  total_abs_cov_draws <- array(0, dim = c(n_draws, n_dims, n_dims))

  for (i in seq_along(dimensions)) {
    d1 <- dimensions[i]
    for (j in seq_along(dimensions)) {
      d2 <- dimensions[j]

      for (comp_idx in seq_along(components)) {
        comp <- components[comp_idx]
        sf <- scale_factors[[comp]]

        # Determine covariance list for off-diagonal elements
        cov_list <- NULL
        if (i != j) {
          pair_name <- paste0(pmin(d1, d2), "_", pmax(d1, d2))
          if (comp == "Residual" || grepl("^Residual", comp)) {
            cov_list <- cov_draws$residual
          } else if (!is.null(cov_draws$random_effect) && comp %in% names(cov_draws$random_effect)) {
            cov_list <- cov_draws$random_effect[[comp]]
          }
        }

        # Universe components: UNSCALED (G-study variance)
        if (comp %in% universe_spec) {
          val_uni <- if (i == j) {
            if (!is.null(vc_draws[[d1]][[comp]])) vc_draws[[d1]][[comp]] else numeric(n_draws)
          } else {
            if (!is.null(cov_list) && pair_name %in% names(cov_list)) {
              cov_list[[pair_name]]
            } else {
              numeric(n_draws)
            }
          }
          uni_cov_draws[, i, j] <- uni_cov_draws[, i, j] + val_uni
        }

        # Observed covariance: SCALED (D-study variance)
        val_obs <- if (i == j) {
          if (!is.null(vc_draws[[d1]][[comp]])) vc_draws[[d1]][[comp]] / sf else numeric(n_draws)
        } else {
          if (!is.null(cov_list) && pair_name %in% names(cov_list)) {
            cov_list[[pair_name]] / sf
          } else {
            numeric(n_draws)
          }
        }

        if (error_info$is_relative_error[comp_idx] || comp %in% universe_spec) {
          total_rel_cov_draws[, i, j] <- total_rel_cov_draws[, i, j] + val_obs
        }
        if (error_info$is_absolute_error[comp_idx] || comp %in% universe_spec) {
          total_abs_cov_draws[, i, j] <- total_abs_cov_draws[, i, j] + val_obs
        }
      }
    }
  }

  # Calculate composite reliability draws from the arrays
  w <- weights[dimensions]
  g_draws <- numeric(n_draws)
  phi_draws <- numeric(n_draws)

  for (draw_idx in 1:n_draws) {
    U <- uni_cov_draws[draw_idx, , ]
    TR <- total_rel_cov_draws[draw_idx, , ]
    TA <- total_abs_cov_draws[draw_idx, , ]

    uni_var_C <- as.numeric(t(w) %*% U %*% w)
    total_rel_var_C <- as.numeric(t(w) %*% TR %*% w)
    total_abs_var_C <- as.numeric(t(w) %*% TA %*% w)

    g_draws[draw_idx] <- if (total_rel_var_C > 0) uni_var_C / total_rel_var_C else NA_real_
    phi_draws[draw_idx] <- if (total_abs_var_C > 0) uni_var_C / total_abs_var_C else NA_real_
  }

  var_result <- compute_var_haberman_draws(
    uni_cov_draws = uni_cov_draws,
    total_rel_cov_draws = total_rel_cov_draws,
    total_abs_cov_draws = total_abs_cov_draws,
    dimensions = dimensions,
    weights = weights
  )

  # Result columns with Composite row
  # prmse_s: univariate reliability (equivalent to G or Phi)
  # prmse_c: prediction using the weighted composite
  # prmse_p: prediction using the whole profile (Haberman's MV framework)

  base_cols <- list(
    dim = c(dimensions, "Composite"),
    prmse_s_rel = c(colMeans(var_result$prmse_s_rel, na.rm = TRUE), mean(g_draws, na.rm = TRUE)),
    prmse_s_abs = c(colMeans(var_result$prmse_s_abs, na.rm = TRUE), mean(phi_draws, na.rm = TRUE)),
    prmse_c_rel = c(colMeans(var_result$prmse_c_rel, na.rm = TRUE), 1.0),
    prmse_c_abs = c(colMeans(var_result$prmse_c_abs, na.rm = TRUE), 1.0),
    prmse_p_rel = c(colMeans(var_result$prmse_p_rel, na.rm = TRUE), mean(var_result$prmse_mv_rel, na.rm = TRUE)),
    prmse_p_abs = c(colMeans(var_result$prmse_p_abs, na.rm = TRUE), mean(var_result$prmse_mv_abs, na.rm = TRUE)),
    var_rel = c(colMeans(var_result$var_rel, na.rm = TRUE), 1.0),
    var_abs = c(colMeans(var_result$var_abs, na.rm = TRUE), 1.0)
  )

  if (!is.null(ci)) {
    if ("prmse" %in% ci) {
      # CIs for S, C, and P metrics
      base_cols$prmse_s_rel_LL <- c(apply(var_result$prmse_s_rel, 2, quantile, probs = probs[1], na.rm = TRUE), quantile(g_draws, probs[1], na.rm = TRUE))
      base_cols$prmse_s_rel_UL <- c(apply(var_result$prmse_s_rel, 2, quantile, probs = probs[2], na.rm = TRUE), quantile(g_draws, probs[2], na.rm = TRUE))
      base_cols$prmse_s_abs_LL <- c(apply(var_result$prmse_s_abs, 2, quantile, probs = probs[1], na.rm = TRUE), quantile(phi_draws, probs[1], na.rm = TRUE))
      base_cols$prmse_s_abs_UL <- c(apply(var_result$prmse_s_abs, 2, quantile, probs = probs[2], na.rm = TRUE), quantile(phi_draws, probs[2], na.rm = TRUE))

      base_cols$prmse_c_rel_LL <- c(apply(var_result$prmse_c_rel, 2, quantile, probs = probs[1], na.rm = TRUE), unname(1.0))
      base_cols$prmse_c_rel_UL <- c(apply(var_result$prmse_c_rel, 2, quantile, probs = probs[2], na.rm = TRUE), unname(1.0))

      base_cols$prmse_p_rel_LL <- c(apply(var_result$prmse_p_rel, 2, quantile, probs = probs[1], na.rm = TRUE), quantile(var_result$prmse_mv_rel, probs[1], na.rm = TRUE))
      base_cols$prmse_p_rel_UL <- c(apply(var_result$prmse_p_rel, 2, quantile, probs = probs[2], na.rm = TRUE), quantile(var_result$prmse_mv_rel, probs[2], na.rm = TRUE))
    }
    if ("var" %in% ci) {
      base_cols$var_rel_LL <- c(apply(var_result$var_rel, 2, quantile, probs = probs[1], na.rm = TRUE), unname(1.0))
      base_cols$var_rel_UL <- c(apply(var_result$var_rel, 2, quantile, probs = probs[2], na.rm = TRUE), unname(1.0))
      base_cols$var_abs_LL <- c(apply(var_result$var_abs, 2, quantile, probs = probs[1], na.rm = TRUE), unname(1.0))
      base_cols$var_abs_UL <- c(apply(var_result$var_abs, 2, quantile, probs = probs[2], na.rm = TRUE), unname(1.0))
    }
  }

  tibble::as_tibble(base_cols)
}

#' Compute Scale Factors for Viable Recalculation
#'
#' Computes scale factors for multiple variance components.
#' Universe components always receive a scale factor of 1 (unscaled).
#' Non-universe components are scaled by the product of non-object facet
#' sample sizes in that component.
#'
#' @param components Character vector of variance component names
#' @param n Named list of sample sizes
#' @param universe_spec Universe components specification
#' @param object_spec Object of measurement specification
#' @keywords internal
compute_scale_factors_for_viable <- function(components, n, universe_spec, object_spec) {
  scale_factors <- list()

  for (comp in components) {
    if (comp %in% universe_spec) {
      scale_factors[[comp]] <- 1
    } else if (comp == "Residual") {
      sf <- 1
      for (facet in names(n)) {
        sf <- sf * n[[facet]]
      }
      scale_factors[[comp]] <- sf
    } else {
      facets <- parse_component_facets(comp)
      sf <- 1
      for (facet in facets) {
        if (facet %in% names(n) && !(facet %in% object_spec)) {
          sf <- sf * n[[facet]]
        }
      }
      scale_factors[[comp]] <- sf
    }
  }

  scale_factors
}

#' Extract Variance Draws from G-Study for Viable Recalculation
#'
#' @param gstudy_obj A gstudy object
#' @param n_draws Number of pseudo-draws to generate for mom backend
#' @keywords internal
extract_variance_draws_from_gstudy <- function(gstudy_obj, n_draws = 1000) {

  if (inherits(gstudy_obj$model, "brmsfit")) {
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required for brms backend.", call. = FALSE)
    }
    draws <- brms::as_draws_matrix(gstudy_obj$model)
    extract_variance_draws(gstudy_obj, draws)
  } else if (inherits(gstudy_obj$model, "momfit")) {
    stop("mom backend should use generate_mom_variance_and_covariance_draws() directly", call. = FALSE)
  } else {
    stop("Unsupported backend. Weight optimization requires brms or mom backend.", call. = FALSE)
  }
}

#' Generate Variance and Covariance Draws from mom Point Estimates
#'
#' Uses variance estimates and SEs to generate normal pseudo-posterior draws,
#' and converts correlation draws to covariance draws.
#'
#' @param gstudy_obj A gstudy object with mom backend
#' @param n_draws Number of pseudo-draws to generate (default 1000)
#'
#' @keywords internal
generate_mom_variance_and_covariance_draws <- function(gstudy_obj, n_draws = 1000) {
  model <- gstudy_obj$model
  vc <- model$variance_components
  dimensions <- gstudy_obj$dimensions
  correlations <- model$correlations

  if (is.null(dimensions) || length(dimensions) == 0) {
    stop("mom gstudy must have dimensions for weight optimization", call. = FALSE)
  }

  components <- unique(vc$component)
  k <- length(dimensions)

  # Compute SE if not present: SE = sqrt(2 * MS^2 / df)
  if (!"se" %in% names(vc)) {
    vc$se <- ifelse(vc$df > 0, sqrt(2 * vc$ms^2 / vc$df), 0)
  }

  vc_draws <- list()
  for (d in dimensions) {
    vc_draws[[d]] <- list()
    vc_d <- vc[vc$dim == d, , drop = FALSE]
    for (comp in components) {
      var_row <- vc_d$var[vc_d$component == comp]
      se_row <- vc_d$se[vc_d$component == comp]
      if (length(var_row) > 0 && !is.na(var_row[1]) && var_row[1] >= 0) {
        var_val <- var_row[1]
        se_val <- if (length(se_row) > 0 && !is.na(se_row[1]) && se_row[1] > 0) se_row[1] else 0
        vc_draws[[d]][[comp]] <- rnorm(n_draws, mean = var_val, sd = se_val)
      } else {
        vc_draws[[d]][[comp]] <- rep(NA_real_, n_draws)
      }
    }
  }

  cov_to_draws <- function(d1, d2, comp, cor_draws) {
    sd_draws_d1 <- vc_draws[[d1]][[comp]]
    sd_draws_d2 <- vc_draws[[d2]][[comp]]

    if (is.null(sd_draws_d1) || is.null(sd_draws_d2)) {
      return(rep(NA_real_, n_draws))
    }

    sd_draws_d1 <- sqrt(pmax(sd_draws_d1, 0))
    sd_draws_d2 <- sqrt(pmax(sd_draws_d2, 0))
    cov_draws <- sd_draws_d1 * sd_draws_d2 * cor_draws
    cov_draws[is.na(sd_draws_d1) | is.na(sd_draws_d2) | is.na(cor_draws)] <- NA_real_
    cov_draws
  }

  cov_draws <- list(residual = list(), random_effect = list())

  if (!is.null(correlations)) {
    if (!is.null(correlations$residual_cor) && nrow(correlations$residual_cor) > 0) {
      for (i in seq_len(nrow(correlations$residual_cor))) {
        row <- correlations$residual_cor[i, ]
        d1 <- row$dim1
        d2 <- row$dim2
        if (!(d1 %in% dimensions) || !(d2 %in% dimensions)) next
        pair_name <- paste0(pmin(d1, d2), "_", pmax(d1, d2))
        cor_est <- row$estimate
        cor_se <- if (!is.null(row$se)) row$se else 0
        cor_draws <- rnorm(n_draws, mean = cor_est, sd = cor_se)
        cov_draws$residual[[pair_name]] <- cov_to_draws(d1, d2, "Residual", cor_draws)
      }
    }

    if (!is.null(correlations$random_effect_cor)) {
      for (facet in names(correlations$random_effect_cor)) {
        re_cor <- correlations$random_effect_cor[[facet]]
        if (is.null(re_cor) || nrow(re_cor) == 0) next
        cov_draws$random_effect[[facet]] <- list()
        for (i in seq_len(nrow(re_cor))) {
          row <- re_cor[i, ]
          d1 <- row$dim1
          d2 <- row$dim2
          if (!(d1 %in% dimensions) || !(d2 %in% dimensions)) next
          pair_name <- paste0(pmin(d1, d2), "_", pmax(d1, d2))
          cor_est <- row$estimate
          cor_se <- if (!is.null(row$se)) row$se else 0
          cor_draws <- rnorm(n_draws, mean = cor_est, sd = cor_se)
          cov_draws$random_effect[[facet]][[pair_name]] <- cov_to_draws(d1, d2, facet, cor_draws)
        }
      }
    }
  }

  list(
    vc_draws = vc_draws,
    cov_draws = cov_draws,
    n_draws = n_draws
  )
}
