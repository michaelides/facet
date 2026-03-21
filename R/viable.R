#' Subscale Viability Analysis
#'
#' Performs Subscale Viability Analysis on a multivariate dstudy object,
#' returning PRMSE and VAR metrics with optional credible intervals.
#'
#' @param dstudy_obj A dstudy object from [dstudy()]
#' @param ci Character vector specifying which metrics to compute CIs for.
#'   Options: "prmse", "var", or both. Default NULL (no CIs).
#' @param probs Numeric vector of length 2 specifying the quantile probabilities
#'   for credible intervals. Default c(0.025, 0.975) for 95% CI.
#' @param weights Numeric vector of weights for computing composite coefficients.
#'   Length must match the number of dimensions. Default NULL uses weights from dstudy.
#'
#' @return A tibble with one row per dimension containing:
#' \describe{
#' \item{dim}{Dimension/subscale name}
#' \item{prmse_c_rel}{PRMSE(C→S_i) (relative) posterior mean}
#' \item{prmse_c_abs}{PRMSE(C→S_i) (absolute) posterior mean}
#' \item{var_rel}{VAR (relative) posterior mean}
#' \item{var_abs}{VAR (absolute) posterior mean}
#' }
#' When `ci` is specified, additional `_LL` and `_UL` columns are included
#' for the corresponding metrics.
#'
#' @details
#' PRMSE(C→S_i) is the proportional reduction in mean squared error when using
#' the composite score to estimate the subscale's universe score.
#'
#' VAR (Value Added Ratio) quantifies whether a subscale's own scores better
#' estimate its universe score than projecting from the composite score.
#' VAR > 1 suggests the subscale adds value and should be reported separately.
#'
#' When alternative weights are provided, all calculations are performed using
#' posterior draws to properly propagate uncertainty. This ensures that changing
#' weights recalculates both the composite coefficients and the resulting
#' VAR/PRMSE values.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Multivariate G-study with brms backend
#' library(brms)
#' g_mv <- gstudy(
#'   bf(Score ~ 0 + Subtest + (0 + Subtest | r | Person)),
#'   data = data, backend = "brms"
#' )
#' d_mv <- dstudy(g_mv, n = list(Person = 5))
#'
#' # Basic viability analysis (no CIs)
#' viable(d_mv)
#'
#' # With 95% CIs for PRMSE only
#' viable(d_mv, ci = "prmse")
#'
#' # With 90% CIs for both metrics
#' viable(d_mv, ci = c("prmse", "var"), probs = c(0.05, 0.95))
#'
#' # With custom weights
#' viable(d_mv, weights = c(Reading = 2, Math = 1, Science = 1))
#' }
viable <- function(dstudy_obj, ci = NULL, probs = c(0.025, 0.975), weights = NULL) {

  if (!inherits(dstudy_obj, "dstudy")) {
    stop("'dstudy_obj' must be a dstudy object", call. = FALSE)
  }

  if (!isTRUE(dstudy_obj$is_multivariate)) {
    stop("No VAR results available. VAR is only computed for multivariate models.", call. = FALSE)
  }

  if (is.null(dstudy_obj$posterior) && is.null(dstudy_obj$composite_posterior) && 
      !inherits(dstudy_obj$gstudy$model, "brmsfit")) {
    stop("viable() requires posterior estimation. Re-run dstudy with a brms backend.", call. = FALSE)
  }

  if (!is.null(ci)) {
    ci <- match.arg(ci, c("prmse", "var"), several.ok = TRUE)
  }

  if (length(probs) != 2) {
    stop("'probs' must have exactly 2 elements", call. = FALSE)
  }
  if (probs[1] >= probs[2]) {
    stop("'probs' must be in increasing order", call. = FALSE)
  }
  if (any(probs < 0) || any(probs > 1)) {
    stop("'probs' must be between 0 and 1", call. = FALSE)
  }

  is_sweep <- isTRUE(dstudy_obj$is_sweep)

  if (is.null(dstudy_obj$var) && !is_sweep) {
    stop("No VAR results available. VAR is only computed for multivariate models.", call. = FALSE)
  }

  if (!is_sweep) {
    dims <- names(dstudy_obj$var$var_rel)
    if (is.null(dims)) {
      dims <- colnames(dstudy_obj$var$var_rel_draws)
    }
  } else {
    dims <- dstudy_obj$gstudy$dimensions
  }
  n_dims <- length(dims)

  if (is.null(weights)) {
    weights <- dstudy_obj$weights
    if (is.null(weights)) {
      weights <- rep(1, n_dims)
      names(weights) <- dims
    }
  } else {
    if (length(weights) != n_dims) {
      stop("weights must have length ", n_dims, " (number of dimensions), got ", length(weights), call. = FALSE)
    }
    names(weights) <- dims
  }

  weights_match <- !is_sweep && isTRUE(all.equal(weights, dstudy_obj$weights))

  if (is_sweep) {
    actual_n <- as.list(dstudy_obj$gstudy$facet_n)
    return(recalculate_var_with_weights(dstudy_obj, weights, dims, ci, probs, actual_n))
  }

  if (weights_match) {
    var <- dstudy_obj$var

    use_stored <- FALSE
    if (!is.null(dstudy_obj$probs) && !is.null(ci)) {
      if (isTRUE(all.equal(probs[1], dstudy_obj$probs[1])) &&
          isTRUE(all.equal(probs[2], dstudy_obj$probs[2]))) {
        if (!is.null(var$var_rel_LL) && !is.null(var$prmse_c_rel_LL)) {
          use_stored <- TRUE
        }
      }
    }

    base_cols <- list(
      dim = dims,
      prmse_c_rel = unname(var$prmse_c_rel),
      prmse_c_abs = unname(var$prmse_c_abs),
      var_rel = unname(var$var_rel),
      var_abs = unname(var$var_abs)
    )

    if (is.null(ci)) {
      return(tibble::as_tibble(base_cols))
    }

    if (use_stored) {
      if ("prmse" %in% ci) {
        base_cols$prmse_c_rel_LL <- unname(var$prmse_c_rel_LL)
        base_cols$prmse_c_rel_UL <- unname(var$prmse_c_rel_UL)
        base_cols$prmse_c_abs_LL <- unname(var$prmse_c_abs_LL)
        base_cols$prmse_c_abs_UL <- unname(var$prmse_c_abs_UL)
      }
      if ("var" %in% ci) {
        base_cols$var_rel_LL <- unname(var$var_rel_LL)
        base_cols$var_rel_UL <- unname(var$var_rel_UL)
        base_cols$var_abs_LL <- unname(var$var_abs_LL)
        base_cols$var_abs_UL <- unname(var$var_abs_UL)
      }
    } else {
      if ("prmse" %in% ci) {
        base_cols$prmse_c_rel_LL <- apply(var$prmse_c_rel_draws, 2, quantile, probs = probs[1], na.rm = TRUE)
        base_cols$prmse_c_rel_UL <- apply(var$prmse_c_rel_draws, 2, quantile, probs = probs[2], na.rm = TRUE)
        base_cols$prmse_c_abs_LL <- apply(var$prmse_c_abs_draws, 2, quantile, probs = probs[1], na.rm = TRUE)
        base_cols$prmse_c_abs_UL <- apply(var$prmse_c_abs_draws, 2, quantile, probs = probs[2], na.rm = TRUE)
      }
      if ("var" %in% ci) {
        base_cols$var_rel_LL <- apply(var$var_rel_draws, 2, quantile, probs = probs[1], na.rm = TRUE)
        base_cols$var_rel_UL <- apply(var$var_rel_draws, 2, quantile, probs = probs[2], na.rm = TRUE)
        base_cols$var_abs_LL <- apply(var$var_abs_draws, 2, quantile, probs = probs[1], na.rm = TRUE)
        base_cols$var_abs_UL <- apply(var$var_abs_draws, 2, quantile, probs = probs[2], na.rm = TRUE)
      }
    }

    return(tibble::as_tibble(base_cols))
  }

  recalculate_var_with_weights(dstudy_obj, weights, dims, ci, probs, n = NULL)
}

#' Recalculate VAR with Alternative Weights
#'
#' @param n Optional named list of sample sizes. If NULL, uses sample sizes
#'   from the dstudy object. For sweep mode, this should be the actual
#'   sample sizes from the G-study.
#' @keywords internal
recalculate_var_with_weights <- function(dstudy_obj, weights, dims, ci, probs, n = NULL) {

  gstudy_obj <- dstudy_obj$gstudy
  dimensions <- dims
  universe_spec <- dstudy_obj$universe
  error_spec <- dstudy_obj$error
  object_spec <- dstudy_obj$object

  if (is.null(n)) {
    n <- dstudy_obj$n
  }

  vc_draws <- extract_variance_draws_from_gstudy(gstudy_obj)
  cov_draws <- gstudy_obj$correlations

  components <- names(vc_draws[[1]])

  components_for_error <- components
  components_for_error <- components_for_error[components_for_error != "Residual" |
                                                  !("Residual" %in% components_for_error)]

  error_info <- identify_error_components_for_draws(
    components = components,
    universe_spec = universe_spec,
    error_spec = error_spec
  )

  n_draws <- length(vc_draws[[1]][[1]])

  scale_factors <- compute_scale_factors_for_viable(components, n, universe_spec, object_spec)

  scaled_vc_draws <- list()
  for (d in dimensions) {
    scaled_vc_draws[[d]] <- list()
    for (comp in components) {
      sf <- scale_factors[[comp]]
      scaled_vc_draws[[d]][[comp]] <- vc_draws[[d]][[comp]] / sf
    }
  }

  composite_variance_draws <- list()
  for (comp in components) {
    comp_draws <- numeric(n_draws)

    for (draw_idx in 1:n_draws) {
      sigma_mat <- matrix(0, length(dimensions), length(dimensions))
      rownames(sigma_mat) <- dimensions
      colnames(sigma_mat) <- dimensions

      for (i in seq_along(dimensions)) {
        d <- dimensions[i]
        sigma_mat[i, i] <- scaled_vc_draws[[d]][[comp]][draw_idx]
      }

      if (!is.null(cov_draws)) {
        cov_list <- NULL
        if (comp == "Residual" || grepl("^Residual", comp)) {
          cov_list <- cov_draws$residual
        } else if (!is.null(cov_draws$random_effect) && comp %in% names(cov_draws$random_effect)) {
          cov_list <- cov_draws$random_effect[[comp]]
        }

        if (!is.null(cov_list)) {
          for (pair_name in names(cov_list)) {
            dims_pair <- strsplit(pair_name, "_")[[1]]
            if (length(dims_pair) == 2 && all(dims_pair %in% dimensions)) {
              idx1 <- match(dims_pair[1], dimensions)
              idx2 <- match(dims_pair[2], dimensions)
              cov_val <- cov_list[[pair_name]][draw_idx] / scale_factors[[comp]]
              sigma_mat[idx1, idx2] <- cov_val
              sigma_mat[idx2, idx1] <- cov_val
            }
          }
        }
      }

      w <- weights[dimensions]
      comp_draws[draw_idx] <- as.numeric(t(w) %*% sigma_mat %*% w)
    }
    composite_variance_draws[[comp]] <- comp_draws
  }

  uni_draws <- rep(0, n_draws)
  for (i in seq_along(components)) {
    comp <- components[i]
    if (comp %in% universe_spec) {
      uni_draws <- uni_draws + composite_variance_draws[[comp]]
    }
  }

  sigma2_delta_draws <- rep(0, n_draws)
  for (i in which(error_info$is_relative_error)) {
    comp <- components[i]
    sigma2_delta_draws <- sigma2_delta_draws + composite_variance_draws[[comp]]
  }

  sigma2_delta_abs_draws <- rep(0, n_draws)
  for (i in which(error_info$is_absolute_error)) {
    comp <- components[i]
    sigma2_delta_abs_draws <- sigma2_delta_abs_draws + composite_variance_draws[[comp]]
  }

  g_draws <- uni_draws / (uni_draws + sigma2_delta_draws)
  g_draws[is.nan(g_draws) | is.infinite(g_draws)] <- NA

  phi_draws <- uni_draws / (uni_draws + sigma2_delta_abs_draws)
  phi_draws[is.nan(phi_draws) | is.infinite(phi_draws)] <- NA

  dim_g_draws <- list()
  dim_phi_draws <- list()

  for (d in dimensions) {
    uni_d <- rep(0, n_draws)
    err_rel_d <- rep(0, n_draws)
    err_abs_d <- rep(0, n_draws)

    for (comp in components) {
      v <- scaled_vc_draws[[d]][[comp]]

      if (comp %in% universe_spec) {
        uni_d <- uni_d + v
      } else {
        err_abs_d <- err_abs_d + v
        is_rel <- comp == "Residual" || grepl("^Residual", comp) || {
          facets <- parse_component_facets(comp)
          any(universe_spec %in% facets)
        }
        if (is_rel) err_rel_d <- err_rel_d + v
      }
    }

    denom_rel <- uni_d + err_rel_d
    g_d <- ifelse(denom_rel > 0, uni_d / denom_rel, NA_real_)
    dim_g_draws[[d]] <- g_d

    denom_abs <- uni_d + err_abs_d
    phi_d <- ifelse(denom_abs > 0, uni_d / denom_abs, NA_real_)
    dim_phi_draws[[d]] <- phi_d
  }

  subscale_g_matrix <- do.call(cbind, lapply(dimensions, function(d) dim_g_draws[[d]]))
  colnames(subscale_g_matrix) <- dimensions

  subscale_phi_matrix <- do.call(cbind, lapply(dimensions, function(d) dim_phi_draws[[d]]))
  colnames(subscale_phi_matrix) <- dimensions

  obs_cov_info <- compute_observed_covariances(
    data = gstudy_obj$data,
    dimensions = dimensions,
    weights = weights,
    dimension_var = gstudy_obj$dimension_var
  )

  var_result <- compute_var_haberman_draws(
    subscale_g_draws = subscale_g_matrix,
    subscale_phi_draws = subscale_phi_matrix,
    composite_g_draws = g_draws,
    composite_phi_draws = phi_draws,
    obs_var = obs_cov_info$obs_var,
    obs_cov = obs_cov_info$obs_cov,
    obs_var_C = obs_cov_info$obs_var_C,
    dimensions = dimensions
  )

  base_cols <- list(
    dim = dimensions,
    prmse_c_rel = colMeans(var_result$prmse_c_rel, na.rm = TRUE),
    prmse_c_abs = colMeans(var_result$prmse_c_abs, na.rm = TRUE),
    var_rel = colMeans(var_result$var_rel, na.rm = TRUE),
    var_abs = colMeans(var_result$var_abs, na.rm = TRUE)
  )

  if (!is.null(ci)) {
    if ("prmse" %in% ci) {
      base_cols$prmse_c_rel_LL <- apply(var_result$prmse_c_rel, 2, quantile, probs = probs[1], na.rm = TRUE)
      base_cols$prmse_c_rel_UL <- apply(var_result$prmse_c_rel, 2, quantile, probs = probs[2], na.rm = TRUE)
      base_cols$prmse_c_abs_LL <- apply(var_result$prmse_c_abs, 2, quantile, probs = probs[1], na.rm = TRUE)
      base_cols$prmse_c_abs_UL <- apply(var_result$prmse_c_abs, 2, quantile, probs = probs[2], na.rm = TRUE)
    }
    if ("var" %in% ci) {
      base_cols$var_rel_LL <- apply(var_result$var_rel, 2, quantile, probs = probs[1], na.rm = TRUE)
      base_cols$var_rel_UL <- apply(var_result$var_rel, 2, quantile, probs = probs[2], na.rm = TRUE)
      base_cols$var_abs_LL <- apply(var_result$var_abs, 2, quantile, probs = probs[1], na.rm = TRUE)
      base_cols$var_abs_UL <- apply(var_result$var_abs, 2, quantile, probs = probs[2], na.rm = TRUE)
    }
  }

  tibble::as_tibble(base_cols)
}

#' Compute Scale Factors for Viable Recalculation
#'
#' @param components Character vector of variance component names
#' @param n Named list of sample sizes
#' @param universe_spec Universe components specification
#' @param object_spec Object of measurement specification
#' @keywords internal
compute_scale_factors_for_viable <- function(components, n, universe_spec, object_spec) {
  scale_factors <- list()

  for (comp in components) {
    if (comp %in% object_spec) {
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
#' @keywords internal
extract_variance_draws_from_gstudy <- function(gstudy_obj) {

  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for posterior recalculation.", call. = FALSE)
  }

  draws <- brms::as_draws_matrix(gstudy_obj$model)

  extract_variance_draws(gstudy_obj, draws)
}
