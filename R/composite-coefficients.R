#' Composite Coefficient Functions
#'
#' Functions for computing composite coefficients in multivariate designs.
#' These functions handle weighted combinations of multiple dimensions.
#' @name composite-coefficients
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect all_of
NULL

#' Build Variance-Covariance Matrix for a Component
#'
#' Constructs a variance-covariance matrix for a specific variance component
#' in a multivariate design. The diagonal elements are the variances for each
#' dimension, and off-diagonal elements are covariances (when available).
#'
#' @param vc A variance components tibble with columns: component, dim, var
#' @param component Character string naming the component (e.g., "Person", "Residual")
#' @param dimensions Character vector of dimension names
#' @param correlations List containing covariance information from the G-study:
#'   - `residual_cov`: Tibble of residual covariances
#'   - `random_effect_cov`: Named list of random effect covariances per facet
#'
#' @return A symmetric variance-covariance matrix with dimensions as row/column names
#'
#' @keywords internal
build_component_covariance_matrix <- function(vc, component, dimensions, correlations = NULL, scale_factor = 1) {
  n_dim <- length(dimensions)
  cov_matrix <- matrix(0, n_dim, n_dim)
  rownames(cov_matrix) <- dimensions
  colnames(cov_matrix) <- dimensions

  for (i in seq_along(dimensions)) {
    dim_var <- vc$var[vc$component == component & vc$dim == dimensions[i]]
    if (length(dim_var) > 0) {
      cov_matrix[i, i] <- dim_var[1] / scale_factor
    }
  }

  if (!is.null(correlations) && component == "Residual") {
    if (!is.null(correlations$residual_cov) && nrow(correlations$residual_cov) > 0) {
      for (i in seq_len(nrow(correlations$residual_cov))) {
        row <- correlations$residual_cov[i, ]
        idx1 <- match(row$dim1, dimensions)
        idx2 <- match(row$dim2, dimensions)
        if (!is.na(idx1) && !is.na(idx2)) {
          # Scale the covariance by the scale_factor
          scaled_cov <- row$estimate / scale_factor
          cov_matrix[idx1, idx2] <- scaled_cov
          cov_matrix[idx2, idx1] <- scaled_cov
        }
      }
    }
  } else if (!is.null(correlations) && !is.null(correlations$random_effect_cov)) {
    facet_name <- component
    if (facet_name %in% names(correlations$random_effect_cov)) {
      facet_cov <- correlations$random_effect_cov[[facet_name]]
      if (!is.null(facet_cov) && nrow(facet_cov) > 0) {
        for (i in seq_len(nrow(facet_cov))) {
          row <- facet_cov[i, ]
          idx1 <- match(row$dim1, dimensions)
          idx2 <- match(row$dim2, dimensions)
          if (!is.na(idx1) && !is.na(idx2)) {
            # Scale the covariance by the scale_factor
            scaled_cov <- row$estimate / scale_factor
            cov_matrix[idx1, idx2] <- scaled_cov
            cov_matrix[idx2, idx1] <- scaled_cov
          }
        }
      }
    }
  }

  cov_matrix
}

#' Compute Weighted Variance from Covariance Matrix
#'
#' Computes the variance of a weighted composite using the quadratic form
#' w' * Sigma * w, where w is the weight vector and Sigma is the
#' variance-covariance matrix.
#'
#' @param cov_matrix A symmetric variance-covariance matrix
#' @param weights Numeric vector of weights (one per dimension)
#'
#' @return The weighted variance as a numeric scalar
#'
#' @keywords internal
compute_weighted_variance <- function(cov_matrix, weights) {
  as.numeric(t(weights) %*% cov_matrix %*% weights)
}

#' Calculate Composite G and Phi Coefficients
#'
#' Computes composite reliability coefficients for multivariate designs by
#' combining variance-covariance information across dimensions using weights.
#'
#' @param vc A variance components tibble (scaled for D-study)
#' @param n Named list of sample sizes for each facet
#' @param weights Named numeric vector of weights (names = dimensions)
#' @param object The object of measurement specification
#' @param error Optional error component specification
#' @param aggregation Optional aggregation specification
#' @param residual_is Optional residual specification
#' @param universe Optional universe specification
#' @param correlations List of covariance information from the G-study
#' @param cut_score Optional cut-score for phi-cut calculation
#' @param mu_y Optional grand mean(s) for phi-cut calculation
#'
#' @return A data.frame with one row containing:
#'   - dim: "Composite"
#'   - uni: Universe score variance (composite)
#'   - sigma2_delta: Relative error variance (composite)
#'   - sigma2_delta_abs: Absolute error variance (composite)
#'   - g: G coefficient (composite)
#'   - phi: Phi coefficient (composite)
#'   - sem_rel: Standard error of measurement (relative)
#'   - sem_abs: Standard error of measurement (absolute)
#'   - phi_cut: Phi-cut coefficient (if cut_score provided)
#'
#' @keywords internal
calculate_composite_coefficients <- function(vc, n, weights, object, error = NULL, aggregation = NULL, residual_is = NULL, universe = NULL, correlations = NULL, cut_score = NULL, mu_y = NULL, gstudy_data = NULL, dimension_var = NULL) {
  dimensions <- names(weights)

  if (is.null(object)) {
    object_spec <- vc$component[1]
  } else {
    object_spec <- parse_specification(object)
  }

  if (is.null(universe) || length(universe) == 0) {
    universe_spec <- object_spec
  } else {
    universe_spec <- parse_specification(universe)
  }

  if (!is.null(error)) {
    error_spec <- parse_specification(error)
  } else {
    error_spec <- NULL
  }

  components <- unique(vc$component)

  # Universe score variance must use G-study (unscaled) variance components.
  # When vc contains d_vc (D-study scaled), we temporarily swap var/unscaled.
  vc_uni <- vc
  if ("var_unscaled" %in% names(vc)) {
    vc_uni$var <- vc$var_unscaled
  }

  uni_composite <- 0
  for (comp in universe_spec) {
    if (comp %in% components) {
      cov_mat <- build_component_covariance_matrix(vc_uni, comp, dimensions, correlations, scale_factor = 1)
      uni_composite <- uni_composite + compute_weighted_variance(cov_mat, weights)
    }
  }

  sigma2_delta_composite <- compute_composite_error_variance(
    vc, dimensions, weights, correlations, universe_spec, error_spec, aggregation, residual_is, "relative", n = n, object_spec = object_spec
  )

  sigma2_delta_abs_composite <- compute_composite_error_variance(
    vc, dimensions, weights, correlations, universe_spec, error_spec, aggregation, residual_is, "absolute", n = n, object_spec = object_spec
  )

  g_composite <- uni_composite / (uni_composite + sigma2_delta_composite)
  phi_composite <- uni_composite / (uni_composite + sigma2_delta_abs_composite)
  sem_rel_composite <- sqrt(sigma2_delta_composite)
  sem_abs_composite <- sqrt(sigma2_delta_abs_composite)

  var_rel <- NA_real_
  var_abs <- NA_real_
  per_subscale_var <- NULL

  if ("dim" %in% names(vc) && length(dimensions) > 1) {
    # Universe score variance uses G-study (unscaled) variance
    uni_var_col <- if ("var_unscaled" %in% names(vc)) "var_unscaled" else "var"
    # Error variance uses D-study (scaled) variance
    err_var_col <- if ("var_scaled" %in% names(vc)) "var_scaled" else "var"

    dim_g <- vapply(dimensions, function(d) {
      vc_d <- vc[vc$dim == d, , drop = FALSE]
      uni_d <- sum(vc_d[[uni_var_col]][vc_d$component %in% universe_spec], na.rm = TRUE)
      rel_comps <- vc_d$component[sapply(vc_d$component, function(c) {
        if (c %in% universe_spec) return(FALSE)
        if (c == "Residual") return(TRUE)
        facets <- parse_component_facets(c)
        any(universe_spec %in% facets)
      })]
      err_rel_d <- sum(vc_d[[err_var_col]][vc_d$component %in% rel_comps], na.rm = TRUE)
      g_d <- if ((uni_d + err_rel_d) > 0) uni_d / (uni_d + err_rel_d) else NA_real_
      g_d
    }, numeric(1))
    names(dim_g) <- dimensions

    dim_phi <- vapply(dimensions, function(d) {
      vc_d <- vc[vc$dim == d, , drop = FALSE]
      uni_d <- sum(vc_d[[uni_var_col]][vc_d$component %in% universe_spec], na.rm = TRUE)
      all_comps <- vc_d$component[!(vc_d$component %in% universe_spec)]
      err_abs_d <- sum(vc_d[[err_var_col]][vc_d$component %in% all_comps], na.rm = TRUE)
      phi_d <- if ((uni_d + err_abs_d) > 0) uni_d / (uni_d + err_abs_d) else NA_real_
      phi_d
    }, numeric(1))
    names(dim_phi) <- dimensions

    # Build model-based covariance matrices from variance components
    # This ensures consistency with the draws-based path and Equation 38's
    # requirement for disattenuated (universe-score) covariances
    k <- length(dimensions)

    # Sigma_tau: universe-score covariance matrix from model VCs
    # Uses G-study (unscaled) variance - vc_uni created earlier has unscaled values
    Sigma_tau <- matrix(0, k, k, dimnames = list(dimensions, dimensions))
    for (comp in universe_spec) {
      if (comp %in% components) {
        cov_mat <- build_component_covariance_matrix(
          vc_uni, comp, dimensions, correlations, scale_factor = 1
        )
        Sigma_tau <- Sigma_tau + cov_mat
      }
    }

    # Sigma_obs_rel: model-based observed covariance matrix (relative decisions)
    # = universe + relative error components (scaled by D-study sample sizes)
    Sigma_obs_rel <- Sigma_tau
    # Sigma_obs_abs: model-based observed covariance matrix (absolute decisions)
    # = universe + all error components (scaled by D-study sample sizes)
    Sigma_obs_abs <- Sigma_tau

    for (comp in components) {
      if (comp %in% universe_spec) next
      # vc$var already contains D-study scaled variance, so use scale_factor = 1
      cov_mat <- build_component_covariance_matrix(
        vc, comp, dimensions, correlations, scale_factor = 1
      )
      # All non-universe components contribute to absolute error
      Sigma_obs_abs <- Sigma_obs_abs + cov_mat
      # Only relative error components contribute to relative error
      is_rel <- comp == "Residual" || {
        facets <- parse_component_facets(comp)
        any(universe_spec %in% facets)
      }
      if (is_rel) {
        Sigma_obs_rel <- Sigma_obs_rel + cov_mat
      }
    }

    # Compute PRMSE_P for all dimensions (profile projection)
    # PRMSE_P[d] = [Sigma_tau * solve(Sigma_obs) * Sigma_tau][d,d] / Sigma_tau[d,d]
    prmse_p_rel_all <- rep(NA_real_, k)
    prmse_p_abs_all <- rep(NA_real_, k)
    names(prmse_p_rel_all) <- dimensions
    names(prmse_p_abs_all) <- dimensions

    S_rel_inv <- tryCatch(solve(Sigma_obs_rel), error = function(e) NULL)
    S_abs_inv <- tryCatch(solve(Sigma_obs_abs), error = function(e) NULL)

    if (!is.null(S_rel_inv)) {
      T_S_rel_T <- Sigma_tau %*% S_rel_inv %*% Sigma_tau
      for (i in seq_along(dimensions)) {
        d <- dimensions[i]
        var_tau_d <- Sigma_tau[i, i]
        if (is.finite(T_S_rel_T[i, i]) && var_tau_d > 0) {
          prmse_p_rel_all[d] <- T_S_rel_T[i, i] / var_tau_d
        }
      }
    }

    if (!is.null(S_abs_inv)) {
      T_S_abs_T <- Sigma_tau %*% S_abs_inv %*% Sigma_tau
      for (i in seq_along(dimensions)) {
        d <- dimensions[i]
        var_tau_d <- Sigma_tau[i, i]
        if (is.finite(T_S_abs_T[i, i]) && var_tau_d > 0) {
          prmse_p_abs_all[d] <- T_S_abs_T[i, i] / var_tau_d
        }
      }
    }

    # Compute PRMSE_C for all dimensions using model-based Sigma_tau %*% w
    # PRMSE_C = [Cov(tau_d, C)]^2 / [Var(tau_d) * Rel(C) * Var(C)]
    # This is equivalent to Equation 38 when expanded, using universe-score
    # (disattenuated) covariances as the paper specifies.
    w <- weights[dimensions]
    cov_tau_C <- as.numeric(Sigma_tau %*% w)  # Cov(tau_d, C) for each d
    var_C_rel <- as.numeric(t(w) %*% Sigma_obs_rel %*% w)
    var_C_abs <- as.numeric(t(w) %*% Sigma_obs_abs %*% w)
    var_tau_C <- as.numeric(t(w) %*% Sigma_tau %*% w)
    Rel_C_rel <- if (var_C_rel > 0) var_tau_C / var_C_rel else NA_real_
    Rel_C_abs <- if (var_C_abs > 0) var_tau_C / var_C_abs else NA_real_

    per_subscale_var <- lapply(dimensions, function(d) {
      d_idx <- which(dimensions == d)
      var_tau_d <- Sigma_tau[d_idx, d_idx]

      # PRMSE_C numerator: (Cov(tau_d, C))^2
      # where Cov(tau_d, C) = [Sigma_tau %*% w][d]
      #   = sigma^2_tau_d * w_d + sum_{j != d} sigma_tau_{d,j} * w_j
      # This matches Equation 38's numerator using disattenuated covariances
      num_c <- cov_tau_C[d_idx]^2

      # PRMSE_C denominator: Var(tau_d) * Rel(C) * Var(C)
      den_rel <- var_tau_d * Rel_C_rel * var_C_rel
      den_abs <- var_tau_d * Rel_C_abs * var_C_abs

      prmse_c_rel <- if (!is.na(den_rel) && is.finite(den_rel) && den_rel > 0) {
        num_c / den_rel
      } else NA_real_

      prmse_c_abs <- if (!is.na(den_abs) && is.finite(den_abs) && den_abs > 0) {
        num_c / den_abs
      } else NA_real_

      var_rel_d <- if (!is.na(prmse_c_rel) && prmse_c_rel > 0) dim_g[d] / prmse_c_rel else NA_real_
      var_abs_d <- if (!is.na(prmse_c_abs) && prmse_c_abs > 0) dim_phi[d] / prmse_c_abs else NA_real_

      list(
        prmse_c_rel = prmse_c_rel,
        prmse_c_abs = prmse_c_abs,
        prmse_p_rel = prmse_p_rel_all[d],
        prmse_p_abs = prmse_p_abs_all[d],
        var_rel = var_rel_d,
        var_abs = var_abs_d
      )
    })
    names(per_subscale_var) <- dimensions
  }

  result <- data.frame(
    dim = "Composite",
    uni = uni_composite,
    sigma2_delta = sigma2_delta_composite,
    sigma2_delta_abs = sigma2_delta_abs_composite,
    g = g_composite,
    phi = phi_composite,
    sem_rel = sem_rel_composite,
    sem_abs = sem_abs_composite,
    stringsAsFactors = FALSE
  )

  if (!is.null(cut_score) && !is.null(mu_y)) {
    mu_y_weighted <- sum(weights * sapply(dimensions, function(d) {
      if (is.list(mu_y)) mu_y[[d]] else mu_y
    }))
    adjustment <- (mu_y_weighted - cut_score)^2
    result$phi_cut <- (uni_composite + adjustment) / (uni_composite + sigma2_delta_abs_composite + adjustment)
  }

  list(
    summary = result,
    var_results = per_subscale_var
  )
}


#' Calculate Composite Variance Draws
#'
#' Computes posterior draws of composite variance for each component type
#' using weighted variance-covariance matrices.
#'
#' @param vc_draws Named list (by dimension) of named lists (by component) of variance draws
#' @param cov_draws Named list from extract_covariance_draws()
#' @param weights Named numeric vector of weights (names = dimensions)
#' @param scale_factor Factor to scale variances and covariances (for D-study scaling)
#'
#' @return Named list (by component) of posterior draws for composite variance
#'
#' @keywords internal
calculate_composite_variance_draws <- function(vc_draws, cov_draws, weights, scale_factor = 1) {
  # Input validation
  if (!is.list(vc_draws) || length(vc_draws) == 0) {
    stop("vc_draws must be a non-empty list", call. = FALSE)
  }
  if (!is.numeric(weights) || length(weights) == 0) {
    stop("weights must be a numeric vector", call. = FALSE)
  }

  dimensions <- names(vc_draws)
  n_dim <- length(dimensions)

  if (!all(names(weights) %in% dimensions)) {
    stop("weights names must match dimensions in vc_draws", call. = FALSE)
  }

  # Check vc_draws structure
  first_dim <- vc_draws[[1]]
  if (!is.list(first_dim) || length(first_dim) == 0) {
    stop("vc_draws must be a nested list structure", call. = FALSE)
  }

  n_draws <- length(first_dim[[1]])

  # Get all component names (from first dimension)
  components <- names(vc_draws[[1]])

  composite_draws <- list()

  for (comp in components) {
    comp_draws <- numeric(n_draws)

    # Get variance draws for each dimension
    var_draws_list <- list()
    for (d in dimensions) {
      var_draws_list[[d]] <- vc_draws[[d]][[comp]] / scale_factor
    }

    # For each posterior draw, compute w' * Sigma * w
    for (draw_idx in 1:n_draws) {
      # Build covariance matrix for this draw
      sigma_mat <- matrix(0, n_dim, n_dim)
      rownames(sigma_mat) <- dimensions
      colnames(sigma_mat) <- dimensions

      for (i in 1:n_dim) {
        sigma_mat[i, i] <- var_draws_list[[dimensions[i]]][draw_idx]
      }

      # Add off-diagonal covariances
      cov_list <- NULL
      if (comp == "Residual") {
        cov_list <- cov_draws$residual
      } else if (comp %in% names(cov_draws$random_effect)) {
        cov_list <- cov_draws$random_effect[[comp]]
      }

      if (!is.null(cov_list)) {
        for (pair_name in names(cov_list)) {
          dims_pair <- strsplit(pair_name, "_")[[1]]
          if (length(dims_pair) == 2 && all(dims_pair %in% dimensions)) {
            idx1 <- match(dims_pair[1], dimensions)
            idx2 <- match(dims_pair[2], dimensions)
            cov_val <- cov_list[[pair_name]][draw_idx] / scale_factor
            sigma_mat[idx1, idx2] <- cov_val
            sigma_mat[idx2, idx1] <- cov_val
          }
        }
      }

      # Compute weighted variance: t(weights) %*% sigma_mat %*% weights
      comp_draws[draw_idx] <- as.numeric(t(weights) %*% sigma_mat %*% weights)
    }

    composite_draws[[comp]] <- comp_draws
  }

  composite_draws
}

#' Calculate D-Study Variance Components with Composite Row
#'
#' Adds composite variance component rows for multivariate models using
#' posterior draws for full uncertainty propagation.
#'
#' @param vc_draws Variance draws (from extract_variance_draws)
#' @param cov_draws Covariance draws (from extract_covariance_draws)
#' @param weights Named vector of weights
#' @param n Sample sizes (named list)
#' @param object Object specification
#' @param n_provided Whether n was provided
#'
#' @return List with:
#'   \describe{
#'     \item{vc_table}{Tibble with variance components including composite rows}
#'     \item{posterior}{List of composite draws}
#'   }
#'
#' @keywords internal
calculate_dstudy_variance_composite <- function(vc_draws, cov_draws, weights, n, object, n_provided) {
  dimensions <- names(weights)
  components <- names(vc_draws[[1]])

  object_spec <- parse_specification(object)

  composite_draws <- list()
  composite_summaries <- list()

  for (comp in components) {
    scale_factor <- compute_component_scale_factor(
      comp, n, object_spec, n_provided
    )

    comp_draws <- calculate_composite_variance_draws(
      vc_draws, cov_draws, weights, scale_factor
    )

    composite_draws[[comp]] <- comp_draws[[comp]]

    var_scaled_mean <- mean(comp_draws[[comp]], na.rm = TRUE)
    var_se <- sd(comp_draws[[comp]], na.rm = TRUE)

    # Unscaled = scaled * scale_factor (reverse the D-study scaling)
    var_unscaled_mean <- var_scaled_mean * scale_factor

    composite_summaries[[comp]] <- list(
      var_scaled = var_scaled_mean,
      var_unscaled = var_unscaled_mean,
      se = var_se
    )
  }

  total_var_scaled <- sum(sapply(composite_summaries, function(x) x$var_scaled))
  total_var_unscaled <- sum(sapply(composite_summaries, function(x) x$var_unscaled))

  composite_rows <- lapply(components, function(comp) {
    data.frame(
      component = comp,
      dim = "Composite",
      var_unscaled = composite_summaries[[comp]]$var_unscaled,
      pct_unscaled = (composite_summaries[[comp]]$var_unscaled / total_var_unscaled) * 100,
      var_scaled = composite_summaries[[comp]]$var_scaled,
      pct_scaled = (composite_summaries[[comp]]$var_scaled / total_var_scaled) * 100,
      stringsAsFactors = FALSE
    )
  })

  list(
    vc_table = do.call(rbind, composite_rows),
    posterior = composite_draws
  )
}

#' Compute Composite Error Variance
#'
#' Calculates the composite relative or absolute error variance by summing
#' weighted variance-covariance matrices for all error components.
#'
#' @param vc A variance components tibble
#' @param dimensions Character vector of dimension names
#' @param weights Named numeric vector of weights
#' @param correlations List of covariance information from the G-study
#' @param universe_spec Character vector of universe components
#' @param error_spec Optional character vector of error components
#' @param aggregation Optional aggregation specification
#' @param residual_is Optional residual specification
#' @param error_type Either "relative" or "absolute"
#' @param scale_factor Factor by which to scale covariances
#'
#' @return The composite error variance as a numeric scalar
#'
#' @keywords internal
compute_composite_error_variance <- function(vc, dimensions, weights, correlations, universe_spec, error_spec, aggregation, residual_is, error_type, n = NULL, object_spec = NULL) {
  components <- unique(vc$component)

  if (!is.null(error_spec)) {
    error_components <- error_spec
  } else {
    if (error_type == "relative") {
      error_components <- vc$component[
        sapply(vc$component, function(comp) {
          if (comp == "Residual") return(TRUE)
          if (comp %in% universe_spec) return(FALSE)
          facets <- parse_component_facets(comp)
          any(universe_spec %in% facets)
        })
      ]
    } else {
      error_components <- setdiff(components, universe_spec)
    }
  }

  total_error <- 0
  for (comp in error_components) {
    if (comp %in% components) {
      # vc$var already contains D-study scaled variance, so use scale_factor = 1
      cov_mat <- build_component_covariance_matrix(vc, comp, dimensions, correlations, scale_factor = 1)
      total_error <- total_error + compute_weighted_variance(cov_mat, weights)
    }
  }

  total_error
}

#' Calculate Composite Coefficients from Posterior Draws
#'
#' Computes composite reliability coefficients from posterior draws of
#' variance components and covariances. This ensures proper uncertainty
#' propagation for multivariate D-studies with brms backend.
#'
#' @param vc_draws Named list (by dimension) of named lists (by component) of variance draws
#' @param cov_draws Named list from extract_covariance_draws()
#' @param weights Named numeric vector of weights (names = dimensions)
#' @param scale_factors Named list of scale factors per component
#' @param object_spec Parsed object specification
#' @param universe_spec Parsed universe specification
#' @param error_spec Parsed error specification (or NULL)
#' @param agg_facets Parsed aggregation specification (or NULL)
#' @param residual_is Residual composition specification
#' @param cut_score Optional cut score for phi-cut
#' @param mu_y Named list of grand means (or single value)
#' @param ci Optional CI specification (character vector)
#' @param probs Probability levels for credible intervals
#'
#' @return List with:
#'   \describe{
#'     \item{summary}{Data frame with coefficient point estimates}
#'     \item{distributions}{Named list of posterior draws for each coefficient}
#'   }
#'
#' @keywords internal
calculate_composite_coefficients_from_draws <- function(vc_draws, cov_draws,
  weights, scale_factors,
  object_spec, universe_spec,
  error_spec, agg_facets,
  residual_is, cut_score = NULL,
  mu_y = NULL, ci = NULL,
  probs = c(0.025, 0.975),
  gstudy_data = NULL,
  dimension_var = NULL) {
  dimensions <- names(weights)
  n_draws <- length(vc_draws[[1]][[1]])
  components <- unique(unlist(lapply(vc_draws, names)))

  error_info <- identify_error_components_for_draws(components, universe_spec, error_spec, agg_facets)

  composite_variance_draws <- list()
  for (comp in components) {
    scale_factor <- scale_factors[[comp]]
    composite_variance_draws[[comp]] <- calculate_composite_variance_draws(
      vc_draws, cov_draws, weights, scale_factor
    )[[comp]]
  }

  # Universe score variance must use G-study (unscaled) variance.
  # composite_variance_draws[[comp]] = t(w) * (vc_draws/sf) * w, already scaled.
  # To get unscaled: multiply by scale_factor to undo the division.
  uni_draws <- rep(0, n_draws)
  for (comp in universe_spec) {
    if (comp %in% components) {
      scale_factor <- scale_factors[[comp]]
      uni_draws <- uni_draws + composite_variance_draws[[comp]] * scale_factor
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
  phi_draws <- uni_draws / (uni_draws + sigma2_delta_abs_draws)
  sem_rel_draws <- sqrt(sigma2_delta_draws)
  sem_abs_draws <- sqrt(sigma2_delta_abs_draws)

  g_draws[is.nan(g_draws) | is.infinite(g_draws)] <- NA
  phi_draws[is.nan(phi_draws) | is.infinite(phi_draws)] <- NA

  # Compute per-dimension G/Phi draws for VAR calculation using Haberman formula
  # Each dim_g_draws[[d]] is a list with g, uni, err_abs
  var_rel_draws <- NULL
  var_abs_draws <- NULL
  var_rel_draws_matrix <- NULL
  var_abs_draws_matrix <- NULL

  if (length(dimensions) > 1 && !is.null(gstudy_data)) {
    dim_g_draws <- lapply(dimensions, function(d) {
      uni_d <- rep(0, n_draws)
      err_rel_d <- rep(0, n_draws)
      err_abs_d <- rep(0, n_draws)
      for (comp in components) {
        sf <- scale_factors[[comp]]
        if (comp %in% universe_spec) {
          # Universe components: use G-study (unscaled) variance
          uni_d <- uni_d + vc_draws[[d]][[comp]]
        } else {
          # Error components: use D-study (scaled) variance
          v <- vc_draws[[d]][[comp]] / sf
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
      denom_abs <- uni_d + err_abs_d
      phi_d <- ifelse(denom_abs > 0, uni_d / denom_abs, NA_real_)
      list(g = g_d, phi = phi_d, uni = uni_d, err_abs = err_abs_d)
    })
    names(dim_g_draws) <- dimensions

    # Construct covariance arrays for Haberman formula
    # uni_cov_draws: universe-score covariance (G-study, UNSCALED)
    # total_rel_cov_draws: observed covariance (D-study, SCALED)
    # total_abs_cov_draws: observed covariance (D-study, SCALED)
    uni_cov_draws <- array(0, dim = c(n_draws, length(dimensions), length(dimensions)))
    total_rel_cov_draws <- array(0, dim = c(n_draws, length(dimensions), length(dimensions)))
    total_abs_cov_draws <- array(0, dim = c(n_draws, length(dimensions), length(dimensions)))

    for (i in seq_along(dimensions)) {
      d1 <- dimensions[i]
      for (j in seq_along(dimensions)) {
        d2 <- dimensions[j]

        for (comp in components) {
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

          if (error_info$is_relative_error[which(components == comp)] || comp %in% universe_spec) {
            total_rel_cov_draws[, i, j] <- total_rel_cov_draws[, i, j] + val_obs
          }
          if (error_info$is_absolute_error[which(components == comp)] || comp %in% universe_spec) {
            total_abs_cov_draws[, i, j] <- total_abs_cov_draws[, i, j] + val_obs
          }
        }
      }
    }

    var_result <- compute_var_haberman_draws(
      uni_cov_draws = uni_cov_draws,
      total_rel_cov_draws = total_rel_cov_draws,
      total_abs_cov_draws = total_abs_cov_draws,
      dimensions = dimensions,
      weights = weights
    )

    var_rel_draws_matrix <- var_result$var_rel
    var_abs_draws_matrix <- var_result$var_abs

    var_rel_draws_matrix <- var_result$var_rel
    var_abs_draws_matrix <- var_result$var_abs

    # Summary values (means) for Composite row - these will be NA since VAR is per-subscale
    var_rel_draws <- NA_real_
    var_abs_draws <- NA_real_
  } else if (length(dimensions) > 1) {
    # No data available - cannot compute VAR
    var_rel_draws <- NA_real_
    var_abs_draws <- NA_real_
  }

  phi_cut_draws <- NULL
  if (!is.null(cut_score) && !is.null(mu_y)) {
    mu_y_weighted <- sum(weights * sapply(dimensions, function(d) {
      if (is.list(mu_y)) mu_y[[d]] else mu_y
    }))
    adjustment <- (mu_y_weighted - cut_score)^2
    phi_cut_draws <- (uni_draws + adjustment) / (uni_draws + sigma2_delta_abs_draws + adjustment)
    phi_cut_draws[is.nan(phi_cut_draws) | is.infinite(phi_cut_draws)] <- NA
  }

  summary_df <- data.frame(
    dim = "Composite",
    uni = mean(uni_draws, na.rm = TRUE),
    sigma2_delta = mean(sigma2_delta_draws, na.rm = TRUE),
    sigma2_delta_abs = mean(sigma2_delta_abs_draws, na.rm = TRUE),
    g = mean(g_draws, na.rm = TRUE),
    phi = mean(phi_draws, na.rm = TRUE),
    sem_rel = mean(sem_rel_draws, na.rm = TRUE),
    sem_abs = mean(sem_abs_draws, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  if (!is.null(cut_score)) {
    summary_df$phi_cut <- mean(phi_cut_draws, na.rm = TRUE)
  }

  if (!is.null(ci)) {
    if ("g" %in% ci) {
      summary_df$g_LL <- quantile(g_draws, probs[1], na.rm = TRUE)
      summary_df$g_UL <- quantile(g_draws, probs[2], na.rm = TRUE)
    }
    if ("phi" %in% ci) {
      summary_df$phi_LL <- quantile(phi_draws, probs[1], na.rm = TRUE)
      summary_df$phi_UL <- quantile(phi_draws, probs[2], na.rm = TRUE)
    }
    if ("phi-cut" %in% ci && !is.null(cut_score)) {
      summary_df$phi_cut_LL <- quantile(phi_cut_draws, probs[1], na.rm = TRUE)
      summary_df$phi_cut_UL <- quantile(phi_cut_draws, probs[2], na.rm = TRUE)
    }
  }

  distributions <- list(
    uni = uni_draws,
    sigma2_delta = sigma2_delta_draws,
    sigma2_delta_abs = sigma2_delta_abs_draws,
    g = g_draws,
    phi = phi_draws,
    sem_rel = sem_rel_draws,
    sem_abs = sem_abs_draws
  )

  if (!is.null(cut_score)) {
    distributions$phi_cut <- phi_cut_draws
  }

  var_results <- NULL
  if (!is.null(var_rel_draws_matrix)) {
    var_results <- list(
      prmse_s_rel_draws = var_result$prmse_s_rel,
      prmse_s_abs_draws = var_result$prmse_s_abs,
      prmse_c_rel_draws = var_result$prmse_c_rel,
      prmse_c_abs_draws = var_result$prmse_c_abs,
      prmse_p_rel_draws = var_result$prmse_p_rel,
      prmse_p_abs_draws = var_result$prmse_p_abs,
      var_rel_draws = var_rel_draws_matrix,
      var_abs_draws = var_abs_draws_matrix,
      prmse_mv_rel_draws = var_result$prmse_mv_rel,
      prmse_mv_abs_draws = var_result$prmse_mv_abs,
      prmse_s_rel = colMeans(var_result$prmse_s_rel, na.rm = TRUE),
      prmse_s_abs = colMeans(var_result$prmse_s_abs, na.rm = TRUE),
      prmse_c_rel = colMeans(var_result$prmse_c_rel, na.rm = TRUE),
      prmse_c_abs = colMeans(var_result$prmse_c_abs, na.rm = TRUE),
      prmse_p_rel = colMeans(var_result$prmse_p_rel, na.rm = TRUE),
      prmse_p_abs = colMeans(var_result$prmse_p_abs, na.rm = TRUE),
      var_rel = colMeans(var_rel_draws_matrix, na.rm = TRUE),
      var_abs = colMeans(var_abs_draws_matrix, na.rm = TRUE),
      prmse_mv_rel = mean(var_result$prmse_mv_rel, na.rm = TRUE),
      prmse_mv_abs = mean(var_result$prmse_mv_abs, na.rm = TRUE),
      prmse_s_rel_LL = apply(var_result$prmse_s_rel, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_s_rel_UL = apply(var_result$prmse_s_rel, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_s_abs_LL = apply(var_result$prmse_s_abs, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_s_abs_UL = apply(var_result$prmse_s_abs, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_c_rel_LL = apply(var_result$prmse_c_rel, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_c_rel_UL = apply(var_result$prmse_c_rel, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_c_abs_LL = apply(var_result$prmse_c_abs, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_c_abs_UL = apply(var_result$prmse_c_abs, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_p_rel_LL = apply(var_result$prmse_p_rel, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_p_rel_UL = apply(var_result$prmse_p_rel, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_p_abs_LL = apply(var_result$prmse_p_abs, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_p_abs_UL = apply(var_result$prmse_p_abs, 2, quantile, probs = probs[2], na.rm = TRUE),
      var_rel_LL = apply(var_rel_draws_matrix, 2, quantile, probs = probs[1], na.rm = TRUE),
      var_rel_UL = apply(var_rel_draws_matrix, 2, quantile, probs = probs[2], na.rm = TRUE),
      var_abs_LL = apply(var_abs_draws_matrix, 2, quantile, probs = probs[1], na.rm = TRUE),
      var_abs_UL = apply(var_abs_draws_matrix, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_mv_rel_LL = quantile(var_result$prmse_mv_rel, probs = probs[1], na.rm = TRUE),
      prmse_mv_rel_UL = quantile(var_result$prmse_mv_rel, probs = probs[2], na.rm = TRUE),
      prmse_mv_abs_LL = quantile(var_result$prmse_mv_abs, probs = probs[1], na.rm = TRUE),
    prmse_mv_abs_UL = quantile(var_result$prmse_mv_abs, probs = probs[2], na.rm = TRUE)
    )
  }

  list(
    summary = summary_df,
    distributions = distributions,
    var_results = var_results
  )
}


#' Compute Observed Score Variances and Covariances for VAR Calculation
#'
#' Computes observed score variances and covariances from the original data,
#' which are required for the Haberman (2008) formula in VAR calculation.
#' These are data quantities (not model estimates) and are constant across
#' posterior draws.
#'
#' @param data A data frame containing the subscale scores.
#' @param dimensions Character vector of dimension/subscale names.
#' @param weights Named numeric vector of weights for each dimension.
#'
#' @return A list with:
#' \describe{
#' \item{obs_var}{Named vector of observed score variances per dimension}
#' \item{obs_cov}{Matrix of observed score covariances between dimensions}
#' \item{obs_var_C}{Observed variance of the weighted composite}
#' }
#'
#' @keywords internal
compute_observed_covariances <- function(data, dimensions, weights, dimension_var = NULL) {
  # Handle long-format data: reshape to wide format
  if (!is.null(dimension_var) && dimension_var %in% names(data)) {
    # Find the score column (first numeric column that isn't dimension_var)
    score_col <- NULL
    for (col in names(data)) {
      if (col != dimension_var && is.numeric(data[[col]])) {
        score_col <- col
        break
      }
    }

    if (!is.null(score_col)) {
      # Reshape from long to wide format
      data <- tidyr::pivot_wider(
        data,
        id_cols = setdiff(names(data), c(dimension_var, score_col)),
        names_from = all_of(dimension_var),
        values_from = all_of(score_col)
      )
    }
  }

  n_dims <- length(dimensions)
  obs_var <- vapply(dimensions, function(d) {
    var(data[[d]], na.rm = TRUE)
  }, numeric(1))
  names(obs_var) <- dimensions

  obs_cov <- matrix(0, nrow = n_dims, ncol = n_dims)
  rownames(obs_cov) <- colnames(obs_cov) <- dimensions
  for (i in seq_len(n_dims)) {
    for (j in seq_len(n_dims)) {
      if (i == j) {
        obs_cov[i, j] <- obs_var[i]
      } else {
        obs_cov[i, j] <- cov(data[[dimensions[i]]], data[[dimensions[j]]], use = "complete.obs")
      }
    }
  }

  w <- weights[dimensions]
  obs_var_C <- as.numeric(t(w) %*% obs_cov %*% w)

  list(
    obs_var = obs_var,
    obs_cov = obs_cov,
    obs_var_C = obs_var_C
  )
}

