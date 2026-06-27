#' Posterior Coefficient Functions
#'
#' Functions for computing coefficients from posterior draws (Bayesian/aov).
#' Used when gstudy estimator is 'brms' or 'mom'.
#' @name coefficients-posterior
#' @keywords internal
NULL

#' Calculate Coefficients Using Posterior Draws
#'
#' Computes generalizability and dependability coefficients using full posterior
#' distributions from a brms model. This provides proper uncertainty quantification
#' for D-study coefficients.
#'
#' @param gstudy_obj A gstudy object fitted with estimator = "brms".
#' @param n Named list of sample sizes for each facet.
#' @param object Specification for object of measurement.
#' @param error Specification for error components.
#'   Note: If the same facet is specified in both `object` and `error`, an error is raised.
#' @param aggregation Character vector of facets to aggregate over.
#'   Note: If the same facet is specified in both `object` and `aggregation`, an error is raised.
#' @param residual_is Character string specifying residual composition.
#' @param is_sweep Logical indicating if this is a sweep over multiple sample sizes.
#' @param n_grid If is_sweep, a data frame with all combinations of sample sizes.
#' @param cut_score Optional cut score for phi-cut calculation.
#' @param mu_y Optional grand mean for phi-cut calculation.
#'
#' @return A list with:
#'   \item{coefficients}{Data frame with coefficient means (and sweep combinations if applicable)}
#'   \item{posterior}{List of posterior distribution vectors (or list of lists for sweep)}
#'
#' @keywords internal
calculate_coefficients_posterior <- function(gstudy_obj, n, object = NULL, universe = NULL, error = NULL,
                                             aggregation = NULL, residual_is = NULL,
                                             is_sweep = FALSE, n_grid = NULL, n_provided = FALSE, use_scaled = TRUE,
                                             cut_score = NULL, mu_y = NULL, ci = NULL, probs = c(0.025, 0.975),
                                             weights = NULL, n_per_dim = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for posterior estimation.", call. = FALSE)
  }

  # Extract posterior draws
  draws <- brms::as_draws_matrix(gstudy_obj$model)

  # Get variance component names and their draws
  vc_draws <- extract_variance_draws(gstudy_obj, draws)

  # Extract covariance draws for multivariate composite calculation
  cov_draws <- NULL
  if (gstudy_obj$is_multivariate && length(gstudy_obj$dimensions) > 1) {
    cov_draws <- extract_covariance_draws(
      gstudy_obj$model,
      gstudy_obj$dimensions
    )
  }

  # Get object of measurement (always from G-study)
  if (is.null(object)) {
    object <- gstudy_obj$object
  }
  object_spec <- parse_specification(object)

  # Process universe specification
  # Default: universe = object only
  if (is.null(universe) || length(universe) == 0) {
    universe_spec <- object_spec
  } else {
    universe_spec <- parse_specification(universe)
  }

  # Parse error specification
  if (!is.null(error)) {
    error_spec <- parse_specification(error)
  } else {
    error_spec <- NULL
  }

  # Parse aggregation specification
  if (!is.null(aggregation)) {
    agg_facets <- parse_specification(aggregation)
  } else {
    agg_facets <- NULL
  }

  # Validate universe vs error
  if (!is.null(error_spec) && length(error_spec) > 0) {
    overlap <- intersect(universe_spec, error_spec)
    if (length(overlap) > 0) {
      stop(
        "Component(s) '", paste(overlap, collapse = "', '"), "' specified as both ",
        "in the universe and in error. ",
        "The same component cannot be in both the universe and error.",
        call. = FALSE
      )
    }
  }

  is_mv <- gstudy_obj$is_multivariate
  mv_dims <- if (is_mv) gstudy_obj$dimensions else character(0)

  if (is_sweep) {
    # Calculate for each sweep combination
    results <- vector("list", nrow(n_grid))
    posterior_list <- vector("list", nrow(n_grid))

    for (i in seq_len(nrow(n_grid))) {
      n_current <- as.list(n_grid[i, , drop = FALSE])
      n_current <- lapply(n_current, as.numeric)

      if (is_mv) {
        # Multivariate: iterate over dimensions
        dim_results <- vector("list", length(mv_dims))
        names(dim_results) <- mv_dims
        dim_posteriors <- vector("list", length(mv_dims))
        names(dim_posteriors) <- mv_dims

 for (d in mv_dims) {
  # Use per-dimension sample sizes when available
  n_dim <- n_current
  if (!is.null(n_per_dim) && !is.null(n_per_dim$n_list[[d]])) {
    dim_n <- n_per_dim$n_list[[d]]
    if (is.data.frame(dim_n)) {
      for (i in seq_len(nrow(dim_n))) {
        f <- dim_n$facet[i]
        if (f %in% names(n_dim)) {
          n_dim[[f]] <- dim_n$n[i]
        }
      }
    } else {
      if (!is.null(dim_n$facet_n)) dim_n <- dim_n$facet_n
      for (facet in names(dim_n)) {
        if (facet %in% names(n_dim)) {
          n_dim[[facet]] <- dim_n[[facet]]
        }
      }
    }
  }

          combo_result <- calculate_single_posterior(
            vc_draws = vc_draws[[d]],
            n = n_dim,
            object_spec = object_spec,
            universe_spec = universe_spec,
            error_spec = error_spec,
            agg_facets = agg_facets,
            residual_is = residual_is,
            gstudy_obj = gstudy_obj,
            n_provided = n_provided,
            use_scaled = use_scaled,
            cut_score = cut_score,
            mu_y = if (is.list(mu_y)) mu_y[[d]] else mu_y,
            ci = ci,
            probs = probs
          )
          dim_results[[d]] <- combo_result$summary
          dim_posteriors[[d]] <- combo_result$distributions
        }

        # Combine results with dimension column
        results[[i]] <- cbind(
          data.frame(n_current, stringsAsFactors = FALSE, check.names = FALSE),
          do.call(rbind, lapply(names(dim_results), function(d) {
            cbind(dim = d, dim_results[[d]])
          }))
        )
        posterior_list[[i]] <- dim_posteriors
      } else {
        # Univariate: single dimension
        combo_result <- calculate_single_posterior(
          vc_draws = vc_draws,
          n = n_current,
          object_spec = object_spec,
          universe_spec = universe_spec,
          error_spec = error_spec,
          agg_facets = agg_facets,
          residual_is = residual_is,
          gstudy_obj = gstudy_obj,
          n_provided = n_provided,
          use_scaled = use_scaled,
          cut_score = cut_score,
          mu_y = mu_y,
          ci = ci,
          probs = probs
        )

        results[[i]] <- cbind(
          data.frame(n_current, stringsAsFactors = FALSE, check.names = FALSE),
          combo_result$summary
        )
        posterior_list[[i]] <- combo_result$distributions
      }
    }

    coefficients <- do.call(rbind, results)
    coefficients <- tibble::as_tibble(coefficients)

    return(list(
      coefficients = coefficients,
      posterior = posterior_list
    ))
  } else {
    # Single sample size
    if (is_mv) {
      dim_results <- vector("list", length(mv_dims))
      names(dim_results) <- mv_dims
      dim_posteriors <- vector("list", length(mv_dims))
      names(dim_posteriors) <- mv_dims

 for (d in mv_dims) {
  # Use per-dimension sample sizes when available
  n_dim <- n
  if (!is.null(n_per_dim) && !is.null(n_per_dim$n_list[[d]])) {
    dim_n <- n_per_dim$n_list[[d]]
    # Handle different structures: sweep returns list with $facet_n,
    # non-sweep returns a data frame or a named list
    if (is.data.frame(dim_n)) {
      for (i in seq_len(nrow(dim_n))) {
        f <- dim_n$facet[i]
        if (f %in% names(n_dim)) {
          n_dim[[f]] <- dim_n$n[i]
        }
      }
    } else {
      if (!is.null(dim_n$facet_n)) dim_n <- dim_n$facet_n
      for (facet in names(dim_n)) {
        if (facet %in% names(n_dim)) {
          n_dim[[facet]] <- dim_n[[facet]]
        }
      }
    }
  }

  result <- calculate_single_posterior(
    vc_draws = vc_draws[[d]],
    n = n_dim,
    object_spec = object_spec,
    universe_spec = universe_spec,
    error_spec = error_spec,
    agg_facets = agg_facets,
    residual_is = residual_is,
    gstudy_obj = gstudy_obj,
    n_provided = n_provided,
    use_scaled = use_scaled,
    cut_score = cut_score,
    mu_y = if (is.list(mu_y)) mu_y[[d]] else mu_y,
    ci = ci,
    probs = probs
  )
  dim_results[[d]] <- cbind(dim = d, result$summary)
        dim_posteriors[[d]] <- result$distributions
      }

      coefficients <- tibble::as_tibble(do.call(rbind, lapply(mv_dims, function(d) dim_results[[d]])))

      composite_result <- NULL
      composite_posterior <- NULL

      if (!is.null(cov_draws) && !is.null(weights)) {
        vc_result <- calculate_dstudy_variance_composite(
          vc_draws = vc_draws,
          cov_draws = cov_draws,
          weights = weights,
          n = n,
          object = object,
          n_provided = n_provided,
          residual_is = residual_is
        )

        components <- names(vc_draws[[1]])
        scale_factors <- list()
        for (comp in components) {
          scale_factors[[comp]] <- compute_component_scale_factor(
            comp, n, object_spec, n_provided, residual_is = residual_is
          )
        }

        coef_result <- calculate_composite_coefficients_from_draws(
          vc_draws = vc_draws,
          cov_draws = cov_draws,
          weights = weights,
          scale_factors = scale_factors,
          object_spec = object_spec,
          universe_spec = universe_spec,
          error_spec = error_spec,
          agg_facets = agg_facets,
          residual_is = residual_is,
          cut_score = cut_score,
          mu_y = mu_y,
          ci = ci,
          probs = probs,
          gstudy_data = gstudy_obj$data,
          dimension_var = gstudy_obj$dimension_var
        )

        coefficients <- dplyr::bind_rows(coefficients, coef_result$summary)

        var_results <- coef_result$var_results

        composite_posterior <- coef_result$distributions

        return(list(
          coefficients = coefficients,
          posterior = dim_posteriors,
          composite_vc = vc_result$vc_table,
          composite_posterior = composite_posterior,
          var_results = var_results
        ))
      }

      return(list(
        coefficients = coefficients,
        posterior = dim_posteriors,
        composite_vc = NULL,
        composite_posterior = NULL,
        var_results = NULL
      ))
    } else {
      result <- calculate_single_posterior(
        vc_draws = vc_draws,
        n = n,
        object_spec = object_spec,
        universe_spec = universe_spec,
        error_spec = error_spec,
        agg_facets = agg_facets,
        residual_is = residual_is,
        gstudy_obj = gstudy_obj,
        n_provided = n_provided,
        use_scaled = use_scaled,
        cut_score = cut_score,
        mu_y = mu_y,
        ci = ci,
        probs = probs
      )

      return(list(
        coefficients = tibble::as_tibble(result$summary),
        posterior = result$distributions
      ))
    }
  }
}

#' Extract Variance Component Draws from brms Model
#'
#' Extracts posterior draws of variance components (SD^2) from a brms model.
#' For multivariate models, returns a nested list structure by dimension.
#'
#' @param gstudy_obj A gstudy object with brms estimator.
#' @param draws A posterior draws_matrix object (e.g., from brms::as_draws_matrix()).
#'
#' @return A named list. For univariate: list of variance draws per component.
#'   For multivariate: nested list with dimensions as names, each containing
#'   variance draws for that dimension.
#'
#' @keywords internal
extract_variance_draws <- function(gstudy_obj, draws) {
  model <- gstudy_obj$model
  vc <- gstudy_obj$variance_components
  is_mv <- gstudy_obj$is_multivariate
  is_long_format <- isTRUE(gstudy_obj$long_format_multivariate)
  dimension_var <- gstudy_obj$dimension_var

  if (is_mv) {
    dims <- unique(vc$dim)
    vc_draws <- vector("list", length(dims))
    names(vc_draws) <- dims

    for (d in dims) {
      vc_dim <- vc[vc$dim == d, ]
      vc_draws[[d]] <- list()

      for (comp in vc_dim$component) {
        sd_draws <- extract_sd_draws_multivariate(
          draws = draws,
          comp = comp,
          dim = d,
          vc_dim = vc_dim,
          is_long_format = is_long_format,
          dimension_var = dimension_var
        )
        vc_draws[[d]][[comp]] <- sd_draws^2
      }
    }
  } else {
    vc_draws <- list()

    for (comp in vc$component) {
      sd_draws <- extract_sd_draws_univariate(
        draws = draws,
        comp = comp,
        vc = vc
      )
      vc_draws[[comp]] <- sd_draws^2
    }
  }

  vc_draws
}

extract_sd_draws_multivariate <- function(draws, comp, dim, vc_dim, is_long_format = FALSE, dimension_var = NULL) {
  param_candidates <- character()

  if (is_long_format && !is.null(dimension_var)) {
    dim_in_params <- paste0(dimension_var, dim)
  } else {
    dim_in_params <- dim
  }

  if (comp == "Residual") {
    param_candidates <- c(
      paste0("sigma_", dim_in_params),
      paste0("sigma_", gsub(" ", "_", dim_in_params)),
      paste0("b_sigma_", dim_in_params),
      paste0("sigma_", dim),
      paste0("b_sigma_", dim)
    )
  } else {
    clean_comp <- gsub(":", "__", comp)
    param_candidates <- c(
      paste0("sd_", clean_comp, "__", dim_in_params, "_Intercept"),
      paste0("sd_", clean_comp, "__", dim_in_params),
      paste0("sd_", comp, "__", dim_in_params, "_Intercept"),
      paste0("sd_", comp, "__", dim_in_params),
      paste0("sd_", clean_comp, "__", dim, "_Intercept"),
      paste0("sd_", clean_comp, "__", dim),
      paste0("sd_", comp, "__", dim, "_Intercept"),
      paste0("sd_", comp, "__", dim)
    )
  }

  for (param_name in param_candidates) {
    if (param_name %in% colnames(draws)) {
      sd_draws <- posterior::extract_variable(draws, param_name)
      if (comp == "Residual" && grepl("^b_sigma_", param_name)) {
        sd_draws <- exp(sd_draws)
      }
      return(sd_draws)
    }
  }

  base_pattern <- if (comp == "Residual") {
    paste0("sigma_", dim_in_params)
  } else {
    paste0("sd_", gsub(":", "__", comp), "__", dim_in_params)
  }
  matching_params <- grep(paste0("^", base_pattern), colnames(draws), value = TRUE)
  if (length(matching_params) > 0) {
    sd_draws <- posterior::extract_variable(draws, matching_params[1])
    return(sd_draws)
  }

  base_pattern_fallback <- if (comp == "Residual") {
    paste0("sigma_", dim)
  } else {
    paste0("sd_", gsub(":", "__", comp), "__", dim)
  }
  matching_params_fallback <- grep(paste0("^", base_pattern_fallback), colnames(draws), value = TRUE)
  if (length(matching_params_fallback) > 0) {
    sd_draws <- posterior::extract_variable(draws, matching_params_fallback[1])
    return(sd_draws)
  }

  var_estimate <- vc_dim$var[vc_dim$component == comp]
  n_draws <- if (is.matrix(draws)) nrow(draws) else length(draws)
  rep(sqrt(var_estimate), n_draws)
}

extract_sd_draws_univariate <- function(draws, comp, vc) {
  if (comp == "Residual") {
    param_name <- "sigma"
    param_name_alt <- "b_sigma"
    matching_pattern <- "^sigma"
  } else {
    clean_comp <- gsub(":", "__", comp)
    param_name <- paste0("sd_", clean_comp, "__Intercept")
    param_name_alt <- paste0("sd_", clean_comp)
    matching_pattern <- paste0("^sd_", clean_comp)
  }

  if (param_name %in% colnames(draws)) {
    sd_draws <- posterior::extract_variable(draws, param_name)
    if (comp == "Residual" && startsWith(param_name, "b_sigma")) {
      sd_draws <- exp(sd_draws)
    }
    return(sd_draws)
  }

  if (param_name_alt %in% colnames(draws)) {
    sd_draws <- posterior::extract_variable(draws, param_name_alt)
    if (comp == "Residual" && startsWith(param_name_alt, "b_sigma")) {
      sd_draws <- exp(sd_draws)
    }
    return(sd_draws)
  }

  matching_params <- grep(matching_pattern, colnames(draws), value = TRUE)
  if (length(matching_params) > 0) {
    sd_draws <- posterior::extract_variable(draws, matching_params[1])
    if (comp == "Residual" && startsWith(matching_params[1], "b_sigma")) {
      sd_draws <- exp(sd_draws)
    }
    return(sd_draws)
  }

  var_estimate <- vc$var[vc$component == comp]
  n_draws <- if (is.matrix(draws)) nrow(draws) else length(draws)
  rep(sqrt(var_estimate), n_draws)
}

#' Calculate Posterior Coefficients for a Single Sample Size
#'
#' @param vc_draws Named list of variance component draws.
#' @param n Named list of sample sizes.
#' @param object_spec Character vector of object components.
#' @param universe_spec Character vector of universe components.
#' @param error_spec Character vector of error components (or NULL).
#' @param agg_facets Character vector of aggregation facets (or NULL).
#' @param residual_is Character string for residual composition.
#' @param gstudy_obj Original gstudy object.
#' @param n_provided Logical indicating if n was explicitly provided.
#' @param use_scaled Logical indicating whether to use scaled variance draws.
#'   When FALSE, uses original variance draws for both universe and error.
#'   When TRUE, uses scaled variance draws for error components only;
#'   universe components always use unscaled (G-study) variance draws.
#' @param cut_score Optional cut score for phi-cut calculation.
#' @param mu_y Optional grand mean for phi-cut calculation.
#'
#' @return List with 'summary' (data frame of means) and 'distributions' (list of vectors).
#'
#' @keywords internal
calculate_single_posterior <- function(vc_draws, n, object_spec, universe_spec,
                                       error_spec, agg_facets, residual_is, gstudy_obj,
                                       n_provided = FALSE, use_scaled = TRUE,
                                       cut_score = NULL, mu_y = NULL,
                                       ci = NULL, probs = c(0.025, 0.975)) {
  n_draws <- length(vc_draws[[1]])

  # Default universe to object if not specified
  if (is.null(universe_spec) || length(universe_spec) == 0) {
    universe_spec <- object_spec
  }

  # Scale variance components for D-study (only if use_scaled)
  if (use_scaled) {
    scaled_draws <- scale_variance_draws(vc_draws, n, object_spec, agg_facets,
      residual_is, gstudy_obj, n_provided)
  } else {
    scaled_draws <- vc_draws
  }

  # Calculate universe score variance (uni)
  # Universe score is always estimated from unscaled (G-study) variance components
  uni <- rep(0, n_draws)
  for (comp in universe_spec) {
    if (comp %in% names(vc_draws)) {
      uni <- uni + vc_draws[[comp]]
    }
  }

  # Calculate relative error variance (sigma2_delta)
  sigma2_delta <- compute_relative_error_draws(scaled_draws, universe_spec,
    error_spec, agg_facets, residual_is)

  # Calculate absolute error variance (sigma2_delta_abs)
  sigma2_delta_abs <- compute_absolute_error_draws(scaled_draws, universe_spec,
    error_spec, agg_facets, residual_is)

  # Calculate coefficients
  g <- uni / (uni + sigma2_delta)
  phi <- uni / (uni + sigma2_delta_abs)

  phi_cut <- NA_real_
  if (!is.null(cut_score) && !is.null(mu_y)) {
    adjustment <- (mu_y - cut_score)^2
    phi_cut <- (uni + adjustment) / (uni + sigma2_delta_abs + adjustment)
  }

  sem_rel <- sqrt(sigma2_delta)
  sem_abs <- sqrt(sigma2_delta_abs)

  # Handle infinite/NaN cases
  g[is.nan(g) | is.infinite(g)] <- NA
  phi[is.nan(phi) | is.infinite(phi)] <- NA
  if (!all(is.na(phi_cut))) {
    phi_cut[is.nan(phi_cut) | is.infinite(phi_cut)] <- NA
  }

  ci_cols <- list()
  if (!is.null(ci)) {
    if ("g" %in% ci) {
      ci_cols$g_LL <- quantile(g, probs[1], na.rm = TRUE)
      ci_cols$g_UL <- quantile(g, probs[2], na.rm = TRUE)
    }
    if ("phi" %in% ci) {
      ci_cols$phi_LL <- quantile(phi, probs[1], na.rm = TRUE)
      ci_cols$phi_UL <- quantile(phi, probs[2], na.rm = TRUE)
    }
    if ("phi-cut" %in% ci && !is.null(cut_score)) {
      ci_cols$phi_cut_LL <- quantile(phi_cut, probs[1], na.rm = TRUE)
      ci_cols$phi_cut_UL <- quantile(phi_cut, probs[2], na.rm = TRUE)
    }
  }

  summary_df <- data.frame(
    uni = mean(uni, na.rm = TRUE),
    sigma2_delta = mean(sigma2_delta, na.rm = TRUE),
    sigma2_delta_abs = mean(sigma2_delta_abs, na.rm = TRUE),
    g = mean(g, na.rm = TRUE),
    phi = mean(phi, na.rm = TRUE),
    sem_rel = mean(sem_rel, na.rm = TRUE),
    sem_abs = mean(sem_abs, na.rm = TRUE)
  )

  if (!is.null(cut_score)) {
    summary_df$phi_cut <- mean(phi_cut, na.rm = TRUE)
  }

  if (length(ci_cols) > 0) {
    summary_df <- cbind(summary_df, as.data.frame(ci_cols))
  }

  distributions <- list(
    uni = uni,
    sigma2_delta = sigma2_delta,
    sigma2_delta_abs = sigma2_delta_abs,
    g = g,
    phi = phi,
    sem_rel = sem_rel,
    sem_abs = sem_abs
  )

  if (!is.null(cut_score)) {
    distributions$phi_cut <- phi_cut
  }

  list(
    summary = summary_df,
    distributions = distributions
  )
}

#' Scale Variance Draws for D-Study
#'
#' Applies D-study scaling to posterior draws of variance components.
#'
#' @param vc_draws Named list of variance component draws.
#' @param n Named list of sample sizes.
#' @param object_spec Character vector of object components.
#' @param agg_facets Character vector of aggregation facets.
#' @param residual_is Character string for residual composition.
#' @param gstudy_obj Original gstudy object.
#' @param n_provided Logical indicating if n was explicitly provided.
#'
#' @return Named list of scaled variance draws.
#'
#' @keywords internal
scale_variance_draws <- function(vc_draws, n, object_spec, agg_facets,
  residual_is, gstudy_obj, n_provided = FALSE) {
  scaled_draws <- vc_draws

  do_rescaling <- !is.null(agg_facets) || (n_provided && length(n) > 0)

  if (!do_rescaling) {
    return(scaled_draws)
  }

  residual_facets <- if (!is.null(residual_is) && residual_is != "") {
    parse_component_facets(residual_is)
  } else {
    NULL
  }

  for (comp in names(scaled_draws)) {
    if (comp %in% object_spec) {
      next
    }

    if (comp == "Residual" || grepl("^Residual", comp)) {
      divisor <- compute_residual_divisor(residual_is, n, object_spec)
      if (!is.null(agg_facets) && length(agg_facets) > 0) {
        non_object_residual_facets <- setdiff(
          if (is.null(residual_facets)) names(n) else residual_facets,
          object_spec
        )
        agg_in_residual <- agg_facets[agg_facets %in% non_object_residual_facets]
        agg_n <- 1
        for (f in agg_in_residual) {
          if (f %in% names(n)) agg_n <- agg_n * n[[f]]
        }
        divisor <- divisor * agg_n
      }
      if (divisor > 1) {
        scaled_draws[[comp]] <- scaled_draws[[comp]] / divisor
      }
      next
    }

    comp_facets <- parse_component_facets(comp)

    has_agg <- !is.null(agg_facets) && any(agg_facets %in% comp_facets)

    if (has_agg) {
      agg_in_comp <- comp_facets[comp_facets %in% agg_facets]
      agg_n <- 1
      for (f in agg_in_comp) {
        if (f %in% names(n)) agg_n <- agg_n * n[[f]]
      }
      if (agg_n > 1) {
        scaled_draws[[comp]] <- scaled_draws[[comp]] / agg_n
      }
    } else {
      scale_factor <- compute_scale_factor_from_facets(
        comp_facets, n, object_spec = object_spec
      )
      if (scale_factor > 1) {
        scaled_draws[[comp]] <- scaled_draws[[comp]] / scale_factor
      }
    }
  }

  if (!is.null(agg_facets) && length(agg_facets) > 0 && !is.null(residual_is)) {
    for (comp in names(scaled_draws)) {
      if (comp == "Residual" || grepl("^Residual", comp)) {
        intersect_facets <- residual_facets[residual_facets %in% agg_facets]
        intersect_facets <- intersect_facets[!(intersect_facets %in% object_spec)]
        if (length(intersect_facets) > 0) {
          agg_factor <- 1
          for (facet in intersect_facets) {
            if (facet %in% names(n)) {
              agg_factor <- agg_factor * n[[facet]]
            }
          }
          if (agg_factor > 1) {
            scaled_draws[[comp]] <- scaled_draws[[comp]] * agg_factor
          }
        }
      }
    }
  }

  scaled_draws
}

#' Compute Relative Error Variance from Draws
#'
#' @param scaled_draws Named list of scaled variance draws.
#' @param universe_spec Character vector of universe components.
#' @param error_spec Character vector of error components (or NULL).
#' @param agg_facets Character vector of aggregation facets.
#'
#' @return Numeric vector of relative error variances.
#'
#' @keywords internal
compute_relative_error_draws <- function(scaled_draws, universe_spec, error_spec,
  agg_facets, residual_is = NULL) {
  n_draws <- length(scaled_draws[[1]])

  if (!is.null(error_spec)) {
    result <- rep(0, n_draws)
    for (comp in error_spec) {
      if (comp %in% names(scaled_draws)) {
        result <- result + scaled_draws[[comp]]
      }
    }
    return(result)
  }

  result <- rep(0, n_draws)

  for (comp in names(scaled_draws)) {
    if (comp %in% universe_spec) {
      next
    }

    if (comp == "Residual" || grepl("^Residual", comp)) {
      result <- result + scaled_draws[[comp]]
      next
    }

    comp_facets <- parse_component_facets(comp)
    is_agg_main <- length(comp_facets) == 1 && comp_facets %in% agg_facets
    if (is_agg_main) {
      next
    }

    if (any(universe_spec %in% comp_facets)) {
      result <- result + scaled_draws[[comp]]
    }
  }

  result
}

#' Compute Absolute Error Variance from Draws
#'
#' @param scaled_draws Named list of scaled variance draws.
#' @param universe_spec Character vector of universe components.
#' @param error_spec Character vector of error components (or NULL).
#' @param agg_facets Character vector of aggregation facets.
#'
#' @return Numeric vector of absolute error variances.
#'
#' @keywords internal
compute_absolute_error_draws <- function(scaled_draws, universe_spec, error_spec,
  agg_facets, residual_is = NULL) {
  n_draws <- length(scaled_draws[[1]])

  if (!is.null(error_spec)) {
    result <- rep(0, n_draws)
    for (comp in error_spec) {
      if (comp %in% names(scaled_draws)) {
        result <- result + scaled_draws[[comp]]
      }
    }
    return(result)
  }

  result <- rep(0, n_draws)

  for (comp in names(scaled_draws)) {
    if (comp %in% universe_spec) {
      next
    }

    if (comp == "Residual" || grepl("^Residual", comp)) {
      result <- result + scaled_draws[[comp]]
      next
    }

    comp_facets <- parse_component_facets(comp)
    is_agg_main <- length(comp_facets) == 1 && comp_facets %in% agg_facets
    if (is_agg_main) {
      next
    }

    result <- result + scaled_draws[[comp]]
  }

  result
}

