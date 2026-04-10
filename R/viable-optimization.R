#' Weight Optimization for Composite Scores
#'
#' Functions for finding optimal weights for composite scores in multivariate
#' generalizability theory. Three optimization methods are provided:
#' - "composite": Maximize composite reliability (PRMSE) via generalized eigenvalue
#' - "subscale": Maximize VAR for subscales (per-subscale and minimax)
#' - "tuning": Grid search guided by VAR diagnostics
#'
#' @name viable-optimization
#' @keywords internal
NULL

#' Extract Variance-Covariance Matrices for Optimization
#'
#' Extracts the universe score and error variance-covariance matrices from
#' a dstudy object for use in weight optimization.
#'
#' @param dstudy_obj A dstudy object
#' @param use_posterior_mean If TRUE, use posterior means; if FALSE, use point estimates
#'
#' @return List with Sigma_tau and Sigma_delta matrices
#'
#' @keywords internal
extract_variance_matrices <- function(dstudy_obj, use_posterior_mean = TRUE) {
  gstudy_obj <- dstudy_obj$gstudy
  dimensions <- gstudy_obj$dimensions
  k <- length(dimensions)

  if (is.null(dimensions) || k < 2) {
    stop("Optimization requires at least 2 dimensions", call. = FALSE)
  }

  is_mom <- inherits(gstudy_obj$model, "momfit")

  if (is_mom) {
    mom_draws <- generate_mom_variance_and_covariance_draws(gstudy_obj)
    vc_draws <- mom_draws$vc_draws
    cov_draws <- mom_draws$cov_draws
    n_draws <- mom_draws$n_draws
  } else {
    vc_draws <- extract_variance_draws_from_gstudy(gstudy_obj)
    cov_draws <- gstudy_obj$correlations
    n_draws <- length(vc_draws[[1]][[1]])
  }

  components <- names(vc_draws[[1]])
  universe_spec <- dstudy_obj$universe
  object_spec <- dstudy_obj$object
  if (isTRUE(dstudy_obj$is_sweep)) {
    n <- as.list(gstudy_obj$facet_n)
  } else {
    n <- dstudy_obj$n
  }

  scale_factors <- compute_scale_factors_for_viable(components, n, universe_spec, object_spec)

  if (use_posterior_mean) {
    vc_means <- lapply(vc_draws, function(dim_draws) {
      sapply(dim_draws, mean, na.rm = TRUE)
    })

    Sigma_tau <- matrix(0, k, k)
    rownames(Sigma_tau) <- colnames(Sigma_tau) <- dimensions

    Sigma_delta_rel <- matrix(0, k, k)
    rownames(Sigma_delta_rel) <- colnames(Sigma_delta_rel) <- dimensions

    Sigma_delta_abs <- matrix(0, k, k)
    rownames(Sigma_delta_abs) <- colnames(Sigma_delta_abs) <- dimensions

    for (i in seq_along(dimensions)) {
      d <- dimensions[i]

      for (comp in components) {
        sf <- scale_factors[[comp]]
        vc_scaled <- vc_means[[d]][comp] / sf

        if (comp %in% universe_spec) {
          Sigma_tau[i, i] <- Sigma_tau[i, i] + vc_scaled
        } else {
          is_rel <- comp == "Residual" || grepl("^Residual", comp) || {
            facets <- parse_component_facets(comp)
            any(universe_spec %in% facets)
          }
          if (is_rel) {
            Sigma_delta_rel[i, i] <- Sigma_delta_rel[i, i] + vc_scaled
          }
          Sigma_delta_abs[i, i] <- Sigma_delta_abs[i, i] + vc_scaled
        }
      }
    }

    if (!is.null(cov_draws)) {
      for (comp in components) {
        if (comp %in% universe_spec) {
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
                cov_mean <- mean(cov_list[[pair_name]], na.rm = TRUE) / scale_factors[[comp]]
                Sigma_tau[idx1, idx2] <- Sigma_tau[idx1, idx2] + cov_mean
                Sigma_tau[idx2, idx1] <- Sigma_tau[idx2, idx1] + cov_mean
              }
            }
          }
        }
      }
    }

    return(list(
      Sigma_tau = Sigma_tau,
      Sigma_delta = Sigma_delta_rel,
      Sigma_delta_abs = Sigma_delta_abs,
      Sigma_total = Sigma_tau + Sigma_delta_rel
    ))
  }

  list(
    vc_draws = vc_draws,
    cov_draws = cov_draws,
    scale_factors = scale_factors,
    n_draws = n_draws
  )
}

#' Solve Generalized Eigenvalue Problem for Composite Optimization
#'
#' Finds optimal weights that maximize composite reliability by solving
#' the generalized eigenvalue problem: Sigma_tau * w = lambda * Sigma_total * w
#'
#' @param Sigma_tau Universe score variance-covariance matrix
#' @param Sigma_total Total variance-covariance matrix (Sigma_tau + Sigma_delta)
#' @param dimensions Character vector of dimension names
#'
#' @return List with optimal weights, eigenvalue, and explained variance
#'
#' @keywords internal
solve_generalized_eigenvalue <- function(Sigma_tau, Sigma_total, dimensions) {
  k <- length(dimensions)

  cond_number <- tryCatch(
    rcond(Sigma_total),
    error = function(e) 1e-10
  )

  if (is.na(cond_number) || cond_number < 1e-10) {
    eps <- 1e-6 * mean(diag(Sigma_total))
    Sigma_total <- Sigma_total + diag(eps, k)
  }

  use_geigen <- requireNamespace("geigen", quietly = TRUE)

  if (use_geigen) {
    eig <- tryCatch(
      geigen::geigen(Sigma_tau, Sigma_total),
      error = function(e) NULL
    )

    if (!is.null(eig) && !is.null(eig$values)) {
      valid_idx <- which(is.finite(eig$values) & eig$values > 0)
      if (length(valid_idx) > 0) {
        idx <- valid_idx[which.max(eig$values[valid_idx])]
        w_raw <- eig$vectors[, idx]
        w_optimal <- abs(w_raw) / sum(abs(w_raw))
        names(w_optimal) <- dimensions
        return(list(
          weights = w_optimal,
          eigenvalue = eig$values[idx],
          method = "geigen"
        ))
      }
    }
  }

  ratio_matrix <- tryCatch(
    solve(Sigma_total) %*% Sigma_tau,
    error = function(e) {
      Sigma_total_inv <- MASS::ginv(Sigma_total)
      Sigma_total_inv %*% Sigma_tau
    }
  )

  eig <- eigen(ratio_matrix, symmetric = FALSE)

  real_eigenvalues <- Re(eig$values)
  valid_idx <- which(is.finite(real_eigenvalues) & real_eigenvalues > 0)

  if (length(valid_idx) == 0) {
    warning("No positive eigenvalues found, using equal weights", call. = FALSE)
    return(list(
      weights = setNames(rep(1/k, k), dimensions),
      eigenvalue = NA_real_,
      method = "fallback"
    ))
  }

  idx <- valid_idx[which.max(real_eigenvalues[valid_idx])]
  w_raw <- Re(eig$vectors[, idx])
  w_optimal <- abs(w_raw) / sum(abs(w_raw))
  names(w_optimal) <- dimensions

  list(
    weights = w_optimal,
    eigenvalue = real_eigenvalues[idx],
    method = "eigen_fallback"
  )
}

#' Compute VAR for Given Weights
#'
#' Computes the Value Added Ratio for all subscales given a weight vector.
#' Used as the objective function for subscale optimization.
#'
#' @param weights Numeric vector of weights (will be normalized)
#' @param dstudy_obj A dstudy object
#' @param type "rel" for var_rel or "abs" for var_abs
#'
#' @return Named vector of VAR values for each subscale
#'
#' @keywords internal
compute_var_for_weights <- function(weights, dstudy_obj, type = "rel") {
  weights <- weights / sum(weights)

  result <- recalculate_var_with_weights(
    dstudy_obj = dstudy_obj,
    weights = weights,
    dims = dstudy_obj$gstudy$dimensions,
    ci = NULL,
    probs = c(0.025, 0.975),
    n = NULL
  )

  var_col <- if (type == "rel") "var_rel" else "var_abs"
  var_values <- result[[var_col]]
  names(var_values) <- result$dim

  var_values
}

#' Softmax Transformation for Unconstrained Optimization
#'
#' Transforms unconstrained parameters to weights that sum to 1.
#'
#' @param u Unconstrained parameter vector
#'
#' @return Weight vector (sums to 1, all positive)
#'
#' @keywords internal
softmax_transform <- function(u) {
  exp_u <- exp(u - max(u))
  exp_u / sum(exp_u)
}

#' Generate Weight Grid for Tuning Method
#'
#' Creates a grid of weight combinations for exhaustive search.
#'
#' @param k Number of dimensions
#' @param resolution Step size for grid (default 0.1)
#' @param min_weight Minimum weight allowed (default 0.01)
#'
#' @return Tibble with all valid weight combinations
#'
#' @keywords internal
generate_weight_grid <- function(k, resolution = 0.1, min_weight = 0.01) {
  if (k == 2) {
    w1_seq <- seq(min_weight, 1 - min_weight, by = resolution)
    grid <- tibble::tibble(w1 = w1_seq)
    grid$w2 <- 1 - grid$w1
    return(grid)
  }

  if (k == 3) {
    max_w1 <- 1 - 2 * min_weight
    max_w2 <- 1 - 2 * min_weight

    grid <- tidyr::crossing(
      w1 = seq(min_weight, max_w1, by = resolution),
      w2 = seq(min_weight, max_w2, by = resolution)
    )
    grid <- grid %>%
      dplyr::filter(w1 + w2 <= 1 - min_weight) %>%
      dplyr::mutate(w3 = 1 - w1 - w2)

    return(grid)
  }

  grid <- expand.grid(
    lapply(1:(k-1), function(i) seq(min_weight, 1 - (k-1)*min_weight, by = resolution))
  )
  names(grid) <- paste0("w", 1:(k-1))

  grid <- grid %>%
    dplyr::rowwise() %>%
    dplyr::filter(sum(dplyr::c_across(dplyr::everything())) <= 1 - min_weight) %>%
    dplyr::ungroup()

  grid[[paste0("w", k)]] <- 1 - rowSums(grid)

  grid
}

#' Optimize Weights for Composite Reliability
#'
#' Main function for "composite" optimization method.
#'
#' @param dstudy_obj A dstudy object
#' @param optimize_target "rel" or "abs"
#'
#' @return List with optimization results
#'
#' @keywords internal
optimize_composite_weights_internal <- function(dstudy_obj, optimize_target = "rel") {
  dimensions <- dstudy_obj$gstudy$dimensions
  k <- length(dimensions)

  matrices <- extract_variance_matrices(dstudy_obj, use_posterior_mean = TRUE)

  if (optimize_target == "abs") {
    Sigma_error <- matrices$Sigma_delta_abs
  } else {
    Sigma_error <- matrices$Sigma_delta
  }
  Sigma_total <- matrices$Sigma_tau + Sigma_error

  result <- solve_generalized_eigenvalue(
    Sigma_tau = matrices$Sigma_tau,
    Sigma_total = Sigma_total,
    dimensions = dimensions
  )

  w_optimal <- result$weights

  metrics <- recalculate_var_with_weights(
    dstudy_obj = dstudy_obj,
    weights = w_optimal,
    dims = dimensions,
    ci = NULL,
    probs = c(0.025, 0.975),
    n = NULL
  )

  var_col <- if (optimize_target == "rel") "var_rel" else "var_abs"

  w_mat <- matrix(w_optimal, nrow = 1) %*% matrices$Sigma_tau %*% matrix(w_optimal, ncol = 1)
  delta_mat <- matrix(w_optimal, nrow = 1) %*% Sigma_error %*% matrix(w_optimal, ncol = 1)

  composite_reliability <- as.numeric(w_mat / (w_mat + delta_mat))

  list(
    method = "composite",
    optimize_target = optimize_target,
    weights = w_optimal,
    metrics = metrics,
    composite_reliability = composite_reliability,
    eigenvalue = result$eigenvalue,
    explained_variance = result$eigenvalue / sum(diag(matrices$Sigma_tau))
  )
}

#' Optimize Weights for Subscale VAR
#'
#' Main function for "subscale" optimization method.
#'
#' @param dstudy_obj A dstudy object
#' @param subscale Target subscale name (NULL for all)
#' @param optimize_target "rel" or "abs"
#'
#' @return List with optimization results
#'
#' @keywords internal
optimize_subscale_weights_internal <- function(dstudy_obj, subscale = NULL, optimize_target = "rel") {
  dimensions <- dstudy_obj$gstudy$dimensions
  k <- length(dimensions)

  per_subscale_results <- list()

  for (target_dim in dimensions) {
    obj_fn <- function(w_unconstrained) {
      w <- softmax_transform(w_unconstrained)
      var_values <- compute_var_for_weights(w, dstudy_obj, optimize_target)
      -var_values[target_dim]
    }

    init_params <- rep(0, k)

    opt_result <- tryCatch(
      optim(
        par = init_params,
        fn = obj_fn,
        method = "BFGS",
        control = list(maxit = 1000)
      ),
      error = function(e) {
        list(par = rep(0, k), value = -1)
      }
    )

    w_optimal <- softmax_transform(opt_result$par)
    names(w_optimal) <- dimensions

    var_values <- compute_var_for_weights(w_optimal, dstudy_obj, optimize_target)

    per_subscale_results[[target_dim]] <- list(
      weights = w_optimal,
      target_var = var_values[target_dim],
      all_var = var_values
    )
  }

  minimax_fn <- function(w_unconstrained) {
    w <- softmax_transform(w_unconstrained)
    var_values <- compute_var_for_weights(w, dstudy_obj, optimize_target)
    -min(var_values)
  }

  init_params <- rep(0, k)

  minimax_result <- tryCatch(
    optim(
      par = init_params,
      fn = minimax_fn,
      method = "BFGS",
      control = list(maxit = 1000)
    ),
    error = function(e) {
      list(par = rep(0, k), value = -1)
    }
  )

  minimax_weights <- softmax_transform(minimax_result$par)
  names(minimax_weights) <- dimensions

  minimax_metrics <- recalculate_var_with_weights(
    dstudy_obj = dstudy_obj,
    weights = minimax_weights,
    dims = dimensions,
    ci = NULL,
    probs = c(0.025, 0.975),
    n = NULL
  )

  minimax_var <- compute_var_for_weights(minimax_weights, dstudy_obj, optimize_target)

  comparison_df <- tibble::tibble(
    method = c("minimax", paste0("per_", dimensions)),
    dplyr::bind_rows(
      tibble::as_tibble_row(minimax_weights),
      dplyr::bind_rows(lapply(per_subscale_results, function(x) {
        tibble::as_tibble_row(x$weights)
      }))
    )
  )

  comparison_df$min_VAR <- c(
    min(minimax_var),
    sapply(per_subscale_results, function(x) min(x$all_var))
  )

  comparison_df$all_viable <- comparison_df$min_VAR > 1

  if (!is.null(subscale)) {
    if (!(subscale %in% dimensions)) {
      stop("subscale must be one of: ", paste(dimensions, collapse = ", "), call. = FALSE)
    }

    target_result <- per_subscale_results[[subscale]]
    target_metrics <- recalculate_var_with_weights(
      dstudy_obj = dstudy_obj,
      weights = target_result$weights,
      dims = dimensions,
      ci = NULL,
      probs = c(0.025, 0.975),
      n = NULL
    )

    return(list(
      method = "subscale",
      target = subscale,
      weights = target_result$weights,
      target_var = target_result$target_var,
      all_var = target_result$all_var,
      metrics = target_metrics
    ))
  }

  list(
    method = "subscale",
    minimax = list(
      weights = minimax_weights,
      metrics = minimax_metrics,
      min_var = min(minimax_var),
      all_viable = all(minimax_var > 1)
    ),
    per_subscale = per_subscale_results,
    comparison = comparison_df
  )
}

#' Optimize Weights via Grid Search (Tuning)
#'
#' Main function for "tuning" optimization method.
#'
#' @param dstudy_obj A dstudy object
#' @param grid_resolution Step size for grid
#' @param optimize_target "rel" or "abs"
#'
#' @return List with optimization results
#'
#' @keywords internal
optimize_tuning_weights_internal <- function(dstudy_obj, grid_resolution = 0.1, optimize_target = "rel") {
  dimensions <- dstudy_obj$gstudy$dimensions
  k <- length(dimensions)

  grid <- generate_weight_grid(k, resolution = grid_resolution)

  if (k <= 3) {
    colnames(grid) <- paste0("w", 1:k)
  }

  var_col <- if (optimize_target == "rel") "var_rel" else "var_abs"

  grid_results <- lapply(1:nrow(grid), function(i) {
    weights <- unlist(grid[i, 1:k])
    weights <- weights / sum(weights)
    names(weights) <- dimensions

    metrics <- tryCatch(
      recalculate_var_with_weights(
        dstudy_obj = dstudy_obj,
        weights = weights,
        dims = dimensions,
        ci = NULL,
        probs = c(0.025, 0.975),
        n = NULL
      ),
      error = function(e) NULL
    )

    if (is.null(metrics)) {
      out <- setNames(as.list(c(rep(NA, k), NA, FALSE)), c(dimensions, "min_VAR", "all_viable"))
      out$all_viable <- FALSE
      return(out)
    }

    var_values <- metrics[[var_col]]
    names(var_values) <- metrics$dim

    out <- as.list(c(var_values, min_VAR = min(var_values)))
    out$all_viable <- all(var_values > 1)
    out
  })

  grid_results_df <- dplyr::bind_rows(grid_results)

  grid <- dplyr::bind_cols(grid, grid_results_df)

  viable_solutions <- grid %>% dplyr::filter(all_viable)

  if (nrow(viable_solutions) > 0) {
    best <- viable_solutions %>%
      dplyr::arrange(dplyr::desc(min_VAR)) %>%
      dplyr::slice(1)
  } else {
    best <- grid %>%
      dplyr::arrange(dplyr::desc(min_VAR)) %>%
      dplyr::slice(1)
  }

  best_weights <- unlist(best[1:k])
  best_weights <- best_weights / sum(best_weights)
  names(best_weights) <- dimensions

  best_metrics <- recalculate_var_with_weights(
    dstudy_obj = dstudy_obj,
    weights = best_weights,
    dims = dimensions,
    ci = NULL,
    probs = c(0.025, 0.975),
    n = NULL
  )

  list(
    method = "tuning",
    best = list(
      weights = best_weights,
      metrics = best_metrics,
      min_var = best$min_VAR,
      all_viable = best$all_viable
    ),
    grid_results = grid,
    resolution = grid_resolution,
    n_grid_points = nrow(grid),
    n_viable = sum(grid$all_viable, na.rm = TRUE)
  )
}
