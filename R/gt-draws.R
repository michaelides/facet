#' Extract Posterior Draws from gstudy or dstudy Objects
#'
#' Extracts pre-calculated posterior draws from gstudy or dstudy objects
#' and returns them as a wide-format tibble where each variance component
#' or coefficient is a separate column.
#'
#' @param object A gstudy or dstudy object
#' @param ... Additional arguments passed to methods
#'
#' @return A tibble or named list of tibbles with posterior draws in wide format
#' @export
gt_draws <- function(object, ...) {
  UseMethod("gt_draws")
}

#' @rdname gt_draws
#' @export
#'
#' @param components Character vector of variance components to extract.
#'   If NULL (default), all components are extracted.
#' @param dims Character vector of dimensions to extract (for multivariate models).
#'   If NULL (default), all dimensions are extracted.
#' @param n Named list of sample sizes for D-study calculation.
#'   If NULL (default), sample sizes are extracted from the G-study data.
#'
#' @return For univariate gstudy objects, a tibble with columns:
#'   \itemize{
#'     \item \code{draw}: Integer draw index
#'     \item Variance component columns (e.g., \code{person}, \code{item}, \code{Residual})
#'     \item Coefficient columns: \code{uni}, \code{sigma2_delta}, \code{sigma2_delta_abs},
#'       \code{g}, \code{phi}, \code{sem_rel}, \code{sem_abs}
#'   }
#'
#'   For multivariate gstudy objects (mgstudy), a named list of tibbles,
#'   one per dimension, each with the same column structure as univariate.
#'
#' @details
#' For gstudy objects, the function extracts:
#' \itemize{
#'   \item Variance component draws for each random effect and residual
#'   \item Derived coefficients (uni, sigma2_delta, g, phi, etc.) calculated
#'     using the sample sizes from the original G-study or provided via \code{n}
#' }
#'
#' @examples
#' \dontrun{
#' # Univariate gstudy
#' g <- gstudy(score ~ (1 | person) + (1 | item), data = my_data, backend = "brms")
#' draws <- gt_draws(g)
#'
#' # Filter specific components
#' draws <- gt_draws(g, components = c("person", "Residual"))
#'
#' # Multivariate gstudy - returns named list
#' mg <- gstudy(cbind(score1, score2) ~ (1 | person) + (1 | item),
#'   data = my_data, backend = "brms")
#' draws <- gt_draws(mg)
#' draws$score1  # Access dimension-specific draws
#' draws$score2
#'
#' # Filter to specific dimensions
#' draws <- gt_draws(mg, dims = "score1")
#' }
gt_draws.gstudy <- function(object, components = NULL, dims = NULL, n = NULL, ...) {
  if (!inherits(object, c("gstudy", "mgstudy"))) {
    stop("object must be a gstudy or mgstudy object", call. = FALSE)
  }

  if (object$backend != "brms") {
    stop("gt_draws() only works with gstudy objects fit with the brms backend",
      call. = FALSE)
  }

  if (!inherits(object$model, "brmsfit")) {
    stop("Model in gstudy object is not a brmsfit object", call. = FALSE)
  }

  draws <- brms::as_draws_matrix(object$model)
  n_draws <- nrow(draws)

  vc_draws <- extract_variance_draws(object, draws)

  if (is.null(n)) {
    n <- extract_sample_sizes(object)
  }

  is_mv <- object$is_multivariate

  if (is_mv) {
    result <- build_gstudy_draws_multivariate(
      object = object,
      vc_draws = vc_draws,
      n = n,
      n_draws = n_draws,
      components = components,
      dims = dims
    )
  } else {
    result <- build_gstudy_draws_univariate(
      object = object,
      vc_draws = vc_draws,
      n = n,
      n_draws = n_draws,
      components = components
    )
  }

  result
}

#' @rdname gt_draws
#' @export
gt_draws.mgstudy <- function(object, components = NULL, dims = NULL, n = NULL, ...) {
  gt_draws.gstudy(object, components = components, dims = dims, n = n, ...)
}

#' @rdname gt_draws
#' @export
#'
#' @param what Character string specifying what to extract:
#'   \describe{
#'     \item{"all"}{Extract all available draws (default)}
#'     \item{"coefficients"}{Extract coefficient draws (uni, sigma2_delta, g, phi, etc.)}
#'     \item{"composite"}{Extract composite posterior draws (multivariate only)}
#'     \item{"var"}{Extract VAR and PRMSE draws (multivariate only)}
#'     \item{"variance"}{Extract variance component draws from underlying gstudy}
#'   }
#' @param coefficients Character vector of coefficients to extract.
#'   If NULL (default), all coefficients are extracted.
#'
#' @return For univariate dstudy objects, a tibble with columns:
#'   \itemize{
#'     \item \code{draw}: Integer draw index
#'     \item Coefficient columns: \code{uni}, \code{sigma2_delta}, \code{sigma2_delta_abs},
#'       \code{g}, \code{phi}, \code{sem_rel}, \code{sem_abs}
#'   }
#'
#'   For multivariate dstudy objects, a named list of tibbles:
#'   \itemize{
#'     \item One tibble per dimension (with coefficients and VAR columns merged)
#'     \item \code{$composite}: Composite posterior draws (if available)
#'   }
#'
#' @details
#' For dstudy objects, the function extracts pre-calculated posterior draws
#' that were computed when \code{estimation = "posterior"} was specified
#' in the \code{\link{dstudy}()} call.
#'
#' For multivariate dstudy objects, VAR/PRMSE columns are merged with coefficient
#' columns in each dimension's tibble.
#'
#' @examples
#' \dontrun{
#' # D-study with posterior estimation
#' g <- gstudy(score ~ (1 | person) + (1 | item), data = my_data, backend = "brms")
#' d <- dstudy(g, n = list(item = 10), estimation = "posterior")
#'
#' # Extract all draws
#' draws <- gt_draws(d)
#'
#' # Extract only coefficient draws
#' draws <- gt_draws(d, what = "coefficients")
#'
#' # Multivariate - returns named list
#' mg <- gstudy(cbind(score1, score2) ~ (1 | person) + (1 | item),
#'   data = my_data, backend = "brms")
#' md <- dstudy(mg, n = list(item = 10), estimation = "posterior")
#' draws <- gt_draws(md)
#' draws$score1  # Per-dimension draws with VAR columns
#' draws$composite  # Composite draws
#' }
gt_draws.dstudy <- function(object, what = c("all", "coefficients", "composite",
                              "var", "variance"),
                            coefficients = NULL, dims = NULL, ...) {
  if (!inherits(object, "dstudy")) {
    stop("object must be a dstudy object", call. = FALSE)
  }

  what <- match.arg(what)

  if (object$estimation != "posterior") {
    stop(
      "gt_draws() only works with dstudy objects created with estimation = 'posterior'. ",
      "Re-run dstudy() with estimation = 'posterior' to obtain posterior draws.",
      call. = FALSE
    )
  }

  is_mv <- object$is_multivariate

  if (is_mv) {
    return(gt_draws_dstudy_multivariate(object, what, coefficients, dims))
  }

  results <- list()

  if (what %in% c("all", "coefficients")) {
    results$coefficients <- extract_dstudy_coefficients(object, coefficients, dims)
  }

  if (what %in% c("all", "variance")) {
    results$variance <- gt_draws.gstudy(object$gstudy, dims = dims)
  }

  if (length(results) == 1) {
    return(results[[1]])
  }

  results
}

#' Handle Multivariate dstudy Draws
#'
#' @keywords internal
gt_draws_dstudy_multivariate <- function(object, what, coefficients, dims) {
  posterior <- object$posterior
  composite_post <- object$composite_posterior
  var <- object$var

  if (is.null(posterior) && is.null(composite_post)) {
    stop("No posterior draws available", call. = FALSE)
  }

  available_dims <- if (!is.null(posterior)) names(posterior) else character()

  if (!is.null(dims)) {
    missing <- setdiff(dims, available_dims)
    if (length(missing) > 0) {
      warning("Dimensions not found: ", paste(missing, collapse = ", "),
        call. = FALSE)
    }
    dims <- intersect(dims, available_dims)
  } else {
    dims <- available_dims
  }

  result <- list()

  if (what %in% c("all", "coefficients") && !is.null(posterior)) {
    coef_draws <- extract_mv_coefficients_as_list(posterior, coefficients, dims)

    if (what == "all" && !is.null(var)) {
      var_draws <- extract_dstudy_var_as_list(object, dims)
      coef_draws <- merge_var_into_coefficients(coef_draws, var_draws)
    }

    result <- c(result, coef_draws)
  }

  if (what %in% c("all", "composite") && !is.null(composite_post)) {
    result$composite <- extract_dstudy_composite(object, coefficients)
  }

  if (what == "var" && !is.null(var)) {
    var_draws <- extract_dstudy_var_as_list(object, dims)
    result <- c(result, var_draws)
  }

  if (what == "variance") {
    gdraws <- gt_draws.gstudy(object$gstudy, dims = dims)
    result <- c(result, gdraws)
  }

  result
}

#' Build Univariate gstudy Draws Tibble
#'
#' @keywords internal
build_gstudy_draws_univariate <- function(object, vc_draws, n, n_draws, components) {
  available_components <- names(vc_draws)

  if (!is.null(components)) {
    missing <- setdiff(components, available_components)
    if (length(missing) > 0) {
      warning("Components not found: ", paste(missing, collapse = ", "),
        call. = FALSE)
    }
    vc_draws <- vc_draws[intersect(components, available_components)]
  }

  coef_result <- calculate_single_posterior(
    vc_draws = vc_draws,
    n = n,
    object_spec = object$object,
    universe_spec = object$object,
    error_spec = NULL,
    agg_facets = NULL,
    residual_is = NULL,
    gstudy_obj = object,
    n_provided = FALSE,
    use_scaled = FALSE
  )

  result <- data.frame(draw = seq_len(n_draws))

  for (comp in names(vc_draws)) {
    result[[comp]] <- vc_draws[[comp]]
  }

  result$uni <- coef_result$distributions$uni
  result$sigma2_delta <- coef_result$distributions$sigma2_delta
  result$sigma2_delta_abs <- coef_result$distributions$sigma2_delta_abs
  result$g <- coef_result$distributions$g
  result$phi <- coef_result$distributions$phi
  result$sem_rel <- coef_result$distributions$sem_rel
  result$sem_abs <- coef_result$distributions$sem_abs

  tibble::as_tibble(result)
}

#' Build Multivariate gstudy Draws as Named List
#'
#' Returns a named list of tibbles, one per dimension.
#' Each tibble has the same structure as univariate.
#'
#' @keywords internal
build_gstudy_draws_multivariate <- function(object, vc_draws, n, n_draws,
                                            components, dims) {
  available_dims <- names(vc_draws)

  if (!is.null(dims)) {
    missing <- setdiff(dims, available_dims)
    if (length(missing) > 0) {
      warning("Dimensions not found: ", paste(missing, collapse = ", "),
        call. = FALSE)
    }
    vc_draws <- vc_draws[intersect(dims, available_dims)]
  }

  result_list <- list()

  for (d in names(vc_draws)) {
    dim_draws <- vc_draws[[d]]

    if (!is.null(components)) {
      available_comps <- names(dim_draws)
      missing <- setdiff(components, available_comps)
      if (length(missing) > 0) {
        warning("Components not found for dimension '", d, "': ",
          paste(missing, collapse = ", "), call. = FALSE)
      }
      dim_draws <- dim_draws[intersect(components, available_comps)]
    }

    coef_result <- calculate_single_posterior(
      vc_draws = dim_draws,
      n = n,
      object_spec = object$object,
      universe_spec = object$object,
      error_spec = NULL,
      agg_facets = NULL,
      residual_is = NULL,
      gstudy_obj = object,
      n_provided = FALSE,
      use_scaled = FALSE
    )

    dim_result <- data.frame(draw = seq_len(n_draws), stringsAsFactors = FALSE)

    for (comp in names(dim_draws)) {
      dim_result[[comp]] <- dim_draws[[comp]]
    }

    dim_result$uni <- coef_result$distributions$uni
    dim_result$sigma2_delta <- coef_result$distributions$sigma2_delta
    dim_result$sigma2_delta_abs <- coef_result$distributions$sigma2_delta_abs
    dim_result$g <- coef_result$distributions$g
    dim_result$phi <- coef_result$distributions$phi
    dim_result$sem_rel <- coef_result$distributions$sem_rel
    dim_result$sem_abs <- coef_result$distributions$sem_abs

    result_list[[d]] <- tibble::as_tibble(dim_result)
  }

  result_list
}

#' Extract Coefficient Draws from dstudy (Univariate)
#'
#' @keywords internal
extract_dstudy_coefficients <- function(object, coefficients, dims) {
  posterior <- object$posterior

  if (is.null(posterior)) {
    return(NULL)
  }

  is_sweep <- object$is_sweep

  if (is_sweep) {
    return(extract_sweep_coefficients(object, coefficients, dims))
  }

  available_coefs <- c("uni", "sigma2_delta", "sigma2_delta_abs",
    "g", "phi", "sem_rel", "sem_abs")

  if (!is.null(coefficients)) {
    missing <- setdiff(coefficients, available_coefs)
    if (length(missing) > 0) {
      warning("Coefficients not found: ", paste(missing, collapse = ", "),
        call. = FALSE)
    }
    coefficients <- intersect(coefficients, available_coefs)
  } else {
    coefficients <- available_coefs
  }

  n_draws <- length(posterior[[1]])

  result <- data.frame(draw = seq_len(n_draws))

  for (coef in coefficients) {
    if (coef %in% names(posterior)) {
      result[[coef]] <- posterior[[coef]]
    }
  }

  tibble::as_tibble(result)
}

#' Extract Multivariate Coefficient Draws as Named List
#'
#' Returns a named list of tibbles (no dim or type columns).
#'
#' @keywords internal
extract_mv_coefficients_as_list <- function(posterior, coefficients, dims) {
  available_dims <- names(posterior)

  if (!is.null(dims)) {
    posterior <- posterior[intersect(dims, available_dims)]
  }

  result_list <- list()

  for (d in names(posterior)) {
    dim_post <- posterior[[d]]

    available_coefs <- c("uni", "sigma2_delta", "sigma2_delta_abs",
      "g", "phi", "sem_rel", "sem_abs")

    if (!is.null(coefficients)) {
      use_coefs <- intersect(coefficients, available_coefs)
    } else {
      use_coefs <- available_coefs
    }

    n_draws <- length(dim_post[[1]])

    dim_result <- data.frame(draw = seq_len(n_draws), stringsAsFactors = FALSE)

    for (coef in use_coefs) {
      if (coef %in% names(dim_post)) {
        dim_result[[coef]] <- dim_post[[coef]]
      }
    }

    result_list[[d]] <- tibble::as_tibble(dim_result)
  }

  result_list
}

#' Extract Sweep Coefficient Draws
#'
#' @keywords internal
extract_sweep_coefficients <- function(object, coefficients, dims) {
  posterior <- object$posterior
  n_grid <- object$n

  result_list <- list()

  for (i in seq_along(posterior)) {
    sweep_post <- posterior[[i]]

    if (is.list(sweep_post) && !is.null(names(sweep_post))) {
      available_coefs <- c("uni", "sigma2_delta", "sigma2_delta_abs",
        "g", "phi", "sem_rel", "sem_abs")

      if (!is.null(coefficients)) {
        use_coefs <- intersect(coefficients, available_coefs)
      } else {
        use_coefs <- available_coefs
      }

      n_draws <- length(sweep_post[[1]])

      sweep_result <- data.frame(
        draw = seq_len(n_draws),
        sweep_index = i,
        stringsAsFactors = FALSE
      )

      for (coef in use_coefs) {
        if (coef %in% names(sweep_post)) {
          sweep_result[[coef]] <- sweep_post[[coef]]
        }
      }

      for (facet in names(n_grid)) {
        sweep_result[[paste0("n_", facet)]] <- n_grid[[facet]][i]
      }

      result_list[[i]] <- tibble::as_tibble(sweep_result)
    }
  }

  dplyr::bind_rows(result_list)
}

#' Extract Composite Posterior Draws from dstudy
#'
#' @keywords internal
extract_dstudy_composite <- function(object, coefficients) {
  composite_post <- object$composite_posterior

  if (is.null(composite_post)) {
    return(NULL)
  }

  available_coefs <- c("uni", "sigma2_delta", "sigma2_delta_abs",
    "g", "phi", "sem_rel", "sem_abs")

  if (!is.null(coefficients)) {
    use_coefs <- intersect(coefficients, available_coefs)
  } else {
    use_coefs <- available_coefs
  }

  n_draws <- length(composite_post[[1]])

  result <- data.frame(draw = seq_len(n_draws), stringsAsFactors = FALSE)

  for (coef in use_coefs) {
    if (coef %in% names(composite_post)) {
      result[[coef]] <- composite_post[[coef]]
    }
  }

  tibble::as_tibble(result)
}

#' Extract VAR Draws as Named List
#'
#' Returns a named list of tibbles, one per dimension.
#'
#' @keywords internal
extract_dstudy_var_as_list <- function(object, dims) {
  var <- object$var

  if (is.null(var)) {
    return(NULL)
  }

  var_matrices <- list(
    prmse_c_rel = var$prmse_c_rel_draws,
    prmse_c_abs = var$prmse_c_abs_draws,
    var_rel = var$var_rel_draws,
    var_abs = var$var_abs_draws
  )

  var_matrices <- var_matrices[!sapply(var_matrices, is.null)]

  if (length(var_matrices) == 0) {
    return(NULL)
  }

  available_dims <- colnames(var_matrices[[1]])

  if (!is.null(dims)) {
    dims_to_use <- intersect(dims, available_dims)
    for (nm in names(var_matrices)) {
      var_matrices[[nm]] <- var_matrices[[nm]][, dims_to_use, drop = FALSE]
    }
  } else {
    dims_to_use <- available_dims
  }

  n_draws <- nrow(var_matrices[[1]])
  result_list <- list()

  for (d_idx in seq_along(dims_to_use)) {
    d <- dims_to_use[d_idx]

    dim_result <- data.frame(draw = seq_len(n_draws), stringsAsFactors = FALSE)

    for (metric in names(var_matrices)) {
      dim_result[[metric]] <- var_matrices[[metric]][, d_idx]
    }

    result_list[[d]] <- tibble::as_tibble(dim_result)
  }

  result_list
}

#' Merge VAR Draws into Coefficient Tibbles
#'
#' Adds VAR columns to each dimension's coefficient tibble.
#'
#' @keywords internal
merge_var_into_coefficients <- function(coef_list, var_list) {
  if (is.null(var_list) || length(var_list) == 0) {
    return(coef_list)
  }

  for (d in names(var_list)) {
    if (d %in% names(coef_list)) {
      var_cols <- setdiff(names(var_list[[d]]), "draw")
      for (col in var_cols) {
        coef_list[[d]][[col]] <- var_list[[d]][[col]]
      }
    }
  }

  coef_list
}
