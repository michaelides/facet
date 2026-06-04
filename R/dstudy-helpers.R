#' D-Study Helper Functions
#'
#' Internal helper functions extracted from dstudy() for maintainability.
#' Each function handles a distinct phase of D-study setup or computation.
#'
#' @name dstudy-helpers
#' @keywords internal
NULL

#' Validate D-Study Inputs and Extract Basic Properties
#'
#' Validates the gstudy object, determines if multivariate,
#' and processes weights.
#'
#' @param gstudy_obj A gstudy or mgstudy object
#' @param weights Optional numeric vector of weights
#' @return A list with is_multivariate, dimensions, weights
#' @keywords internal
validate_dstudy_inputs <- function(gstudy_obj, weights = NULL) {
  if (!inherits(gstudy_obj, "gstudy") && !inherits(gstudy_obj, "mgstudy")) {
    stop("'gstudy_obj' must be an object of class 'gstudy' or 'mgstudy'", call. = FALSE)
  }

  is_multivariate <- inherits(gstudy_obj, "mgstudy")

  if (is_multivariate) {
    dimensions <- gstudy_obj$dimensions
    if (is.null(weights)) {
      weights <- rep(1, length(dimensions))
    } else {
      if (length(weights) != length(dimensions)) {
        stop(
          "weights must have length ", length(dimensions),
          " (number of dimensions), got ", length(weights),
          call. = FALSE
        )
      }
    }
    names(weights) <- dimensions
  } else {
    weights <- NULL
  }

  list(
    is_multivariate = is_multivariate,
    dimensions = if (is_multivariate) dimensions else NULL,
    weights = weights
  )
}

#' Resolve Estimation Method
#'
#' Determines the appropriate estimation method based on the backend.
#' Warns and refits if posterior requested with non-brms backend.
#' Warns and overrides if simple requested with brms backend.
#'
#' @param gstudy_obj A gstudy or mgstudy object
#' @param estimation Requested estimation method (NULL, "simple", or "posterior")
#' @return A list with gstudy_obj (possibly refit) and estimation (resolved)
#' @keywords internal
resolve_estimation_method <- function(gstudy_obj, estimation = NULL) {
  check_estimation_issues(gstudy_obj)

  if (is.null(estimation)) {
    if (gstudy_obj$backend == "brms") {
      estimation <- "posterior"
    } else {
      estimation <- "simple"
    }
  } else {
    estimation <- match.arg(estimation, c("simple", "posterior"))
  }

  if (estimation == "posterior" && gstudy_obj$backend != "brms") {
    warning(
      "estimation = 'posterior' requires backend = 'brms'. ",
      "Refitting gstudy model with backend = 'brms'.",
      call. = FALSE
    )
    gstudy_obj <- gstudy(
      formula = gstudy_obj$formula,
      data = gstudy_obj$data,
      backend = "brms"
    )
  }

  if (estimation == "simple" && gstudy_obj$backend == "brms") {
    warning(
      "estimation = 'simple' is not recommended for brms backend. ",
      "The variance estimates displayed in the variance components table ",
      "are computed from squared posterior draws (mean(SD^2)), which properly ",
      "accounts for uncertainty. Using estimation = 'posterior' ensures ",
      "consistency between variance components and coefficient calculations.",
      call. = FALSE
    )
    estimation <- "posterior"
  }

  list(
    gstudy_obj = gstudy_obj,
    estimation = estimation
  )
}

#' Validate Cut Score and Credible Interval Parameters
#'
#' Extracts grand mean if cut_score provided, validates ci parameter
#' against backend, validates probs, and warns if phi-cut CI requested
#' without cut_score.
#'
#' @param gstudy_obj A gstudy or mgstudy object
#' @param cut_score Optional cutoff score
#' @param ci Credible interval specification
#' @param probs Probability levels for credible intervals
#' @return A list with mu_y, ci (possibly modified), probs
#' @keywords internal
validate_cut_score_ci <- function(gstudy_obj, cut_score = NULL, ci = NULL,
  probs = c(0.025, 0.975)) {
  mu_y <- NULL
  if (!is.null(cut_score)) {
    mu_y <- extract_grand_mean(gstudy_obj)
  }

  if (!is.null(ci)) {
    ci <- match.arg(ci, c("g", "phi", "phi-cut"), several.ok = TRUE)
    if (gstudy_obj$backend != "brms") {
      warning(
        "Credible intervals for 'mom' and 'lme4' backends are not yet implemented. ",
        "Consider using backend = 'brms' for uncertainty quantification.",
        call. = FALSE
      )
      ci <- NULL
    }
  }

  if (!is.null(ci)) {
    if (length(probs) != 2) {
      stop("'probs' must have exactly 2 elements", call. = FALSE)
    }
    if (probs[1] >= probs[2]) {
      stop("'probs' must be in increasing order", call. = FALSE)
    }
    if (any(probs < 0) || any(probs > 1)) {
      stop("'probs' must be between 0 and 1", call. = FALSE)
    }

    if ("phi-cut" %in% ci && is.null(cut_score)) {
      warning(
        "Credible intervals for 'phi-cut' require a 'cut_score' to be specified. ",
        "The phi-cut coefficient is used for criterion-referenced (absolute) decisions. ",
        "Specify cut_score to compute phi-cut credible intervals.",
        call. = FALSE
      )
      ci <- setdiff(ci, "phi-cut")
      if (length(ci) == 0) {
        ci <- NULL
      }
    }
  }

  list(
    mu_y = mu_y,
    ci = ci,
    probs = probs
  )
}

#' Parse D-Study Specifications
#'
#' Extracts variance components, object of measurement, and parses
#' universe/error/aggregation specifications. Validates that universe
#' and error don't overlap and warns about non-interacting universe components.
#'
#' @param gstudy_obj A gstudy or mgstudy object
#' @param universe Universe specification
#' @param error Error specification
#' @param aggregation Aggregation specification
#' @return A list with vc, object, object_spec, universe_spec, error_spec, agg_spec
#' @keywords internal
parse_dstudy_specifications <- function(gstudy_obj, universe = NULL,
  error = NULL, aggregation = NULL) {
  vc <- gstudy_obj$variance_components
  object <- gstudy_obj$object
  object_spec <- parse_specification(object)

  universe_spec <- parse_specification(universe)

  if (is.null(universe) || length(universe_spec) == 0) {
    universe_spec <- object_spec
  } else {
    if (!all(object_spec %in% universe_spec)) {
      warning(
        "The specification of the universe did not include the object of measurement '",
        paste(object_spec, collapse = ", "),
        "'. It has been added automatically.",
        call. = FALSE
      )
      universe_spec <- unique(c(object_spec, universe_spec))
    }

    for (comp in setdiff(universe_spec, object_spec)) {
      comp_facets <- parse_component_facets(comp)
      if (!any(object_spec %in% comp_facets)) {
        potential_interaction <- paste0(object, ":", comp)
        alternative <- paste0(comp, ":", object)
        warning(
          "Component '", comp, "' in universe does not interact with the object '", object, "'. ",
          "The universe score will include '", comp, "' variance scaled by n_", comp, ". ",
          "If you meant to include the interaction '", potential_interaction, "' or '", alternative, "', specify it explicitly.",
          call. = FALSE
        )
      }
    }
  }

  error_spec <- parse_specification(error)
  agg_spec <- parse_specification(aggregation)

  if (!is.null(error) && length(error_spec) > 0) {
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

  list(
    vc = vc,
    object = object,
    object_spec = object_spec,
    universe_spec = universe_spec,
    error_spec = error_spec,
    agg_spec = agg_spec
  )
}

#' Resolve D-Study Sample Sizes
#'
#' Determines sample sizes from multiple possible sources: user-provided n,
#' per-dimension sample sizes, or extracted from the G-study. Also determines
#' whether this is a sweep (multiple sample sizes per facet).
#'
#' @param gstudy_obj A gstudy or mgstudy object
#' @param n User-provided sample sizes
#' @param is_multivariate Whether this is a multivariate model
#' @return A list with n, n_provided, n_per_dim, n_tibble, is_sweep
#' @keywords internal
resolve_dstudy_sample_sizes <- function(gstudy_obj, n = list(),
  is_multivariate = FALSE) {
  n_provided <- length(n) > 0
  n_per_dim <- NULL
  n_tibble <- NULL

  needs_per_dim_n <- is_multivariate && (
    !is.null(gstudy_obj$sample_size_info_per_dim) ||
      !is.null(gstudy_obj$long_format_multivariate) ||
      (!is.null(gstudy_obj$is_unbalanced) && gstudy_obj$is_unbalanced)
  )

  if (length(n) == 0) {
    if (needs_per_dim_n) {
      n_tibble <- extract_sample_sizes_per_dim(gstudy_obj)
      if (!is.null(n_tibble) && nrow(n_tibble) > 0) {
        n_per_dim <- expand_n_per_dim(n_tibble, sweep = FALSE)
        n <- extract_sample_sizes(gstudy_obj)
        message(
          "No sample sizes provided in 'n'. Using per-dimension G-study sample sizes."
        )
      } else {
        n <- extract_sample_sizes(gstudy_obj)
        message(
          "No sample sizes provided in 'n'. Using G-study sample sizes: ",
          paste(names(n), n, sep = " = ", collapse = ", ")
        )
      }
    } else {
      n <- extract_sample_sizes(gstudy_obj)
      message(
        "No sample sizes provided in 'n'. Using G-study sample sizes: ",
        paste(names(n), n, sep = " = ", collapse = ", ")
      )
    }
  } else if (is.data.frame(n)) {
    validation <- validate_n_tibble(n, gstudy_obj$dimensions, gstudy_obj$facets)
    if (!validation$valid) {
      stop(validation$error, call. = FALSE)
    }
    n_tibble <- n
    n_per_dim <- expand_n_per_dim(n_tibble, sweep = TRUE)
    n_provided <- TRUE
    n <- extract_sample_sizes(gstudy_obj)
  } else if (is.list(n) && needs_per_dim_n) {
    n_tibble <- create_n_tibble_from_list(n, gstudy_obj$dimensions)
    if (!is.null(n_tibble)) {
      warning(
        "For multivariate models, providing n as a simple list will apply ",
        "the same sample sizes to all dimensions. Consider providing a tibble ",
        "with columns: dim, facet, n for per-dimension control.",
        call. = FALSE
      )
    }
  }

  is_sweep <- if (!is.null(n_per_dim)) {
    n_per_dim$is_sweep
  } else {
    any(sapply(n, length) > 1)
  }

  list(
    n = n,
    n_provided = n_provided,
    n_per_dim = n_per_dim,
    n_tibble = n_tibble,
    is_sweep = is_sweep
  )
}

#' Resolve Per-Dimension Sample Size for One Dimension
#'
#' Extracts a per-dim named list of facet -> n from the `n_per_dim` structure
#' built by [resolve_dstudy_sample_sizes()]. Falls back to the global `n` for
#' facets not specified per-dim. Returns a list keyed by dim when iterating
#' in a loop. Mirrors the per-dim override pattern in
#' [calculate_coefficients_posterior()] so non-posterior and posterior paths
#' agree on what "per-dim n" means.
#'
#' @param dim_name Character: the dimension to extract.
#' @param n_per_dim The `n_per_dim` slot of [resolve_dstudy_sample_sizes()] output.
#' @param n_global The flat `n` list (per-facet, possibly scalar) used as fallback.
#' @return A named list of facet -> n. Returns `n_global` unchanged when no
#'   per-dim data is available for `dim_name`.
#' @keywords internal
resolve_dim_n <- function(dim_name, n_per_dim, n_global) {
  if (is.null(n_per_dim) || is.null(n_per_dim$n_list)) {
    return(n_global)
  }
  dim_n_raw <- n_per_dim$n_list[[dim_name]]
  if (is.null(dim_n_raw)) {
    return(n_global)
  }

  n_dim <- n_global

  if (is.data.frame(dim_n_raw)) {
    for (i in seq_len(nrow(dim_n_raw))) {
      f <- dim_n_raw$facet[i]
      if (f %in% names(n_dim)) {
        n_dim[[f]] <- dim_n_raw$n[i]
      }
    }
  } else {
    if (!is.null(dim_n_raw$facet_n)) {
      dim_n_raw <- dim_n_raw$facet_n
    }
    for (facet in names(dim_n_raw)) {
      if (facet %in% names(n_dim)) {
        n_dim[[facet]] <- dim_n_raw[[facet]]
      }
    }
  }

  n_dim
}

#' Should the D-Study Path Iterate Per Dimension?
#'
#' Returns TRUE when the mgstudy object has per-dim sample-size information
#' and the user-provided `n` is a per-dim tibble (i.e., we should call the
#' variance functions separately for each dim rather than once with scalar n).
#'
#' @param is_multivariate Logical: whether the model is multivariate.
#' @param n_per_dim The `n_per_dim` slot from [resolve_dstudy_sample_sizes()].
#' @param n_tibble The `n_tibble` slot from [resolve_dstudy_sample_sizes()].
#' @return Logical scalar.
#' @keywords internal
needs_per_dim_vc_loop <- function(is_multivariate, n_per_dim, n_tibble) {
  isTRUE(is_multivariate) &&
    !is.null(n_per_dim) &&
    !is.null(n_tibble) &&
    nrow(n_tibble) > 0
}

#' Apply a Variance Function Per Dimension
#'
#' When [needs_per_dim_vc_loop()] is TRUE, splits the variance components
#' tibble by `dim` and calls `fn(vc_dim, n_dim, ...)` for each dimension,
#' then row-binds the results. Otherwise calls `fn(vc, n, ...)` once with
#' the full vc and global n. Used to give the non-posterior D-study paths
#' per-dim awareness without changing their signatures.
#'
#' @param vc Variance components tibble with a `dim` column.
#' @param n Named list of facet -> n (global n; per-dim n overrides apply).
#' @param n_per_dim The `n_per_dim` slot from [resolve_dstudy_sample_sizes()].
#' @param n_tibble The `n_tibble` slot from [resolve_dstudy_sample_sizes()].
#' @param fn A variance function such as [calculate_divided_variance()] or
#'   [calculate_dstudy_variance()].
#' @param ... Additional arguments passed to `fn`.
#' @return Tibble of variance components (per-dim, row-bound).
#' @keywords internal
apply_per_dim <- function(vc, n, n_per_dim, n_tibble, fn, ...) {
  if (!needs_per_dim_vc_loop(TRUE, n_per_dim, n_tibble)) {
    return(fn(vc, n, ...))
  }
  dims <- unique(vc$dim)
  parts <- lapply(dims, function(d) {
    vc_dim <- vc[vc$dim == d, , drop = FALSE]
    n_dim <- resolve_dim_n(d, n_per_dim, n)
    fn(vc_dim, n_dim, ...)
  })
  do.call(rbind, parts)
}

#' Process Residual and Error Specification Overlap
#'
#' Determines residual composition from formula/data, handles overlap
#' between error and aggregation specifications, and removes residual
#' from error spec if double-counting would occur.
#'
#' @param gstudy_obj A gstudy or mgstudy object
#' @param error Error specification
#' @param aggregation Aggregation specification
#' @param residual_is User-specified residual composition
#' @return A list with residual_composition, residual_is_effective, error (cleaned)
#' @keywords internal
process_residual_and_error <- function(gstudy_obj, error = NULL,
  aggregation = NULL, residual_is = NULL) {
  residual_composition <- parse_residual_facets(
    gstudy_obj$formula,
    gstudy_obj$data
  )

  residual_is_effective <- if (is.null(residual_is)) residual_composition else residual_is

  if (!is.null(error)) {
    error_spec <- parse_specification(error)

    if (!is.null(aggregation)) {
      agg_facets <- parse_specification(aggregation)

      overlap <- intersect(error_spec, agg_facets)
      if (length(overlap) > 0) {
        warning(
          "Component(s) '", paste(overlap, collapse = "', '"), "' specified for both ",
          "aggregation and error. Removing from error specification. ",
          "Note: interaction terms containing aggregation facets (e.g., 'p:", overlap[1], "') ",
          "are NOT removed.",
          call. = FALSE
        )

        error_spec <- setdiff(error_spec, agg_facets)

        if (length(error_spec) > 0) {
          error <- error_spec
        } else {
          error <- NULL
        }
      }
    }

    error_spec <- parse_specification(error)
    if (residual_composition %in% error_spec) {
      warning(
        "The residual component '", residual_composition, "' is already included ",
        "in the model variance components. Removing '", residual_composition,
        "' from error specification to avoid double-counting.",
        call. = FALSE
      )

      error_spec_filtered <- error_spec[error_spec != residual_composition]

      if (length(error_spec_filtered) > 0) {
        error <- error_spec_filtered
      } else {
        error <- NULL
      }
    }
  }

  list(
    residual_composition = residual_composition,
    residual_is_effective = residual_is_effective,
    error = error
  )
}

#' Append Composite Coefficients if Multivariate
#'
#' Conditionally appends composite coefficient rows to the coefficients tibble
#' when the model is multivariate with more than one dimension.
#'
#' @param coefficients Existing coefficients tibble
#' @param vc Variance components tibble
#' @param n Sample sizes
#' @param weights Dimension weights
#' @param object Object of measurement name
#' @param error Error specification
#' @param aggregation Aggregation specification
#' @param residual_is_effective Effective residual specification
#' @param universe_spec Universe specification
#' @param gstudy_obj The G-study object (for correlations, data, dimension_var)
#' @param cut_score Optional cut score
#' @param mu_y Grand mean
#' @param estimate_label Optional label for the estimate column
#' @return A list with coefficients (possibly with composite row) and var_results
#' @keywords internal
maybe_append_composite <- function(coefficients, vc, n, weights, object,
  error, aggregation, residual_is_effective, universe_spec,
  gstudy_obj, cut_score = NULL, mu_y = NULL, estimate_label = NULL) {
  is_multivariate <- inherits(gstudy_obj, "mgstudy")
  dimensions <- if (is_multivariate) gstudy_obj$dimensions else NULL

  if (!is_multivariate || length(dimensions) <= 1) {
    return(list(coefficients = coefficients, var_results = NULL))
  }

  composite_coefs <- calculate_composite_coefficients(
    vc = vc,
    n = n,
    weights = weights,
    object = object,
    error = error,
    aggregation = aggregation,
    residual_is = residual_is_effective,
    universe = universe_spec,
    correlations = gstudy_obj$correlations,
    cut_score = cut_score,
    mu_y = mu_y,
    gstudy_data = gstudy_obj$data,
    dimension_var = gstudy_obj$dimension_var
  )

  composite_summary <- composite_coefs$summary
  var_results <- composite_coefs$var_results

  if (!is.null(estimate_label)) {
    composite_summary$estimate <- estimate_label
  } else if ("estimate" %in% names(coefficients)) {
    composite_summary$estimate <- coefficients$estimate[1]
  }

  coefficients <- dplyr::bind_rows(coefficients, composite_summary)

  list(coefficients = coefficients, var_results = var_results)
}

#' Build D-Study Result Object
#'
#' Constructs the final dstudy result list from computed components.
#'
#' @param gstudy_obj The G-study object
#' @param d_vc D-study variance components
#' @param coefficients Coefficients tibble
#' @param n Sample sizes
#' @param n_tibble Per-dimension sample sizes tibble
#' @param n_per_dim Per-dimension sample sizes list
#' @param object Object of measurement
#' @param universe_spec Universe specification
#' @param error Error specification
#' @param aggregation Aggregation specification
#' @param residual_is Original residual_is parameter
#' @param residual_composition Residual composition
#' @param is_sweep Whether this is a sweep
#' @param estimation Estimation method used
#' @param posterior Posterior draws (NULL if simple)
#' @param composite_post Composite posterior draws
#' @param var_results VAR results
#' @param is_multivariate Whether multivariate
#' @param cut_score Cut score (NULL if not provided)
#' @param mu_y Grand mean (NULL if no cut_score)
#' @param ci CI specification
#' @param probs Probability levels
#' @param weights Dimension weights
#' @return A dstudy object
#' @keywords internal
build_dstudy_result <- function(gstudy_obj, d_vc, coefficients, n, n_tibble,
  n_per_dim, object, universe_spec, error, aggregation, residual_is,
  residual_composition, is_sweep, estimation, posterior, composite_post,
  var_results, is_multivariate, cut_score, mu_y, ci, probs, weights) {
  result <- list(
    gstudy = gstudy_obj,
    variance_components = d_vc,
    coefficients = coefficients,
    n = n,
    n_tibble = n_tibble,
    n_per_dim = n_per_dim,
    object = object,
    universe = universe_spec,
    error = error,
    aggregation = aggregation,
    residual_is = residual_is,
    residual_composition = residual_composition,
    is_sweep = is_sweep,
    estimation = estimation,
    posterior = if (estimation == "posterior") posterior else NULL,
    composite_posterior = composite_post,
    var = var_results,
    is_multivariate = is_multivariate,
    cut_score = cut_score,
    mu_y = mu_y,
    ci = ci,
    probs = if (!is.null(ci)) probs else NULL,
    weights = weights,
    long_format_multivariate = gstudy_obj$long_format_multivariate,
    dimension_var = gstudy_obj$dimension_var
  )

  class(result) <- "dstudy"
  result
}