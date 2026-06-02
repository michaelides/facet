#' Conduct a Generalizability (G) Study
#'
#' Estimate variance components for a generalizability theory analysis using
#' mixed effects models. Supports frequentist (lme4), Bayesian (brms), and
#' method of moments (mom) backends.
#'
#' @param formula A formula specifying the G-study model. The left-hand side
#'   should be the outcome variable, and the right-hand side should specify
#'   the variance components using lme4-style syntax (e.g., `y ~ (1|person) + (1|item)`).
#'   For multivariate models, use brms syntax: `mvbind(y1, y2) ~ (1|person)`.
#' @param data A data frame containing the variables in the formula.
#' @param backend Character string specifying the backend to use. One of
#'   "auto" (default, chooses based on formula type), "lme4", "brms", or "mom".
#'   "mom" uses the method of moments (ANOVA-based) for variance component estimation.
#' @param facets Character vector of facet names (optional, auto-detected from formula if NULL).
#' @param nested Optional named list specifying nesting relationships.
#' Names are nested facets and values are nesting facets.
#' For example, `list(task = "rater")` means task is nested within rater.
#' If NULL (default), nesting is auto-detected from the data structure.
#' @param ci_method Character string specifying the method for confidence intervals
#' for the lme4 backend. One of "none" (default, no CIs), "profile" (more accurate,
#' slower), or "boot" (bootstrap, most accurate, slowest).
#' Only applicable for lme4 backend; brms and mom provide CIs automatically.
#' @param nsim Integer: number of bootstrap simulations (only for ci_method = "boot").
#'   Default is 1000.
#' @param boot.type Character: bootstrap type, "perc" (percentile), "basic",
#' or "norm" (normal-theory) (only for ci_method = "boot"). Default is "perc".
#' @param prior A brmsprior object or list of priors created by [set_prior()]
#' or related functions. Only applicable when using brms backend.
#' Use [default_prior()] to see available parameters for priors.
#' @family model fitting
#' @param unbalanced Logical indicating whether to enable unbalanced multivariate
#' estimation. Default is FALSE. When TRUE:
#' \itemize{
#' \item For \strong{mom backend}: Each dimension is analyzed with its available
#' data using Henderson's Method III for variance component estimation.
#' Correlations are computed using pairwise complete cases.
#' \item For \strong{brms backend}: Not implemented. Use long-format specification
#' with \code{bf(Score ~ 0 + Dimension + (0+Dimension|Facet), sigma ~ 0 + Dimension)}
#' for sparse multivariate data.
#' \item For \strong{lme4 backend}: Not applicable (lme4 does not support
#' multivariate models).
#' }
#' @param ... Additional arguments passed to the backend fitting function
#' (e.g., `lmer()` or `brm()`) and to `confint.merMod` for confidence intervals.
#'
#' @return An object of class "gstudy" containing:
#' \item{model}{The fitted model object from the backend}
#' \item{variance_components}{A tibble of estimated variance components}
#' \item{facets}{Character vector of facet names}
#' \item{facet_n}{Named numeric vector of sample sizes for each main effect facet}
#' \item{sample_size_info}{Comprehensive sample size information including main effects, interactions, residual, and nested effects}
#' \item{backend}{The backend used for fitting}
#' \item{is_multivariate}{Logical indicating if the model is multivariate}
#' \item{formula}{The formula used}
#' \item{data}{The original data}
#' \item{n_obs}{Number of observations}
#'
#' @seealso [dstudy()] for conducting D-studies
#'
#' @export
#'
#' @examples
#' # Basic univariate G-study with lme4 (default)
#' g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#'   (1 | Person:Task) + (1 | Person:Rater) + (1 | Task:Rater),
#' data = brennan
#' )
#'
#' # G-study with profile confidence intervals
#' \dontrun{
#' g_prof <- gstudy(Score ~ (1 | Person) + (1 | Task),
#'   data = brennan,
#'   ci_method = "profile"
#' )
#' }
#'
#' # Method of moments G-study (ANOVA-based)
#' g_mom <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#'   (1 | Person:Task) + (1 | Person:Rater) + (1 | Task:Rater),
#' data = brennan,
#' backend = "mom"
#' )
#'
#' \donttest{
#' # Bayesian G-study with brms
#' g_bayes <- gstudy(Score ~ (1 | Person) + (1 | Task),
#'   data = brennan,
#'   backend = "brms"
#' )
#'
#' # Bayesian G-study with custom priors
#' my_prior <- set_prior("normal(0, 1)", class = "sd", group = "Person")
#' g_bayes_prior <- gstudy(Score ~ (1 | Person) + (1 | Task),
#'   data = brennan,
#'   prior = my_prior,
#'   backend = "brms"
#' )
#'
#' # Multivariate G-study (automatically uses brms)
#' # Formatting data for multivariate example
#' b_wide <- tidyr::pivot_wider(brennan, names_from = Task, values_from = Score,
#'   names_prefix = "Task")
#' g_multi <- gstudy(mvbind(Task1, Task2) ~ (1 | Person) + (1 | Rater),
#'   data = b_wide
#' )
#' }
gstudy <- function(formula, data, backend = c("auto", "lme4", "brms", "mom"),
                   facets = NULL, nested = NULL, unbalanced = FALSE,
                   ci_method = c("none", "profile", "boot"),
                   nsim = 1000, boot.type = c("perc", "basic", "norm"),
                   prior = NULL, ...) {
  # 1. Match backend argument
  backend <- match.arg(backend)

  # 2. Match ci_method argument
  ci_method <- match.arg(ci_method)

  # 2.5 Match boot.type argument
  boot.type <- match.arg(boot.type)

  # 3. Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }

  # 4. Validate formula
  validate_formula(formula, backend)

  # 5. Detect if multivariate
  is_mv <- is_multivariate(formula)

  # 5.5. Detect if long-format multivariate
  long_format_info <- is_long_format_multivariate(formula)

  # If long-format detected, set is_mv to TRUE
  if (long_format_info$is_long) {
    is_mv <- TRUE
  }

  # 5.6. Check multivariate balance for wide-format models
  if (is_mv && !long_format_info$is_long) {
    balance_info <- check_multivariate_balance(formula, data, FALSE, unbalanced = unbalanced)

    # Only show warning if NOT in unbalanced mode (user explicitly opted in)
    if (!unbalanced && !is.null(balance_info$warning_message)) {
      warning(balance_info$warning_message, call. = FALSE)
    }

    if (!is.null(balance_info$info_message)) {
      message(balance_info$info_message)
    }
  }

  # 6. Select backend
  selected_backend <- select_backend(formula, backend)

  # 6.1. Validate unbalanced parameter for selected backend
  if (unbalanced) {
    if (selected_backend == "lme4") {
      if (is_mv) {
        stop(
          "Multivariate models are not supported by lme4 backend. ",
          "Use backend = 'brms' or backend = 'mom'.",
          call. = FALSE
        )
      } else {
        warning(
          "The 'unbalanced' parameter is only applicable to multivariate models. ",
          "It will be ignored for this univariate model.",
          call. = FALSE
        )
      }
    } else if (selected_backend == "brms") {
      warning(
        "unbalanced = TRUE is not implemented for brms backend. ",
        "For sparse/unbalanced multivariate data, use long-format specification:\n",
        "  bf(Score ~ 0 + Dimension + (0+Dimension|Facet), sigma ~ 0 + Dimension)\n",
        "See ?gstudy for long-format multivariate examples.",
        call. = FALSE
      )
    }
  }

  # 6.5. Validate interaction levels for lme4 backend
  validate_interaction_levels(formula, data, selected_backend)

  # 7. Warn if ci_method specified for brms backend
  if (ci_method != "none" && selected_backend == "brms") {
    warning(
      "ci_method is only applicable for lme4 backend. ",
      "brms provides CIs automatically from posterior samples.",
      call. = FALSE
    )
  }

  # 7b. Warn if ci_method specified for mom backend (has its own CIs)
  if (ci_method != "none" && selected_backend == "mom") {
    warning(
      "ci_method is only applicable for lme4 backend. ",
      "Method of moments provides CIs automatically using asymptotic approximations.",
      call. = FALSE
    )
  }

  # 7d. Validate prior is only used with brms backend
  if (!is.null(prior) && selected_backend != "brms") {
    stop(
      "prior is only supported with brms backend. ",
      "Use backend = 'brms' to use custom priors.",
      call. = FALSE
    )
  }

  # 8. Fit model
  model <- if (selected_backend == "lme4") {
    fit_lme4(formula, data, ...)
  } else if (selected_backend == "mom") {
    fit_mom(formula, data, unbalanced = unbalanced, ...)
  } else {
    fit_brms(formula, data, prior = prior, ...)
  }

  # 8.5. Handle long-format multivariate extraction
  if (long_format_info$is_long && selected_backend == "brms") {
    dimension_var <- long_format_info$dimension_var
    dimensions <- unique(data[[dimension_var]])

    # Calculate sample sizes per dimension
    sample_size_info_per_dim <- calculate_sample_size_info_per_dimension(
      formula, data, dimension_var
    )
    sample_size_tibble <- calculate_sample_size_tibble_per_dim(sample_size_info_per_dim)

    # Extract variance components with long-format handling
    vc <- extract_vc_brms_long_format(
      model,
      formula = formula,
      dimension_var = dimension_var,
      data = data
    )

    # Extract facet specifications FIRST (preserves formula order)
    facet_specs <- extract_facet_specs(formula)

    # Determine facets from variance components
    facets <- unique(vc$component[vc$component != "Residual"])

    # Determine object of measurement from formula order (first random effect)
    # For long-format models, facet_specs preserves the order from the formula
    object <- if (length(facet_specs) > 0) facet_specs[1] else if (length(facets) > 0) facets[1] else NULL

    # Extract covariances for long-format multivariate models
    # This is needed for composite coefficient calculation
    covariances <- extract_covariances_brms_long_format(model, dimension_var, data, vc)

    # For backward compatibility, also store in correlations slot
    # The build_component_covariance_matrix function uses the covariances structure
    correlations <- covariances

    # Detect estimation issues
    estimation_issues <- detect_brms_issues(model, vc)

    # Create mgstudy object
    result <- list(
      model = model,
      variance_components = vc,
      facets = facets,
      facet_specs = facet_specs,
      sample_size_info = calculate_sample_size_info(formula, data),
      sample_size_info_per_dim = sample_size_info_per_dim,
      sample_size_tibble = sample_size_tibble,
      backend = selected_backend,
      is_multivariate = TRUE,
      long_format_multivariate = TRUE,
      dimension_var = dimension_var,
      dimensions = dimensions,
      correlations = correlations,
      formula = formula,
      data = data,
      n_obs = nrow(data),
      object = object,
      estimation_issues = estimation_issues
    )
    class(result) <- "mgstudy"
    return(result)
  }

  # 9. Extract variance components
  vc <- extract_variance_components(
    model,
    selected_backend,
    ci_method = ci_method,
    nsim = nsim,
    boot.type = boot.type,
    formula = formula,
    ...
  )

  # 9.5. Extract facet specifications preserving user's original order
  facet_specs <- extract_facet_specs(formula)

  # 9.6. Reorder variance components to match formula specification
  # This ensures consistent ordering across different backends
  if (length(facet_specs) > 0) {
    vc <- reorder_variance_components(vc, facet_specs)
  }

  # 10. Parse facets if not provided
  if (is.null(facets)) {
    facets <- extract_facets(formula, vc)
  }

  # 10.6. Calculate sample sizes for each facet (main effects only)
  facet_n <- sapply(facets, function(f) {
    if (f %in% names(data)) {
      length(unique(data[[f]]))
    } else {
      NA
    }
  }, USE.NAMES = TRUE)

  # 10.7. Calculate comprehensive sample size information
  sample_size_info <- calculate_sample_size_info(formula, data, nested)

  # 10.8. Convert to tibble format for display
  sample_size_tibble <- calculate_single_sample_size_tibble(sample_size_info)

  # 11. Determine object of measurement (always first facet)
  object <- if (length(facets) > 0) facets[1] else NULL

  # 12. Extract dimension names
  dimensions <- unique(vc$dim)

  # 13. Extract correlations and covariances for multivariate models
  # Covariances are needed for composite coefficient calculation
  # Correlations are needed for display in summary()
  # Also propagate per-dimension sample-size metadata from mom fits so that
  # unbalanced designs can be honored by downstream D-study paths.
  correlations <- NULL
  is_unbalanced <- FALSE
  n_per_dim <- NULL
  sample_size_info_per_dim_mom <- NULL
  if (is_mv) {
    if (selected_backend == "brms") {
      cov_result <- extract_covariances_brms(model)
      cor_result <- extract_correlations_brms(model)

      correlations <- list(
        residual_cov = cov_result$residual_cov,
        residual_cov_matrix = cov_result$residual_cov_matrix,
        random_effect_cov = cov_result$random_effect_cov,
        random_effect_cov_matrix = cov_result$random_effect_cov_matrix,
        residual_cor = cor_result$residual_cor,
        residual_cor_matrix = cor_result$residual_cor_matrix,
        random_effect_cor = cor_result$random_effect_cor,
        random_effect_cor_matrix = cor_result$random_effect_cor_matrix
      )
    } else if (selected_backend == "mom") {
      if (!is.null(model$correlations)) {
        correlations <- model$correlations
      }
      if (isTRUE(model$is_unbalanced)) {
        is_unbalanced <- TRUE
        n_per_dim <- model$n_per_dim
        sample_size_info_per_dim_mom <- build_sample_size_info_per_dim_from_mom(model)
      }
    }
  }

  # 13.5. Detect and store estimation issues
  estimation_issues <- NULL
  if (selected_backend == "lme4") {
    estimation_issues <- detect_lme4_issues(model)
  } else if (selected_backend == "mom") {
    estimation_issues <- detect_mom_issues(model)
  } else if (selected_backend == "brms") {
    estimation_issues <- detect_brms_issues(model, vc)
  }

  # 14. Create gstudy or mgstudy object
  result <- list(
    model = model,
    variance_components = vc,
    facets = facets,
    facet_specs = facet_specs,
    facet_n = facet_n,
    sample_size_info = sample_size_info,
    sample_size_tibble = sample_size_tibble,
    sample_size_info_per_dim = sample_size_info_per_dim_mom,
    object = object,
    backend = selected_backend,
    is_multivariate = is_mv,
    is_unbalanced = is_unbalanced,
    n_per_dim = n_per_dim,
    formula = formula,
    data = data,
    n_obs = nrow(data),
    dimensions = dimensions,
    correlations = correlations,
    estimation_issues = estimation_issues
  )

  # Assign class based on whether multivariate
  if (is_mv) {
    class(result) <- "mgstudy"
  } else {
    class(result) <- "gstudy"
  }

  result
}
