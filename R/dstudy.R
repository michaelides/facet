#' Conduct a Decision (D) Study
#'
#' Compute generalizability and dependability coefficients based on the
#' results of a G-study. Allows for exploring different measurement designs
#' by modifying the number of levels for each facet.
#'
#' @param gstudy_obj An object of class "gstudy" from [gstudy()].
#' @param n A named list specifying the number of levels for each facet
#' in the D-study design (e.g., `list(items = 10, raters = 3)`).
#' Can also be a list of vectors to explore multiple sample sizes.
#' @param universe Specification for components that contribute to the universe score.
#' The object of measurement (first component from the G-study formula) is always
#' included in the universe. Can be:
#' - NULL (default): universe includes only the object of measurement
#' - A character string: "person:item" (object is automatically added)
#' - A character vector: c("person", "person:item") (object should be included)
#' - A formula: ~ person + person:item
#' Non-object components in the universe are scaled by their corresponding sample sizes.
#' For example, if the object is "p" and universe includes "p:o", the universe score
#' variance is computed as: var(p) + var(p:o) / n_o.
#' @param error Specification for error components. Can be:
#' - NULL (default): all components not in universe (and not in aggregation)
#' - A character string: "person:item"
#' - A character vector: c("person:item", "person:rater")
#' - A formula: ~ person:item + person:rater
#' Note: If a component is specified in both `universe` and `error`, an error is raised.
#' Note: If a component is specified in both `aggregation` and `error`, a warning
#' is issued and the component is removed from error (interaction terms are preserved).
#' @param aggregation Character vector of facets to aggregate over.
#' Components containing these facets will be divided by the sample size.
#' For example, `aggregation = "item"` divides item and person:item by n_item.
#' When aggregation is specified, the main effects of aggregation facets are
#' automatically excluded from error variance (but interaction terms are kept).
#' Default is NULL (no aggregation).
#' @param residual_is Character string specifying which facets make up the residual.
#' For example, "person:item:rater" means the residual represents the
#' three-way interaction. Required when aggregation is specified and you want
#' the residual to be rescaled. Default is NULL (residual is not rescaled).
#' @param estimation Character string specifying how to calculate coefficients:
#' - "simple": uses point estimates from variance components (for lme4/mom backends)
#' - "posterior": uses full posterior distributions (for brms backend)
#' For brms backend, posterior estimation is always used to ensure consistency
#' between variance component estimates and coefficient calculations. This avoids
#' Jensen's inequality bias that would occur if variance estimates were computed
#' as mean(SD)^2 rather than mean(SD^2).
#' When estimation = "posterior" is requested with non-brms backend,
#' a warning is issued and the gstudy model is refit with backend = "brms".
#' @param cut_score Optional numeric value specifying a cutoff score for criterion-referenced
#' decisions. When provided, calculates phi-cut coefficient (phi_cut) in addition to standard
#' phi coefficient. For multivariate models, can be a single value applied to all dimensions.
#' Default is NULL (no phi-cut calculation).
#' @param ci Character vector specifying which coefficients to compute credible intervals for.
#'   Options: "g", "phi", "phi-cut". Can specify multiple: `ci = c("g", "phi")`.
#'   Credible intervals are only available when using the brms backend.
#'   Default is NULL (no credible intervals computed).
#' @param probs Numeric vector of length 2 specifying the quantiles for credible interval
#' calculation. Default is `c(0.025, 0.975)` for a 95% credible interval.
#' Works like the `quantile()` function. Only used when `ci` is specified.
#' @param weights Numeric vector of weights for computing composite coefficients.
#' Length must match the number of dimensions in the multivariate design.
#' Default is NULL, which uses equal weights (1 for each dimension).
#' Only applicable for multivariate G-studies (mgstudy objects).
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class "dstudy" containing:
#' \item{gstudy}{The original G-study object}
#' \item{variance_components}{A tibble of variance components for the D-study with columns:
#' \itemize{
#' \item component: Name of variance component
#' \item var_unscaled: Unscaled variance estimate (from G-study)
#' \item pct_unscaled: Percentage of total unscaled variance
#' \item var_scaled: Scaled variance (divided by D-study sample sizes)
#' \item pct_scaled: Percentage of total scaled variance
#' \item dim: Dimension/response variable (for multivariate models)
#' }}
#' \item{coefficients}{A tibble with G and D coefficients. For multivariate designs,
#' the coefficients tibble includes an additional row with dim = "Composite" showing
#' the weighted composite G and Phi coefficients.}
#' \item{n}{The number of levels for each facet}
#' \item{object}{The object of measurement (first component from G-study)}
#' \item{universe}{The universe components specification}
#' \item{error}{The error specification (if provided)}
#' \item{aggregation}{The aggregation specification (if provided)}
#' \item{residual_is}{The residual specification (if provided)}
#' \item{residual_composition}{The facets that make up the residual (from parse_residual_facets)}
#' \item{estimation}{The estimation method used ("simple" or "posterior")}
#' \item{posterior}{List of posterior distributions (only when estimation = "posterior")}
#' \item{cut_score}{The cutoff score used for phi-cut calculation (if provided)}
#' \item{mu_y}{The grand mean(s) used for phi-cut calculation (if cut_score provided)}
#' \item{ci}{The credible interval specification (if provided)}
#' \item{probs}{The probability levels used for credible interval calculation (if ci provided)}
#'
#' @seealso [gstudy()] for conducting G-studies
#'
#' @export
#'
#' @details
#' ## Universe, Error, and Aggregation
#'
#' By default, the *universe score* variance contains only the object of measurement
#' (e.g., `Person`), and all remaining components contribute to *error* variance.
#'
#' - Use `universe` to include additional components in the universe score
#' (e.g., a facet that is fixed in your applied setting).
#' - Use `error` to restrict which components count as error
#' (e.g., when some facets are considered fixed and not a source of
#' measurement error for your inference).
#' - Use `aggregation` when the measurement procedure involves averaging over
#' a facet (e.g., averaging across multiple raters), so that the variance
#' components for that facet are divided by their respective sample sizes
#' before computing coefficients.
#'
#' ## Composite Coefficients for Multivariate Designs
#'
#' For multivariate G-studies (mgstudy objects), composite coefficients provide
#' an overall reliability estimate for a weighted sum of dimension scores.
#' The composite score is:
#'
#' `Y_composite = sum(w_d * Y_d)`
#'
#' The variance of this composite for any variance component is computed as:
#'
#' `sigma^2_composite = w' * Sigma * w`
#'
#' Where:
#' - `w` is the vector of weights (one per dimension)
#' - `Sigma` is the variance-covariance matrix for that component
#'
#' Expanded, this equals:
#'
#' `sigma^2_composite = sum(w_d^2 * sigma^2_d) + 2 * sum(w_d * w_d' * sigma_dd')`
#'
#' The covariances (off-diagonal elements) contribute when dimensions are correlated.
#' This means composite reliability often exceeds the average of dimension-specific
#' reliabilities, as the composite "borrows strength" from correlations.
#'
#' The `weights` parameter controls the contribution of each dimension:
#' - Default: `NULL` uses equal weights (1 for each dimension)
#' - Custom: Provide a numeric vector with length matching the number of dimensions
#'
#' For posterior estimation with brms backend, composite coefficients are computed
#' for each posterior draw, properly propagating uncertainty to the final estimates.
#'
#' @examples
#' # First conduct a G-study using the brennan dataset
#' # (Person crossed with Task and Rater)
#' g <- gstudy(Score ~ (1|Person) + (1|Task) + (1|Rater) +
#'             (1|Person:Task) + (1|Person:Rater) + (1|Task:Rater),
#'             data = brennan)
#'
#' # D-study with specific sample sizes
#' d <- dstudy(g, n = list(Task = 3, Rater = 4))
#' print(d)
#'
#' # D-study exploring multiple sample sizes (sweep)
#' d_sweep <- dstudy(g, n = list(Task = c(3, 5, 10), Rater = c(2, 4, 8)))
#' print(d_sweep)
#'
#' # D-study with aggregation (averaging over Raters)
#' d_agg <- dstudy(g, n = list(Task = 3, Rater = 4),
#'   aggregation = "Rater",
#'   residual_is = "Person:Task:Rater")
#'
#' \donttest{
#' # D-study with posterior estimation (requires brms backend)
#' g_brms <- gstudy(Score ~ (1|Person) + (1|Task) + (1|Rater),
#'                  data = brennan, backend = "brms")
#' d_post <- dstudy(g_brms, n = list(Task = 3, Rater = 4))
#'
#' # D-study with credible intervals (requires brms backend)
#' d_ci <- dstudy(g_brms, n = list(Task = 3, Rater = 4), ci = c("g", "phi"))
#' print(d_ci)
#'
#' # Custom probability levels (90% credible interval)
#' d_ci_90 <- dstudy(g_brms, n = list(Task = 3, Rater = 4),
#'                   ci = "g", probs = c(0.05, 0.95))
#' }
dstudy <- function(gstudy_obj, n = list(), universe = NULL,
  error = NULL, aggregation = NULL, residual_is = NULL,
  estimation = NULL, cut_score = NULL, ci = NULL,
  probs = c(0.025, 0.975), weights = NULL, ...) {
  # 1. Validate input
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

  # 1.1. Check for estimation issues and warn
  check_estimation_issues(gstudy_obj)

# 1.5. Set default estimation based on backend
  # For brms backend, always use posterior draws-based estimation
  # For lme4/mom backends, use simple estimation
  if (is.null(estimation)) {
    if (gstudy_obj$backend == "brms") {
      estimation <- "posterior"
    } else {
      estimation <- "simple"
    }
  } else {
    estimation <- match.arg(estimation, c("simple", "posterior"))
  }

  # 1.6. If posterior estimation requested but not brms backend, warn and refit
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

  # 1.7. Warn if simple estimation requested for brms backend
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

  # 1.8. Extract grand mean for phi-cut calculation if cut_score provided
  mu_y <- NULL
  if (!is.null(cut_score)) {
    mu_y <- extract_grand_mean(gstudy_obj)
  }

  # 1.9. Validate ci parameter
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

# 1.10. Validate probs parameter
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
  
  # 1.11. Warn if phi-cut CI requested but no cut_score provided
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

  # 2. Get variance components from G-study
vc <- gstudy_obj$variance_components

# 3. Get object of measurement (always first component from G-study)
object <- gstudy_obj$object
object_spec <- parse_specification(object)

# 3.5. Process universe specification
universe_spec <- parse_specification(universe)

# Default: universe = object only
if (is.null(universe) || length(universe_spec) == 0) {
universe_spec <- object_spec
} else {
# Ensure object is in universe
if (!all(object_spec %in% universe_spec)) {
warning(
"The specification of the universe did not include the object of measurement '",
paste(object_spec, collapse = ", "),
"'. It has been added automatically.",
call. = FALSE
)
universe_spec <- unique(c(object_spec, universe_spec))
}

# Check if non-object universe components interact with the object
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

# 3.6. Validate that universe doesn't overlap with error
error_spec <- parse_specification(error)
agg_spec <- parse_specification(aggregation)

# Check universe vs error (exact specification match)
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

# 4. Calculate n from G-study data if not provided
n_provided <- length(n) > 0
if (length(n) == 0) {
  n <- extract_sample_sizes(gstudy_obj)
  n_provided <- TRUE  # n was extracted, so treat as provided
  message(
    "No sample sizes provided in 'n'. Using G-study sample sizes: ",
    paste(names(n), n, sep = " = ", collapse = ", ")
  )
}

# 5. Check if n contains vectors (for sweeping)
is_sweep <- any(sapply(n, length) > 1)

# 8. Determine residual composition from formula and data
residual_composition <- parse_residual_facets(
  gstudy_obj$formula,
  gstudy_obj$data
)

  # 9. Use residual_composition as default for residual_is if not specified
  residual_is_effective <- if (is.null(residual_is)) residual_composition else residual_is

  # 10. Check if residual matches any error component specified by user
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

if (estimation == "posterior") {
# Use posterior-based calculation
if (is_sweep) {
n_grid <- expand.grid(n, stringsAsFactors = FALSE)

# Only calculate both unscaled and scaled estimates when n was not provided by user
if (!n_provided) {
# Calculate unscaled estimates using posterior
posterior_results_unscaled <- calculate_coefficients_posterior(
    gstudy_obj = gstudy_obj,
    n = n,
    object = object,
    universe = universe_spec,
    error = error,
    aggregation = aggregation,
    residual_is = residual_is_effective,
    is_sweep = TRUE,
    n_grid = n_grid,
    n_provided = n_provided,
    use_scaled = FALSE,
    cut_score = cut_score,
    mu_y = mu_y,
    ci = ci,
    probs = probs,
    weights = weights
  )
    unscaled_coefs <- posterior_results_unscaled$coefficients
    unscaled_coefs$estimate <- "unscaled"

# Calculate scaled estimates using posterior
  posterior_results_scaled <- calculate_coefficients_posterior(
    gstudy_obj = gstudy_obj,
    n = n,
    object = object,
    universe = universe_spec,
    error = error,
    aggregation = aggregation,
    residual_is = residual_is_effective,
    is_sweep = TRUE,
    n_grid = n_grid,
    n_provided = n_provided,
    use_scaled = TRUE,
    cut_score = cut_score,
    mu_y = mu_y,
    ci = ci,
    probs = probs,
    weights = weights
  )
scaled_coefs <- posterior_results_scaled$coefficients
scaled_coefs$estimate <- "scaled"

# Combine - need to align columns
coef_cols <- intersect(names(unscaled_coefs), names(scaled_coefs))
n_cols <- setdiff(names(unscaled_coefs), coef_cols)
n_cols <- n_cols[!n_cols %in% c("estimate")]

coefficients <- rbind(
unscaled_coefs[, c(n_cols, "estimate", setdiff(coef_cols, c(n_cols, "estimate")))],
scaled_coefs[, c(n_cols, "estimate", setdiff(coef_cols, c(n_cols, "estimate")))]
)
coefficients <- tibble::as_tibble(coefficients)
# Reorder columns to put estimate first
coefficients <- coefficients[, c("estimate", setdiff(names(coefficients), "estimate"))]

posterior <- posterior_results_scaled$posterior
} else {
  # When n is provided, only return scaled estimates (no estimate column)
  posterior_results <- calculate_coefficients_posterior(
    gstudy_obj = gstudy_obj,
    n = n,
    object = object,
    universe = universe_spec,
    error = error,
    aggregation = aggregation,
    residual_is = residual_is_effective,
    is_sweep = TRUE,
    n_grid = n_grid,
    n_provided = n_provided,
    use_scaled = TRUE,
    cut_score = cut_score,
    mu_y = mu_y,
    ci = ci,
    probs = probs,
    weights = weights
  )
  coefficients <- posterior_results$coefficients
  posterior <- posterior_results$posterior
}
  } else {
    # Single sample size
# Only calculate both unscaled and scaled estimates when n was not provided by user
if (!n_provided) {
  # Calculate unscaled estimates using posterior
  posterior_results_unscaled <- calculate_coefficients_posterior(
    gstudy_obj = gstudy_obj,
    n = n,
    object = object,
    universe = universe_spec,
    error = error,
    aggregation = aggregation,
    residual_is = residual_is_effective,
    is_sweep = FALSE,
    n_provided = n_provided,
    use_scaled = FALSE,
    cut_score = cut_score,
    mu_y = mu_y,
    ci = ci,
    probs = probs,
    weights = weights
  )
  unscaled_coefs <- posterior_results_unscaled$coefficients
  unscaled_coefs$estimate <- "unscaled"

  # Calculate scaled estimates using posterior
  posterior_results_scaled <- calculate_coefficients_posterior(
    gstudy_obj = gstudy_obj,
      n = n,
      object = object,
      universe = universe_spec,
      error = error,
      aggregation = aggregation,
      residual_is = residual_is_effective,
      is_sweep = FALSE,
      n_provided = n_provided,
      use_scaled = TRUE,
      cut_score = cut_score,
      mu_y = mu_y,
      ci = ci,
      probs = probs
    )
    scaled_coefs <- posterior_results_scaled$coefficients
    scaled_coefs$estimate <- "scaled"

has_dim <- "dim" %in% names(unscaled_coefs)
# Dynamically include all coefficient columns (including CI columns like g_LL, g_UL, etc.)
n_cols_pattern <- "^n_"
coef_cols <- c("estimate", "dim", "uni", "sigma2_delta", "sigma2_delta_abs",
"g", "phi", "phi_cut", "sem_rel", "sem_abs")
# Add any CI columns that might be present
ci_cols <- grep("_LL|_UL$", names(unscaled_coefs), value = TRUE)
coef_cols <- c(coef_cols, ci_cols)
coef_cols <- intersect(coef_cols, names(unscaled_coefs))

coefficients <- tibble::as_tibble(rbind(
unscaled_coefs[, coef_cols, drop = FALSE],
scaled_coefs[, coef_cols, drop = FALSE]
))

posterior <- posterior_results_scaled$posterior
} else {
    # When n is provided, only return scaled estimates (no estimate column)
    posterior_results <- calculate_coefficients_posterior(
gstudy_obj = gstudy_obj,
    n = n,
    object = object,
    universe = universe_spec,
    error = error,
    aggregation = aggregation,
    residual_is = residual_is_effective,
    is_sweep = FALSE,
    n_provided = n_provided,
    use_scaled = TRUE,
    cut_score = cut_score,
    mu_y = mu_y,
    ci = ci,
    probs = probs,
    weights = weights
  )
      coefficients <- posterior_results$coefficients
      posterior <- posterior_results$posterior
    }
  }

  # Initialize composite variables for all paths
  composite_vc <- NULL
  composite_post <- NULL

  # Extract composite variance components from posterior_results (single sample size path)
  if (estimation == "posterior" && !is_sweep && n_provided && is_multivariate && 
      length(dimensions) > 1 && exists("posterior_results") && !is.null(posterior_results)) {
    composite_vc <- posterior_results$composite_vc
    composite_post <- posterior_results$composite_posterior
  }

  # Calculate composite coefficients for posterior estimation path
  if (estimation == "posterior" && is_multivariate && length(dimensions) > 1) {
  if (is_sweep) {
    # Sweep with posterior estimation
    composite_posterior_results <- lapply(seq_len(nrow(n_grid)), function(i) {
      n_current <- as.list(n_grid[i, , drop = FALSE])
      n_current <- lapply(n_current, function(x) as.numeric(x))
      
      # Calculate composite coefficients using the original variance components
      # and scaling appropriately for each sample size combination
      composite_coefs <- calculate_composite_coefficients(
        vc = vc,  # Use original unscaled variance components
        n = n_current,
        weights = weights,
        object = object,
        error = error,
        aggregation = aggregation,
        residual_is = residual_is_effective,
        universe = universe_spec,
        correlations = gstudy_obj$correlations,
        cut_score = cut_score,
        mu_y = mu_y
      )
      
      cbind(data.frame(n_current, stringsAsFactors = FALSE), composite_coefs)
    })
    
    composite_combined <- do.call(rbind, composite_posterior_results)
    coefficients <- dplyr::bind_rows(coefficients, composite_combined)
  } else {
    # Single sample size with posterior estimation
    if (!n_provided) {
      # Unscaled composite
      composite_unscaled <- calculate_composite_coefficients(
        vc = vc,
        n = NULL,
        weights = weights,
        object = object,
        error = error,
        aggregation = aggregation,
        residual_is = residual_is_effective,
        universe = universe_spec,
        correlations = gstudy_obj$correlations,
        cut_score = cut_score,
        mu_y = mu_y
      )
      composite_unscaled$estimate <- "unscaled"
      
      # Scaled composite
      composite_scaled <- calculate_composite_coefficients(
        vc = calculate_dstudy_variance(vc, n, object, aggregation, TRUE, residual_is_effective, facet_n = gstudy_obj$facet_n),
        n = n,
        weights = weights,
        object = object,
        error = error,
        aggregation = aggregation,
        residual_is = residual_is_effective,
        universe = universe_spec,
        correlations = gstudy_obj$correlations,
        cut_score = cut_score,
        mu_y = mu_y
      )
      composite_scaled$estimate <- "scaled"
      
      coefficients <- dplyr::bind_rows(coefficients, composite_unscaled, composite_scaled)
    } else {
      # n was provided
      composite_coefs <- calculate_composite_coefficients(
        vc = calculate_dstudy_variance(vc, n, object, aggregation, n_provided, residual_is_effective, facet_n = gstudy_obj$facet_n),
        n = n,
        weights = weights,
        object = object,
        error = error,
        aggregation = aggregation,
        residual_is = residual_is_effective,
        universe = universe_spec,
        correlations = gstudy_obj$correlations,
        cut_score = cut_score,
        mu_y = mu_y
      )
      coefficients <- dplyr::bind_rows(coefficients, composite_coefs)
    }
  }
}

  # For posterior, calculate variance components with both scaled and unscaled
  # Remove diagnostic columns from brms backend
  # But only if NOT in sweep mode (sweep mode uses base variance components)
  if (!is_sweep) {
    d_vc <- calculate_dstudy_variance(vc, n, object, aggregation, TRUE, residual_is_effective, facet_n = gstudy_obj$facet_n)
  } else {
    # For sweep mode, use base variance components (sweep handles scaling separately)
    d_vc <- vc
  }
  
  # Append composite variance components for multivariate posterior estimation
  if (!is_sweep && is_multivariate && length(dimensions) > 1 && !is.null(composite_vc)) {
    d_vc <- dplyr::bind_rows(d_vc, composite_vc)
  }
  } else if (is_sweep) {
    # Create grid of all sample size combinations
    n_grid <- expand.grid(n, stringsAsFactors = FALSE)

    # Calculate coefficients for each combination
    results <- lapply(seq_len(nrow(n_grid)), function(i) {
      n_current <- as.list(n_grid[i, , drop = FALSE])
      n_current <- lapply(n_current, function(x) as.numeric(x))

      if (!n_provided) {
        # When n not provided: use divided variance for scaled estimates
        # Scaled: divide variance components by sample sizes of non-object facets
        vc_scaled <- calculate_divided_variance(vc, n_current, object, residual_is_effective)
scaled_coefs <- calculate_divided_coefficients(vc_scaled, object, error,
aggregation, residual_is_effective, universe_spec, cut_score, mu_y)

# Unscaled: use original variance components directly
unscaled_coefs <- calculate_coefficients(vc, n_current, object, error,
aggregation, residual_is_effective, universe_spec, cut_score, mu_y)

        # Add estimate type and combine
        unscaled_coefs$estimate <- "unscaled"
        scaled_coefs$estimate <- "scaled"

        combined <- rbind(unscaled_coefs, scaled_coefs)

        # Reorder columns to put estimate first (after n columns if sweep)
        combined <- combined[, c("estimate", setdiff(names(combined), "estimate"))]

        # Combine with sample sizes
        cbind(data.frame(n_current, stringsAsFactors = FALSE), combined)
      } else {
        # When n is provided: use standard D-study variance calculation
        d_vc <- calculate_dstudy_variance(vc, n_current, object, aggregation, n_provided, residual_is_effective, facet_n = gstudy_obj$facet_n)

        # Calculate coefficients (scaled only, no estimate column)
scaled_coefs <- calculate_coefficients(d_vc, n_current, object, error,
aggregation, residual_is_effective, universe_spec, cut_score, mu_y)

        cbind(data.frame(n_current, stringsAsFactors = FALSE), scaled_coefs)
      }
    })

coefficients <- do.call(rbind, results)
  coefficients <- tibble::as_tibble(coefficients)

  # Calculate composite coefficients for multivariate designs (sweep path)
  if (is_multivariate && length(dimensions) > 1) {
    composite_results <- lapply(seq_len(nrow(n_grid)), function(i) {
      n_current <- as.list(n_grid[i, , drop = FALSE])
      n_current <- lapply(n_current, function(x) as.numeric(x))

      # Get variance components for this iteration
      if (!n_provided) {
        d_vc_iter <- calculate_dstudy_variance(vc, n_current, object, aggregation, FALSE, residual_is_effective, facet_n = gstudy_obj$facet_n)
      } else {
        d_vc_iter <- calculate_dstudy_variance(vc, n_current, object, aggregation, n_provided, residual_is_effective, facet_n = gstudy_obj$facet_n)
      }

      composite_coefs <- calculate_composite_coefficients(
        vc = d_vc_iter,
        n = n_current,
        weights = weights,
        object = object,
        error = error,
        aggregation = aggregation,
        residual_is = residual_is_effective,
        universe = universe_spec,
        correlations = gstudy_obj$correlations,
        cut_score = cut_score,
        mu_y = mu_y
      )

      if ("estimate" %in% names(coefficients)) {
        est_type <- coefficients$estimate[1]
        composite_coefs$estimate <- est_type
      }

      cbind(data.frame(n_current, stringsAsFactors = FALSE), composite_coefs)
    })

    composite_combined <- do.call(rbind, composite_results)
    coefficients <- dplyr::bind_rows(coefficients, composite_combined)
  }

  # For sweeping, variance_components is the base (unscaled) version (sweep handles scaling separately)
 # Remove diagnostic columns if present (brms, mom backends)
 d_vc <- vc %>%
 select(-any_of(c('error', 'se', 'lower', 'upper', 'sd', 'Rhat', 'Bulk_ESS', 'Tail_ESS')))
 } else {
    # 11. Non-sweep path: calculate coefficients
    if (!n_provided) {
      # When n not provided: use divided variance for scaled estimates
      # Scaled: divide variance components by sample sizes of non-object facets
      vc_scaled <- calculate_divided_variance(vc, n, object, residual_is_effective)
scaled_coefs <- calculate_divided_coefficients(vc_scaled, object, error,
aggregation, residual_is_effective, universe_spec, cut_score, mu_y)

# Unscaled: use original variance components directly
# Pass n = NULL so that universe components are NOT scaled
unscaled_coefs <- calculate_coefficients(vc, n = NULL, object, error,
aggregation, residual_is_effective, universe_spec, cut_score, mu_y)

      # Combine into long format
      unscaled_coefs$estimate <- "unscaled"
      scaled_coefs$estimate <- "scaled"

      has_dim <- "dim" %in% names(scaled_coefs)
      if (has_dim) {
        coef_cols <- c("estimate", "dim", "uni", "sigma2_delta", "sigma2_delta_abs",
                       "g", "phi", "phi_cut", "sem_rel", "sem_abs")
      } else {
        coef_cols <- c("estimate", "uni", "sigma2_delta", "sigma2_delta_abs",
                       "g", "phi", "phi_cut", "sem_rel", "sem_abs")
      }
      coef_cols <- intersect(coef_cols, names(scaled_coefs))

        coefficients <- tibble::as_tibble(rbind(
      unscaled_coefs[, coef_cols, drop = FALSE],
      scaled_coefs[, coef_cols, drop = FALSE]
    ))

    # For non-sweep with n not provided, still calculate both scaled and unscaled
    d_vc <- calculate_dstudy_variance(vc, n, object, aggregation, TRUE, residual_is_effective, facet_n = gstudy_obj$facet_n)

    # Calculate composite coefficients for multivariate designs
    if (is_multivariate && length(dimensions) > 1) {
      # Need to compute composite for both unscaled and scaled
      # Unscaled composite
      composite_unscaled <- calculate_composite_coefficients(
        vc = vc,
        n = NULL,
        weights = weights,
        object = object,
        error = error,
        aggregation = aggregation,
        residual_is = residual_is_effective,
        universe = universe_spec,
        correlations = gstudy_obj$correlations,
        cut_score = cut_score,
        mu_y = mu_y
      )
      composite_unscaled$estimate <- "unscaled"

      # Scaled composite
      composite_scaled <- calculate_composite_coefficients(
        vc = d_vc,
        n = n,
        weights = weights,
        object = object,
        error = error,
        aggregation = aggregation,
        residual_is = residual_is_effective,
        universe = universe_spec,
        correlations = gstudy_obj$correlations,
        cut_score = cut_score,
        mu_y = mu_y
      )
      composite_scaled$estimate <- "scaled"

      coefficients <- dplyr::bind_rows(coefficients, composite_unscaled, composite_scaled)
    }
    } else {
# When n IS provided: use standard D-study variance calculation
  d_vc <- calculate_dstudy_variance(vc, n, object, aggregation, n_provided, residual_is_effective, facet_n = gstudy_obj$facet_n)

  # Calculate coefficients (scaled only, no estimate column)
  coefficients <- calculate_coefficients(d_vc, n, object, error, aggregation, residual_is_effective, universe_spec, cut_score, mu_y)

  # Calculate composite coefficients for multivariate designs
  if (is_multivariate && length(dimensions) > 1) {
    composite_coefs <- calculate_composite_coefficients(
      vc = d_vc,
      n = n,
      weights = weights,
      object = object,
      error = error,
      aggregation = aggregation,
      residual_is = residual_is_effective,
      universe = universe_spec,
      correlations = gstudy_obj$correlations,
      cut_score = cut_score,
      mu_y = mu_y
    )

    # Add estimate column if present in other coefficients
    if ("estimate" %in% names(coefficients)) {
      composite_coefs$estimate <- coefficients$estimate[1]
    }

    # Append composite row to coefficients
    coefficients <- dplyr::bind_rows(coefficients, composite_coefs)
  }
}
}

# 13. Create dstudy object
result <- list(
  gstudy = gstudy_obj,
  variance_components = d_vc,
  coefficients = coefficients,
  n = n,
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
  is_multivariate = is_multivariate,
  cut_score = cut_score,
  mu_y = mu_y,
  ci = ci,
  probs = if (!is.null(ci)) probs else NULL,
  weights = weights
)

  class(result) <- "dstudy"
  result
}

extract_sample_sizes <- function(gstudy_obj) {
  facet_n <- gstudy_obj$facet_n
  
  if (is.null(facet_n) || length(facet_n) == 0) {
    data <- gstudy_obj$data
    facets <- gstudy_obj$facets
    
    if (is.null(data) || is.null(facets)) {
      return(list())
    }
    
    n_list <- list()
    for (facet in facets) {
      if (facet %in% names(data)) {
        n_list[[facet]] <- length(unique(data[[facet]]))
      }
    }
    
    n_list
  } else {
    as.list(facet_n)
  }
}