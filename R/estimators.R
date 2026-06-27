#' Estimator Management for G-Studies
#'
#' Functions to manage and interact with different fitting estimators (lme4, brms).
#' These are internal functions used by [gstudy()] and [dstudy()].
#'
#' @name estimators
#' @keywords internal
NULL

#' Select the Appropriate Estimator
#'
#' Determines which estimator to use based on formula characteristics and user preference.
#' If estimator is "auto", selects brms for multivariate formulas and lme4 for univariate.
#'
#' @param formula The model formula.
#' @param estimator Character string specifying the desired estimator: "auto", "lme4", or "brms".
#' @return Character string indicating the estimator to use.
#'
#' @keywords internal
select_estimator <- function(formula, estimator = c("auto", "lme4", "brms", "aov")) {
  estimator <- match.arg(estimator)

  # Check if multivariate
  is_mv <- is_multivariate(formula)

  if (estimator == "auto") {
    # Use brms for multivariate, lme4 for univariate (mom is not default)
    if (is_mv) {
      # Prefer brms for multivariate but allow mom
      if (requireNamespace("brms", quietly = TRUE)) {
        return("brms")
      } else {
        # Fall back to mom if brms not available
        return("aov")
      }
    } else {
      # Default to lme4 for univariate (faster)
      if (requireNamespace("lme4", quietly = TRUE)) {
        return("lme4")
      } else if (requireNamespace("brms", quietly = TRUE)) {
        return("brms")
      } else {
        return("aov")
      }
    }
  }

  # Check if requested estimator is available
  if (estimator == "lme4" && !requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for the lme4 estimator.", call. = FALSE)
  }
  if (estimator == "brms" && !requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for the brms estimator.", call. = FALSE)
  }

  # aov estimator is always available (uses base R)
  # No additional checks needed

  # Check compatibility
  if (is_mv && estimator == "lme4") {
    stop(
      "Multivariate formulas (containing mvbind or set_rescor) require brms or aov estimator.\n",
      "Use: estimator = 'brms' or estimator = 'aov'",
      call. = FALSE
    )
  }

  estimator
}

#' Fit a Model Using lme4
#'
#' Fits a mixed effects model using lme4::lmer(). Converts brms-style formulas
#' to lme4 format if necessary.
#'
#' @param formula A formula for the model (lme4 style).
#' @param data A data frame.
#' @param ... Additional arguments passed to [lme4::lmer()].
#' @return A fitted lmerMod object.
#'
#' @keywords internal
fit_lme4 <- function(formula, data, ...) {
  # Ensure lme4 is available
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required.", call. = FALSE)
  }

  # Convert formula if it's a brmsformula
  if (inherits(formula, "brmsformula")) {
    formula <- formula$formula
  }

  # Check for multivariate and stop if present
  if (is_multivariate(formula)) {
    stop(
      "Multivariate formulas cannot be fit with lme4. Use brms estimator instead.",
      call. = FALSE
    )
  }

  # Fit the model
  model <- tryCatch(
    {
      lme4::lmer(formula = formula, data = data, ...)
    },
    error = function(e) {
      # Provide more informative error message
      if (grepl("convergence", e$message, ignore.case = TRUE)) {
        warning(
          "Model convergence issue. Consider scaling variables or simplifying the model.",
          call. = FALSE
        )
      }
      stop("Error fitting model with lme4: ", e$message, call. = FALSE)
    }
  )

  model
}

#' Fit a Model Using brms
#'
#' Fits a Bayesian mixed effects model using brms::brm(). Handles both
#' univariate and multivariate formulas.
#'
#' @param formula A formula for the model (can be brmsformula).
#' @param data A data frame.
#' @param prior A brmsprior object or list of priors as created by [set_prior()]
#'   or related functions. If NULL, default priors are used.
#' @param ... Additional arguments passed to [brms::brm()].
#' @return A fitted brmsfit object.
#'
#' @keywords internal
fit_brms <- function(formula, data, prior = NULL, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required.", call. = FALSE)
  }

  formula_char <- paste(deparse(formula), collapse = " ")
  has_set_rescor_in_formula <- grepl("set_rescor\\s*\\(", formula_char)

  if (has_set_rescor_in_formula && !inherits(formula, "brmsformula")) {
    stop(
      "'set_rescor()' should not be part of the formula right-hand side.\n",
      "Use bf() to create a brmsformula and add set_rescor() with '+':\n",
      "  bf(mvbind(y1, y2) ~ (1|person)) + set_rescor(TRUE)\n",
      "Or use brms::bf() if bf is not available:\n",
      "  brms::bf(mvbind(y1, y2) ~ (1|person)) + brms::set_rescor(TRUE)",
      call. = FALSE
    )
  }

  if (is_multivariate(formula) && !inherits(formula, "brmsformula") && !inherits(formula, "mvbrmsformula")) {
    formula <- brms::bf(formula)
  }

  args <- list(...)

  if (is.null(args$chains)) {
    args$chains <- 4
  }
  if (is.null(args$iter)) {
    args$iter <- 2000
  }
  if (is.null(args$refresh)) {
    args$refresh <- args$iter / 10
  }

  args$formula <- formula
  args$data <- data
  if (!is.null(prior)) {
    args$prior <- prior
  }

  model <- tryCatch(
    {
      do.call(brms::brm, args)
    },
    error = function(e) {
      stop("Error fitting model with brms: ", e$message, call. = FALSE)
    }
  )

  model
}

#' Extract Variance Components from a Fitted Model
#'
#' Generic function that dispatches to the appropriate estimator-specific extractor.
#'
#' @param model A fitted model object.
#' @param estimator Character string indicating the estimator used.
#' @param ci_method Character string specifying the CI method for lme4 estimator.
#'   One of "none", "profile", or "boot". Default is "none".
#' @param nsim Integer: number of bootstrap simulations (only for ci_method = "boot").
#'   Default is 1000.
#' @param boot.type Character: bootstrap type, "perc" or "basic" (only for ci_method = "boot").
#'   Default is "perc".
#' @param ... Additional arguments passed to confint.merMod.
#' @return A tibble of variance components with columns:
#'   \item{component}{Name of the variance component}
#'   \item{facet}{Associated facet name}
#'   \item{type}{Type: "main", "interaction", or "residual"}
#'   \item{var}{Point estimate of variance}
#'   \item{pct}{Percentage of total variance}
#'
#'   For lme4 with ci_method != "none" or brms estimator, also includes:
#'   \item{lower}{Lower confidence interval bound}
#'   \item{upper}{Upper confidence interval bound}
#'
#'   For brms estimator, also includes:
#' \item{error}{Standard error of the estimate}
#' \item{sd}{Standard deviation}
#' \item{Rhat}{Convergence diagnostic}
#' \item{ESS}{Effective sample size}
#'
#' @keywords internal
extract_variance_components <- function(model, estimator, ci_method = "none",
                                        nsim = 1000, boot.type = "perc",
                                        formula = NULL, ...) {
  if (estimator == "lme4") {
    extract_vc_lme4(model,
      ci_method = ci_method, nsim = nsim,
      boot.type = boot.type, formula = formula, ...
    )
  } else if (estimator == "brms") {
    extract_vc_brms(model, formula = formula)
  } else if (estimator == "aov") {
    extract_vc_aov(model, formula = formula)
  } else {
    stop("Unknown estimator: ", estimator, call. = FALSE)
  }
}
