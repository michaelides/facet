#' Estimation Issue Detection Functions
#'
#' Functions to detect estimation issues from different backends (lme4, mom, brms)
#' and issue appropriate warnings in dstudy.
#'
#' @name estimation-warnings
#' @keywords internal
NULL

#' Detect lme4 Estimation Issues
#'
#' Checks for convergence issues and singularity in lme4 models.
#'
#' @param model An lmerMod object from lme4::lmer().
#' @return A list with components:
#'   \item{convergence}{Logical indicating if convergence issues were detected}
#'   \item{singularity}{Logical indicating if the model is singular}
#'
#' @keywords internal
detect_lme4_issues <- function(model) {
  issues <- list()

  # Check convergence - lme4 stores convergence info in optinfo
  if (!is.null(model@optinfo$conv$lme4)) {
    conv_code <- model@optinfo$conv$lme4$code
    if (!is.null(conv_code) && conv_code != 0) {
      issues$convergence <- TRUE
      # Store the convergence message if available
      if (!is.null(model@optinfo$conv$lme4$messages)) {
        issues$convergence_messages <- model@optinfo$conv$lme4$messages
      }
    }
  }

  # Check singularity
  if (requireNamespace("lme4", quietly = TRUE)) {
    if (lme4::isSingular(model)) {
      issues$singularity <- TRUE
    }
  }

  issues
}

#' Detect Method of Moments Estimation Issues
#'
#' Checks for truncated variance estimates (when MS < residual MS, indicating
#' negative variance estimates that were truncated to zero).
#'
#' @param model A momfit object from fit_mom().
#' @return A list with components:
#'   \item{truncated_variance}{Logical indicating if variance was truncated}
#'   \item{truncated_components}{Character vector of components with truncated variance}
#'
#' @keywords internal
detect_mom_issues <- function(model) {
  issues <- list()

  # Get variance components from the model
  vc <- model$variance_components

  if (!is.null(vc) && "ms" %in% names(vc) && "resid_ms" %in% names(vc)) {
    # Find components where MS < residual MS (excluding Residual itself)
    truncated <- vc$ms < vc$resid_ms & vc$component != "Residual"
    truncated[is.na(truncated)] <- FALSE

    if (any(truncated)) {
      issues$truncated_variance <- TRUE
      issues$truncated_components <- vc$component[truncated]
    }
  }

  issues
}

#' Detect brms Estimation Issues
#'
#' Checks for low ESS, high R-hat, and treedepth issues in brms models.
#' ESS threshold is 100 per chain.
#'
#' @param model A brmsfit object from brms::brm().
#' @param vc Optional variance components tibble (to avoid re-extraction).
#' @return A list with components:
#'   \item{high_rhat}{Logical indicating if R-hat > 1.05 was detected}
#'   \item{high_rhat_components}{Character vector of components with high R-hat}
#'   \item{low_bulk_ess}{Logical indicating if low Bulk ESS was detected}
#'   \item{low_bulk_ess_components}{Character vector of components with low Bulk ESS}
#'   \item{low_tail_ess}{Logical indicating if low Tail ESS was detected}
#'   \item{low_tail_ess_components}{Character vector of components with low Tail ESS}
#'   \item{treedepth_exceeded}{Logical indicating if max treedepth was exceeded}
#'   \item{treedepth_count}{Number of iterations that hit max treedepth}
#'
#' @keywords internal
detect_brms_issues <- function(model, vc = NULL) {
  issues <- list()

  # Get number of chains
  n_chains <- 4 # default
  if (!is.null(model$fit) && !is.null(model$fit@sim)) {
    n_chains <- model$fit@sim$chains
  }
  ess_threshold <- 100 * n_chains

  # Extract variance components if not provided
  if (is.null(vc)) {
    vc <- extract_vc_brms(model)
  }

  # Check Rhat > 1.05
  if ("Rhat" %in% names(vc)) {
    high_rhat_idx <- !is.na(vc$Rhat) & vc$Rhat > 1.05
    if (any(high_rhat_idx)) {
      issues$high_rhat <- TRUE
      issues$high_rhat_components <- vc$component[high_rhat_idx]
    }
  }

  # Check low Bulk ESS
  if ("Bulk_ESS" %in% names(vc)) {
    low_bulk_ess_idx <- !is.na(vc$Bulk_ESS) & vc$Bulk_ESS < ess_threshold
    if (any(low_bulk_ess_idx)) {
      issues$low_bulk_ess <- TRUE
      issues$low_bulk_ess_components <- vc$component[low_bulk_ess_idx]
    }
  }

  # Check low Tail ESS
  if ("Tail_ESS" %in% names(vc)) {
    low_tail_ess_idx <- !is.na(vc$Tail_ESS) & vc$Tail_ESS < ess_threshold
    if (any(low_tail_ess_idx)) {
      issues$low_tail_ess <- TRUE
      issues$low_tail_ess_components <- vc$component[low_tail_ess_idx]
    }
  }

  # Check treedepth exceeded using brms::check_hmc_diagnostics
  if (requireNamespace("brms", quietly = TRUE)) {
    tryCatch(
      {
        diag <- brms::check_hmc_diagnostics(model)
        if (!is.null(diag) && !is.null(diag$max_treedepth)) {
          treedepth_count <- sum(diag$max_treedepth)
          if (treedepth_count > 0) {
            issues$treedepth_exceeded <- TRUE
            issues$treedepth_count <- treedepth_count
          }
        }
      },
      error = function(e) {
        # Silently ignore errors in diagnostic checking
      })
  }

  issues
}

#' Check and Warn About Estimation Issues
#'
#' Called by dstudy() to check for estimation issues in the gstudy object
#' and issue appropriate warnings.
#'
#' @param gstudy_obj A gstudy or mgstudy object.
#'
#' @keywords internal
check_estimation_issues <- function(gstudy_obj) {
  # Handle old gstudy objects without estimation_issues field
  issues <- gstudy_obj$estimation_issues
  if (is.null(issues) || length(issues) == 0) {
    return(invisible(NULL))
  }

  backend <- gstudy_obj$backend

  if (backend == "lme4") {
    # Check convergence
    if (!is.null(issues$convergence) && isTRUE(issues$convergence)) {
      warning(
        "Model convergence issue detected in gstudy. ",
        "The g and phi coefficients cannot be trusted.",
        call. = FALSE
      )
    }

    # Check singularity
    if (!is.null(issues$singularity) && isTRUE(issues$singularity)) {
      warning(
        "Model is singular (some variance components are zero). ",
        "The g and phi coefficients cannot be trusted.",
        call. = FALSE
      )
    }

  } else if (backend == "mom") {
    # Check truncated variance
    if (!is.null(issues$truncated_variance) && isTRUE(issues$truncated_variance)) {
      components_str <- paste(issues$truncated_components, collapse = ", ")
      warning(
        "Negative variance estimates truncated to zero for: ", components_str, ". ",
        "The g and phi coefficients cannot be trusted.",
        call. = FALSE
      )
    }

  } else if (backend == "brms") {
    # Get threshold for ESS warnings
    n_chains <- 4
    if (!is.null(gstudy_obj$model$fit) && !is.null(gstudy_obj$model$fit@sim)) {
      n_chains <- gstudy_obj$model$fit@sim$chains
    }
    ess_threshold <- 100 * n_chains

    # Check high R-hat
    if (!is.null(issues$high_rhat) && isTRUE(issues$high_rhat)) {
      components_str <- paste(issues$high_rhat_components, collapse = ", ")
      warning(
        "R-hat > 1.05 detected for: ", components_str, ". ",
        "This indicates non-convergence. ",
        "The g and phi coefficients cannot be trusted.",
        call. = FALSE
      )
    }

    # Check low Bulk ESS
    if (!is.null(issues$low_bulk_ess) && isTRUE(issues$low_bulk_ess)) {
      components_str <- paste(issues$low_bulk_ess_components, collapse = ", ")
      warning(
        "Low bulk effective sample size (ESS < ", ess_threshold, ") detected for: ",
        components_str, ". ",
        "The g and phi coefficients cannot be trusted.",
        call. = FALSE
      )
    }

    # Check low Tail ESS
    if (!is.null(issues$low_tail_ess) && isTRUE(issues$low_tail_ess)) {
      components_str <- paste(issues$low_tail_ess_components, collapse = ", ")
      warning(
        "Low tail effective sample size (ESS < ", ess_threshold, ") detected for: ",
        components_str, ". ",
        "The g and phi coefficients cannot be trusted.",
        call. = FALSE
      )
    }

    # Check treedepth exceeded
    if (!is.null(issues$treedepth_exceeded) && isTRUE(issues$treedepth_exceeded)) {
      warning(
        "Maximum treedepth exceeded during sampling for ", issues$treedepth_count,
        " iterations. Consider increasing adapt_delta or reparameterizing. ",
        "The g and phi coefficients may be unreliable.",
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}
