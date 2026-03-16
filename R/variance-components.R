#' Variance Component Extraction and Manipulation
#'
#' Functions for extracting, summarizing, and manipulating variance components
#' from fitted models in generalizability theory analyses.
#'
#' @name variance-components
#' @keywords internal
NULL

#' Extract Variance Components from a Model
#'
#' Generic function for extracting variance components from different
#' model types.
#'
#' @param model A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @return A tibble with columns: component, variance, percent.
#'
#' @keywords internal
extract_vc <- function(model, ...) {
  UseMethod("extract_vc")
}

#' Extract Variance Components from lme4 Model
#'
#' Extracts variance components from an lme4 model using VarCorr.
#' Calculates confidence intervals using confint.merMod.
#'
#' @param model An lmerMod object from lme4::lmer().
#' @param ci_method Method for confidence intervals: "none", "profile", or "boot".
#'   - "none": No confidence intervals (default).
#'   - "profile": Profile likelihood (more accurate, slower).
#'   - "boot": Parametric bootstrap (most accurate, slowest).
#' @param conf_level Confidence level for intervals (default 0.95).
#' @param nsim Number of bootstrap simulations (if ci_method = "boot", default 1000).
#' @param boot.type Bootstrap type: "perc" (percentile), "basic", or "norm" (normal-theory) (default "perc").
#' @param formula The original formula (used to extract response name).
#' @param ... Additional arguments passed to confint.merMod.
#' @return A tibble with columns: component, dim, type, var, pct.
#'   If ci_method != "none", also includes: lower, upper.
#'
#' @keywords internal
extract_vc_lme4 <- function(model, ci_method = "none", conf_level = 0.95,
                            nsim = 1000, boot.type = "perc", formula = NULL, ...) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required.", call. = FALSE)
  }

  ci_method <- match.arg(ci_method, c("none", "profile", "boot"))

  vc <- lme4::VarCorr(model)

  component_names <- names(vc)

  dim_name <- "response"
  if (!is.null(formula)) {
    dim_name <- extract_response_names(formula)
    if (length(dim_name) == 0) dim_name <- "response"
  }

  results <- lapply(component_names, function(nm) {
    vc_comp <- vc[[nm]]

    variance <- as.numeric(attr(vc_comp, "stddev"))^2

    type <- if (grepl(":", nm)) "interaction" else "main"

    data.frame(
      component = nm,
      dim = dim_name,
      type = type,
      var = variance,
      stringsAsFactors = FALSE
    )
  })

  residual_var <- attr(vc, "sc")^2
  results[[length(results) + 1]] <- data.frame(
    component = "Residual",
    dim = dim_name,
    type = "residual",
    var = residual_var,
    stringsAsFactors = FALSE
  )

  vc_df <- do.call(rbind, results)
  vc_tibble <- tibble::as_tibble(vc_df)

  total_var <- sum(vc_tibble$var)
  vc_tibble$pct <- (vc_tibble$var / total_var) * 100

  if (ci_method != "none") {
    vc_tibble$lower <- NA_real_
    vc_tibble$upper <- NA_real_

    tryCatch(
      {
        ci_args <- list(
          object = model,
          method = ci_method
        )

        if (ci_method == "boot") {
          ci_args$nsim <- nsim
          ci_args$boot.type <- boot.type
          ci_args$.progress <- "txt"
        }

        ci_args <- modifyList(ci_args, list(...))

        ci_result <- do.call(lme4::confint.merMod, ci_args)

        ci_rows <- grep("^\\.sig|^\\.sigma", rownames(ci_result))

        if (length(ci_rows) > 0) {
          ci_sd <- ci_result[ci_rows, , drop = FALSE]

          n_components <- length(component_names)

          for (i in seq_along(component_names)) {
            if (i <= nrow(ci_sd) - 1) {
              lower_sd <- ci_sd[i, 1]
              upper_sd <- ci_sd[i, 2]

              if (!is.na(lower_sd) && !is.na(upper_sd)) {
                lower_sd <- max(0, lower_sd)
                vc_tibble$lower[i] <- lower_sd^2
                vc_tibble$upper[i] <- upper_sd^2
              }
            }
          }

          residual_row <- nrow(ci_sd)
          if (residual_row > 0) {
            lower_sd <- ci_sd[residual_row, 1]
            upper_sd <- ci_sd[residual_row, 2]

            if (!is.na(lower_sd) && !is.na(upper_sd)) {
              lower_sd <- max(0, lower_sd)
              residual_idx <- nrow(vc_tibble)
              vc_tibble$lower[residual_idx] <- lower_sd^2
              vc_tibble$upper[residual_idx] <- upper_sd^2
            }
          }
        }
      },
      error = function(e) {
        warning(
          "Confidence intervals could not be computed using method '",
          ci_method, "': ", e$message,
          "\nReturning variance components without CIs.",
          call. = FALSE
        )
      },
      warning = function(w) {
        if (grepl("convergence|singular", w$message, ignore.case = TRUE)) {
          warning(
            "Model may be singular or have convergence issues. ",
            "Confidence intervals may be unreliable.\n",
            "Original warning: ", w$message,
            call. = FALSE
          )
        }
      }
    )
  }

  vc_tibble
}

#' Extract Variance Components from brms Model
#'
#' Extracts variance components from a brms model using posterior samples.
#' Computes variance estimates as mean(SD^2) from posterior draws, ensuring
#' proper uncertainty quantification and avoiding Jensen's inequality bias
#' that would occur with mean(SD)^2.
#'
#' @param model A brmsfit object from brms::brm().
#' @param conf_level Credible interval level (default 0.95).
#' @param formula The original formula (used to extract response names).
#' @return A tibble with columns:
#' \itemize{
#' \item component: Name of the variance component
#' \item dim: Response variable name (for multivariate models)
#' \item type: "main", "interaction", or "residual"
#' \item var: Mean of variance draws (mean(SD^2))
#' \item error: Standard error of variance (sd(variance_draws))
#' \item lower: 2.5th percentile of variance draws
#' \item upper: 97.5th percentile of variance draws
#' \item sd: Mean of SD draws (for display purposes)
#' \item Rhat: Convergence diagnostic
#' \item Bulk_ESS: Effective sample size (bulk)
#' \item Tail_ESS: Effective sample size (tail)
#' \item pct: Percentage of total variance
#' }
#'
#' @keywords internal
extract_vc_brms <- function(model, conf_level = 0.95, formula = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required.", call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required.", call. = FALSE)
  }

  model_summary <- suppressWarnings(summary(model))
  random_summary <- model_summary$random
  spec_pars <- model_summary$spec_pars

  vc <- brms::VarCorr(model)

  draws <- brms::as_draws_matrix(model)

  resp_names <- NULL

  if (!is.null(model$formula) && inherits(model$formula, "mvbrmsformula")) {
    resp_names <- model$formula$responses
  }

  if (is.null(resp_names) && !is.null(model$formula) && !is.null(model$formula$formula)) {
    resp_names <- tryCatch(
      {
        bform <- model$formula$formula
        if (inherits(bform, "list") && length(bform) > 0) {
          names(bform)
        } else {
          NULL
        }
      },
      error = function(e) NULL
    )
  }

  if (is.null(resp_names) && !is.null(formula)) {
    resp_names <- extract_response_names(formula)
  }

  is_mv <- !is.null(resp_names) && length(resp_names) > 1

  if (!is_mv) {
    if (is.null(resp_names) || length(resp_names) == 0) {
      resp_names <- "response"
    }
  }

  normalize_group_name <- function(name) {
    name <- gsub("__", ":", name)
    name
  }

  vc_groups <- names(vc)
  normalized_vc_groups <- vapply(vc_groups, normalize_group_name, character(1))
  names(normalized_vc_groups) <- vc_groups

  expected_facet_specs <- character()
  if (!is.null(formula)) {
    expected_facet_specs <- extract_facet_specs(formula)
  }

  results <- list()

  for (grp in names(vc)) {
    normalized_grp <- normalized_vc_groups[grp]
    type <- if (grp == "residual__") "residual" else if (grepl(":", normalized_grp)) "interaction" else "main"

    if (is_mv) {
      for (resp in resp_names) {
        result <- extract_single_variance_from_draws(
          draws = draws,
          grp = grp,
          resp = resp,
          type = type,
          random_summary = random_summary,
          spec_pars = spec_pars,
          is_mv = TRUE
        )
        if (!is.null(result)) {
          if (grp != "residual__" && grp != normalized_grp) {
            result$component <- normalized_grp
          }
          results[[length(results) + 1]] <- result
        }
      }
    } else {
      result <- extract_single_variance_from_draws(
        draws = draws,
        grp = grp,
        resp = if (length(resp_names) == 1) resp_names[1] else "response",
        type = type,
        random_summary = random_summary,
        spec_pars = spec_pars,
        is_mv = FALSE
      )
      if (!is.null(result)) {
        if (grp != "residual__" && grp != normalized_grp) {
          result$component <- normalized_grp
        }
        results[[length(results) + 1]] <- result
      }
    }
  }

  if (length(expected_facet_specs) > 0) {
    processed_components <- vapply(results, function(r) r$component[1], character(1))
    
    for (spec in expected_facet_specs) {
      if (!(spec %in% processed_components)) {
        grp_to_try <- spec
        
        type <- if (grepl(":", spec)) "interaction" else "main"
        
        if (is_mv) {
          for (resp in resp_names) {
            result <- extract_single_variance_from_draws(
              draws = draws,
              grp = grp_to_try,
              resp = resp,
              type = type,
              random_summary = random_summary,
              spec_pars = spec_pars,
              is_mv = TRUE
            )
            if (!is.null(result)) {
              result$component <- spec
              results[[length(results) + 1]] <- result
            }
          }
        } else {
          result <- extract_single_variance_from_draws(
            draws = draws,
            grp = grp_to_try,
            resp = if (length(resp_names) == 1) resp_names[1] else "response",
            type = type,
            random_summary = random_summary,
            spec_pars = spec_pars,
            is_mv = FALSE
          )
          if (!is.null(result)) {
            result$component <- spec
            results[[length(results) + 1]] <- result
          }
        }
      }
    }
  }

  vc_df <- do.call(rbind, results)
  vc_tibble <- tibble::as_tibble(vc_df)

  total_var <- sum(vc_tibble$var, na.rm = TRUE)
  vc_tibble$pct <- (vc_tibble$var / total_var) * 100

  vc_tibble
}

extract_single_variance_from_draws <- function(draws, grp, resp, type,
                                                random_summary, spec_pars, is_mv) {
  if (grp == "residual__") {
    if (is_mv && resp != "response") {
      param_name <- paste0("sigma_", resp)
    } else {
      param_name <- "sigma"
    }
    component <- "Residual"
  } else {
    if (is_mv && resp != "response") {
      param_name <- paste0("sd_", grp, "__", resp, "_Intercept")
    } else {
      param_name <- paste0("sd_", grp, "__Intercept")
    }
    component <- grp
  }

  if (!param_name %in% colnames(draws)) {
    param_candidates <- character()
    
    if (grp != "residual__") {
      grp_double_underscore <- gsub(":", "__", grp)
      if (is_mv && resp != "response") {
        param_candidates <- c(
          paste0("sd_", grp_double_underscore, "__", resp, "_Intercept"),
          paste0("sd_", grp, "__", resp, "_Intercept")
        )
      } else {
        param_candidates <- c(
          paste0("sd_", grp_double_underscore, "__Intercept"),
          paste0("sd_", grp, "__Intercept")
        )
      }
      param_candidates <- unique(param_candidates)
    } else {
      param_candidates <- if (is_mv && resp != "response") paste0("sigma_", resp) else "sigma"
    }
    
    found <- FALSE
    for (cand in param_candidates) {
      if (cand %in% colnames(draws)) {
        param_name <- cand
        found <- TRUE
        break
      }
    }
    
    if (!found) {
      base_name <- if (grp == "residual__") {
        if (is_mv && resp != "response") paste0("sigma_", resp) else "sigma"
      } else {
        paste0("sd_", gsub(":", "__", grp))
      }
      matching_params <- grep(paste0("^", base_name), colnames(draws), value = TRUE)
      if (length(matching_params) > 0) {
        param_name <- matching_params[1]
      } else {
        return(NULL)
      }
    }
  }

  sd_draws <- posterior::extract_variable(draws, param_name)

  var_draws <- sd_draws^2

  estimate <- mean(var_draws)
  lower <- as.numeric(quantile(var_draws, 0.025))
  upper <- as.numeric(quantile(var_draws, 0.975))
  se <- sd(var_draws)
  sd_mean <- mean(sd_draws)

  rhat_val <- NA_real_
  bulk_ess_val <- NA_real_
  tail_ess_val <- NA_real_

  if (grp == "residual__") {
    if (!is.null(spec_pars)) {
      sigma_name <- if (is_mv && resp != "response") paste0("sigma_", resp) else "sigma"
      if (sigma_name %in% rownames(spec_pars)) {
        rhat_val <- spec_pars[sigma_name, "Rhat"]
        bulk_ess_val <- spec_pars[sigma_name, "Bulk_ESS"]
        tail_ess_val <- spec_pars[sigma_name, "Tail_ESS"]
      }
    }
    if (is.na(rhat_val) && !is.null(random_summary) && "residual__" %in% names(random_summary)) {
      resid_rand <- random_summary[["residual__"]]
      resp_row <- if (is_mv && resp != "response") resp else "sigma"
      if (resp_row %in% rownames(resid_rand)) {
        rhat_val <- resid_rand[resp_row, "Rhat"]
        bulk_ess_val <- resid_rand[resp_row, "Bulk_ESS"]
        tail_ess_val <- resid_rand[resp_row, "Tail_ESS"]
      }
    }
  } else {
    if (!is.null(random_summary) && grp %in% names(random_summary)) {
      grp_rand <- random_summary[[grp]]
      if (param_name %in% rownames(grp_rand)) {
        rhat_val <- grp_rand[param_name, "Rhat"]
        bulk_ess_val <- grp_rand[param_name, "Bulk_ESS"]
        tail_ess_val <- grp_rand[param_name, "Tail_ESS"]
        } else {
          alt_row <- gsub("^sd_", "sd(", param_name)
          alt_row <- paste0(alt_row, ")")
          if (alt_row %in% rownames(grp_rand)) {
            rhat_val <- grp_rand[alt_row, "Rhat"]
            bulk_ess_val <- grp_rand[alt_row, "Bulk_ESS"]
            tail_ess_val <- grp_rand[alt_row, "Tail_ESS"]
          } else {
            mv_row <- gsub("^sd_[^_]+__", "sd(", param_name)
            mv_row <- paste0(mv_row, ")")
            if (mv_row %in% rownames(grp_rand)) {
              rhat_val <- grp_rand[mv_row, "Rhat"]
              bulk_ess_val <- grp_rand[mv_row, "Bulk_ESS"]
              tail_ess_val <- grp_rand[mv_row, "Tail_ESS"]
            }
          }
        }
    }
  }

  data.frame(
    component = component,
    dim = resp,
    type = type,
    var = estimate,
    error = se,
    lower = lower,
    upper = upper,
    sd = sd_mean,
    Rhat = rhat_val,
    Bulk_ESS = bulk_ess_val,
    Tail_ESS = tail_ess_val,
    stringsAsFactors = FALSE
  )
}

#' @rdname extract_vc
#' @keywords internal
extract_vc.merMod <- function(model, ...) {
  extract_vc_lme4(model, ...)
}

#' @rdname extract_vc
#' @keywords internal
extract_vc.brmsfit <- function(model, ...) {
  extract_vc_brms(model, ...)
}

#' Compute Variance Component Percentages
#'
#' Calculates the percentage of total variance accounted for by each component.
#'
#' @param vc A tibble of variance components with a "var" column.
#' @return The input tibble with an added "pct" column.
#'
#' @keywords internal
compute_vc_percentages <- function(vc) {
  if (!"var" %in% names(vc)) {
    stop("Variance components tibble must have a 'var' column", call. = FALSE)
  }

  total_var <- sum(vc$var, na.rm = TRUE)
  vc$pct <- (vc$var / total_var) * 100

  vc
}

#' Summarize Variance Components for dstudy Objects
#'
#' Creates a summary table of variance components for dstudy objects,
#' showing both unscaled and scaled estimates.
#'
#' @param vc A tibble of variance components from a dstudy object.
#' @param digits Number of decimal places for rounding.
#' @return A formatted tibble.
#'
#' @keywords internal
summarize_vc_dstudy <- function(vc, digits = 3) {
  required_cols <- c("component", "var_unscaled", "pct_unscaled",
                     "var_scaled", "pct_scaled")
  missing_cols <- setdiff(required_cols, names(vc))
  if (length(missing_cols) > 0) {
    stop("dstudy variance components tibble is missing columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  summary_df <- vc[, c("component", "var_unscaled", "pct_unscaled",
                       "var_scaled", "pct_scaled")]

  summary_df$var_unscaled <- round(summary_df$var_unscaled, digits)
  summary_df$pct_unscaled <- round(summary_df$pct_unscaled, digits)
  summary_df$var_scaled <- round(summary_df$var_scaled, digits)
  summary_df$pct_scaled <- round(summary_df$pct_scaled, digits)

  summary_df
}

#' Summarize Variance Components
#'
#' Creates a summary table of variance components suitable for reporting.
#'
#' @param vc A tibble of variance components.
#' @param digits Number of decimal places for rounding.
#' @param scale Scale for displaying results: "variance" (default) or "sd".
#' @return A formatted tibble.
#'
#' @keywords internal
summarize_vc <- function(vc, digits = 3, scale = c("variance", "sd")) {
  # Check if this is a dstudy variance components tibble
  if (all(c("var_unscaled", "var_scaled", "pct_unscaled", "pct_scaled") %in% names(vc))) {
    return(summarize_vc_dstudy(vc, digits = digits))
  }

  scale <- match.arg(scale)

  required_cols <- c("component", "var", "pct")
  missing_cols <- setdiff(required_cols, names(vc))
  if (length(missing_cols) > 0) {
    stop("Variance components tibble is missing columns: ",
         paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  summary_df <- vc[, c("component", "var", "pct")]

  if (scale == "sd") {
    if ("sd" %in% names(vc)) {
      summary_df$estimate <- round(vc$sd, digits)
    } else {
      summary_df$estimate <- round(sqrt(vc$var), digits)
    }

    if (all(c("lower", "upper") %in% names(vc))) {
      if ("sd" %in% names(vc)) {
        summary_df$`2.5% CI` <- round(sqrt(vc$lower), digits)
        summary_df$`97.5% CI` <- round(sqrt(vc$upper), digits)
      } else {
        summary_df$`2.5% CI` <- round(sqrt(vc$lower), digits)
        summary_df$`97.5% CI` <- round(sqrt(vc$upper), digits)
      }
    }
  } else {
    summary_df$estimate <- round(vc$var, digits)

    if (all(c("lower", "upper") %in% names(vc))) {
      summary_df$`2.5% CI` <- round(vc$lower, digits)
      summary_df$`97.5% CI` <- round(vc$upper, digits)
    }
  }

  if ("error" %in% names(vc)) {
    summary_df$SE <- round(vc$error, digits)
  }

  if ("Rhat" %in% names(vc)) {
    summary_df$Rhat <- round(vc$Rhat, digits)
  }

  if ("Bulk_ESS" %in% names(vc)) {
    summary_df$Bulk_ESS <- vc$Bulk_ESS
  }

  if ("Tail_ESS" %in% names(vc)) {
    summary_df$Tail_ESS <- vc$Tail_ESS
  }

  summary_df$pct <- round(summary_df$pct, digits)

  col_order <- c("component", "dim", "estimate", "pct")
  ci_cols <- c("2.5% CI", "97.5% CI")
  diag_cols <- c("SE", "Rhat", "Bulk_ESS", "Tail_ESS")

  for (col in ci_cols) {
    if (col %in% names(summary_df)) {
      col_order <- c(col_order, col)
    }
  }
  for (col in diag_cols) {
    if (col %in% names(summary_df)) {
      col_order <- c(col_order, col)
    }
  }

  final_cols <- col_order[col_order %in% names(summary_df)]
  summary_df <- summary_df[, final_cols, drop = FALSE]

  summary_df
}

#' Summarize Correlation Tibble for Display
#'
#' Formats a correlation tibble for printing, similar to summarize_vc.
#'
#' @param cor_tibble A tibble with columns: dim1, dim2, estimate, se, lower, upper, Rhat, Bulk_ESS, Tail_ESS
#' @param digits Number of digits for rounding (default 3).
#' @return A tibble formatted for display.
#'
#' @keywords internal
summarize_cor <- function(cor_tibble, digits = 3) {
  if (is.null(cor_tibble) || nrow(cor_tibble) == 0) {
    return(cor_tibble)
  }

  summary_df <- cor_tibble

  summary_df$estimate <- round(summary_df$estimate, digits)
  summary_df$se <- round(summary_df$se, digits)
  summary_df$`2.5% CI` <- round(summary_df$lower, digits)
  summary_df$`97.5% CI` <- round(summary_df$upper, digits)

  if ("Rhat" %in% names(summary_df)) {
    summary_df$Rhat <- round(summary_df$Rhat, digits)
  }

  col_order <- c("dim1", "dim2", "estimate", "se", "2.5% CI", "97.5% CI", "Rhat", "Bulk_ESS", "Tail_ESS")
  final_cols <- col_order[col_order %in% names(summary_df)]
  summary_df <- summary_df[, final_cols, drop = FALSE]

  summary_df
}

#' Summarize Covariance Tibble for Display
#'
#' Formats a covariance tibble for printing with rounded values.
#'
#' @param cov_tibble A tibble with covariance estimates.
#' @param digits Number of digits to round to.
#'
#' @return A formatted data frame for printing.
#'
#' @keywords internal
summarize_cov <- function(cov_tibble, digits = 3) {
  if (is.null(cov_tibble) || nrow(cov_tibble) == 0) {
    return(cov_tibble)
  }

  summary_df <- cov_tibble

  summary_df$estimate <- round(summary_df$estimate, digits)
  summary_df$se <- round(summary_df$se, digits)
  summary_df$`2.5% CI` <- round(summary_df$lower, digits)
  summary_df$`97.5% CI` <- round(summary_df$upper, digits)

  if ("Rhat" %in% names(summary_df)) {
    summary_df$Rhat <- round(summary_df$Rhat, digits)
  }

  col_order <- c("dim1", "dim2", "estimate", "se", "2.5% CI", "97.5% CI", "Rhat", "Bulk_ESS", "Tail_ESS")
  final_cols <- col_order[col_order %in% names(summary_df)]
  summary_df <- summary_df[, final_cols, drop = FALSE]

  summary_df
}

#' @rdname extract_vc
#' @keywords internal
extract_vc.brmsfit <- function(model, ...) {
  extract_vc_brms(model, ...)
}

#' Compute Variance Component Percentages
#'
#' Calculates the percentage of total variance accounted for by each component.
#'
#' @param vc A tibble of variance components with a "var" column.
#' @return The input tibble with an added "pct" column.
#'
#' @keywords internal
compute_vc_percentages <- function(vc) {
  if (!"var" %in% names(vc)) {
    stop("Variance components tibble must have a 'var' column",
      call. = FALSE
    )
  }

  total_var <- sum(vc$var, na.rm = TRUE)
  vc$pct <- (vc$var / total_var) * 100

  vc
}

#' Reorder Variance Components to Match Formula Order
#'
#' Reorders the variance components tibble to follow the user-specified
#' formula order, with Residual always last. Also renames interaction
#' components to match the user's specification (e.g., "Rater:Task" vs "Task:Rater").
#'
#' @param vc A tibble of variance components
#' @param facet_specs Character vector of facet specifications in user order
#' @return The vc tibble reordered to match facet_specs, with Residual last
#'
#' @keywords internal
reorder_variance_components <- function(vc, facet_specs) {
  if (is.null(vc) || nrow(vc) == 0 || length(facet_specs) == 0) {
    return(vc)
  }

  is_residual <- grepl("^Residual", vc$component)
  residual_rows <- vc[is_residual, , drop = FALSE]
  non_residual <- vc[!is_residual, , drop = FALSE]

  if (nrow(non_residual) == 0) {
    return(rbind(non_residual, residual_rows))
  }

  name_mapping <- list()
  for (spec in facet_specs) {
    if (grepl(":", spec)) {
      parts <- strsplit(spec, ":")[[1]]
      parts <- trimws(parts)
      alt_ordering <- paste(rev(parts), collapse = ":")
      if (alt_ordering %in% unique(non_residual$component) &&
          !(spec %in% unique(non_residual$component))) {
        name_mapping[[alt_ordering]] <- spec
      }
    }
  }

  if (length(name_mapping) > 0) {
    for (old_name in names(name_mapping)) {
      new_name <- name_mapping[[old_name]]
      idx <- non_residual$component == old_name
      if (any(idx)) {
        non_residual$component[idx] <- new_name
      }
    }
  }

  unique_components <- unique(non_residual$component)

  component_order <- character()
  for (spec in facet_specs) {
    if (spec %in% unique_components) {
      component_order <- c(component_order, spec)
      next
    }

    if (grepl(":", spec)) {
      parts <- strsplit(spec, ":")[[1]]
      parts <- trimws(parts)
      alt_ordering <- paste(rev(parts), collapse = ":")
      if (alt_ordering %in% unique_components) {
        component_order <- c(component_order, alt_ordering)
        next
      }
    }
  }

  idx <- match(component_order, unique_components)
  idx <- idx[!is.na(idx)]

  missing <- setdiff(seq_along(unique_components), idx)
  if (length(missing) > 0) {
    idx <- c(idx, missing)
  }

  ordered_non_residual <- non_residual[non_residual$component %in% unique_components[idx], ]
  ordered_non_residual$component <- factor(ordered_non_residual$component, levels = unique_components[idx])
  ordered_non_residual <- ordered_non_residual[order(ordered_non_residual$component), ]
  ordered_non_residual$component <- as.character(ordered_non_residual$component)

  rbind(ordered_non_residual, residual_rows)
}

#' Extract Correlation Matrices from Multivariate brms Models
#'
#' Extracts residual correlations and correlated random effects from
#' multivariate brms models.
#'
#' @param model A brmsfit object from brms::brm().
#' @return A list containing:
#' \item{residual_cor}{Tibble with columns: dim1, dim2, estimate, se, lower, upper, Rhat, Bulk_ESS, Tail_ESS}
#' \item{random_effect_cor}{Named list of tibbles for random effects (same columns as residual_cor)}
#' \item{residual_cor_matrix}{Residual correlation matrix (for backward compatibility)}
#' \item{random_effect_cor_matrix}{Named list of correlation matrices (for backward compatibility)}
#'
#' @keywords internal
extract_correlations_brms <- function(model) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required.", call. = FALSE)
  }

  vc <- brms::VarCorr(model)
  model_summary <- suppressWarnings(summary(model))
  random_summary <- model_summary$random
  spec_pars <- model_summary$spec_pars
  rescor_pars <- model_summary$rescor_pars

  resp_names <- NULL
  if (!is.null(model$formula) && inherits(model$formula, "mvbrmsformula")) {
    resp_names <- model$formula$responses
  }

  if (is.null(resp_names) && !is.null(model$formula) && !is.null(model$formula$formula)) {
    bform <- model$formula$formula
    if (inherits(bform, "list") && length(bform) > 0) {
      resp_names <- names(bform)
    }
  }

  result <- list(
    residual_cor = NULL,
    random_effect_cor = list(),
    residual_cor_matrix = NULL,
    random_effect_cor_matrix = list()
  )

  if (is.null(resp_names) || length(resp_names) <= 1) {
    return(result)
  }

  cor_array_to_tibble <- function(cor_array, resp_names, grp_name = NULL, random_summary, rescor_pars) {
    n <- length(resp_names)
    rows <- list()

    for (i in 2:n) {
      for (j in 1:(i - 1)) {
        estimate <- cor_array[i, "Estimate", j]
        se <- cor_array[i, "Est.Error", j]
        lower <- cor_array[i, "Q2.5", j]
        upper <- cor_array[i, "Q97.5", j]

        rhat_val <- NA_real_
        bulk_ess_val <- NA_real_
        tail_ess_val <- NA_real_

        if (is.null(grp_name) || grp_name == "residual__") {
          if (!is.null(rescor_pars) && nrow(rescor_pars) > 0) {
            possible_names <- c(
              paste0("rescor(", resp_names[j], ",", resp_names[i], ")"),
              paste0("rescor(", resp_names[i], ",", resp_names[j], ")")
            )
            for (pname in possible_names) {
              if (pname %in% rownames(rescor_pars)) {
                rhat_val <- rescor_pars[pname, "Rhat"]
                bulk_ess_val <- rescor_pars[pname, "Bulk_ESS"]
                tail_ess_val <- rescor_pars[pname, "Tail_ESS"]
                break
              }
            }
          }
        } else {
          if (!is.null(grp_name) && grp_name %in% names(random_summary)) {
            grp_summary <- random_summary[[grp_name]]
            if (!is.null(grp_summary) && nrow(grp_summary) > 0) {
              possible_names <- c(
                paste0("cor(", resp_names[j], "_Intercept,", resp_names[i], "_Intercept)"),
                paste0("cor(", resp_names[i], "_Intercept,", resp_names[j], "_Intercept)")
              )
              for (pname in possible_names) {
                if (pname %in% rownames(grp_summary)) {
                  rhat_val <- grp_summary[pname, "Rhat"]
                  bulk_ess_val <- grp_summary[pname, "Bulk_ESS"]
                  tail_ess_val <- grp_summary[pname, "Tail_ESS"]
                  break
                }
              }
            }
          }
        }

        rows[[length(rows) + 1]] <- data.frame(
          dim1 = resp_names[j],
          dim2 = resp_names[i],
          estimate = estimate,
          se = se,
          lower = lower,
          upper = upper,
          Rhat = rhat_val,
          Bulk_ESS = bulk_ess_val,
          Tail_ESS = tail_ess_val,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(rows) > 0) {
      do.call(rbind, rows)
    } else {
      data.frame(
        dim1 = character(),
        dim2 = character(),
        estimate = numeric(),
        se = numeric(),
        lower = numeric(),
        upper = numeric(),
        Rhat = numeric(),
        Bulk_ESS = numeric(),
        Tail_ESS = numeric(),
        stringsAsFactors = FALSE
      )
    }
  }

  if ("residual__" %in% names(vc) && "cor" %in% names(vc[["residual__"]])) {
    cor_array <- vc[["residual__"]][["cor"]]
    if (!is.null(cor_array) && length(dim(cor_array)) == 3) {
      n_resp <- length(resp_names)
      result$residual_cor <- cor_array_to_tibble(cor_array, resp_names, "residual__", random_summary, rescor_pars)

      estimate_mat <- cor_array[, "Estimate", ]
      if (n_resp == 2 && all(dim(estimate_mat) == c(2, 2))) {
        estimate_mat <- t(estimate_mat)
      }
      colnames(estimate_mat) <- resp_names
      rownames(estimate_mat) <- resp_names
      result$residual_cor_matrix <- estimate_mat
    }
  }

  for (grp in names(vc)) {
    if (grp == "residual__") next
    if ("cor" %in% names(vc[[grp]])) {
      cor_array <- vc[[grp]][["cor"]]
      if (!is.null(cor_array) && length(dim(cor_array)) == 3) {
        n_resp <- length(resp_names)
        result$random_effect_cor[[grp]] <- cor_array_to_tibble(cor_array, resp_names, grp, random_summary, rescor_pars)

        estimate_mat <- cor_array[, "Estimate", ]
        if (n_resp == 2 && all(dim(estimate_mat) == c(2, 2))) {
          estimate_mat <- t(estimate_mat)
        }
        colnames(estimate_mat) <- resp_names
        rownames(estimate_mat) <- resp_names
        result$random_effect_cor_matrix[[grp]] <- estimate_mat
      }
    }
  }

  if (length(result$random_effect_cor) == 0) {
    result$random_effect_cor <- NULL
  }
  if (length(result$random_effect_cor_matrix) == 0) {
    result$random_effect_cor_matrix <- NULL
  }

  result
}

#' Extract Covariances from brms Multivariate Model
#'
#' Extracts residual covariances and random effect covariances from a
#' multivariate brms model by computing them from the full posterior.
#'
#' @param model A brmsfit object from a multivariate model.
#'
#' @return A list with:
#'   \item{residual_cov}{Tibble with residual covariances (dim1, dim2, estimate, se, lower, upper, Rhat, Bulk_ESS, Tail_ESS)}
#'   \item{random_effect_cov}{Named list of tibbles for each facet with correlated random effects}
#'   \item{residual_cov_matrix}{Matrix of residual covariance point estimates}
#'   \item{random_effect_cov_matrix}{Named list of matrices for random effect covariances}
#'
#' @keywords internal
extract_covariances_brms <- function(model) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required.", call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required.", call. = FALSE)
  }

  draws <- brms::as_draws_matrix(model)
  model_summary <- suppressWarnings(summary(model))
  random_summary <- model_summary$random
  spec_pars <- model_summary$spec_pars
  rescor_pars <- model_summary$rescor_pars

  resp_names <- NULL
  if (!is.null(model$formula) && inherits(model$formula, "mvbrmsformula")) {
    resp_names <- model$formula$responses
  }

  if (is.null(resp_names) && !is.null(model$formula) && !is.null(model$formula$formula)) {
    bform <- model$formula$formula
    if (inherits(bform, "list") && length(bform) > 0) {
      resp_names <- names(bform)
    }
  }

  result <- list(
    residual_cov = NULL,
    random_effect_cov = list(),
    residual_cov_matrix = NULL,
    random_effect_cov_matrix = list()
  )

  if (is.null(resp_names) || length(resp_names) <= 1) {
    return(result)
  }

  n_resp <- length(resp_names)

  cov_array_to_tibble <- function(cov_draws_list, resp_names, grp_name = NULL, random_summary, rescor_pars) {
    n <- length(resp_names)
    rows <- list()

    for (i in 2:n) {
      for (j in 1:(i - 1)) {
        cov_draws <- cov_draws_list[[paste0(resp_names[j], "_", resp_names[i])]]
        if (is.null(cov_draws)) {
          cov_draws <- cov_draws_list[[paste0(resp_names[i], "_", resp_names[j])]]
        }

        if (is.null(cov_draws) || length(cov_draws) == 0) {
          next
        }

        estimate <- mean(cov_draws, na.rm = TRUE)
        se <- sd(cov_draws, na.rm = TRUE)
        lower <- quantile(cov_draws, 0.025, na.rm = TRUE)
        upper <- quantile(cov_draws, 0.975, na.rm = TRUE)

        rhat_val <- NA_real_
        bulk_ess_val <- NA_real_
        tail_ess_val <- NA_real_

        if (is.null(grp_name) || grp_name == "residual__") {
          if (!is.null(rescor_pars) && nrow(rescor_pars) > 0) {
            possible_names <- c(
              paste0("rescor(", resp_names[j], ",", resp_names[i], ")"),
              paste0("rescor(", resp_names[i], ",", resp_names[j], ")")
            )
            for (pname in possible_names) {
              if (pname %in% rownames(rescor_pars)) {
                rhat_val <- rescor_pars[pname, "Rhat"]
                bulk_ess_val <- rescor_pars[pname, "Bulk_ESS"]
                tail_ess_val <- rescor_pars[pname, "Tail_ESS"]
                break
              }
            }
          }
        } else {
          if (!is.null(grp_name) && grp_name %in% names(random_summary)) {
            grp_summary <- random_summary[[grp_name]]
            if (!is.null(grp_summary) && nrow(grp_summary) > 0) {
              possible_names <- c(
                paste0("cor(", resp_names[j], "_Intercept,", resp_names[i], "_Intercept)"),
                paste0("cor(", resp_names[i], "_Intercept,", resp_names[j], "_Intercept)")
              )
              for (pname in possible_names) {
                if (pname %in% rownames(grp_summary)) {
                  rhat_val <- grp_summary[pname, "Rhat"]
                  bulk_ess_val <- grp_summary[pname, "Bulk_ESS"]
                  tail_ess_val <- grp_summary[pname, "Tail_ESS"]
                  break
                }
              }
            }
          }
        }

        rows[[length(rows) + 1]] <- data.frame(
          dim1 = resp_names[j],
          dim2 = resp_names[i],
          estimate = estimate,
          se = se,
          lower = lower,
          upper = upper,
          Rhat = rhat_val,
          Bulk_ESS = bulk_ess_val,
          Tail_ESS = tail_ess_val,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(rows) > 0) {
      do.call(rbind, rows)
    } else {
      data.frame(
        dim1 = character(),
        dim2 = character(),
        estimate = numeric(),
        se = numeric(),
        lower = numeric(),
        upper = numeric(),
        Rhat = numeric(),
        Bulk_ESS = numeric(),
        Tail_ESS = numeric(),
        stringsAsFactors = FALSE
      )
    }
  }

  extract_sd_draws <- function(draws, param_base, resp) {
    possible_names <- c(
      paste0(param_base, "_", resp, "_Intercept"),
      paste0(param_base, "__", resp, "_Intercept"),
      paste0(param_base, "_", resp),
      paste0(param_base, "__", resp)
    )
    for (pname in possible_names) {
      if (pname %in% colnames(draws)) {
        return(posterior::extract_variable(draws, pname))
      }
    }
    NULL
  }

  extract_cor_draws <- function(draws, param_base, resp1, resp2) {
    possible_names <- c(
      paste0(param_base, "(", resp1, "_Intercept,", resp2, "_Intercept)"),
      paste0(param_base, "(", resp2, "_Intercept,", resp1, "_Intercept)")
    )
    for (pname in possible_names) {
      if (pname %in% colnames(draws)) {
        return(posterior::extract_variable(draws, pname))
      }
    }
    NULL
  }

  sigma_draws <- list()
  for (resp in resp_names) {
    sigma_draws[[resp]] <- posterior::extract_variable(draws, paste0("sigma_", resp))
  }

  if ("rescor" %in% colnames(draws) || any(grepl("^rescor[(]", colnames(draws))) || any(grepl("^rescor_", colnames(draws)))) {
    cov_draws_list <- list()
    cov_mat <- matrix(NA, n_resp, n_resp)
    rownames(cov_mat) <- resp_names
    colnames(cov_mat) <- resp_names

for (i in 2:n_resp) {
  for (j in 1:(i - 1)) {
    cor_draws <- NULL
    possible_names <- c(
      paste0("rescor(", resp_names[j], ",", resp_names[i], ")"),
      paste0("rescor(", resp_names[i], ",", resp_names[j], ")"),
      paste0("rescor__", resp_names[j], "__", resp_names[i]),
      paste0("rescor__", resp_names[i], "__", resp_names[j])
    )
    for (pname in possible_names) {
      if (pname %in% colnames(draws)) {
        cor_draws <- posterior::extract_variable(draws, pname)
        break
      }
    }

    if (!is.null(cor_draws) && length(cor_draws) > 0) {
      sd_i <- sigma_draws[[resp_names[i]]]
      sd_j <- sigma_draws[[resp_names[j]]]
      cov_draws <- cor_draws * sd_i * sd_j
      cov_draws_list[[paste0(resp_names[j], "_", resp_names[i])]] <- cov_draws
      cov_mat[resp_names[j], resp_names[i]] <- mean(cov_draws, na.rm = TRUE)
      cov_mat[resp_names[i], resp_names[j]] <- mean(cov_draws, na.rm = TRUE)
    }
  }
}

    for (i in seq_len(n_resp)) {
      var_draws <- sigma_draws[[resp_names[i]]]^2
      cov_mat[resp_names[i], resp_names[i]] <- mean(var_draws, na.rm = TRUE)
    }

    if (length(cov_draws_list) > 0) {
      result$residual_cov <- cov_array_to_tibble(cov_draws_list, resp_names, "residual__", random_summary, rescor_pars)
      result$residual_cov_matrix <- cov_mat
    }
  }

  vc <- brms::VarCorr(model)
  for (grp in names(vc)) {
    if (grp == "residual__") next

    sd_param_base <- paste0("sd_", grp)
    cor_param_base <- paste0("cor_", grp)

    has_cor <- FALSE
    for (i in 2:n_resp) {
      for (j in 1:(i - 1)) {
        possible_names <- c(
          paste0(cor_param_base, "(", resp_names[j], "_Intercept,", resp_names[i], "_Intercept)"),
          paste0(cor_param_base, "(", resp_names[i], "_Intercept,", resp_names[j], "_Intercept)")
        )
        for (pname in possible_names) {
          if (pname %in% colnames(draws)) {
            has_cor <- TRUE
            break
          }
        }
        if (has_cor) break
      }
      if (has_cor) break
    }

    if (!has_cor) next

    cov_draws_list <- list()
    cov_mat <- matrix(NA, n_resp, n_resp)
    rownames(cov_mat) <- resp_names
    colnames(cov_mat) <- resp_names

    sd_draws_list <- list()
    for (resp in resp_names) {
      sd_draws_list[[resp]] <- extract_sd_draws(draws, sd_param_base, resp)
    }

    for (i in 2:n_resp) {
      for (j in 1:(i - 1)) {
        cor_draws <- extract_cor_draws(draws, cor_param_base, resp_names[j], resp_names[i])

        if (!is.null(cor_draws) && length(cor_draws) > 0 &&
            !is.null(sd_draws_list[[resp_names[i]]]) &&
            !is.null(sd_draws_list[[resp_names[j]]])) {
          cov_draws <- cor_draws * sd_draws_list[[resp_names[i]]] * sd_draws_list[[resp_names[j]]]
          cov_draws_list[[paste0(resp_names[j], "_", resp_names[i])]] <- cov_draws
          cov_mat[resp_names[j], resp_names[i]] <- mean(cov_draws, na.rm = TRUE)
          cov_mat[resp_names[i], resp_names[j]] <- mean(cov_draws, na.rm = TRUE)
        }
      }
    }

    for (i in seq_len(n_resp)) {
      if (!is.null(sd_draws_list[[resp_names[i]]])) {
        var_draws <- sd_draws_list[[resp_names[i]]]^2
        cov_mat[resp_names[i], resp_names[i]] <- mean(var_draws, na.rm = TRUE)
      }
    }

    if (length(cov_draws_list) > 0) {
      result$random_effect_cov[[grp]] <- cov_array_to_tibble(cov_draws_list, resp_names, grp, random_summary, rescor_pars)
      result$random_effect_cov_matrix[[grp]] <- cov_mat
    }
  }

  if (length(result$random_effect_cov) == 0) {
    result$random_effect_cov <- NULL
  }
  if (length(result$random_effect_cov_matrix) == 0) {
    result$random_effect_cov_matrix <- NULL
  }

  result
}

#' Extract Covariance Draws from brms Model
#'
#' Extracts posterior draws of covariances between dimensions for
#' multivariate models. Used for computing composite variance components.
#'
#' @param model A fitted brms model object
#' @param dimensions Character vector of dimension/response names
#'
#' @return Named list with:
#' \describe{
#' \item{residual}{List of residual covariance draws (named by dimension pairs)}
#' \item{random_effect}{Named list of covariances per facet}
#' }
#'
#' @keywords internal
extract_covariance_draws <- function(model, dimensions) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required.", call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required.", call. = FALSE)
  }
  
  if (!inherits(model, "brmsfit")) {
    stop("model must be a brmsfit object", call. = FALSE)
  }

  draws <- brms::as_draws_matrix(model)
  n_dim <- length(dimensions)

  result <- list(
    residual = list(),
    random_effect = list()
  )

  if (n_dim <= 1) {
    return(result)
  }

  sigma_draws <- list()
  for (d in dimensions) {
    param_name <- paste0("sigma_", d)
    if (param_name %in% colnames(draws)) {
      sigma_draws[[d]] <- posterior::extract_variable(draws, param_name)
    }
  }

  for (i in 2:n_dim) {
    for (j in 1:(i - 1)) {
      param_name <- paste0("rescor(", dimensions[j], ",", dimensions[i], ")")
      alt_name <- paste0("rescor(", dimensions[i], ",", dimensions[j], ")")

      cor_draws <- NULL
      if (param_name %in% colnames(draws)) {
        cor_draws <- posterior::extract_variable(draws, param_name)
      } else if (alt_name %in% colnames(draws)) {
        cor_draws <- posterior::extract_variable(draws, alt_name)
      }

      if (!is.null(cor_draws) && !is.null(sigma_draws[[dimensions[i]]]) && !is.null(sigma_draws[[dimensions[j]]])) {
        cov_draws <- sigma_draws[[dimensions[j]]] * sigma_draws[[dimensions[i]]] * cor_draws
        result$residual[[paste0(dimensions[j], "_", dimensions[i])]] <- cov_draws
      }
    }
  }

  vc_groups <- names(brms::VarCorr(model))
  vc_groups <- setdiff(vc_groups, "residual__")

  for (grp in vc_groups) {
    clean_grp <- gsub("__", ":", grp)
    result$random_effect[[clean_grp]] <- list()

    for (i in 2:n_dim) {
      for (j in 1:(i - 1)) {
        param_name <- paste0("cor_", grp, "(", dimensions[j], "_Intercept,", dimensions[i], "_Intercept)")
        alt_name <- paste0("cor_", grp, "(", dimensions[i], "_Intercept,", dimensions[j], "_Intercept)")

        cor_draws <- NULL
        if (param_name %in% colnames(draws)) {
          cor_draws <- posterior::extract_variable(draws, param_name)
        } else if (alt_name %in% colnames(draws)) {
          cor_draws <- posterior::extract_variable(draws, alt_name)
        }

        if (!is.null(cor_draws)) {
          sd_param1 <- paste0("sd_", grp, "__", dimensions[i], "_Intercept")
          sd_param2 <- paste0("sd_", grp, "__", dimensions[j], "_Intercept")

          sd_draws1 <- NULL
          sd_draws2 <- NULL

          if (sd_param1 %in% colnames(draws)) {
            sd_draws1 <- posterior::extract_variable(draws, sd_param1)
          }
          if (sd_param2 %in% colnames(draws)) {
            sd_draws2 <- posterior::extract_variable(draws, sd_param2)
          }

          if (!is.null(sd_draws1) && !is.null(sd_draws2)) {
            cov_draws <- sd_draws2 * sd_draws1 * cor_draws
            result$random_effect[[clean_grp]][[paste0(dimensions[j], "_", dimensions[i])]] <- cov_draws
          }
        }
      }
    }

    if (length(result$random_effect[[clean_grp]]) == 0) {
      result$random_effect[[clean_grp]] <- NULL
    }
  }

  if (length(result$random_effect) == 0 || all(sapply(result$random_effect, is.null))) {
    result$random_effect <- NULL
  }

  result
}

#' Extract Covariances from MOM Multivariate Model
#'
#' Extracts residual covariances and random effect covariances from a
#' multivariate method-of-moments model.
#'
#' @param model A momfit object from a multivariate model.
#'
#' @return A list with:
#'   \item{residual_cov}{Tibble with residual covariances (dim1, dim2, estimate, se, lower, upper)}
#'   \item{random_effect_cov}{Named list of tibbles for each facet with correlated random effects}
#'   \item{residual_cov_matrix}{Matrix of residual covariance point estimates}
#'   \item{random_effect_cov_matrix}{Named list of matrices for random effect covariances}
#'
#' @keywords internal
extract_covariances_mom <- function(model) {
  if (!inherits(model, "momfit")) {
    stop("model must be a momfit object", call. = FALSE)
  }

  if (!isTRUE(model$is_multivariate) || length(model$responses) <= 1) {
    return(list(
      residual_cov = NULL,
      random_effect_cov = list(),
      residual_cov_matrix = NULL,
      random_effect_cov_matrix = list()
    ))
  }

  responses <- model$responses
  n_resp <- length(responses)

  result <- list(
    residual_cov = NULL,
    random_effect_cov = list(),
    residual_cov_matrix = NULL,
    random_effect_cov_matrix = list()
  )

  vc <- model$variance_components
  correlations <- model$correlations

  if (is.null(correlations)) {
    return(result)
  }

  residual_sd <- numeric(n_resp)
  names(residual_sd) <- responses
  for (i in seq_along(responses)) {
    resp <- responses[i]
    resid_vc <- vc[vc$dim == resp & vc$component == "Residual", ]
    if (nrow(resid_vc) > 0) {
      residual_sd[i] <- sqrt(max(0, resid_vc$var[1]))
    }
  }

  if (!is.null(correlations$residual_cor) && nrow(correlations$residual_cor) > 0) {
    cov_list <- list()
    cov_mat <- matrix(NA, n_resp, n_resp)
    rownames(cov_mat) <- responses
    colnames(cov_mat) <- responses

    for (i in seq_len(nrow(correlations$residual_cor))) {
      cor_row <- correlations$residual_cor[i, ]
      dim1 <- cor_row$dim1
      dim2 <- cor_row$dim2
      cor_est <- cor_row$estimate

      sd1 <- residual_sd[dim1]
      sd2 <- residual_sd[dim2]

      cov_est <- cor_est * sd1 * sd2

      se_cov <- NA_real_
      if ("se" %in% names(cor_row) && !is.na(cor_row$se)) {
        se_cov <- cor_row$se * sd1 * sd2
      }

      lower_cov <- NA_real_
      upper_cov <- NA_real_
      if ("lower" %in% names(cor_row) && !is.na(cor_row$lower)) {
        lower_cov <- cor_row$lower * sd1 * sd2
      }
      if ("upper" %in% names(cor_row) && !is.na(cor_row$upper)) {
        upper_cov <- cor_row$upper * sd1 * sd2
      }

      cov_list[[length(cov_list) + 1]] <- data.frame(
        dim1 = dim1,
        dim2 = dim2,
        estimate = cov_est,
        se = se_cov,
        lower = lower_cov,
        upper = upper_cov,
        stringsAsFactors = FALSE
      )

      cov_mat[dim1, dim2] <- cov_est
      cov_mat[dim2, dim1] <- cov_est
    }

    for (i in seq_len(n_resp)) {
      cov_mat[responses[i], responses[i]] <- residual_sd[i]^2
    }

    if (length(cov_list) > 0) {
      result$residual_cov <- do.call(rbind, cov_list)
      result$residual_cov_matrix <- cov_mat
    }
  }

  if (!is.null(correlations$random_effect_cor) && length(correlations$random_effect_cor) > 0) {
    for (facet in names(correlations$random_effect_cor)) {
      cor_tibble <- correlations$random_effect_cor[[facet]]
      if (is.null(cor_tibble) || nrow(cor_tibble) == 0) next

      facet_sd <- numeric(n_resp)
      names(facet_sd) <- responses
      for (i in seq_along(responses)) {
        resp <- responses[i]
        facet_vc <- vc[vc$dim == resp & vc$component == facet, ]
        if (nrow(facet_vc) > 0) {
          facet_sd[i] <- sqrt(max(0, facet_vc$var[1]))
        }
      }

      cov_list <- list()
      cov_mat <- matrix(NA, n_resp, n_resp)
      rownames(cov_mat) <- responses
      colnames(cov_mat) <- responses

      for (i in seq_len(nrow(cor_tibble))) {
        cor_row <- cor_tibble[i, ]
        dim1 <- cor_row$dim1
        dim2 <- cor_row$dim2
        cor_est <- cor_row$estimate

        sd1 <- facet_sd[dim1]
        sd2 <- facet_sd[dim2]

        cov_est <- cor_est * sd1 * sd2

        se_cov <- NA_real_
        if ("se" %in% names(cor_row) && !is.na(cor_row$se)) {
          se_cov <- cor_row$se * sd1 * sd2
        }

        lower_cov <- NA_real_
        upper_cov <- NA_real_
        if ("lower" %in% names(cor_row) && !is.na(cor_row$lower)) {
          lower_cov <- cor_row$lower * sd1 * sd2
        }
        if ("upper" %in% names(cor_row) && !is.na(cor_row$upper)) {
          upper_cov <- cor_row$upper * sd1 * sd2
        }

        cov_list[[length(cov_list) + 1]] <- data.frame(
          dim1 = dim1,
          dim2 = dim2,
          estimate = cov_est,
          se = se_cov,
          lower = lower_cov,
          upper = upper_cov,
          stringsAsFactors = FALSE
        )

        cov_mat[dim1, dim2] <- cov_est
        cov_mat[dim2, dim1] <- cov_est
      }

      for (i in seq_len(n_resp)) {
        cov_mat[responses[i], responses[i]] <- facet_sd[i]^2
      }

      if (length(cov_list) > 0) {
        result$random_effect_cov[[facet]] <- do.call(rbind, cov_list)
        result$random_effect_cov_matrix[[facet]] <- cov_mat
      }
    }
  }

  if (length(result$random_effect_cov) == 0) {
    result$random_effect_cov <- NULL
  }
  if (length(result$random_effect_cov_matrix) == 0) {
    result$random_effect_cov_matrix <- NULL
  }

  result
}

#' @rdname extract_vc
#' @keywords internal
extract_vc.brmsfit <- function(model, ...) {
  extract_vc_brms(model, ...)
}

#' Extract Variance Components from brms Long-Format Model
#'
#' Extracts variance components for long-format multivariate models
#' where dimensions are defined by a dimension variable.
#'
#' @param model A brmsfit object.
#' @param conf_level Credible interval level (default 0.95).
#' @param formula The original formula.
#' @param dimension_var Name of the dimension variable.
#' @param data The original data.
#'
#' @return A tibble with variance components per dimension.
#'
#' @keywords internal
extract_vc_brms_long_format <- function(model, conf_level = 0.95, formula = NULL, dimension_var = NULL, data = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required.", call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required.", call. = FALSE)
  }

  draws <- brms::as_draws_matrix(model)
  model_summary <- suppressWarnings(summary(model))
  random_summary <- model_summary$random
  spec_pars <- model_summary$spec_pars

  # Get dimension levels from data
  dimension_levels <- unique(data[[dimension_var]])
  dimension_levels <- sort(as.character(dimension_levels))

  # Detect if sigma uses log link
  is_log_link <- detect_sigma_log_link(model)

  # Get facet names from formula
  vc <- brms::VarCorr(model)
  facet_names <- names(vc)
  facet_names <- facet_names[facet_names != "residual__"]

  # Normalize facet names (replace __ with :)
  normalize_group_name <- function(name) {
    gsub("__", ":", name)
  }

  results <- list()

  for (dim in dimension_levels) {
    # Extract random effect variances for this dimension
    for (facet in facet_names) {
      normalized_facet <- normalize_group_name(facet)
      type <- if (grepl(":", normalized_facet)) "interaction" else "main"

      result <- extract_single_variance_from_draws_long_format(
        draws = draws,
        grp = facet,
        resp = dim,
        type = type,
        random_summary = random_summary,
        spec_pars = spec_pars,
        is_mv = TRUE,
        is_log_link = FALSE
      )

      if (!is.null(result)) {
        result$component <- normalized_facet
        results[[length(results) + 1]] <- result
      }
    }

    # Extract residual (sigma) variance for this dimension
    result <- extract_single_variance_from_draws_long_format(
      draws = draws,
      grp = "residual__",
      resp = dim,
      type = "residual",
      random_summary = random_summary,
      spec_pars = spec_pars,
      is_mv = TRUE,
      is_log_link = is_log_link
    )

    if (!is.null(result)) {
      results[[length(results) + 1]] <- result
    }
  }

  # Combine results
  vc_df <- do.call(rbind, results)
  vc_tibble <- tibble::as_tibble(vc_df)

  # Calculate percentages
  vc_tibble <- vc_tibble %>%
    dplyr::group_by(dim) %>%
    dplyr::mutate(pct = (var / sum(var, na.rm = TRUE)) * 100) %>%
    dplyr::ungroup()

  vc_tibble
}

#' Detect if Sigma Uses Log Link in brms Model
#'
#' Checks whether the sigma parameters in a brms model use a log link
#' (indicated by 'b_sigma_' prefix in parameter names).
#'
#' @param model A brmsfit object.
#' @return Logical; TRUE if sigma uses log link, FALSE otherwise.
#'
#' @keywords internal
detect_sigma_log_link <- function(model) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    return(FALSE)
  }

  draws <- brms::as_draws_matrix(model)
  param_names <- colnames(draws)

  has_b_sigma <- any(grepl("^b_sigma_", param_names))
  has_sigma <- any(grepl("^sigma_", param_names))

  has_b_sigma && !has_sigma
}

#' Extract Single Variance Component from Long-Format brms Model
#'
#' Helper function to extract a single variance component from posterior draws
#' for long-format multivariate models.
#'
#' @param draws Posterior draws matrix.
#' @param grp Group name (facet name or "residual__").
#' @param resp Response dimension name.
#' @param type Component type: "main", "interaction", or "residual".
#' @param random_summary Random effects summary from brms.
#' @param spec_pars Special parameters summary from brms.
#' @param is_mv Whether this is a multivariate model.
#' @param is_log_link Whether sigma uses log link.
#'
#' @return A data.frame with variance component estimates.
#'
#' @keywords internal
extract_single_variance_from_draws_long_format <- function(draws, grp, resp, type,
                                                            random_summary, spec_pars, is_mv,
                                                            is_log_link = FALSE) {
  if (grp == "residual__") {
    if (is_log_link) {
      param_name <- paste0("b_sigma_", resp)
    } else {
      param_name <- paste0("sigma_", resp)
    }
    component <- "Residual"
  } else {
    param_name <- paste0("sd_", grp, "__", resp, "_Intercept")
    component <- grp
  }

  if (!param_name %in% colnames(draws)) {
    param_candidates <- character()

    if (grp != "residual__") {
      grp_double_underscore <- gsub(":", "__", grp)
      param_candidates <- c(
        paste0("sd_", grp_double_underscore, "__", resp, "_Intercept"),
        paste0("sd_", grp, "__", resp, "_Intercept")
      )
      param_candidates <- unique(param_candidates)
    } else {
      if (is_log_link) {
        param_candidates <- c(
          paste0("b_sigma_", resp),
          paste0("sigma_", resp)
        )
      } else {
        param_candidates <- paste0("sigma_", resp)
      }
    }

    found <- FALSE
    for (cand in param_candidates) {
      if (cand %in% colnames(draws)) {
        param_name <- cand
        found <- TRUE
        break
      }
    }

    if (!found) {
      base_name <- if (grp == "residual__") {
        if (is_log_link) paste0("b_sigma_", resp) else paste0("sigma_", resp)
      } else {
        paste0("sd_", gsub(":", "__", grp), "__", resp)
      }
      matching_params <- grep(paste0("^", base_name), colnames(draws), value = TRUE)
      if (length(matching_params) > 0) {
        param_name <- matching_params[1]
      } else {
        return(NULL)
      }
    }
  }

  draws_values <- posterior::extract_variable(draws, param_name)

  if (grp == "residual__" && is_log_link) {
    sd_draws <- exp(draws_values)
  } else {
    sd_draws <- draws_values
  }

  var_draws <- sd_draws^2

  estimate <- mean(var_draws)
  lower <- as.numeric(quantile(var_draws, 0.025))
  upper <- as.numeric(quantile(var_draws, 0.975))
  se <- sd(var_draws)
  sd_mean <- mean(sd_draws)

  rhat_val <- NA_real_
  bulk_ess_val <- NA_real_
  tail_ess_val <- NA_real_

  if (grp == "residual__") {
    if (!is.null(spec_pars)) {
      sigma_name <- paste0("sigma_", resp)
      if (sigma_name %in% rownames(spec_pars)) {
        rhat_val <- spec_pars[sigma_name, "Rhat"]
        bulk_ess_val <- spec_pars[sigma_name, "Bulk_ESS"]
        tail_ess_val <- spec_pars[sigma_name, "Tail_ESS"]
      }
    }
  } else {
    if (!is.null(random_summary) && grp %in% names(random_summary)) {
      grp_rand <- random_summary[[grp]]
        if (param_name %in% rownames(grp_rand)) {
          rhat_val <- grp_rand[param_name, "Rhat"]
          bulk_ess_val <- grp_rand[param_name, "Bulk_ESS"]
          tail_ess_val <- grp_rand[param_name, "Tail_ESS"]
        }
      }
    }

    data.frame(
      component = component,
      dim = resp,
      type = type,
      var = estimate,
      error = se,
      lower = lower,
      upper = upper,
      sd = sd_mean,
      Rhat = rhat_val,
      Bulk_ESS = bulk_ess_val,
      Tail_ESS = tail_ess_val,
      stringsAsFactors = FALSE
    )
  }

  extract_correlations_brms_long_format <- function(model, dimension_var, data) {
    if (!requireNamespace("brms", quietly = TRUE)) {
      return(NULL)
    }

    draws <- tryCatch(
      brms::as_draws_matrix(model),
      error = function(e) NULL
    )

    if (is.null(draws)) {
      return(list(random_effect_cor = list()))
    }

    dimension_levels <- unique(data[[dimension_var]])
    dimension_levels <- sort(as.character(dimension_levels))

    cor_params <- grep("^cor_", colnames(draws), value = TRUE)

    if (length(cor_params) == 0) {
      return(list(random_effect_cor = list()))
    }

    cor_by_facet <- list()

    for (param in cor_params) {
      match_result <- regmatches(param, regexec("cor_([^\\[]+)\\[([^,]+),([^\\]]+)\\]", param))[[1]]

      if (length(match_result) == 4) {
        facet <- match_result[2]
        dim1 <- match_result[3]
        dim2 <- match_result[4]

        if (is.null(cor_by_facet[[facet]])) {
          cor_by_facet[[facet]] <- list()
        }

        cor_draws <- tryCatch(
          posterior::extract_variable(draws, param),
          error = function(e) NULL
        )

        if (!is.null(cor_draws)) {
          cor_by_facet[[facet]] <- c(
            cor_by_facet[[facet]],
            list(list(
              dim1 = dim1,
              dim2 = dim2,
              estimate = mean(cor_draws),
              se = sd(cor_draws),
              lower = as.numeric(quantile(cor_draws, 0.025)),
              upper = as.numeric(quantile(cor_draws, 0.975))
            ))
          )
        }
      }
    }

    random_effect_cor <- list()

    for (facet in names(cor_by_facet)) {
      cor_list <- cor_by_facet[[facet]]

      cor_tibble <- dplyr::bind_rows(lapply(cor_list, function(x) {
        tibble::tibble(
          dim1 = x$dim1,
          dim2 = x$dim2,
          estimate = x$estimate,
          se = x$se,
          lower = x$lower,
          upper = x$upper
        )
      }))

      random_effect_cor[[facet]] <- cor_tibble
    }

  list(
    random_effect_cor = random_effect_cor,
    residual_cor = NULL
  )
}

#' Extract Covariances from brms Long-Format Multivariate Model
#'
#' Extracts residual covariances and random effect covariances from a
#' long-format multivariate brms model. Converts correlations to covariances
#' using Cov = Cor * SD_i * SD_j.
#'
#' @param model A fitted brms model
#' @param dimension_var Name of the dimension variable in the data
#' @param data The original data used to fit the model
#' @param vc Optional variance components tibble (for extracting SDs)
#'
#' @return A list with:
#'   - residual_cov: Tibble with residual covariances (dim1, dim2, estimate, se, lower, upper)
#'   - random_effect_cov: Named list of tibbles for each facet with covariances
#'   - residual_cov_matrix: Matrix of residual covariance point estimates
#'   - random_effect_cov_matrix: Named list of matrices for random effect covariances
#'
#' @keywords internal
extract_covariances_brms_long_format <- function(model, dimension_var, data, vc = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    return(list(
      residual_cov = NULL,
      random_effect_cov = list(),
      residual_cov_matrix = NULL,
      random_effect_cov_matrix = list()
    ))
  }

  draws <- tryCatch(
    brms::as_draws_matrix(model),
    error = function(e) NULL
  )

  if (is.null(draws)) {
    return(list(
      residual_cov = NULL,
      random_effect_cov = list(),
      residual_cov_matrix = NULL,
      random_effect_cov_matrix = list()
    ))
  }

  dimension_levels <- unique(data[[dimension_var]])
  dimension_levels <- sort(as.character(dimension_levels))
  n_dim <- length(dimension_levels)

  result <- list(
    residual_cov = NULL,
    random_effect_cov = list(),
    residual_cov_matrix = NULL,
    random_effect_cov_matrix = list()
  )

  # Extract sigma (residual SD) parameters for each dimension
  sigma_draws_list <- list()
  for (dim in dimension_levels) {
    possible_names <- c(
      paste0("sigma_", dim),
      paste0("sigma_", dim, "_Intercept"),
      paste0("sigma_", gsub(" ", "_", dim))
    )
    for (pname in possible_names) {
      if (pname %in% colnames(draws)) {
        sigma_draws_list[[dim]] <- posterior::extract_variable(draws, pname)
        break
      }
    }
  }

  # Extract residual correlations (rescor) if present
  rescor_draws_list <- list()
  rescor_params <- grep("^rescor\\(", colnames(draws), value = TRUE)

  for (param in rescor_params) {
    match_result <- regmatches(param, regexec("rescor\\(([^,]+),([^)]+)\\)", param))[[1]]
    if (length(match_result) == 3) {
      dim1 <- match_result[1]
      dim2 <- match_result[2]
      key <- paste0(dim1, "_", dim2)
      rescor_draws_list[[key]] <- posterior::extract_variable(draws, param)
    }
  }

  # Compute residual covariances from sigma and rescor
  if (length(sigma_draws_list) > 0 && length(rescor_draws_list) > 0) {
    residual_cov_rows <- list()
    residual_cov_matrix <- matrix(0, n_dim, n_dim)
    rownames(residual_cov_matrix) <- dimension_levels
    colnames(residual_cov_matrix) <- dimension_levels

    # Diagonal: variance
    for (dim in dimension_levels) {
      if (!is.null(sigma_draws_list[[dim]])) {
        var_draws <- sigma_draws_list[[dim]]^2
        residual_cov_matrix[dim, dim] <- mean(var_draws, na.rm = TRUE)
      }
    }

    # Off-diagonal: covariance from correlation
    for (i in 2:n_dim) {
      for (j in 1:(i - 1)) {
        dim1 <- dimension_levels[j]
        dim2 <- dimension_levels[i]

        key1 <- paste0(dim1, "_", dim2)
        key2 <- paste0(dim2, "_", dim1)
        cor_draws <- rescor_draws_list[[key1]]
        if (is.null(cor_draws)) cor_draws <- rescor_draws_list[[key2]]

        if (!is.null(cor_draws) &&
            !is.null(sigma_draws_list[[dim1]]) &&
            !is.null(sigma_draws_list[[dim2]])) {
          cov_draws <- cor_draws * sigma_draws_list[[dim1]] * sigma_draws_list[[dim2]]

          estimate <- mean(cov_draws, na.rm = TRUE)
          se <- sd(cov_draws, na.rm = TRUE)
          lower <- as.numeric(quantile(cov_draws, 0.025, na.rm = TRUE))
          upper <- as.numeric(quantile(cov_draws, 0.975, na.rm = TRUE))

          residual_cov_rows[[length(residual_cov_rows) + 1]] <- tibble::tibble(
            dim1 = dim1,
            dim2 = dim2,
            estimate = estimate,
            se = se,
            lower = lower,
            upper = upper
          )

          residual_cov_matrix[dim1, dim2] <- estimate
          residual_cov_matrix[dim2, dim1] <- estimate
        }
      }
    }

    if (length(residual_cov_rows) > 0) {
      result$residual_cov <- dplyr::bind_rows(residual_cov_rows)
      result$residual_cov_matrix <- residual_cov_matrix
    }
  }

  # Extract random effect correlations and convert to covariances
  cor_params <- grep("^cor_", colnames(draws), value = TRUE)

  if (length(cor_params) > 0) {
    # Group correlation parameters by facet
    cor_by_facet <- list()
    sd_by_facet <- list()

    for (param in cor_params) {
      match_result <- regmatches(param, regexec("cor_([^\\[]+)\\[([^,]+),([^\\]]+)\\]", param))[[1]]

      if (length(match_result) == 4) {
        facet <- match_result[2]
        dim1 <- match_result[3]
        dim2 <- match_result[4]

        if (is.null(cor_by_facet[[facet]])) {
          cor_by_facet[[facet]] <- list()
        }

        cor_draws <- tryCatch(
          posterior::extract_variable(draws, param),
          error = function(e) NULL
        )

        if (!is.null(cor_draws)) {
          cor_by_facet[[facet]][[paste0(dim1, "_", dim2)]] <- list(
            dim1 = dim1,
            dim2 = dim2,
            cor_draws = cor_draws
          )
        }
      }
    }

    # Extract SD parameters for each facet and dimension
    sd_params <- grep("^sd_", colnames(draws), value = TRUE)

    for (param in sd_params) {
      match_result <- regmatches(param, regexec("sd_([^\\[]+)\\[([^\\]]+)\\]", param))[[1]]

      if (length(match_result) == 3) {
        facet <- match_result[2]
        dim_combo <- match_result[3]

        if (is.null(sd_by_facet[[facet]])) {
          sd_by_facet[[facet]] <- list()
        }

        sd_draws <- tryCatch(
          posterior::extract_variable(draws, param),
          error = function(e) NULL
        )

        if (!is.null(sd_draws)) {
          sd_by_facet[[facet]][[dim_combo]] <- sd_draws
        }
      }
    }

    # Convert correlations to covariances for each facet
    for (facet in names(cor_by_facet)) {
      cov_rows <- list()
      n_facet_dim <- length(unique(unlist(lapply(cor_by_facet[[facet]], function(x) c(x$dim1, x$dim2)))))

      # Get all dimensions for this facet
      facet_dims <- unique(unlist(lapply(cor_by_facet[[facet]], function(x) c(x$dim1, x$dim2))))

      # Build covariance matrix for this facet
      cov_matrix <- matrix(0, length(facet_dims), length(facet_dims))
      rownames(cov_matrix) <- facet_dims
      colnames(cov_matrix) <- facet_dims

      # Fill diagonal with variances
      for (dim in facet_dims) {
        sd_param_name <- paste0(dim, "_Intercept")
        if (is.null(sd_by_facet[[facet]][[sd_param_name]])) {
          sd_param_name <- dim
        }

        if (!is.null(sd_by_facet[[facet]][[sd_param_name]])) {
          sd_draws <- sd_by_facet[[facet]][[sd_param_name]]
          var_draws <- sd_draws^2
          cov_matrix[dim, dim] <- mean(var_draws, na.rm = TRUE)
        }
      }

      # Fill off-diagonal with covariances
      for (key in names(cor_by_facet[[facet]])) {
        cor_info <- cor_by_facet[[facet]][[key]]
        dim1 <- cor_info$dim1
        dim2 <- cor_info$dim2
        cor_draws <- cor_info$cor_draws

        # Get SDs for both dimensions
        sd1_name <- paste0(dim1, "_Intercept")
        if (is.null(sd_by_facet[[facet]][[sd1_name]])) sd1_name <- dim1
        sd2_name <- paste0(dim2, "_Intercept")
        if (is.null(sd_by_facet[[facet]][[sd2_name]])) sd2_name <- dim2

        sd1_draws <- sd_by_facet[[facet]][[sd1_name]]
        sd2_draws <- sd_by_facet[[facet]][[sd2_name]]

        if (!is.null(sd1_draws) && !is.null(sd2_draws)) {
          cov_draws <- cor_draws * sd1_draws * sd2_draws

          estimate <- mean(cov_draws, na.rm = TRUE)
          se <- sd(cov_draws, na.rm = TRUE)
          lower <- as.numeric(quantile(cov_draws, 0.025, na.rm = TRUE))
          upper <- as.numeric(quantile(cov_draws, 0.975, na.rm = TRUE))

          cov_rows[[length(cov_rows) + 1]] <- tibble::tibble(
            dim1 = dim1,
            dim2 = dim2,
            estimate = estimate,
            se = se,
            lower = lower,
            upper = upper
          )

          cov_matrix[dim1, dim2] <- estimate
          cov_matrix[dim2, dim1] <- estimate
        }
      }

      if (length(cov_rows) > 0) {
        result$random_effect_cov[[facet]] <- dplyr::bind_rows(cov_rows)
        result$random_effect_cov_matrix[[facet]] <- cov_matrix
      }
    }
  }

  result
}
