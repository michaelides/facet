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

  model_summary <- summary(model)
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
  model_summary <- summary(model)
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

#' @rdname extract_vc
#' @keywords internal
extract_vc.brmsfit <- function(model, ...) {
  extract_vc_brms(model, ...)
}
