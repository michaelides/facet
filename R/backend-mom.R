#' Method of Moments Backend for G-Studies
#'
#' Functions to fit G-study models using the method of moments (ANOVA-based
#' variance component estimation). This is an internal backend used by [gstudy()].
#'
#' @name mom-backend
#' @keywords internal
NULL

#' Method of Moments Fit Object
#'
#' An object representing a G-study model fitted using the method of moments.
#'
#' @format A list with class "momfit" containing:
#' \describe{
#'   \item{formula}{The model formula}
#'   \item{data}{The data frame}
#'   \item{aov_model}{The aov model object}
#'   \item{anova_results}{List of ANOVA results per stratum}
#'   \item{variance_components}{Tibble of estimated variance components}
#'   \item{response}{Name of the response variable}
#'   \item{random_facets}{Character vector of facet names}
#' }
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   score = rnorm(100),
#'   person = factor(rep(1:20, each = 5)),
#'   item = factor(rep(1:5, times = 20))
#' )
#' fit <- fit_mom(score ~ (1 | person) + (1 | item), data)
#' summary(fit)
#' }
#'
#' @name momfit
#' @rdname momfit
NULL

#' Fit a Model Using Method of Moments
#'
#' Fits a mixed effects model using the method of moments (ANOVA-based)
#' variance component estimation. Works best with balanced designs.
#'
#' @param formula A formula for the model (lme4 style), e.g., `y ~ (1 | facet1) + (1 | facet2)`.
#' @param data A data frame containing the variables in the formula.
#' @param ... Additional arguments (currently unused).
#' @return An object of class "momfit" containing:
#'   \item{formula}{The formula used}
#'   \item{data}{The data}
#'   \item{aov_model}{The aov model object}
#'   \item{anova_results}{The ANOVA results for each random effect}
#'   \item{variance_components}{The estimated variance components as a tibble}
#'   \item{response}{Name of the response variable}
#'   \item{random_facets}{Character vector of random effect facet names}
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   score = rnorm(100),
#'   person = factor(rep(1:20, each = 5)),
#'   item = factor(rep(1:5, times = 20))
#' )
#' fit_mom(score ~ (1 | person) + (1 | item), data)
#' }
#'
#' @param unbalanced Logical. If TRUE, enables unbalanced multivariate estimation
#' using Henderson's Method III. Each dimension is analyzed with its available data.
#' Default is FALSE.
#'
#' @export
#' @rdname fit_mom
fit_mom <- function(formula, data, unbalanced = FALSE, ...) {
  # Check if formula is multivariate
  if (is_multivariate(formula)) {
    return(fit_mom_multivariate(formula, data, unbalanced = unbalanced, ...))
  }

  # Parse the formula
  parsed <- parse_g_formula(formula)
  response <- parsed$response
  random_terms <- parsed$random
  random_facets <- parsed$random_facets
  random_facet_specs <- parsed$random_facet_specs

  # Check that we have random effects
  if (length(random_terms) == 0) {
    stop("Method of moments requires at least one random effect term.",
      call. = FALSE
    )
  }

  # Fit using the aovlist approach to get all variance components
  # Build formula with all random effects as Error strata
  aov_formula <- build_aov_formula(formula, data)

  # Fit the ANOVA model
  aov_model <- tryCatch(
    {
      aov(as.formula(aov_formula), data = data)
    },
    error = function(e) {
      stop("Error fitting model with method of moments: ", e$message,
        call. = FALSE
      )
    }
  )

  # Extract variance components from the aovlist object
  anova_results <- extract_anova_from_aovlist(aov_model, random_facets)

  # Compute variance components from ANOVA results
  vc <- compute_mom_vc_from_results(anova_results, data, response, random_facets, random_facet_specs)

  # Create the momfit object
  result <- structure(
    list(
      formula = formula,
      data = data,
      aov_model = aov_model,
      anova_results = anova_results,
      variance_components = vc,
      response = response,
      random_facets = random_facets,
      random_facet_specs = random_facet_specs
    ),
    class = "momfit"
  )

  result
}

#' Build aov Formula from lme4-style Formula
#'
#' @keywords internal
build_aov_formula <- function(formula, data = NULL) {
  formula_char <- paste(deparse(formula), collapse = " ")

  # Check for multivariate formula
  is_mv <- is_multivariate(formula)

  # Extract response - handle multivariate case
  if (is_mv) {
    # For multivariate, we just return the formula char as-is
    # The caller will replace the response
    resp_match <- regmatches(formula_char, regexpr("mvbind\\s*\\([^)]+\\)", formula_char))
    if (length(resp_match) > 0) {
      # Keep the mvbind part, will be replaced by caller
      response <- "RESPONSE_PLACEHOLDER"
    } else {
      response <- sub("~.+", "", formula_char)
    }
  } else {
    response <- sub("~.+", "", formula_char)
  }

  response <- sub("\\s+$", "", response)

  # Extract random effect terms
  # Pattern: (1|facet) or (1|cor_name|facet)
  # Must start with 1 followed by |
  random_pattern <- "\\(\\s*1\\s*\\|[^)]+\\)"
  random_terms <- regmatches(formula_char, gregexpr(random_pattern, formula_char))[[1]]

  if (length(random_terms) == 0) {
    return(formula_char)
  }

  # Extract facet names
  # Handle both (1|facet) and (1|cor_name|facet) patterns
  facets <- character()
  for (term in random_terms) {
    # Split by | and get the last part
    parts <- strsplit(term, "\\|")[[1]]
    if (length(parts) >= 2) {
      last_part <- trimws(parts[length(parts)])
      # Remove trailing ) if present
      last_part <- sub("\\)$", "", last_part)
      facets <- c(facets, last_part)
    }
  }

  # Build aov formula: y ~ facet1 + facet2 + ... + Error(strata)
  # For mixed nested/crossed designs, use additive Error terms
  # The crossing detection is handled in compute_n_effective, not here
  if (length(facets) == 1) {
    aov_formula <- paste(response, "~", facets[1], "+ Error(", facets[1], ")")
  } else {
    # Use all facets as additive Error terms to handle mixed designs
    error_strata <- paste(facets, collapse = " + ")
    aov_formula <- paste(
      response, "~", paste(facets, collapse = " + "),
      "+ Error(", error_strata, ")"
    )
  }

  aov_formula
}

#' Extract ANOVA Results from aovlist Object
#'
#' @keywords internal
extract_anova_from_aovlist <- function(aov_model, random_facets) {
  # aov with Error() returns an aovlist object
  # Each element corresponds to an error stratum

  results <- list()
  aov_summary <- summary(aov_model)

  # The summary is a list of summary.aovlist objects
  # Each contains a list of anova data.frames
  for (i in seq_along(aov_summary)) {
    # Each element is a list containing anova data.frames
    table_list <- aov_summary[[i]]

    # Extract each data.frame from the list
    for (j in seq_along(table_list)) {
      table_df <- table_list[[j]]

      if (is.data.frame(table_df)) {
        for (k in seq_len(nrow(table_df))) {
          results[[length(results) + 1]] <- list(
            source = trimws(rownames(table_df)[k]),
            df = as.numeric(table_df$Df[k]),
            ss = as.numeric(table_df$"Sum Sq"[k]),
            ms = as.numeric(table_df$"Mean Sq"[k]),
            stratum = paste0("stratum_", i)
          )
        }
      }
    }
  }

  results
}

#' Detect Crossing Patterns Between Facets in Data
#'
#' Analyzes the data to determine if pairs of facets are nested, crossed, or partially crossed.
#' - Nested: each level of A appears with exactly 1 level of B
#' - Crossed: each level of A appears with all levels of B
#' - Partial: some levels appear with 1, others with multiple
#'
#' @param data The data frame
#' @param facet_specs Character vector of facet specifications (e.g., c("Person", "Task", "Rater:Task"))
#' @return List containing:
#'   - pairwise: list of crossing patterns for each facet pair
#'   - facets: list with crossing info for each facet
#'   - pair_counts: number of unique observed pairs for each combination
#'
#' @keywords internal
detect_crossing_patterns_from_data <- function(data, facet_specs) {
  # Extract unique base facets
  base_facets <- unique(unlist(lapply(facet_specs, function(x) {
    strsplit(x, ":")[[1]]
  })))
  base_facets <- base_facets[base_facets %in% names(data)]

  result <- list(
    pairwise = list(),
    facets = list(),
    pair_counts = list()
  )

  # For each pair of facets, analyze crossing pattern
  for (i in seq_along(base_facets)) {
    for (j in seq_along(base_facets)) {
      if (i == j) next

      f_i <- base_facets[i]
      f_j <- base_facets[j]

      if (!(f_i %in% names(data) && f_j %in% names(data))) next

      # Get unique levels
      levels_i <- unique(data[[f_i]])
      levels_j <- unique(data[[f_j]])
      n_i <- length(levels_i)
      n_j <- length(levels_j)

      # Count how many levels of j each level of i appears with
      levels_j_per_i <- sapply(levels_i, function(lvl_i) {
        subset_data <- data[data[[f_i]] == lvl_i, ]
        length(unique(subset_data[[f_j]]))
      })

      # Count unique pairs
      pairs <- paste(data[[f_i]], data[[f_j]], sep = ":")
      n_unique_pairs <- length(unique(pairs))

      # Determine pattern
      min_levels <- min(levels_j_per_i)
      max_levels <- max(levels_j_per_i)

      if (min_levels == 1 && max_levels == 1) {
        pattern <- "nested"
      } else if (min_levels == n_j && max_levels == n_j) {
        pattern <- "crossed"
      } else {
        pattern <- "partial"
      }

      # Store pairwise info
      pair_name <- paste(sort(c(f_i, f_j)), collapse = ":")
      result$pairwise[[pair_name]] <- list(
        facet_a = f_i,
        facet_b = f_j,
        pattern = pattern,
        min_levels = min_levels,
        max_levels = max_levels,
        n_a = n_i,
        n_b = n_j,
        n_unique_pairs = n_unique_pairs
      )

      result$pair_counts[[pair_name]] <- n_unique_pairs
    }
  }

  # Store per-facet info
  for (f in base_facets) {
    result$facets[[f]] <- list(
      nested_in = character(),
      crossed_with = character(),
      partial_with = character()
    )
  }

  # Fill in facet-level info from pairwise
  for (pair_name in names(result$pairwise)) {
    info <- result$pairwise[[pair_name]]
    f_a <- info$facet_a
    f_b <- info$facet_b

    if (info$pattern == "nested") {
      result$facets[[f_a]]$nested_in <- c(result$facets[[f_a]]$nested_in, f_b)
    } else if (info$pattern == "crossed") {
      result$facets[[f_a]]$crossed_with <- c(result$facets[[f_a]]$crossed_with, f_b)
    } else {
      result$facets[[f_a]]$partial_with <- c(result$facets[[f_a]]$partial_with, f_b)
    }
  }

  result
}

#' Compute Effective Sample Size for a Variance Component
#'
#' Calculates the effective n for a given effect, accounting for nested vs crossed designs.
#'
#' @keywords internal
compute_n_effective <- function(source, n_total, facet_n, random_facets, crossing_info) {
  effect_facets <- unlist(strsplit(source, ":"))

  # For 3-way or higher interactions, treat as crossed
  if (length(effect_facets) >= 3) {
    n_effective <- n_total
    for (f in effect_facets) {
      if (f %in% names(facet_n)) {
        n_effective <- n_effective / facet_n[f]
      }
    }
    return(n_effective)
  }

  if (length(effect_facets) == 1 && effect_facets[1] %in% names(facet_n)) {
    # Main effect: n = total / (levels of this facet)
    main_facet <- effect_facets[1]
    n_effective <- n_total / facet_n[main_facet]
    return(n_effective)
  } else if (length(effect_facets) == 2) {
    # 2-way interaction: need to determine if nested, crossed, or partial
    f_a <- effect_facets[1]
    f_b <- effect_facets[2]

    # Look up crossing pattern
    pair_name <- paste(sort(c(f_a, f_b)), collapse = ":")

    if (!is.null(crossing_info$pairwise[[pair_name]])) {
      pattern_info <- crossing_info$pairwise[[pair_name]]
      pattern <- pattern_info$pattern
      n_unique_pairs <- pattern_info$n_unique_pairs

      if (pattern == "nested") {
        # For nested: n_eff = n_unique_pairs (number of "container" levels)
        # Actually for nested, the n_eff should be the number of observations
        # at the "between" level, which is n_total / n_unique_pairs
        n_effective <- n_total / n_unique_pairs
        return(n_effective)
      } else if (pattern == "crossed") {
        # For crossed: n_eff = n_total / (n_A * n_B)
        n_effective <- n_total / (facet_n[f_a] * facet_n[f_b])
        return(n_effective)
      } else if (pattern == "partial") {
        # For partially crossed: n_eff = n_total / n_unique_pairs
        n_effective <- n_total / n_unique_pairs
        return(n_effective)
      }
    }

    # Fallback: treat as crossed
    n_effective <- n_total
    for (f in effect_facets) {
      if (f %in% names(facet_n)) {
        n_effective <- n_effective / facet_n[f]
      }
    }
    return(n_effective)
  }

  1
}

#' Compute Variance Components from ANOVA Results
#'
#' @keywords internal
compute_mom_vc_from_results <- function(anova_results, data, response, random_facets, random_facet_specs = NULL) {
  if (length(anova_results) == 0) {
    stop("No ANOVA results to extract variance components from", call. = FALSE)
  }

  results <- list()

  # Get total observations
  n_total <- nrow(data)

  # If random_facet_specs is NULL, use random_facets
  if (is.null(random_facet_specs)) {
    random_facet_specs <- random_facets
  }

  # Get sample sizes for each facet
  # Use facet_specs to get the actual names as user specified
  all_facet_names <- unique(c(random_facets, unlist(strsplit(random_facet_specs, ":"))))
  facet_n <- sapply(all_facet_names, function(f) {
    if (f %in% names(data)) {
      length(unique(data[[f]]))
    } else {
      NA
    }
  })
  names(facet_n) <- all_facet_names

  # Build mapping from ANOVA source names to user-specified names
  # ANOVA may alphabetize interaction terms (e.g., "Task:Rater" instead of "Rater:Task")
  facet_name_mapping <- build_facet_name_mapping(random_facet_specs, anova_results)

  # Detect crossing patterns from the data
  crossing_info <- detect_crossing_patterns_from_data(data, random_facet_specs)

  # Build residual MS/DF mapping using stratum information
  resid_map <- build_source_resid_map(anova_results)

  # Process each result
  for (i in seq_along(anova_results)) {
    item <- anova_results[[i]]
    source <- item$source
    df <- item$df
    ms <- item$ms

    # Skip Intercept and residual within the same stratum
    if (source == "Intercept" || source == "Residuals") {
      next
    }

    # Determine type
    type <- if (grepl(":", source)) "interaction" else "main"

    # Get the user-specified name for this source
    user_source <- if (!is.null(facet_name_mapping[[source]])) {
      facet_name_mapping[[source]]
    } else {
      source
    }

    # Find the effective n for this effect
    n_effective <- compute_n_effective(source, n_total, facet_n, random_facets, crossing_info)

    # Method of moments estimate: var = (MS_effect - MS_residual) / n_effective
    # Use stratum-aware residual MS and df from the mapping
    resid_ms <- if (!is.null(resid_map$ms[[source]])) resid_map$ms[[source]] else 0
    resid_df <- if (!is.null(resid_map$df[[source]])) resid_map$df[[source]] else df

    var_est <- max(0, (ms - resid_ms) / n_effective)

    # Satterthwaite approximate degrees of freedom for the variance component
    if (var_est > 0 && ms > resid_ms) {
      num <- (ms - resid_ms)^2
      denom <- (ms^2 / df) + (resid_ms^2 / max(resid_df, 1))
      approx_df <- max(num / denom, 1)
    } else {
      approx_df <- df
    }

    results[[length(results) + 1]] <- data.frame(
      component = user_source,
      facet = user_source,
      type = type,
      var = var_est,
      df = df,
      ms = ms,
      resid_ms = resid_ms,
      n_eff = n_effective,
      approx_df = approx_df,
      stringsAsFactors = FALSE
    )
  }

  # Add residual variance
  resid_result <- find_residual_item(anova_results)
  if (!is.null(resid_result)) {
    results[[length(results) + 1]] <- data.frame(
      component = "Residual",
      facet = "Residual",
      type = "residual",
      var = resid_result$ms,
      df = resid_result$df,
      ms = resid_result$ms,
      resid_ms = resid_result$ms,
      n_eff = 1,
      approx_df = resid_result$df,
      stringsAsFactors = FALSE
    )
  }

  vc_df <- do.call(rbind, results)
  vc_tibble <- tibble::as_tibble(vc_df)

  # Calculate percentage of total variance
  total_var <- sum(vc_tibble$var, na.rm = TRUE)
  vc_tibble$pct <- (vc_tibble$var / total_var) * 100

  vc_tibble
}

#' Build Mapping from ANOVA Source Names to User-Specified Names
#'
#' @keywords internal
build_facet_name_mapping <- function(random_facet_specs, anova_results) {
  mapping <- list()

  # Get all source names from ANOVA results
  anova_sources <- sapply(anova_results, function(x) x$source)

  # For each user-specified facet spec, try to match to ANOVA source
  for (spec in random_facet_specs) {
    # Check if spec contains ":" (interaction)
    if (grepl(":", spec)) {
      # Create all possible orderings
      facets <- strsplit(spec, ":")[[1]]
      facets <- trimws(facets)

      if (length(facets) == 2) {
        # Try both orderings
        possible_names <- c(
          paste(facets[1], facets[2], sep = ":"),
          paste(facets[2], facets[1], sep = ":")
        )

        # Find which one exists in ANOVA results
        for (name in possible_names) {
          if (name %in% anova_sources) {
            mapping[[name]] <- spec
            break
          }
        }
      }
    }
  }

  mapping
}

#' Build Residual MS/DF Mapping for Each Source
#'
#' Builds a mapping from each ANOVA source to the appropriate residual MS and df,
#' using stratum information to select the correct residual when multiple strata
#' contain "Residuals" entries (common in nested and crossed designs with Error() terms).
#'
#' @keywords internal
build_source_resid_map <- function(anova_results) {
  sources <- sapply(anova_results, function(x) x$source)
  strata <- sapply(anova_results, function(x) x$stratum)
  unique_strata <- unique(strata[!is.na(strata)])

  # Get all non-Intercept, non-Residuals sources
  effect_idx <- which(sources != "Intercept" & sources != "Residuals")
  effect_sources <- sources[effect_idx]

  # Collect Residuals entries with their strata
  resid_idx <- which(sources == "Residuals")
  resid_strata <- strata[resid_idx]
  resid_ms_vals <- sapply(resid_idx, function(i) anova_results[[i]]$ms)
  resid_df_vals <- sapply(resid_idx, function(i) anova_results[[i]]$df)

  # Build a list of interaction sources (non-main effects) with their MS and strata
  interaction_sources <- effect_sources[grepl(":", effect_sources)]

  resid_map <- list()
  resid_df_map <- list()

  for (src in effect_sources) {
    src_stratum <- strata[which(sources == src)[1]]
    src_facets <- unlist(strsplit(src, ":"))

    # Step 1: look for a higher-order interaction containing this source
    best_match <- NULL
    best_len <- Inf

    for (other in interaction_sources) {
      if (other == src) next
      other_facets <- unlist(strsplit(other, ":"))
      if (all(src_facets %in% other_facets) && length(other_facets) > length(src_facets)) {
        if (length(other_facets) < best_len) {
          best_match <- other
          best_len <- length(other_facets)
        }
      }
    }

    if (!is.null(best_match)) {
      for (item in anova_results) {
        if (item$source == best_match) {
          resid_map[[src]] <- item$ms
          resid_df_map[[src]] <- item$df
          break
        }
      }
    }

    if (!is.null(resid_map[[src]])) next

    # Step 2: no interaction found; pick Residuals from the appropriate stratum.
    # The correct residual is from a stratum "below" this source's stratum.
    # If only one Residuals exists, use it. If multiple exist, prefer the one
    # from the last (innermost) stratum, which represents the within-cell error.
    if (length(resid_idx) == 0) {
      resid_map[[src]] <- 0
      resid_df_map[[src]] <- 1
    } else if (length(resid_idx) == 1) {
      resid_map[[src]] <- resid_ms_vals[1]
      resid_df_map[[src]] <- resid_df_vals[1]
    } else {
      # Multiple Residuals from different strata (nested/crossed design).
      # Use the last Residuals (innermost stratum = within-cell error).
      resid_map[[src]] <- resid_ms_vals[length(resid_ms_vals)]
      resid_df_map[[src]] <- resid_df_vals[length(resid_df_vals)]
    }
  }

  list(ms = resid_map, df = resid_df_map)
}

#' Find Residual MS for a Given Source
#'
#' For nested designs, the residual for a given effect is the next interaction
#' that contains all the factors in the effect plus one more factor.
#'
#' @keywords internal
find_residual_ms <- function(anova_results, source) {
  map <- build_source_resid_map(anova_results)
  if (!is.null(map$ms[[source]])) {
    return(map$ms[[source]])
  }
  0
}

#' Find Residual Item
#'
#' Returns the Residuals entry from the innermost (last) stratum.
#' This is the within-cell error, which is the correct residual variance
#' for G-study estimates.
#'
#' @keywords internal
find_residual_item <- function(anova_results) {
  # In designs with Error() terms, there can be multiple "Residuals" entries
  # from different strata. The last one is from the innermost (within-cell) stratum.
  last_resid <- NULL
  for (i in seq_along(anova_results)) {
    item <- anova_results[[i]]
    if (item$source == "Residuals") {
      last_resid <- item
    }
  }
  last_resid
}

#' Fit Multivariate Model Using Method of Moments
#'
#' @param formula A formula for the model
#' @param data A data frame containing the variables
#' @param unbalanced Logical. If TRUE, uses Henderson's Method III for unbalanced designs
#' @param ... Additional arguments
#' @keywords internal
fit_mom_multivariate <- function(formula, data, unbalanced = FALSE, ...) {
  if (unbalanced) {
    return(fit_mom_multivariate_unbalanced(formula, data, ...))
  }

  # Parse the formula to get responses
  formula_char <- deparse(formula)

  # Extract response variables from mvbind
  resp_match <- regmatches(formula_char, regexpr("mvbind\\s*\\([^)]+\\)", formula_char))
  if (length(resp_match) == 0 || resp_match == "") {
    stop("Could not parse multivariate formula for method of moments",
      call. = FALSE
    )
  }

  # Extract variable names from mvbind(...)
  resp_str <- sub("mvbind\\s*\\(", "", resp_match)
  resp_str <- sub("\\)\\s*$", "", resp_str)
  responses <- unlist(strsplit(resp_str, "\\s*,\\s*"))
  responses <- trimws(responses)

  # Parse random facets from formula
  parsed <- parse_g_formula(formula)
  random_facets <- parsed$random_facets
  random_facet_specs <- parsed$random_facet_specs

  # Extract set_rescor setting from formula (default FALSE)
  rescor_setting <- extract_rescor_setting(formula)

  # Extract random effect correlation specifications
  # Pattern: (1|cor_name|facet) means correlate random effects for 'facet'
  random_effect_cors <- extract_random_effect_cor(formula)

  # Check that we have random effects
  if (length(random_facets) == 0) {
    stop("Method of moments requires at least one random effect term.",
      call. = FALSE
    )
  }

  # Build aov formula using the same approach as univariate
  aov_formula <- build_aov_formula(formula, data)

  # We'll fit each response separately using aov
  # and then combine the results
  all_vc <- list()
  aov_models <- list()
  # Store formulas for each response to use when re-fitting
  # (formula() doesn't work on aovlist objects with Error terms)
  formulas_by_response <- list()

  for (resp in responses) {
    # Build formula for this response: y ~ facets + Error(facets)
    # The aov_formula already has the response, we need to swap it
    rhs_formula <- sub("^[^~]+~", "", aov_formula)
    rhs_formula <- sub("^\\s+", "", rhs_formula) # Remove leading whitespace
    resp_aov_formula <- paste(resp, "~", rhs_formula)

    # Fit ANOVA model for this response
    aov_model <- tryCatch(
      {
        aov(as.formula(resp_aov_formula), data = data)
      },
      error = function(e) {
        stop("Error fitting model with method of moments for response ", resp, ": ",
          e$message,
          call. = FALSE
        )
      }
    )

    aov_models[[resp]] <- aov_model
    # Store the formula string for this model to use later when re-fitting
    # (formula() doesn't work on aovlist objects with Error terms)
    formulas_by_response[[resp]] <- resp_aov_formula

    # Extract ANOVA results
    anova_results <- extract_anova_from_aovlist(aov_model, random_facets)

    # Compute variance components
    vc <- compute_mom_vc_from_results(
      anova_results, data, resp, random_facets, random_facet_specs
    )

    # Add dimension column
    vc$dim <- resp
    all_vc[[resp]] <- vc
  }

  # Combine variance components
  combined_vc <- do.call(rbind, all_vc)

  # Compute correlations
  # - If rescor_setting is TRUE: compute residual correlations
  # - If random_effect_cors has entries: compute random effect correlations for those facets
  # These are independent - can have both, either, or neither
  correlations <- NULL

  # Always compute correlations (let the function handle what to compute based on flags)
  # But pass flags to control what gets computed
  compute_rescor <- rescor_setting # residual correlations
  compute_re_cors <- length(random_effect_cors) > 0 # random effect correlations

  if (compute_rescor || compute_re_cors) {
    correlations <- compute_mom_multivariate_correlations(
      data, responses, random_facets, aov_models, formulas_by_response, random_effect_cors,
      compute_rescor = compute_rescor,
      compute_re_cors = compute_re_cors
    )
  }

  # Store rescor setting in result
  rescor_used <- rescor_setting

  # Create result object
  result <- structure(
    list(
      formula = formula,
      data = data,
      aov_models = aov_models,
      variance_components = combined_vc,
      responses = responses,
      random_facets = random_facets,
      random_facet_specs = random_facet_specs,
      correlations = correlations,
      rescor_used = rescor_used,
      random_effect_cors = random_effect_cors,
      is_multivariate = TRUE
    ),
    class = "momfit"
  )

  result
}

#' Fit Unbalanced Multivariate Model Using Method of Moments (Henderson's Method III)
#'
#' Implements Henderson's Method III for variance component estimation in
#' unbalanced multivariate designs. Each dimension is analyzed separately
#' using its available data, and correlations are computed using pairwise
#' complete cases.
#'
#' @param formula A formula for the model
#' @param data A data frame containing the variables
#' @param ... Additional arguments
#' @return A momfit object with unbalanced estimation
#' @keywords internal
fit_mom_multivariate_unbalanced <- function(formula, data, ...) {
  # Parse the formula to get responses
  formula_char <- deparse(formula)

  # Extract response variables from mvbind
  resp_match <- regmatches(formula_char, regexpr("mvbind\\s*\\([^)]+\\)", formula_char))
  if (length(resp_match) == 0 || resp_match == "") {
    stop("Could not parse multivariate formula for method of moments", call. = FALSE)
  }

  # Extract variable names from mvbind(...)
  resp_str <- sub("mvbind\\s*\\(", "", resp_match)
  resp_str <- sub("\\)\\s*$", "", resp_str)
  responses <- unlist(strsplit(resp_str, "\\s*,\\s*"))
  responses <- trimws(responses)

  # Parse random facets from formula
  parsed <- parse_g_formula(formula)
  random_facets <- parsed$random_facets
  random_facet_specs <- parsed$random_facet_specs

  # Extract set_rescor setting from formula (default FALSE)
  rescor_setting <- extract_rescor_setting(formula)

  # Extract random effect correlation specifications
  random_effect_cors <- extract_random_effect_cor(formula)

  # Check that we have random effects
  if (length(random_facets) == 0) {
    stop("Method of moments requires at least one random effect term.", call. = FALSE)
  }

  # Build aov formula template
  aov_formula_template <- build_aov_formula(formula, data)
  rhs_formula <- sub("^[^~]+~", "", aov_formula_template)
  rhs_formula <- sub("^\\s+", "", rhs_formula)

  # Track observations per dimension
  n_per_dim <- list()

  # Fit each response separately using its available data (Henderson's Method III)
  all_vc <- list()
  aov_models <- list()
  data_per_dim <- list()
  formulas_by_response <- list()

  for (resp in responses) {
    # Subset to non-missing for THIS response only
    dim_data <- data[!is.na(data[[resp]]), , drop = FALSE]
    n_per_dim[[resp]] <- nrow(dim_data)
    data_per_dim[[resp]] <- dim_data

    # Build formula for this response
    resp_aov_formula <- paste(resp, "~", rhs_formula)

    # Fit ANOVA model for this response using its available data
    aov_model <- tryCatch(
      {
        aov(as.formula(resp_aov_formula), data = dim_data)
      },
      error = function(e) {
        stop("Error fitting model with method of moments for response ", resp, ": ",
          e$message, call. = FALSE)
      }
    )

    aov_models[[resp]] <- aov_model
    formulas_by_response[[resp]] <- resp_aov_formula

    # Extract ANOVA results
    anova_results <- extract_anova_from_aovlist(aov_model, random_facets)

    # Compute variance components using this dimension's data
    vc <- compute_mom_vc_from_results(
      anova_results, dim_data, resp, random_facets, random_facet_specs
    )

    # Add dimension column
    vc$dim <- resp
    all_vc[[resp]] <- vc
  }

  # Combine variance components
  combined_vc <- do.call(rbind, all_vc)

  # Compute correlations using pairwise complete cases
  correlations <- NULL
  compute_rescor <- rescor_setting
  compute_re_cors <- length(random_effect_cors) > 0

  if (compute_rescor || compute_re_cors) {
    correlations <- compute_mom_correlations_unbalanced(
      data, responses, random_facets, aov_models, data_per_dim,
      formulas_by_response, random_effect_cors, compute_rescor, compute_re_cors
    )
  }

  # Create result object
  result <- structure(
    list(
      formula = formula,
      data = data,
      aov_models = aov_models,
      variance_components = combined_vc,
      responses = responses,
      random_facets = random_facets,
      random_facet_specs = random_facet_specs,
      correlations = correlations,
      rescor_used = rescor_setting,
      random_effect_cors = random_effect_cors,
      n_per_dim = n_per_dim,
      is_unbalanced = TRUE,
      is_multivariate = TRUE
    ),
    class = "momfit"
  )

  result
}

#' Compute Correlations for Unbalanced Multivariate Mom Model
#'
#' Computes residual and random effect correlations using pairwise complete cases
#' for unbalanced multivariate designs.
#'
#' @param data Original data frame
#' @param responses Character vector of response variable names
#' @param random_facets Character vector of random effect facet names
#' @param aov_models List of fitted aov models (one per response variable)
#' @param data_per_dim List of data frames used for each dimension
#' @param aov_formulas List of aov formula strings per response (with Error terms)
#' @param random_effect_cors List of random effect correlation specifications
#' @param compute_rescor Whether to compute residual correlations
#' @param compute_re_cors Whether to compute random effect correlations
#' @return List with residual and random effect correlations
#' @keywords internal
compute_mom_correlations_unbalanced <- function(data, responses, random_facets,
                                                aov_models, data_per_dim,
                                                aov_formulas = list(),
                                                random_effect_cors = list(),
                                                compute_rescor = FALSE,
                                                compute_re_cors = FALSE) {
  if (length(responses) < 2) {
    return(NULL)
  }

  result <- list(
    residual_cor = NULL,
    random_effect_cor = list(),
    residual_cor_matrix = NULL,
    random_effect_cor_matrix = list(),
    n_pairwise = list()
  )

  # Compute residual correlations using pairwise complete cases
  if (compute_rescor) {
    cor_list <- list()
    cor_matrix <- matrix(NA, length(responses), length(responses))
    rownames(cor_matrix) <- responses
    colnames(cor_matrix) <- responses
    diag(cor_matrix) <- 1

    for (i in seq_along(responses)) {
      for (j in seq_along(responses)) {
        if (i < j) {
          resp_i <- responses[i]
          resp_j <- responses[j]

          # Use pairwise complete cases for this pair of responses
          pairwise_idx <- !is.na(data[[resp_i]]) & !is.na(data[[resp_j]])
          pairwise_data <- data[pairwise_idx, , drop = FALSE]
          n_pairwise <- nrow(pairwise_data)

          if (n_pairwise < 3) {
            warning("Not enough pairwise complete cases for correlation between ",
              resp_i, " and ", resp_j, " (n = ", n_pairwise, ")", call. = FALSE)
            next
          }

          # Get residuals from models (need to refit on pairwise data if necessary)
          resp_formula_i <- if (!is.null(aov_formulas[[resp_i]])) aov_formulas[[resp_i]] else NULL
          resp_formula_j <- if (!is.null(aov_formulas[[resp_j]])) aov_formulas[[resp_j]] else NULL
          resids_i <- get_aov_residuals_for_pair(aov_models[[resp_i]], pairwise_data, resp_i, data_per_dim[[resp_i]], resp_formula_i)
          resids_j <- get_aov_residuals_for_pair(aov_models[[resp_j]], pairwise_data, resp_j, data_per_dim[[resp_j]], resp_formula_j)

          # Compute correlation
          cor_ij <- cor(resids_i, resids_j, use = "complete.obs")

          if (!is.na(cor_ij)) {
            # Standard error and CI using Fisher's z transformation
            se_ij <- sqrt((1 - cor_ij^2) / (n_pairwise - 2))
            z_ij <- atanh(cor_ij)
            se_z <- 1 / sqrt(n_pairwise - 3)
            z_lower <- z_ij - 1.96 * se_z
            z_upper <- z_ij + 1.96 * se_z
            lower <- tanh(z_lower)
            upper <- tanh(z_upper)

            cor_list[[length(cor_list) + 1]] <- data.frame(
              dim1 = resp_i,
              dim2 = resp_j,
              estimate = cor_ij,
              se = se_ij,
              lower = lower,
              upper = upper,
              n_pairwise = n_pairwise,
              stringsAsFactors = FALSE
            )

            cor_matrix[i, j] <- cor_ij
            cor_matrix[j, i] <- cor_ij
          }
        }
      }
    }

    if (length(cor_list) > 0) {
      result$residual_cor <- do.call(rbind, cor_list)
      result$residual_cor_matrix <- cor_matrix
    }
  }

  # Compute random effect correlations using pairwise complete cases per facet
  if (compute_re_cors) {
    for (facet in names(random_effect_cors)) {
      if (!facet %in% random_facets) next
      if (!facet %in% names(data)) next

      facet_cor_list <- list()

      for (i in seq_along(responses)) {
        for (j in seq_along(responses)) {
          if (i < j) {
            resp_i <- responses[i]
            resp_j <- responses[j]

            # Use pairwise complete cases
            pairwise_idx <- !is.na(data[[resp_i]]) & !is.na(data[[resp_j]])
            pairwise_data <- data[pairwise_idx, , drop = FALSE]

            if (!facet %in% names(pairwise_data)) next

            # Compute facet means for each response using pairwise data
            facet_means_i <- aggregate(
              pairwise_data[[resp_i]],
              by = list(facet = pairwise_data[[facet]]),
              FUN = mean,
              na.rm = TRUE
            )
            facet_means_j <- aggregate(
              pairwise_data[[resp_j]],
              by = list(facet = pairwise_data[[facet]]),
              FUN = mean,
              na.rm = TRUE
            )

            # Match facet levels
            common_facets <- intersect(facet_means_i$facet, facet_means_j$facet)
            if (length(common_facets) < 3) next

            y_i <- facet_means_i$x[facet_means_i$facet %in% common_facets]
            y_j <- facet_means_j$x[facet_means_j$facet %in% common_facets]

            n_facet <- length(common_facets)
            cor_ij <- cor(y_i, y_j)

            if (!is.na(cor_ij)) {
              se_ij <- sqrt((1 - cor_ij^2) / (n_facet - 2))
              z_ij <- atanh(cor_ij)
              se_z <- 1 / sqrt(n_facet - 3)
              z_lower <- z_ij - 1.96 * se_z
              z_upper <- z_ij + 1.96 * se_z
              lower <- tanh(z_lower)
              upper <- tanh(z_upper)

              facet_cor_list[[length(facet_cor_list) + 1]] <- data.frame(
                dim1 = resp_i,
                dim2 = resp_j,
                estimate = cor_ij,
                se = se_ij,
                lower = lower,
                upper = upper,
                n_facet = n_facet,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }

      if (length(facet_cor_list) > 0) {
        result$random_effect_cor[[facet]] <- do.call(rbind, facet_cor_list)

        # Build matrix
        n_resp <- length(responses)
        cor_matrix <- matrix(NA, n_resp, n_resp)
        rownames(cor_matrix) <- responses
        colnames(cor_matrix) <- responses
        for (k in seq_len(nrow(result$random_effect_cor[[facet]]))) {
          r <- result$random_effect_cor[[facet]][k, ]
          cor_matrix[r$dim1, r$dim2] <- r$estimate
          cor_matrix[r$dim2, r$dim1] <- r$estimate
        }
        diag(cor_matrix) <- 1
        result$random_effect_cor_matrix[[facet]] <- cor_matrix
      }
    }
  }

  result
}

#' Get AOv Residuals for Pairwise Correlation
#'
#' Extracts residuals from an aov model, handling the case where the model
#' was fit on a different subset of data than the pairwise data.
#'
#' @param aov_model Fitted aov model
#' @param pairwise_data Data frame with pairwise complete cases
#' @param resp Response variable name
#' @param orig_data Original data used to fit the model
#' @param resp_aov_formula The original aov formula string (with Error terms)
#'   used to fit the model. When provided, the model is refit on pairwise data
#'   using this formula to preserve the Error() structure.
#' @return Vector of residuals aligned with pairwise_data
#' @keywords internal
get_aov_residuals_for_pair <- function(aov_model, pairwise_data, resp, orig_data,
                                       resp_aov_formula = NULL) {
  resids <- NULL

  if (inherits(aov_model, "aovlist")) {
    if (!is.null(aov_model[["Within"]])) {
      resids_raw <- aov_model[["Within"]]$residuals
      if (is.matrix(resids_raw)) {
        resids <- as.vector(resids_raw)
      } else {
        resids <- resids_raw
      }
    }
  } else {
    resids <- residuals(aov_model)
  }

  if (is.null(resids)) {
    return(rep(NA, nrow(pairwise_data)))
  }

  # For unbalanced designs, we need to refit on pairwise data
  # because residuals length won't match pairwise_data
  if (length(resids) != nrow(pairwise_data)) {
    # Need to refit model on pairwise data
    # Use the original formula (with Error terms) if available,
    # otherwise fall back to a simple model
    tryCatch(
      {
        if (!is.null(resp_aov_formula)) {
          # Refit with the original formula (preserves Error() structure)
          refit_model <- aov(as.formula(resp_aov_formula), data = pairwise_data)
          if (inherits(refit_model, "aovlist") && !is.null(refit_model[["Within"]])) {
            resids_raw <- refit_model[["Within"]]$residuals
            resids <- if (is.matrix(resids_raw)) as.vector(resids_raw) else resids_raw
          } else {
            resids <- residuals(refit_model)
          }
        } else {
          # Fallback: simple model without Error() terms
          formula_str <- paste(resp, "~ .")
          refit_model <- aov(as.formula(formula_str), data = pairwise_data)
          resids <- residuals(refit_model)
        }
      },
      error = function(e) {
        resids <<- rep(NA, nrow(pairwise_data))
      })
  }

  resids
}

#' Compute Correlations for Multivariate Mom Model
#'
#' @param data The data frame
#' @param responses Character vector of response variable names
#' @param random_facets Character vector of random effect facet names
#' @param aov_models List of fitted aov models (one per response variable)
#' @param random_effect_cors List of random effect correlation specifications
#' @param compute_rescor Whether to compute residual correlations
#' @param compute_re_cors Whether to compute random effect correlations
#' @keywords internal
compute_mom_multivariate_correlations <- function(data, responses, random_facets, aov_models, aov_formulas, random_effect_cors = list(),
                                                  compute_rescor = FALSE, compute_re_cors = FALSE) {
  if (length(responses) < 2) {
    return(NULL)
  }

  result <- list(
    residual_cor = NULL,
    random_effect_cor = list(),
    residual_cor_matrix = NULL,
    random_effect_cor_matrix = list()
  )

  n <- nrow(data)

  # Handle missing values - remove rows with any NA in response variables
  # This is needed for random effect correlation computation
  complete_rows <- complete.cases(data[, responses, drop = FALSE])
  data_clean <- data[complete_rows, , drop = FALSE]
  n_clean_orig <- nrow(data_clean)

  if (n_clean_orig < n) {
    warning(
      "Cases with missing values in response variables have been removed. ",
      "Computing correlations on ", n_clean_orig, " complete cases out of ", n, " original cases.",
      call. = FALSE
    )
  }

  # Re-fit aov models on complete data to ensure residuals align properly
  # This is necessary because aov with Error() terms may handle missing values differently
  aov_models_clean <- list()
  for (resp in responses) {
    if (!is.null(aov_models[[resp]])) {
      aov_model <- aov_models[[resp]]
      # Re-fit on clean data if original model had different number of observations
      # This ensures residuals will have correct length
      if (inherits(aov_model, "aovlist")) {
        # Get the original model residuals length
        if (!is.null(aov_model[["Within"]])) {
          orig_resids <- aov_model[["Within"]]$residuals
          orig_len <- if (is.matrix(orig_resids)) nrow(orig_resids) else length(orig_resids)

          if (orig_len != n_clean_orig) {
            # Need to re-fit on clean data
            # Use the stored formula string instead of trying to extract from aovlist
            # (formula() doesn't work on aovlist objects with Error terms)
            if (!is.null(aov_formulas[[resp]])) {
              resp_aov_formula <- aov_formulas[[resp]]

              aov_model_clean <- tryCatch(
                aov(as.formula(resp_aov_formula), data = data_clean),
                error = function(e) {
                  NULL
                }
              )
              if (!is.null(aov_model_clean)) {
                aov_models_clean[[resp]] <- aov_model_clean
              }
            } else {
              # Fallback: skip re-fitting if formula not available
              aov_models_clean[[resp]] <- aov_model
            }
          } else {
            aov_models_clean[[resp]] <- aov_model
          }
        } else {
          aov_models_clean[[resp]] <- aov_model
        }
      } else {
        aov_models_clean[[resp]] <- aov_model
      }
    }
  }

  # Compute residual correlations ONLY if compute_rescor is TRUE
  if (compute_rescor) {
    # Extract residuals from each aov model
    # For aovlist objects (from aov with Error() terms), residuals() returns NULL
    # We need to extract residuals from the "Within" stratum manually
    # The residuals from each response model should be aligned by position
    resid_matrix <- NULL

    for (resp in responses) {
      # Use the clean models (re-fitted on complete cases if needed)
      aov_model <- aov_models_clean[[resp]]

      if (!is.null(aov_model)) {
        # Try to get residuals from the aovlist
        resids <- NULL

        # Check if it's an aovlist
        if (inherits(aov_model, "aovlist")) {
          # For aovlist, residuals are in the "Within" stratum
          # This represents the residual variance after accounting for all random effects
          if (!is.null(aov_model[["Within"]])) {
            resids_raw <- aov_model[["Within"]]$residuals
            # Extract as vector if it's a matrix
            if (!is.null(resids_raw) && is.matrix(resids_raw)) {
              resids <- as.vector(resids_raw)
            } else if (!is.null(resids_raw)) {
              resids <- resids_raw
            }
          }
        } else {
          # For regular aov objects
          resids <- residuals(aov_model)
        }

        # Only proceed if we got valid residuals
        if (!is.null(resids) && length(resids) > 0) {
          # Handle case where residuals length doesn't match clean data
          # This can happen when aov with Error() terms has different structure
          if (length(resids) != n_clean_orig) {
            if (length(resids) > n_clean_orig) {
              # Subset residuals to match clean data length
              # This handles the case where model was fit on more data
              resids <- resids[1:n_clean_orig]
            } else {
              # Pad with NA if residuals are fewer than expected
              resids <- c(resids, rep(NA, n_clean_orig - length(resids)))
            }
          }

          # Create a data frame - residuals are aligned by position
          if (is.null(resid_matrix)) {
            resid_matrix <- data.frame(resid = resids)
            names(resid_matrix)[1] <- resp
          } else {
            resid_matrix[[resp]] <- resids
          }
        }
      }
    }

    # Use residuals directly - they are already aligned with clean data from the first filtering step
    # The models were re-fitted on complete cases (data_clean) in the earlier step
    resid_clean <- resid_matrix
    n_clean <- n_clean_orig

    if (n_clean < 2 || is.null(resid_clean) || ncol(resid_clean) < 2) {
      warning("Not enough complete cases to compute residual correlations", call. = FALSE)
    } else {
      # Use the actual response columns that were successfully extracted
      valid_responses <- names(resid_clean)

      # Compute covariance matrix of residuals
      cov_matrix <- cov(resid_clean, use = "complete.obs")

      # Convert covariance matrix to correlation matrix
      # r_ij = sigma_ij / sqrt(sigma_ii * sigma_jj)
      sd_vec <- sqrt(diag(cov_matrix))
      cor_matrix <- cov_matrix
      for (i in seq_along(valid_responses)) {
        for (j in seq_along(valid_responses)) {
          if (i != j && !is.na(cov_matrix[i, j])) {
            cor_matrix[i, j] <- cov_matrix[i, j] / (sd_vec[i] * sd_vec[j])
          }
        }
      }
      diag(cor_matrix) <- 1

      # Store as data frame for compatibility
      resid_cor_list <- list()
      for (i in seq_along(valid_responses)) {
        for (j in seq_along(valid_responses)) {
          if (i < j) {
            resp_i <- valid_responses[i]
            resp_j <- valid_responses[j]
            cor_ij <- cor_matrix[i, j]

            if (!is.na(cor_ij)) {
              # Compute standard error and confidence intervals
              se_ij <- sqrt((1 - cor_ij^2) / (n_clean - 2))

              # Fisher's z transformation for CI
              z_ij <- atanh(cor_ij)
              se_z <- 1 / sqrt(n_clean - 3)
              z_lower <- z_ij - 1.96 * se_z
              z_upper <- z_ij + 1.96 * se_z
              lower <- tanh(z_lower)
              upper <- tanh(z_upper)

              resid_cor_list[[length(resid_cor_list) + 1]] <- data.frame(
                dim1 = resp_i,
                dim2 = resp_j,
                estimate = cor_ij,
                se = se_ij,
                lower = lower,
                upper = upper,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }

      # Set up matrix format
      rownames(cor_matrix) <- valid_responses
      colnames(cor_matrix) <- valid_responses
      result$residual_cor_matrix <- cor_matrix

      # Also store as data frame
      if (length(resid_cor_list) > 0) {
        result$residual_cor <- do.call(rbind, resid_cor_list)
      }
    }
  }

  # Compute random effect correlations ONLY if compute_re_cors is TRUE
  if (compute_re_cors) {
    # Only compute for facets that have correlation specifications
    for (facet in random_facets) {
      # Skip if this facet is not in random_effect_cors
      if (!(facet %in% names(random_effect_cors))) {
        next
      }

      if (!facet %in% names(data_clean)) next

      facet_means <- aggregate(
        data_clean[, responses, drop = FALSE],
        by = list(facet = data_clean[[facet]]),
        FUN = mean
      )

      if (nrow(facet_means) > 2) {
        facet_cor_list <- list()

        for (i in seq_along(responses)) {
          for (j in seq_along(responses)) {
            if (i < j) {
              resp_i <- responses[i]
              resp_j <- responses[j]

              # Compute correlation of facet means
              y_i <- facet_means[[resp_i]]
              y_j <- facet_means[[resp_j]]
              cor_ij <- cor(y_i, y_j)

              n_facet <- nrow(facet_means)
              se_ij <- sqrt((1 - cor_ij^2) / (n_facet - 2))

              z_ij <- atanh(cor_ij)
              se_z <- 1 / sqrt(n_facet - 3)
              z_lower <- z_ij - 1.96 * se_z
              z_upper <- z_ij + 1.96 * se_z
              lower <- tanh(z_lower)
              upper <- tanh(z_upper)

              facet_cor_list[[length(facet_cor_list) + 1]] <- data.frame(
                dim1 = resp_i,
                dim2 = resp_j,
                estimate = cor_ij,
                se = se_ij,
                lower = lower,
                upper = upper,
                stringsAsFactors = FALSE
              )
            }
          }
        }

        if (length(facet_cor_list) > 0) {
          result$random_effect_cor[[facet]] <- do.call(rbind, facet_cor_list)

          # Build matrix
          n_resp <- length(responses)
          cor_matrix <- matrix(NA, n_resp, n_resp)
          rownames(cor_matrix) <- responses
          colnames(cor_matrix) <- responses
          for (i in seq_len(nrow(result$random_effect_cor[[facet]]))) {
            r <- result$random_effect_cor[[facet]][i, ]
            cor_matrix[r$dim1, r$dim2] <- r$estimate
            cor_matrix[r$dim2, r$dim1] <- r$estimate
          }
          diag(cor_matrix) <- 1
          result$random_effect_cor_matrix[[facet]] <- cor_matrix
        }
      }
    }
  }
  result
}

#' Build aov Formula from lme4-style Formula
#'
#' @keywords internal
extract_design_info <- function(anova_table, random_facets) {
  # Get row names from ANOVA table (excluding Intercept and Error)
  row_names <- rownames(anova_table)

  # Identify the sources
  sources <- row_names[!row_names %in% c("Intercept", "Residuals")]

  # Determine which are facets vs interactions
  info <- list(
    sources = sources,
    n_sources = length(sources),
    df = anova_table$Df,
    ms = anova_table$`Mean Sq`
  )

  info
}

#' Compute Variance Components from ANOVA Table using Method of Moments
#'
#' @keywords internal
compute_mom_variance_components <- function(anova_table, design_info) {
  # Get sources (excluding intercept and residuals)
  sources <- rownames(anova_table)
  sources <- sources[!sources %in% c("Intercept", "Residuals")]

  if (length(sources) == 0) {
    stop("No variance components to extract from ANOVA table", call. = FALSE)
  }

  # Extract df and MS for each source
  df_vec <- anova_table$Df[sources]
  ms_vec <- anova_table$`Mean Sq`[sources]

  # Get residual df and ms
  resid_idx <- which(rownames(anova_table) == "Residuals")
  if (length(resid_idx) == 0) {
    stop("No residual term found in ANOVA table", call. = FALSE)
  }
  resid_df <- anova_table$Df[resid_idx]
  resid_ms <- anova_table$`Mean Sq`[resid_idx]

  # For method of moments, the variance component is estimated as:
  # sigma^2 = MS_between - MS_within / n
  # where n is the number of observations per group

  # Build a simple model: each source's variance component
  # This is a simplified version - full implementation would build EMS matrix

  results <- list()

  for (i in seq_along(sources)) {
    source <- sources[i]
    df <- df_vec[i]
    ms <- ms_vec[i]

    # Determine if this is a main effect or interaction
    type <- if (grepl(":", source)) "interaction" else "main"

    # Simple MoM estimate: var = (MS - residual_MS) / n_effective
    # For a balanced design, n_effective = total_n / (df + 1)
    # This is simplified - proper implementation needs EMS system

    # Use residual MS as estimate if MS < resid_MS (can happen with MoM)
    var_est <- max(0, (ms - resid_ms))

    # Compute approximate degrees of freedom using Satterthwaite
    # df approx = (MS - resid_MS)^2 / (MS^2/df + resid_MS^2/resid_df)
    if (var_est > 0 && ms > resid_ms) {
      num <- (ms - resid_ms)^2
      denom <- (ms^2 / df) + (resid_ms^2 / resid_df)
      approx_df <- num / denom
    } else {
      approx_df <- df
    }

    results[[length(results) + 1]] <- data.frame(
      component = source,
      facet = source,
      type = type,
      var = var_est,
      df = df,
      ms = ms,
      approx_df = approx_df,
      stringsAsFactors = FALSE
    )
  }

  # Add residual variance
  results[[length(results) + 1]] <- data.frame(
    component = "Residual",
    facet = "Residual",
    type = "residual",
    var = resid_ms,
    df = resid_df,
    ms = resid_ms,
    approx_df = resid_df,
    stringsAsFactors = FALSE
  )

  vc_df <- do.call(rbind, results)
  vc_tibble <- tibble::as_tibble(vc_df)

  # Calculate percentage of total variance
  total_var <- sum(vc_tibble$var, na.rm = TRUE)
  vc_tibble$pct <- (vc_tibble$var / total_var) * 100

  vc_tibble
}

#' Extract Variance Components from Method of Moments Fit
#'
#' Extracts variance components from a MoM fit with confidence intervals
#' using asymptotic approximations.
#'
#' @param model A momfit object from [fit_mom()].
#' @param conf_level Confidence level for intervals (default 0.95).
#' @param formula The original formula (used to extract response name). Optional.
#' @return A tibble with columns:
#'   \item{component}{Name of the variance component}
#'   \item{facet}{Associated facet name}
#'   \item{type}{Type: "main", "interaction", or "residual"}
#'   \item{var}{Point estimate of variance}
#'   \item{pct}{Percentage of total variance}
#'   \item{lower}{Lower confidence interval bound}
#'   \item{upper}{Upper confidence interval bound}
#'   \item{se}{Standard error of the variance estimate}
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   score = rnorm(100),
#'   person = factor(rep(1:20, each = 5)),
#'   item = factor(rep(1:5, times = 20))
#' )
#' model <- fit_mom(score ~ (1 | person) + (1 | item), data)
#' extract_vc_mom(model)
#' }
#'
#' @export
#' @rdname extract_vc_mom
extract_vc_mom <- function(model, conf_level = 0.95, formula = NULL) {
  if (!inherits(model, "momfit")) {
    stop("model must be a momfit object", call. = FALSE)
  }

  vc <- model$variance_components

  alpha <- 1 - conf_level
  z_crit <- qnorm(1 - alpha / 2)

  vc$lower <- NA_real_
  vc$upper <- NA_real_
  vc$se <- NA_real_

  for (i in seq_len(nrow(vc))) {
    var_est <- vc$var[i]

    # Use Satterthwaite approximate df if available, otherwise fall back to df
    approx_df <- if ("approx_df" %in% names(vc)) vc$approx_df[i] else vc$df[i]

    # SE of variance component using Satterthwaite approximation
    # Var(sigma^2) ≈ 2 * sigma^4 / df_approx
    if (!is.na(var_est) && var_est > 0 && !is.na(approx_df) && approx_df > 0) {
      se_var <- sqrt(2 * var_est^2 / approx_df)
    } else {
      se_var <- NA_real_
    }

    vc$se[i] <- se_var

    if (!is.na(se_var) && se_var > 0) {
      vc$lower[i] <- max(0, var_est - z_crit * se_var)
      vc$upper[i] <- var_est + z_crit * se_var
    }
  }

  # The 'facet' column from compute_mom_vc_from_results is not needed
  # and should not be used to set 'dim'. The dim should be the response variable name.
  # Remove facet column if present (it duplicates component information)
  if ("facet" %in% names(vc)) {
    vc$facet <- NULL
  }

  if (!"dim" %in% names(vc)) {
    dim_name <- "response"
    if (!is.null(formula)) {
      dim_name <- extract_response_names(formula)
      if (length(dim_name) == 0) dim_name <- "response"
    }
    vc$dim <- dim_name
  }

  output_cols <- c("component", "dim", "type", "var", "pct", "lower", "upper", "se")
  output_cols <- output_cols[output_cols %in% names(vc)]

  vc_out <- vc[, output_cols, drop = FALSE]

  # Round numeric columns to 3 decimal places for consistent display
  numeric_cols <- c("var", "pct", "lower", "upper", "se")
  for (col in numeric_cols) {
    if (col %in% names(vc_out)) {
      vc_out[[col]] <- round(vc_out[[col]], 3)
    }
  }

  vc_out
}

#' Summary Method for momfit Objects
#'
#' Prints a summary of a momfit object including variance components.
#'
#' @param object A momfit object from [fit_mom()].
#' @param ... Additional arguments passed to other methods.
#' @return Side effect: prints summary to console. Returns the object invisibly.
#'
#' @method summary momfit
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   score = rnorm(100),
#'   person = factor(rep(1:20, each = 5)),
#'   item = factor(rep(1:5, times = 20))
#' )
#' model <- fit_mom(score ~ (1 | person) + (1 | item), data)
#' summary(model)
#' }
#'
#' @export
summary.momfit <- function(object, ...) {
  cat("Method of Moments Fit\n")
  cat("====================\n\n")
  cat("Response:", object$response, "\n")
  cat("Random facets:", paste(object$random_facets, collapse = ", "), "\n")
  cat("\nVariance Components:\n")
  print(object$variance_components[, c("component", "var", "pct")])
  invisible(object)
}
