#' Formula Helper Functions
#'
#' Utilities for parsing and manipulating formulas in the context of
#' generalizability theory analyses.
#'
#' @name formula-helpers
#' @keywords internal
NULL

#' Detect if Formula is Multivariate
#'
#' Checks whether a formula represents a multivariate model by looking for
#' brms-specific multivariate syntax (mvbind or set_rescor).
#'
#' @param formula A formula or brmsformula object.
#' @return Logical indicating whether the formula is multivariate.
#'
#' @keywords internal
is_multivariate <- function(formula) {
  if (inherits(formula, "mvbrmsformula")) {
    return(TRUE)
  }

  if (inherits(formula, "bform") && !is.null(formula$responses) && length(formula$responses) > 1) {
    return(TRUE)
  }

  if (!is.null(formula$responses) && length(formula$responses) > 1) {
    return(TRUE)
  }

  formula_char <- deparse(formula)

  patterns <- c(
    "mvbind\\s*\\(",
    "set_rescor\\s*\\(",
    "bf\\s*\\("
  )

  any(sapply(patterns, function(p) any(grepl(p, formula_char))))
}

#' Extract set_rescor Setting from Formula
#'
#' Parses the formula to detect if set_rescor() is specified and returns its value.
#' Returns FALSE by default (no residual correlations) if not specified.
#'
#' @param formula A formula object.
#' @return Logical indicating whether to estimate residual correlations (default FALSE).
#'
#' @keywords internal
extract_rescor_setting <- function(formula) {
  # Check if it's a brmsformula with rescor set
  if (inherits(formula, "brmsformula") || inherits(formula, "mvbrmsformula")) {
    # Try to get rescor from the formula object
    if (!is.null(formula$rescor)) {
      return(formula$rescor)
    }
  }
  
  # Parse from character representation
  # Collapse multiple lines into single string for pattern matching
  formula_char <- paste(deparse(formula), collapse = " ")

  # Look for set_rescor() in the formula
  if (grepl("set_rescor\\s*\\(\\s*TRUE\\s*\\)", formula_char, ignore.case = TRUE)) {
    return(TRUE)
  }
  
  if (grepl("set_rescor\\s*\\(\\s*FALSE\\s*\\)", formula_char, ignore.case = TRUE)) {
    return(FALSE)
  }
  
  # Default: no residual correlations
  FALSE
}

#' Extract Random Effect Correlation Specifications from Formula
#'
#' Parses the formula to detect random effects with correlation structures.
#' In brms, this is specified as (1|cor_name|facet) where cor_name is the
#' name of a correlation matrix.
#'
#' @param formula A formula object.
#' @return Named list where names are facet names and values are correlation matrix names.
#'   For example: list(person = "cor_p", item = "cor_i")
#'
#' @keywords internal
extract_random_effect_cor <- function(formula) {
  # Collapse multiple lines into single string for pattern matching
  formula_char <- paste(deparse(formula), collapse = " ")
  
  # Split formula by + to get individual terms
  # This prevents matching across multiple terms
  rhs <- sub(".*~", "", formula_char)
  terms <- strsplit(rhs, "\\+")[[1]]
  
  result <- list()
  
  # Pattern to match (1|cor_name|facet) - random effect with correlation
  # This requires exactly 2 pipe characters in the term
  pattern <- "^\\s*\\(\\s*1\\s*\\|\\s*(.+)\\s*\\|\\s*(.+)\\s*\\)\\s*$"
  
  for (term in terms) {
    term <- trimws(term)
    
    # Count pipes - correlation syntax has exactly 2 pipes
    pipe_count <- length(gregexpr("\\|", term)[[1]])
    
    if (pipe_count == 2) {
      # Extract the two parts
      # Remove leading ( and trailing )
      term_clean <- sub("^\\s*\\(\\s*1\\s*\\|\\s*", "", term)
      term_clean <- sub("\\s*\\)\\s*$", "", term_clean)
      
      parts <- strsplit(term_clean, "\\s*\\|\\s*")[[1]]
      if (length(parts) == 2) {
        cor_name <- trimws(parts[1])
        facet_name <- trimws(parts[2])
        result[[facet_name]] <- cor_name
      }
    }
  }
  
  result
}

#' Parse Formula to Extract Facets
#'
#' Extracts information about the response variable, fixed effects, and
#' random effects from a G-study formula.
#'
#' @param formula A formula object.
#' @return A list with components:
#'   \item{response}{Character string naming the response variable}
#'   \item{fixed}{Character vector of fixed effect terms}
#'   \item{random}{Character vector of random effect terms (as parsed)}
#'   \item{random_facets}{Character vector of individual facet names from random effects}
#'   \item{random_facet_specs}{Character vector of facet specifications as user specified (e.g., "Rater:Task")}
#'
#' @keywords internal
parse_g_formula <- function(formula) {
  if (inherits(formula, "mvbrmsformula") || inherits(formula, "bform")) {
    if (!is.null(formula$forms) && is.list(formula$forms) && length(formula$forms) > 0) {
      first_formula <- formula$forms[[1]]
      formula_terms <- terms(first_formula$formula)
      response <- names(formula$forms)[1]
      formula_for_findbars <- first_formula$formula
    } else {
      stop("Cannot parse multivariate brmsformula: no forms found", call. = FALSE)
    }
  } else {
    formula_terms <- terms(formula)
    response <- all.vars(formula, max.names = 1)[1]
    formula_for_findbars <- formula
  }

  all_vars <- all.vars(formula)

  random_terms <- NULL
  random_facets <- character()
  random_facet_specs <- character()

  if (requireNamespace("reformulas", quietly = TRUE)) {
    random_terms <- reformulas::findbars(formula_for_findbars)

    if (!is.null(random_terms) && length(random_terms) > 0) {
      random_facet_specs <- unlist(lapply(random_terms, function(term) {
        term_char <- deparse(term)
        matches <- regmatches(term_char, regexpr("\\|[^)]+", term_char))
        if (length(matches) > 0) {
          facet_part <- sub("^\\|", "", matches)
          if (grepl("\\|", facet_part)) {
            parts <- strsplit(facet_part, "\\|")[[1]]
            facet_part <- trimws(parts[length(parts)])
          }
          return(trimws(facet_part))
        }
        character()
      }))

      random_facets <- unlist(lapply(random_terms, function(term) {
        term_char <- deparse(term)
        matches <- regmatches(term_char, regexpr("\\|[^)]+", term_char))
        if (length(matches) > 0) {
          facet_part <- sub("^\\|", "", matches)
          if (grepl("\\|", facet_part)) {
            parts <- strsplit(facet_part, "\\|")[[1]]
            facet_part <- trimws(parts[length(parts)])
          }
          facets <- strsplit(facet_part, ":")[[1]]
          facets <- trimws(facets)
          return(facets)
        }
        character()
      }))
    }
  }

  fixed_terms <- attr(formula_terms, "term.labels")
  if (!is.null(random_terms) && length(random_terms) > 0) {
    random_pattern <- sapply(random_terms, function(x) {
      paste0("\\Q", deparse(x), "\\E")
    })
    for (rp in random_pattern) {
      fixed_terms <- fixed_terms[!grepl(rp, fixed_terms, fixed = TRUE)]
    }
  }

  list(
    response = response,
    fixed = fixed_terms,
    random = if (!is.null(random_terms)) sapply(random_terms, deparse) else character(),
    random_facets = random_facets,
    random_facet_specs = random_facet_specs
  )
}

#' Validate Formula for G-Study
#'
#' Checks that a formula is valid for a G-study analysis. The formula should
#' have a response variable and at least one random effect term.
#'
#' @param formula A formula object.
#' @param backend Character string indicating the backend to use ("auto", "lme4", "brms", or "mom").
#' @return TRUE if valid, otherwise raises an error.
#'
#' @keywords internal
validate_formula <- function(formula, backend = "auto") {
  valid_classes <- c("formula", "brmsformula", "mvbrmsformula", "bform")
  if (!any(sapply(valid_classes, function(cls) inherits(formula, cls)))) {
    stop("formula must be a formula or brmsformula object", call. = FALSE)
  }

  parsed <- parse_g_formula(formula)

  if (length(parsed$response) == 0 || parsed$response == "") {
    stop("Formula must have a response variable on the left-hand side", call. = FALSE)
  }

  if (length(parsed$random) == 0) {
    stop(
      "Formula must contain at least one random effect term.\n",
      "Use syntax like (1 | facet) to specify random effects.",
      call. = FALSE
    )
  }

  is_mv <- is_multivariate(formula)
  if (is_mv && backend == "lme4") {
    stop(
      "Multivariate formulas (containing mvbind or set_rescor) require brms backend.\n",
      "Use: backend = 'brms'",
      call. = FALSE
    )
  }

  TRUE
}

#' Detect Facets from a Formula
#'
#' Identifies which variables in the formula are facets (sources of variance)
#' in the generalizability theory sense.
#'
#' @param formula A formula object.
#' @param data A data frame (optional, used for validation).
#' @return A character vector of facet names.
#'
#' @keywords internal
detect_facets <- function(formula, data = NULL) {
  parsed <- parse_g_formula(formula)
  
  # Get unique facet names from random effects
  facets <- unique(parsed$random_facets)
  
  # Validate against data if provided
  if (!is.null(data)) {
    missing_facets <- setdiff(facets, names(data))
    if (length(missing_facets) > 0) {
      stop(
        "Facet(s) not found in data: ", 
        paste(missing_facets, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  facets
}

#' Convert Formula to Backend-Specific Format
#'
#' Converts a G-study formula to the format required by a specific backend.
#' For lme4, this removes any brms-specific syntax.
#'
#' @param formula A formula object.
#' @param backend Character string indicating the target backend.
#' @return A formula object in the backend-specific format.
#'
#' @keywords internal
convert_formula <- function(formula, backend) {
  if (backend == "lme4") {
    # Check if formula is a brmsformula
    if (inherits(formula, "brmsformula")) {
      # Extract the base formula
      formula <- formula$formula
    }
    
    # Check for multivariate syntax and warn
    if (is_multivariate(formula)) {
      warning(
        "Multivariate formulas cannot be converted to lme4 format. ",
        "Using brms backend instead.",
        call. = FALSE
      )
    }
  }
  
  formula
}

#' Extract Facet Names from Formula or Variance Components
#'
#' Extracts the names of facets from either a formula or variance components tibble.
#'
#' @param formula Model formula (optional).
#' @param vc Variance components tibble (optional).
#' @return Character vector of facet names.
#'
#' @keywords internal
extract_facets <- function(formula = NULL, vc = NULL) {
  facets <- character()
  
  # Try to extract from formula first
  if (!is.null(formula)) {
    parsed <- parse_g_formula(formula)
    facets <- parsed$random_facets
  }
  
  # If no facets from formula, try variance components
  if (length(facets) == 0 && !is.null(vc)) {
    if ("facet" %in% names(vc)) {
      facets <- unique(vc$facet)
      facets <- facets[facets != "Residual" & facets != ""]
    } else if ("component" %in% names(vc)) {
      facets <- unique(vc$component)
      facets <- facets[facets != "Residual" & facets != ""]
    }
  }
  
  unique(facets)
}

#' Extract Facet Specifications Preserving User's Order
#'
#' Extracts facet specifications as the user specified them (e.g., "Rater:Task"
#' vs "Task:Rater"), preserving the original order.
#'
#' @param formula A formula object.
#' @return Character vector of facet specifications in user-specified order.
#'
#' @keywords internal
extract_facet_specs <- function(formula = NULL) {
  if (is.null(formula)) {
    return(character())
  }

  parsed <- parse_g_formula(formula)
  parsed$random_facet_specs
}

#' Extract Response Variable Names from Formula
#'
#' Extracts the names of response (dependent) variables from a formula.
#' Works for both univariate and multivariate (brms mvbind) formulas.
#'
#' @param formula A formula or brmsformula object.
#' @return Character vector of response variable names. For univariate formulas,
#'   returns a single name. For multivariate formulas, returns all response names.
#'
#' @keywords internal
extract_response_names <- function(formula) {
  formula_char <- deparse(formula)

  if (is_multivariate(formula)) {
    resp_names <- character()

    mvbind_match <- regmatches(formula_char, regexpr("mvbind\\s*\\([^)]+\\)", formula_char))
    if (length(mvbind_match) > 0) {
      resp_str <- sub("mvbind\\s*\\(", "", mvbind_match)
      resp_str <- sub("\\)\\s*$", "", resp_str)
      resp_names <- unlist(strsplit(resp_str, "\\s*,\\s*"))
      resp_names <- trimws(resp_names)
    }

    if (length(resp_names) == 0 && inherits(formula, "brmsformula")) {
      if (requireNamespace("brms", quietly = TRUE)) {
        tryCatch({
          resp_names <- brms::response_names(formula)
        }, error = function(e) {
          resp_names <- character()
        })
      }
    }

    if (length(resp_names) == 0) {
      resp_names <- "response"
    }

    return(resp_names)
  } else {
    if (inherits(formula, "brmsformula") && requireNamespace("brms", quietly = TRUE)) {
      tryCatch({
        return(brms::response_names(formula))
      }, error = function(e) {})
    }

    all_vars <- all.vars(formula)
    if (length(all_vars) > 0) {
      return(all_vars[1])
    }
    return("response")
  }
}

#' Parse Residual Facets from Formula
#'
#' Determines which facets make up the residual (lowest-level) variance component
#' based on the formula and data structure. The residual is typically the interaction
#' of all non-nested facets.
#'
#' @param formula A formula object with lme4-style random effects (e.g., `score ~ (1|p) + (1|i:d)`).
#' @param data Optional data frame to detect actual nesting patterns from the data.
#'   If NULL, assumes all facets are crossed (most conservative).
#'
#' @return A character string representing the facets that make up the residual,
#'   with facets separated by colons (e.g., `"p:i:d"`).
#'
#' @details
#' The function determines the residual composition by:
#' \enumerate{
#'   \item Parsing the formula to extract all unique facets from random effects
#'   \item If `data` is provided, detecting which facets are nested within others
#'   \item Returning the interaction of all non-nested (crossed) facets
#' }
#'
#' ## Nesting Detection
#' When `data` is provided, the function analyzes the actual data structure to
#' determine if facets appearing in interaction terms (like `i:d`) are nested
#' or crossed:
#' \itemize{
#'   \item **Nested**: If each level of `i` appears with only one level of `d`,
#'     then `i` is nested in `d`. The nested facet is absorbed and the residual
#'     excludes the nesting facet (e.g., residual = `"p:i"` instead of `"p:i:d"`).
#'   \item **Crossed**: If levels of `i` appear with multiple levels of `d`,
#'     the facets are crossed. The residual includes all facets (e.g., `"p:i:d"`).
#' }
#'
#' @examples
#' # Fully crossed design: p, i, d all crossed with each other
#' f1 <- score ~ (1|p) + (1|i) + (1|d)
#' parse_residual_facets(f1) # Returns "p:i:d"
#'
#' # With nested data: i nested in d
#' nested_data <- data.frame(
#'   score = rnorm(100),
#'   p = factor(rep(1:10, 10)),
#'   i = factor(rep(1:5, each = 20)), # i levels repeat within d
#'   d = factor(rep(1:2, each = 50))
#' )
#' f2 <- score ~ (1|p) + (1|i:d) + (1|d)
#' parse_residual_facets(f2, nested_data) # May return "p:i" if i nested in d
#'
#' @seealso
#' \code{\link{detect_crossing_patterns_from_data}} for the underlying nesting detection
#' \code{\link{parse_g_formula}} for parsing formula structure
#'
#' @keywords internal
parse_residual_facets <- function(formula, data = NULL) {
  parsed <- parse_g_formula(formula)
  facet_specs <- parsed$random_facet_specs

  if (length(facet_specs) == 0) {
    return("")
  }

  base_facets <- unique(parsed$random_facets)

  if (length(base_facets) == 1) {
    return(base_facets[1])
  }

  if (is.null(data)) {
    return(paste(base_facets, collapse = ":"))
  }

  base_facets <- base_facets[base_facets %in% names(data)]

  nesting_facets <- character()
  for (i in seq_along(base_facets)) {
    for (j in seq_along(base_facets)) {
      if (i == j) next

      f_a <- base_facets[i]
      f_b <- base_facets[j]

      levels_a <- unique(data[[f_a]])

      levels_b_per_a <- sapply(levels_a, function(lvl_a) {
        subset_data <- data[data[[f_a]] == lvl_a, ]
        length(unique(subset_data[[f_b]]))
      })

      min_levels <- min(levels_b_per_a)
      max_levels <- max(levels_b_per_a)

      if (min_levels == 1 && max_levels == 1) {
        nesting_facets <- c(nesting_facets, f_b)
      }
    }
  }

  nesting_facets <- unique(nesting_facets)
  residual_facets <- base_facets[!base_facets %in% nesting_facets]

  if (length(residual_facets) == 0) {
    return(paste(base_facets, collapse = ":"))
  }

  paste(residual_facets, collapse = ":")
}
