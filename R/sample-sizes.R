#' Sample Size Calculations for G-Studies
#'
#' Functions for calculating sample sizes for main effects, interactions,
#' residual, and nested effects in generalizability theory analyses.
#'
#' @name sample-sizes
#' @keywords internal
NULL

#' Calculate Comprehensive Sample Size Information
#'
#' Calculates sample sizes for all variance components in a G-study design,
#' including main effects, interactions, residual, and nested effects.
#'
#' @param formula A formula object with lme4-style random effects.
#' @param data A data frame containing the variables in the formula.
#' @param nested Optional list specifying nesting relationships.
#'   Named list where names are nested facets and values are nesting facets.
#'   For example, `list(task = "rater")` means task is nested within rater.
#'   If NULL, nesting is auto-detected from the data.
#'
#' @return A list with components:
#'   \item{main}{Named vector of sample sizes for main effects}
#'   \item{interactions}{Named vector of sample sizes for interaction terms specified in formula}
#'   \item{residual}{List with facets (character) and n (numeric) for the residual}
#'   \item{nested}{List of nested effect information, or NULL if none detected}
#'
#' @keywords internal
calculate_sample_size_info <- function(formula, data, nested = NULL) {
  if (is.null(data) || !is.data.frame(data)) {
    return(list(
      main = integer(0),
      interactions = integer(0),
      residual = list(facets = "", n = NA_integer_),
      nested = NULL
    ))
  }

  parsed <- parse_g_formula(formula)
  facet_specs <- parsed$random_facet_specs

  if (length(facet_specs) == 0) {
    return(list(
      main = integer(0),
      interactions = integer(0),
      residual = list(facets = "", n = NA_integer_),
      nested = NULL
    ))
  }

  main_n <- integer(0)
  interactions_n <- integer(0)

  for (spec in facet_specs) {
    if (grepl(":", spec)) {
      facets_in_interaction <- strsplit(spec, ":")[[1]]
      facets_in_interaction <- trimws(facets_in_interaction)

      if (all(facets_in_interaction %in% names(data))) {
        unique_combos <- unique(data[, facets_in_interaction, drop = FALSE])
        n_combos <- nrow(unique_combos)
        interactions_n <- c(interactions_n, n_combos)
        names(interactions_n)[length(interactions_n)] <- spec
      } else {
        interactions_n <- c(interactions_n, NA_integer_)
        names(interactions_n)[length(interactions_n)] <- spec
      }
    } else {
      facet_name <- trimws(spec)
      if (facet_name %in% names(data)) {
        n_levels <- length(unique(data[[facet_name]]))
        main_n <- c(main_n, n_levels)
        names(main_n)[length(main_n)] <- facet_name
      } else {
        main_n <- c(main_n, NA_integer_)
        names(main_n)[length(main_n)] <- facet_name
      }
    }
  }

  residual_facets <- parse_residual_facets(formula, data)
  residual_n <- calculate_residual_sample_size(data, residual_facets)

  all_facets <- unique(parsed$random_facets)
  nesting_info <- detect_nesting_patterns(data, all_facets)

  if (!is.null(nested) && length(nested) > 0) {
    nesting_info <- override_nesting_detection(nesting_info, nested, all_facets)
  }

  nested_info <- NULL
  if (length(nesting_info) > 0) {
    nested_info <- calculate_nested_sample_sizes(data, nesting_info)
  }

  list(
    main = main_n,
    interactions = interactions_n,
    residual = list(
      facets = residual_facets,
      n = residual_n
    ),
    nested = nested_info
  )
}

#' Calculate Sample Size for Residual
#'
#' Calculates the number of unique observations for the residual component.
#'
#' @param data A data frame.
#' @param residual_facets Character string of facets making up the residual,
#'   separated by colons (e.g., "person:rater:item").
#'
#' @return Integer number of unique combinations.
#'
#' @keywords internal
calculate_residual_sample_size <- function(data, residual_facets) {
  if (is.null(residual_facets) || residual_facets == "") {
    return(NA_integer_)
  }

  facets <- strsplit(residual_facets, ":")[[1]]
  facets <- trimws(facets)

  if (!all(facets %in% names(data))) {
    return(NA_integer_)
  }

  unique_combos <- unique(data[, facets, drop = FALSE])
  nrow(unique_combos)
}

#' Detect Nesting Patterns from Data
#'
#' Analyzes the data structure to determine if facets are nested.
#' A facet A is considered nested in facet B if each level of A appears
#' with only one level of B.
#'
#' @param data A data frame.
#' @param facets Character vector of facet names to analyze.
#'
#' @return A list where each element describes a nesting relationship:
#'   \item{nested_facet}{Name of the nested facet}
#'   \item{nesting_facet}{Name of the facet that the nested facet is within}
#'
#' @keywords internal
detect_nesting_patterns <- function(data, facets) {
  if (is.null(data) || length(facets) < 2) {
    return(list())
  }

  facets <- facets[facets %in% names(data)]

  nesting_relationships <- list()

  for (i in seq_along(facets)) {
    for (j in seq_along(facets)) {
      if (i == j) next

      f_nested <- facets[i]
      f_nesting <- facets[j]

      key <- paste0(f_nested, "_in_", f_nesting)
      if (key %in% names(nesting_relationships)) next

      levels_nested <- unique(data[[f_nested]])

      levels_nesting_per_nested <- sapply(levels_nested, function(lvl) {
        subset_data <- data[data[[f_nested]] == lvl, ]
        length(unique(subset_data[[f_nesting]]))
      })

      min_levels <- min(levels_nesting_per_nested)
      max_levels <- max(levels_nesting_per_nested)

      if (min_levels == 1 && max_levels == 1) {
        nesting_relationships[[key]] <- list(
          nested_facet = f_nested,
          nesting_facet = f_nesting
        )
      }
    }
  }

  nesting_relationships
}

#' Override Nesting Detection with User Specification
#'
#' Merges auto-detected nesting with user-specified nesting relationships.
#'
#' @param detected List of auto-detected nesting relationships.
#' @param user_nested Named list of user-specified nesting (nested = nesting).
#' @param all_facets Character vector of all facet names.
#'
#' @return Updated list of nesting relationships.
#'
#' @keywords internal
override_nesting_detection <- function(detected, user_nested, all_facets) {
  for (nested_facet in names(user_nested)) {
    nesting_facet <- user_nested[[nested_facet]]

    if (!nested_facet %in% all_facets) {
      warning(
        "Nested facet '", nested_facet, "' not found in formula facets. ",
        "Ignoring this nesting specification.",
        call. = FALSE
      )
      next
    }

    if (!nesting_facet %in% all_facets) {
      warning(
        "Nesting facet '", nesting_facet, "' not found in formula facets. ",
        "Ignoring this nesting specification for '", nested_facet, "'.",
        call. = FALSE
      )
      next
    }

    key <- paste0(nested_facet, "_in_", nesting_facet)

    conflicting_key <- paste0(nested_facet, "_in_")
    conflicting <- names(detected)[startsWith(names(detected), conflicting_key)]

    if (length(conflicting) > 0) {
      for (ck in conflicting) {
        if (ck != key) {
          detected[[ck]] <- NULL
        }
      }
    }

    detected[[key]] <- list(
      nested_facet = nested_facet,
      nesting_facet = nesting_facet
    )
  }

  detected
}

#' Calculate Sample Sizes for Nested Effects
#'
#' For each nesting relationship, calculates the number of nested units
#' per nesting unit, including mean and harmonic mean for unbalanced designs.
#'
#' @param data A data frame.
#' @param nesting_info List of nesting relationships from detect_nesting_patterns().
#'
#' @return A list where each element contains:
#'   \item{nested_facet}{Name of the nested facet}
#'   \item{nesting_facet}{Name of the nesting facet}
#'   \item{n_groups}{Number of groups (levels of nesting facet)}
#'   \item{mean_per_group}{Mean number of nested units per group}
#'   \item{harmonic_mean_per_group}{Harmonic mean of nested units per group}
#'   \item{counts_per_group}{Named vector of counts for each group}
#'
#' @keywords internal
calculate_nested_sample_sizes <- function(data, nesting_info) {
  if (is.null(data) || length(nesting_info) == 0) {
    return(NULL)
  }

  result <- list()

  for (key in names(nesting_info)) {
    rel <- nesting_info[[key]]

    nested_facet <- rel$nested_facet
    nesting_facet <- rel$nesting_facet

    if (!nested_facet %in% names(data) || !nesting_facet %in% names(data)) {
      next
    }

    nesting_levels <- unique(data[[nesting_facet]])

    counts_per_group <- sapply(nesting_levels, function(lvl) {
      subset_data <- data[data[[nesting_facet]] == lvl, ]
      length(unique(subset_data[[nested_facet]]))
    })
    names(counts_per_group) <- nesting_levels

    n_groups <- length(nesting_levels)
    mean_per_group <- mean(counts_per_group)
    harmonic_mean_per_group <- length(counts_per_group) / sum(1 / counts_per_group)

    result[[key]] <- list(
      nested_facet = nested_facet,
      nesting_facet = nesting_facet,
      n_groups = n_groups,
      mean_per_group = mean_per_group,
      harmonic_mean_per_group = harmonic_mean_per_group,
      counts_per_group = counts_per_group
    )
  }

  if (length(result) == 0) {
    return(NULL)
  }

  result
}

#' Format Sample Size Info for Display
#'
#' Creates a formatted string representation of sample size information
#' for use in print and summary methods.
#'
#' @param sample_size_info List from calculate_sample_size_info().
#' @param indent Number of spaces for indentation (default 1).
#'
#' @return Character vector of formatted lines.
#'
#' @keywords internal
format_sample_size_info <- function(sample_size_info, indent = 1) {
  if (is.null(sample_size_info)) {
    return(character(0))
  }

  indent_str <- paste(rep(" ", indent), collapse = "")

  lines <- character(0)

  if (length(sample_size_info$main) > 0) {
    for (i in seq_along(sample_size_info$main)) {
      n_val <- sample_size_info$main[i]
      if (!is.na(n_val)) {
        lines <- c(lines, paste0(indent_str, names(sample_size_info$main)[i], ": ", n_val))
      }
    }
  }

  if (length(sample_size_info$interactions) > 0) {
    for (i in seq_along(sample_size_info$interactions)) {
      n_val <- sample_size_info$interactions[i]
      if (!is.na(n_val)) {
        lines <- c(lines, paste0(indent_str, names(sample_size_info$interactions)[i], ": ", n_val))
      }
    }
  }

  if (!is.null(sample_size_info$residual) &&
      !is.null(sample_size_info$residual$facets) &&
      sample_size_info$residual$facets != "") {
    residual_text <- paste0("Residual (", sample_size_info$residual$facets, ")")
    if (!is.na(sample_size_info$residual$n)) {
      lines <- c(lines, paste0(indent_str, residual_text, ": ", sample_size_info$residual$n))
    }
  }

  if (!is.null(sample_size_info$nested) && length(sample_size_info$nested) > 0) {
    nested_lines <- format_nested_info(sample_size_info$nested, indent + 1)
    lines <- c(lines, nested_lines)
  }

  lines
}

#' Format Nested Info for Display
#'
#' Creates a formatted string representation of nested effect information.
#'
#' @param nested_info List from calculate_nested_sample_sizes().
#' @param indent Number of spaces for indentation (default 1).
#'
#' @return Character vector of formatted lines.
#'
#' @keywords internal
format_nested_info <- function(nested_info, indent = 1) {
  if (is.null(nested_info) || length(nested_info) == 0) {
    return(character(0))
  }

  indent_str <- paste(rep(" ", indent), collapse = "")
  lines <- character(0)

  for (key in names(nested_info)) {
    info <- nested_info[[key]]

    header <- paste0(info$nested_facet, " nested in ", info$nesting_facet)
    lines <- c(lines, paste0(indent_str, header, ":"))

    lines <- c(lines, paste0(indent_str, "  ", info$nesting_facet, "s: ", info$n_groups))

    if (abs(info$mean_per_group - info$harmonic_mean_per_group) > 0.001) {
      lines <- c(lines, paste0(
        indent_str, "  ", info$nested_facet, "s per ", info$nesting_facet,
        ": mean=", round(info$mean_per_group, 1),
        ", harmonic_mean=", round(info$harmonic_mean_per_group, 1)
      ))
    } else {
      lines <- c(lines, paste0(
        indent_str, "  ", info$nested_facet, "s per ", info$nesting_facet,
        ": ", round(info$mean_per_group, 1)
      ))
    }
  }

  lines
}

#' Convert Sample Size Info to Tibble Format
#'
#' Converts the list-based sample_size_info into a compact tibble format
#' suitable for display in print and summary methods.
#'
#' @param sample_size_info List from calculate_sample_size_info().
#'
#' @return A tibble with columns:
#'   \item{effect}{Name of the effect (facet, interaction, or nested description)}
#'   \item{type}{Type of effect: "main", "interaction", "residual", or "nested"}
#'   \item{n}{Sample size for the effect}
#'
#' @keywords internal
calculate_single_sample_size_tibble <- function(sample_size_info) {
  if (is.null(sample_size_info)) {
    return(tibble::tibble(
      effect = character(),
      type = character(),
      n = numeric()
    ))
  }
  
  rows <- list()
  
  # Main effects
  if (length(sample_size_info$main) > 0) {
    for (i in seq_along(sample_size_info$main)) {
      n_val <- sample_size_info$main[i]
      if (!is.na(n_val)) {
        rows[[length(rows) + 1]] <- tibble::tibble(
          effect = names(sample_size_info$main)[i],
          type = "main",
          n = n_val
        )
      }
    }
  }
  
  # Interactions
  if (length(sample_size_info$interactions) > 0) {
    for (i in seq_along(sample_size_info$interactions)) {
      n_val <- sample_size_info$interactions[i]
      if (!is.na(n_val)) {
        rows[[length(rows) + 1]] <- tibble::tibble(
          effect = names(sample_size_info$interactions)[i],
          type = "interaction",
          n = n_val
        )
      }
    }
  }
  
  # Residual
  if (!is.null(sample_size_info$residual) && 
      sample_size_info$residual$facets != "" && 
      !is.na(sample_size_info$residual$n)) {
    rows[[length(rows) + 1]] <- tibble::tibble(
      effect = paste0(sample_size_info$residual$facets, " (res)"),
      type = "residual",
      n = sample_size_info$residual$n
    )
  }
  
  # Nested effects
  if (!is.null(sample_size_info$nested) && length(sample_size_info$nested) > 0) {
    for (key in names(sample_size_info$nested)) {
      info <- sample_size_info$nested[[key]]
      compound_name <- paste0(info$nested_facet, " per ", info$nesting_facet)
      rows[[length(rows) + 1]] <- tibble::tibble(
        effect = compound_name,
        type = "nested",
        n = info$mean_per_group
      )
    }
  }
  
  if (length(rows) == 0) {
    return(tibble::tibble(
      effect = character(),
      type = character(),
      n = numeric()
    ))
  }
  
  dplyr::bind_rows(rows)
}

#' Format Nested Details as Footnote
#'
#' Creates formatted footnote text for nested effects, showing group counts,
#' means, and harmonic means when applicable.
#'
#' @param nested_info List of nested effect information from calculate_nested_sample_sizes().
#'
#' @return Character vector of formatted footnote lines, one per nested effect.
#'   Returns character(0) if nested_info is NULL or empty.
#'
#' @keywords internal
format_nested_footnote <- function(nested_info) {
  if (is.null(nested_info) || length(nested_info) == 0) {
    return(character(0))
  }
  
  lines <- character(0)
  
  for (key in names(nested_info)) {
    info <- nested_info[[key]]
    
    nested_facet <- info$nested_facet
    nesting_facet <- info$nesting_facet
    n_groups <- info$n_groups
    mean_per_group <- info$mean_per_group
    harmonic_mean <- info$harmonic_mean_per_group
    
    # Capitalize first letter of facet names for display
    nested_cap <- paste0(toupper(substring(nested_facet, 1, 1)), 
                         substring(nested_facet, 2))
    nesting_cap <- paste0(toupper(substring(nesting_facet, 1, 1)), 
                          substring(nesting_facet, 2))
    
    # Build the main text
    main_text <- sprintf("%s per %s: %d %ss",
                         nested_cap,
                         nesting_cap,
                         n_groups,
                         nesting_facet)
    
    # Add mean info
    mean_text <- sprintf(", %.1f %ss per %s",
                         mean_per_group,
                         nested_facet,
                         nesting_facet)
    
    # Add harmonic mean if unbalanced
    if (abs(mean_per_group - harmonic_mean) > 0.01) {
      harmonic_text <- sprintf(" (harmonic: %.1f)", harmonic_mean)
      line <- paste0(main_text, mean_text, harmonic_text)
    } else {
      line <- paste0(main_text, mean_text)
    }
    
    lines <- c(lines, line)
  }
  
  lines
}
