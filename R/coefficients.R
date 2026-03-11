#' Generalizability and Dependability Coefficients
#'
#' Functions for computing G (generalizability) and D (dependability) coefficients
#' from variance components in generalizability theory.
#'
#' @name coefficients
#' @keywords internal
#' @importFrom magrittr %>%
NULL

#' Parse Facet Names from Component String
#'
#' Extracts facet names from a variance component string.
#' Handles interactions (e.g., "p:r:s"), main effects, and residual.
#'
#' @param component Character string naming the variance component.
#' @return A character vector of facet names.
#'
#' @keywords internal
parse_component_facets <- function(component) {
  if (identical(component, "Residual")) {
    return("Residual")
  }
  strsplit(component, ":")[[1]]
}

#' Check if Component is an Interaction
#'
#' Determines whether a variance component represents an interaction
#' (contains ":" in the name).
#'
#' @param component Character string naming the variance component.
#' @return TRUE if the component is an interaction, FALSE otherwise.
#'
#' @keywords internal
is_interaction <- function(component) {
grepl(":", component)
}

#' Get Universe Scale Factor
#'
#' Calculates the scale factor for a non-object universe component.
#' The scale factor is the product of sample sizes for all facets in the component
#' except the object facet (which is not scaled).
#'
#' @param component Character string naming the variance component (e.g., "p:o").
#' @param n Named list of sample sizes for each facet.
#' @param agg_facets Character vector of aggregation facets (optional).
#' @return The scale factor (product of relevant sample sizes).
#'
#' @keywords internal
get_universe_scale_factor <- function(component, n, agg_facets = NULL) {
if (is.null(n) || length(n) == 0) {
return(1)
}

# Extract facets from component (e.g., "p:o" -> c("p", "o"))
comp_facets <- parse_component_facets(component)

# Build scale factor from sample sizes of all facets in the component
scale_factor <- 1
for (facet in comp_facets) {
if (facet %in% names(n)) {
scale_factor <- scale_factor * n[[facet]]
}
}

scale_factor
}

#' Get Universe Scale Factor for a Component
#'
#' Calculates the scale factor for a non-object universe component,
#' handling the special case where the component might be mapped to Residual.
#'
#' @param comp User-specified component name (e.g., "person:item").
#' @param n Named list of sample sizes for each facet.
#' @param agg_facets Character vector of aggregation facets (optional).
#' @param residual_is Character string specifying what the residual represents.
#' @param actual_comp The actual component name in variance components (e.g., "Residual").
#' @return The scale factor (product of relevant sample sizes).
#'
#' @keywords internal
get_universe_scale_factor_for_component <- function(comp, n, agg_facets = NULL, residual_is = NULL, actual_comp = NULL) {
if (is.null(n) || length(n) == 0) {
return(1)
}

# Determine which component to use for scale factor calculation
component_for_scale <- if (!is.null(actual_comp)) actual_comp else comp

# Special handling for Residual
if (component_for_scale == "Residual") {
# Use residual_is to determine the facets
if (!is.null(residual_is) && residual_is != "") {
comp_facets <- parse_component_facets(residual_is)
} else {
# Can't determine facets, return 1
return(1)
}
} else {
comp_facets <- parse_component_facets(component_for_scale)
}

# Build scale factor from sample sizes of all facets in the component
scale_factor <- 1
for (facet in comp_facets) {
if (facet %in% names(n)) {
scale_factor <- scale_factor * n[[facet]]
}
}

scale_factor
}

#' Check if Component is the Object of Measurement
#'
#' Determines whether a variance component represents the main effect
#' of the object of measurement (not an interaction).
#'
#' @param component Character string naming the variance component.
#' @param object Character string naming the object of measurement.
#' @return TRUE if the component is the object main effect.
#'
#' @keywords internal
is_object_effect <- function(component, object) {
  !grepl(":", component) && identical(component, object)
}

#' Check if Object is in Component
#'
#' Determines whether the object of measurement appears in a variance component
#' (as main effect or part of an interaction).
#'
#' @param component Character string naming the variance component.
#' @param object Character string naming the object of measurement.
#' @return TRUE if the object appears in the component.
#'
#' @keywords internal
object_in_component <- function(component, object) {
  if (identical(component, "Residual")) {
    return(FALSE)
  }
  facets <- parse_component_facets(component)
  object %in% facets
}

#' Compute Scale Factor for D-Study Variance Component
#'
#' Calculates the scaling factor for a variance component based on
#' the D-study sample sizes. For main effects, divides by n_facet.
#' For interactions, divides by the product of n for each facet.
#'
#' @param facets Character vector of facet names in the component.
#' @param n Named list of sample sizes for each facet.
#' @return The scale factor (numeric).
#'
#' @keywords internal
compute_scale_factor <- function(facets, n) {
  if (length(facets) == 0 || identical(facets, "Residual")) {
    return(1)
  }
  
  scale_factor <- 1
  for (facet in facets) {
    if (facet %in% names(n)) {
      scale_factor <- scale_factor * n[[facet]]
    }
  }
  scale_factor
}

#' Compute Aggregation Factor for Variance Component
#'
#' Calculates the aggregation scaling factor for a variance component.
#' When aggregating over a facet, components containing that facet are
#' divided by the sample size for that facet.
#'
#' @param component Character string naming the variance component.
#' @param aggregation_facets Character vector of facets to aggregate over.
#' @param n Named list of sample sizes for each facet.
#' @return The aggregation scale factor (numeric).
#'
#' @keywords internal
compute_aggregation_factor <- function(component, aggregation_facets, n) {
  if (is.null(aggregation_facets) || length(aggregation_facets) == 0) {
    return(1)
  }
  
  if (identical(component, "Residual")) {
    return(1)
  }
  
  component_facets <- parse_component_facets(component)
  intersect_facets <- component_facets[component_facets %in% aggregation_facets]
  
  if (length(intersect_facets) == 0) {
    return(1)
  }
  
  factor <- 1
  for (facet in intersect_facets) {
    if (facet %in% names(n)) {
      factor <- factor * n[[facet]]
    }
  }
  factor
}

#' Compute Residual Aggregation Factor
#'
#' Calculates the aggregation scaling factor for the residual component.
#' The residual is scaled based on the intersection of aggregation facets
#' and the facets specified in residual_is.
#'
#' @param residual_is Character string specifying which facets make up the residual.
#' @param aggregation_facets Character vector of facets to aggregate over.
#' @param n Named list of sample sizes for each facet.
#' @return The aggregation scale factor for residual (numeric).
#'
#' @keywords internal
compute_residual_aggregation_factor <- function(residual_is, aggregation_facets, n) {
  if (is.null(residual_is) || identical(residual_is, "")) {
    return(1)
  }
  
  if (is.null(aggregation_facets) || length(aggregation_facets) == 0) {
    return(1)
  }
  
  residual_facets <- parse_component_facets(residual_is)
  intersect_facets <- aggregation_facets[aggregation_facets %in% residual_facets]
  
  if (length(intersect_facets) == 0) {
    return(1)
  }
  
  factor <- 1
  for (facet in intersect_facets) {
    if (facet %in% names(n)) {
      factor <- factor * n[[facet]]
    }
  }
  factor
}

#' Parse Specification for Object or Error Components
#'
#' Parses a specification for object of measurement or error components.
#' Accepts:
#' - A single character string: "p"
#' - A character vector: c("p", "p:d")
#' - A formula: obj ~ p + p:d (LHS is ignored for error spec, or can be used for object name)
#' - A one-sided formula: ~ p + p:d (extracts variables from RHS)
#'
#' @param spec Specification (character, character vector, or formula).
#' @return A character vector of component names.
#'
#' @keywords internal
extract_facet_names <- function(spec) {
  if (is.null(spec) || length(spec) == 0) return(character(0))
  unique(unlist(strsplit(spec, ":")))
}

parse_specification <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  
  if (inherits(x, "formula")) {
    if (length(x) == 2) {
      rhs <- x[[2]]
    } else {
      rhs <- x[[3]]
    }
    
    components <- character(0)
    
    if (is.call(rhs)) {
      parse_formula_terms <- function(expr) {
        result <- character(0)
        if (is.call(expr)) {
          if (expr[[1]] == as.symbol("+")) {
            result <- c(result, parse_formula_terms(expr[[2]]))
            result <- c(result, parse_formula_terms(expr[[3]]))
          } else if (expr[[1]] == as.symbol(":")) {
            result <- c(result, deparse(expr))
          } else {
            result <- c(result, parse_formula_terms(expr[[2]]))
          }
        } else if (is.symbol(expr)) {
          result <- c(result, as.character(expr))
        }
        result
      }
      components <- parse_formula_terms(rhs)
    } else if (is.symbol(rhs)) {
      components <- as.character(rhs)
    }
    
    components <- components[components != "."]
    return(components)
  }
  
  if (is.character(x)) {
    if (length(x) == 1) {
      return(x)
    } else {
      return(x)
    }
  }
  
  stop("Specification must be a character string, character vector, or formula",
       call. = FALSE)
}

#' Check if Component is Part of Object Specification
#'
#' Determines whether a variance component is part of the object of measurement
#' specification (which can include multiple components).
#'
#' @param component Character string naming the variance component.
#' @param object_spec Character vector of object component names.
#' @return TRUE if the component is part of the object specification.
#'
#' @keywords internal
is_object_component <- function(component, object_spec) {
  component %in% object_spec
}

#' Check if Component is Part of Error Specification
#'
#' Determines whether a variance component is part of the error specification.
#'
#' @param component Character string naming the variance component.
#' @param error_spec Character vector of error component names, or NULL for default.
#' @param object_spec Character vector of object component names.
#' @param vc_all Variance components tibble (for default error calculation).
#' @return TRUE if the component is part of the error specification.
#'
#' @keywords internal
is_error_component <- function(component, error_spec, object_spec, vc_all) {
  if (!is.null(error_spec)) {
    return(component %in% error_spec)
  }
  
  !is_object_component(component, object_spec) && component != "Residual"
}

#' Calculate D-Study Variance Components
#'
#' Scales variance components according to D-study sample sizes.
#' For each facet, the variance is divided by the number of levels in the D-study.
#'
#' @param vc G-study variance components tibble.
#' @param n Named list of sample sizes for each facet.
#' @param object Character vector naming the object(s) of measurement.
#' @param aggregation Character vector of facets to aggregate over. When specified,
#' components containing aggregation facets are divided only by the aggregation
#' facet's n (not by the full scale factor).
#' @param n_provided Logical indicating if n was explicitly provided by the user.
#' When FALSE (default), no rescaling is applied.
#' @param residual_is Character string specifying which facets make up the residual.
#' @param facet_n Named vector of facet sample sizes from G-study. Used to fill
#' in missing aggregation facets when they are not specified in `n`.
#' @return A tibble of D-study variance components with columns:
#' \itemize{
#' \item component: Name of variance component
#' \item var_unscaled: Unscaled variance estimate (from G-study)
#' \item pct_unscaled: Percentage of total unscaled variance
#' \item var_scaled: Scaled variance (divided by D-study sample sizes)
#' \item pct_scaled: Percentage of total scaled variance
#' }
#' When aggregation is specified, returns both scaled and unscaled estimates.
#'
#' @keywords internal
calculate_dstudy_variance <- function(vc, n, object, aggregation = NULL, n_provided = FALSE, residual_is = NULL, facet_n = NULL) {
  d_vc <- vc

  object_spec <- parse_specification(object)

  agg_facets <- if (!is.null(aggregation)) parse_specification(aggregation) else NULL
  residual_facets <- if (!is.null(residual_is) && residual_is != "") parse_component_facets(residual_is) else NULL

  has_aggregation <- !is.null(aggregation) && length(aggregation) > 0

  # Store original n for scale_factor computation
  # scale_factor should only include facets that were originally in n (user-specified)
  n_original <- n

  # Fill in missing aggregation facets from facet_n (G-study sample sizes)
  if (has_aggregation && !is.null(facet_n)) {
    for (f in agg_facets) {
      if (!(f %in% names(n)) && f %in% names(facet_n)) {
        n[[f]] <- facet_n[[f]]
      }
    }
  }

  n_facets <- if (n_provided && length(n) > 0) names(n) else character(0)
  n_original_facets <- if (n_provided && length(n_original) > 0) names(n_original) else character(0)

  if (has_aggregation) {
    effective_agg_facets <- unique(c(agg_facets, n_facets))
    additional_agg_facets <- setdiff(agg_facets, n_original_facets)
  } else {
    effective_agg_facets <- n_facets
    additional_agg_facets <- character(0)
  }

  d_vc <- d_vc %>%
    mutate(
      facets_list = purrr::map(component, function(c) {
        if (c == "Residual" && !is.null(residual_facets)) {
          return(residual_facets)
        }
        parse_component_facets(c)
      }),
      is_object = component %in% object_spec,
      is_residual = (component == "Residual"),
      is_interaction = purrr::map_lgl(facets_list, ~ length(.x) > 1),
      scale_factor = purrr::map_dbl(facets_list, compute_scale_factor, n = n_original),
      n_facet = purrr::map_dbl(facets_list, function(flist) {
        if (length(flist) == 1 && flist %in% names(n)) {
          n[[flist]]
        } else {
          1
        }
      }),
      has_additional_agg = purrr::map_lgl(component, function(c) {
        if (length(additional_agg_facets) == 0) return(FALSE)
        if (c == "Residual" && !is.null(residual_facets)) {
          return(any(additional_agg_facets %in% residual_facets))
        }
        facets <- parse_component_facets(c)
        any(facets %in% additional_agg_facets)
      }),
      additional_agg_factor = purrr::map_dbl(component, function(c) {
        if (length(additional_agg_facets) == 0) return(1)
        if (c == "Residual" && !is.null(residual_facets)) {
          agg_in_comp <- residual_facets[residual_facets %in% additional_agg_facets]
        } else {
          facets <- parse_component_facets(c)
          agg_in_comp <- facets[facets %in% additional_agg_facets]
        }
        if (length(agg_in_comp) > 0) {
          agg_factor <- 1
          for (f in agg_in_comp) {
            if (f %in% names(n)) agg_factor <- agg_factor * n[[f]]
          }
          return(agg_factor)
        }
        return(1)
      }),
      has_agg = purrr::map_lgl(component, function(c) {
        if (is.null(effective_agg_facets) || length(effective_agg_facets) == 0) return(FALSE)
        if (c == "Residual" && !is.null(residual_facets)) {
          return(any(effective_agg_facets %in% residual_facets))
        }
        facets <- parse_component_facets(c)
        any(facets %in% effective_agg_facets)
      }),
      agg_n = purrr::map_dbl(component, function(c) {
        if (is.null(effective_agg_facets) || length(effective_agg_facets) == 0) return(1)
        if (c == "Residual" && !is.null(residual_facets)) {
          agg_in_comp <- residual_facets[residual_facets %in% effective_agg_facets]
        } else {
          facets <- parse_component_facets(c)
          agg_in_comp <- facets[facets %in% effective_agg_facets]
        }
        if (length(agg_in_comp) > 0) {
          agg_factor <- 1
          for (f in agg_in_comp) {
            if (f %in% names(n)) agg_factor <- agg_factor * n[[f]]
          }
          return(agg_factor)
        }
        return(1)
      })
    )

 do_rescaling <- !is.null(aggregation) || (n_provided && length(n) > 0)

 if (!do_rescaling) {
 d_vc <- d_vc %>%
 mutate(
 var_unscaled = .data$var,
 var_scaled = .data$var,
 pct_unscaled = (.data$var / sum(.data$var, na.rm = TRUE)) * 100,
 pct_scaled = (.data$var / sum(.data$var, na.rm = TRUE)) * 100
 ) %>%
 select(-any_of(c('var', 'pct', 'facets_list', 'is_object', 'is_residual', 'is_interaction',
 'scale_factor', 'n_facet', 'has_additional_agg', 'additional_agg_factor',
 'has_agg', 'agg_n', 'error', 'se', 'lower', 'upper', 'sd', 'Rhat', 'Bulk_ESS', 'Tail_ESS')))
 } else {
 total_n <- if (length(n) > 0) prod(unlist(n)) else 1

 if (has_aggregation) {
 d_vc <- d_vc %>%
 mutate(
 var_unscaled = .data$var,
 var = case_when(
 is_object ~ .data$var,
 has_additional_agg ~ .data$var / (scale_factor * additional_agg_factor),
 is_residual & has_agg ~ .data$var / (scale_factor * additional_agg_factor),
 is_residual ~ .data$var / total_n,
 has_agg & !is_interaction ~ .data$var / (scale_factor * additional_agg_factor),
 has_agg & is_interaction ~ .data$var / (scale_factor * additional_agg_factor),
 is_interaction ~ .data$var / scale_factor,
n_facet > 1 ~ .data$var / n_facet,
			TRUE ~ .data$var
		),
		var_scaled = .data$var,
		pct_scaled = (.data$var / sum(.data$var, na.rm = TRUE)) * 100,
		pct_unscaled = (var_unscaled / sum(var_unscaled, na.rm = TRUE)) * 100
	) %>%
		select(-any_of(c('pct', 'facets_list', 'is_object', 'is_residual', 'is_interaction',
			'scale_factor', 'n_facet', 'has_additional_agg', 'additional_agg_factor',
			'has_agg', 'agg_n', 'error', 'se', 'lower', 'upper', 'sd', 'Rhat', 'Bulk_ESS', 'Tail_ESS')))
	} else {
 d_vc <- d_vc %>%
 mutate(
 var_unscaled = .data$var,
 var = case_when(
is_object ~ .data$var,
		has_agg ~ .data$var / scale_factor,
		is_residual ~ .data$var / total_n,
		is_interaction ~ .data$var / scale_factor,
		n_facet > 1 ~ .data$var / n_facet,
		TRUE ~ .data$var
	),
	var_scaled = .data$var,
	pct_unscaled = (var_unscaled / sum(var_unscaled, na.rm = TRUE)) * 100,
	pct_scaled = (.data$var / sum(.data$var, na.rm = TRUE)) * 100
) %>%
	select(-any_of(c('pct', 'facets_list', 'is_object', 'is_residual', 'is_interaction',
		'scale_factor', 'n_facet', 'has_additional_agg', 'additional_agg_factor',
		'has_agg', 'agg_n', 'error', 'se', 'lower', 'upper', 'sd', 'Rhat', 'Bulk_ESS', 'Tail_ESS')))
}
}

 d_vc
}

#' Calculate Divided Variance Components
#'
#' Divides each variance component by the sample sizes of non-object facets
#' that appear in that component. This provides an alternative scaling where:
#' - Object component is NOT scaled (same as standard D-study)
#' - For non-object components, divides only by n for non-object facets in component
#' - For Residual, divides by product of all non-object facet sample sizes
#'
#' This differs from standard D-study variance scaling in how interactions are handled.
#' Standard D-study scales interactions by all facets in the component, while this
#' method scales only by the non-object facets.
#'
#' @param vc Variance components tibble (unscaled, from G-study)
#' @param n Named list of sample sizes for each facet
#' @param object Specification for object of measurement
#' @param residual_is Character string specifying which facets make up the residual
#' @return Tibble with additional 'var_divided' column
#'
#' @keywords internal
calculate_divided_variance <- function(vc, n, object, residual_is = NULL) {
  if (is.null(n) || length(n) == 0) {
    d_vc <- vc %>%
      mutate(var_divided = var)
    return(d_vc)
  }

  object_spec <- parse_specification(object)

  all_facets <- names(n)
  non_object_facets <- setdiff(all_facets, object_spec)

  residual_facets <- if (!is.null(residual_is) && residual_is != "") {
    parse_component_facets(residual_is)
  } else {
    NULL
  }

  d_vc <- vc %>%
    mutate(
      comp_facets = purrr::map(component, function(c) {
        if (c == "Residual") {
          if (!is.null(residual_facets)) {
            residual_facets
          } else {
            non_object_facets
          }
        } else {
          parse_component_facets(c)
        }
      }),
      is_object = component %in% object_spec,
      divisor = purrr::map_dbl(comp_facets, function(facets) {
        if (length(facets) == 0) return(1)

        facets_to_use <- facets[facets %in% non_object_facets]
        if (length(facets_to_use) == 0) return(1)

        div <- 1
        for (f in facets_to_use) {
          if (f %in% names(n)) {
            div <- div * n[[f]]
          }
        }
        div
      }),
      var_divided = ifelse(is_object, var, var / divisor)
    ) %>%
    select(-comp_facets, -is_object, -divisor)

  d_vc
}

#' Calculate Coefficients Using Divided Variance Components
#'
#' Computes coefficients using the "divided" variance estimates where each
#' component is divided only by the sample sizes of non-object facets.
#'
#' @param vc_divided Variance components tibble with 'var_divided' column
#' @param object Specification for object of measurement
#' @param error Specification for error components (optional)
#' @param aggregation Aggregation facets (optional)
#' @param residual_is Residual composition specification
#' @param universe Universe specification (optional)
#' @return Data frame with coefficient estimates
#'
#' @keywords internal
calculate_divided_coefficients <- function(vc_divided, object, error = NULL,
aggregation = NULL, residual_is = NULL, universe = NULL) {
vc_for_calc <- vc_divided %>%
mutate(var = var_divided) %>%
select(-var_divided)

calculate_coefficients(vc_for_calc, n = NULL, object, error, aggregation, residual_is, universe)
}

#' Calculate G and D Coefficients
#'
#' Computes generalizability (G) and dependability (D/Phi) coefficients
#' from D-study variance components.
#'
#' @param vc D-study variance components tibble.
#' @param n Named list of sample sizes for each facet (required for aggregation).
#' @param object Specification for object of measurement. Can be:
#'   - NULL (default): uses first component in vc
#'   - A character string: "p"
#'   - A character vector: c("p", "p:d")
#'   - A formula: obj ~ p + p:d
#' @param error Specification for error components. Can be:
#'   - NULL (default): all components not in object specification
#'   - A character string: "p:i"
#'   - A character vector: c("p:i", "p:r")
#'   - A formula: ~ p:i + p:r
#'   Note: If the same facet is specified in both `object` and `error`, an error is raised.
#' @param aggregation Character vector of facets to aggregate over.
#'   Components containing these facets will be divided by the sample size.
#'   Default is NULL (no aggregation).
#'   Note: If the same facet is specified in both `object` and `aggregation`, an error is raised.
#' @param residual_is Character string specifying which facets make up the residual.
#'   For example, "p:i:r" means the residual represents the p:i:r interaction.
#'   Default is NULL (residual is not rescaled).
#' @return A data frame with columns:
#' \item{uni}{Universe score variance. Represents the true individual differences
#' in the trait being measured. It reflects the systematic variance attributable
#' to the object of measurement (e.g., persons) after accounting for measurement
#' error. In classical test theory terms, this corresponds to the "true score"
#' variance. A larger universe score variance indicates greater heterogeneity
#' among objects of measurement, which generally improves the precision of
#' relative comparisons.}
#' \item{sigma2_delta}{Relative error variance. Error variance affecting relative
#' decisions (rankings). Includes only variance components that interact with the
#' object of measurement (i.e., components that contain the object facet). For
#' example, in a person x rater design with person as object, this includes
#' person:rater and residual components.}
#' \item{sigma2_delta_abs}{Absolute error variance. Error variance affecting absolute
#' decisions. Includes ALL variance components except the object of measurement.
#' For example, in a person x rater design with person as object, this includes
#' rater, person:rater, and residual components.}
#' \item{g}{Generalizability coefficient. Quantifies the reliability of
#' measurements for relative decisions—situations where the focus is on
#' rank-ordering or comparing objects of measurement (e.g., identifying high
#' vs. low performers). Calculated as: g = uni / (uni + sigma2_delta).
#' The G coefficient ranges from 0 to 1, with higher values indicating better
#' reliability. Values >= 0.70 are typically considered acceptable for group-level
#' decisions, while >= 0.80 is preferred for individual-level decisions.}
#' \item{phi}{Dependability coefficient (Phi coefficient). Quantifies the
#' reliability of measurements for absolute decisions—situations where the
#' focus is on the absolute level of a measurement (e.g., determining if a
#' score exceeds a cutoff). Calculated as: phi = uni / (uni + sigma2_delta_abs).
#' The Phi coefficient is always <= G coefficient because absolute error includes
#' all error sources, not just those affecting relative rankings.}
#' \item{sem_rel}{Standard error of measurement for relative decisions.
#' Calculated as sqrt(sigma2_delta). Indicates the typical error when
#' interpreting relative standing or rank-order among objects of measurement.}
#' \item{sem_abs}{Standard error of measurement for absolute decisions.
#' Calculated as sqrt(sigma2_delta_abs). Indicates the typical error when
#' interpreting absolute score levels against criterion-referenced standards.}
#'
#' @keywords internal
calculate_coefficients <- function(vc, n = NULL, object = NULL, error = NULL,
aggregation = NULL, residual_is = NULL, universe = NULL) {
if (is.null(object) && is.null(error)) {
default_object <- vc$component[1]
warning(
"No object specified. Using first component '", default_object,
"' as the object of measurement.",
call. = FALSE
)
object_spec <- default_object
} else if (is.null(object)) {
object_spec <- character(0)
} else {
object_spec <- parse_specification(object)
}

# Process universe specification
# Default: universe = object only
if (is.null(universe) || length(universe) == 0) {
universe_spec <- object_spec
} else {
universe_spec <- parse_specification(universe)
}

if (!is.null(error)) {
error_spec <- parse_specification(error)
} else {
error_spec <- NULL
}

if (!is.null(aggregation)) {
agg_facets <- parse_specification(aggregation)
} else {
agg_facets <- NULL
}

if (!is.null(error_spec) && length(error_spec) > 0) {
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

has_dim <- "dim" %in% names(vc)
dims <- if (has_dim) unique(vc$dim) else "response"

# Check if we have unscaled variance (from calculate_dstudy_variance with aggregation)
has_var_unscaled <- "var_unscaled" %in% names(vc)

results <- lapply(dims, function(d) {
if (has_dim) {
vc_dim <- vc[vc$dim == d, , drop = FALSE]
} else {
vc_dim <- vc
}

vc_modified <- vc_dim

# Calculate universe score variance
# Use unscaled variance if available for proper scaling
uni <- 0
for (comp in universe_spec) {
# Check if component exists (might be named differently)
actual_comp <- comp
if (!(comp %in% vc_modified$component)) {
# Check if this is a residual specification
if (!is.null(residual_is) && comp == residual_is && "Residual" %in% vc_modified$component) {
actual_comp <- "Residual"
}
}

if (actual_comp %in% vc_modified$component) {
# Use unscaled variance if available, otherwise use var
if (has_var_unscaled) {
comp_var <- vc_modified$var_unscaled[vc_modified$component == actual_comp]
} else {
comp_var <- vc_modified$var[vc_modified$component == actual_comp]
}

if (comp %in% object_spec) {
# Object component: NOT scaled
uni <- uni + comp_var
} else {
# Non-object universe component: scaled by sample sizes
# Need to determine scale factor based on component facets
scale_factor <- get_universe_scale_factor_for_component(comp, n, agg_facets, residual_is, actual_comp)
uni <- uni + comp_var / scale_factor
}
}
}

# For relative error, exclude universe components (not just object)
sigma2_delta <- compute_relative_error(vc_modified, universe_spec, error_spec, agg_facets, residual_is)

# For absolute error, exclude universe components (not just object)
sigma2_delta_abs <- compute_absolute_error(vc_modified, universe_spec, error_spec, agg_facets, residual_is)

g <- uni / (uni + sigma2_delta)

phi <- uni / (uni + sigma2_delta_abs)

sem_rel <- sqrt(sigma2_delta)
sem_abs <- sqrt(sigma2_delta_abs)

result_df <- data.frame(
uni = uni,
sigma2_delta = sigma2_delta,
sigma2_delta_abs = sigma2_delta_abs,
g = g,
phi = phi,
sem_rel = sem_rel,
sem_abs = sem_abs
)

if (has_dim) {
result_df$dim <- d
}

result_df
})

do.call(rbind, results)
}

#' Compute Relative Error Variance
#'
#' Calculates the relative error variance for a D-study design.
#' Relative error includes variance components that affect relative rankings
#' (interactions with the universe components).
#'
#' @param vc Variance components tibble.
#' @param universe_spec Character vector of universe component names.
#' @param error_spec Character vector of error component names, or NULL for default.
#' @return The relative error variance.
#'
#' @keywords internal
compute_relative_error <- function(vc, universe_spec, error_spec = NULL, agg_facets = NULL, residual_is = NULL) {
if (!is.null(error_spec)) {
vc_filtered <- vc %>%
filter(component %in% error_spec)

sum(vc_filtered$var, na.rm = TRUE)
} else {
residual_facets <- if (!is.null(residual_is)) parse_component_facets(residual_is) else NULL

vc_filtered <- vc %>%
mutate(
is_universe = component %in% universe_spec,
is_residual = (component == "Residual"),
is_interaction = purrr::map_lgl(component, is_interaction),
has_universe = purrr::map_lgl(component, function(c) {
facets <- parse_component_facets(c)
any(universe_spec %in% facets)
}),
is_agg_main = if (!is.null(agg_facets)) {
purrr::map_lgl(component, function(c) {
if (c == "Residual") return(FALSE)
facets <- parse_component_facets(c)
length(facets) == 1 && facets %in% agg_facets
})
} else {
FALSE
}
) %>%
filter(
!is_universe,
!is_agg_main,
is_residual | has_universe
)

sum(vc_filtered$var, na.rm = TRUE)
}
}

#' Compute Absolute Error Variance
#'
#' Calculates the absolute error variance for a D-study design.
#' Absolute error includes all variance components that affect absolute scores
#' (all components except the universe components).
#'
#' @param vc Variance components tibble.
#' @param universe_spec Character vector of universe component names.
#' @param error_spec Character vector of error component names, or NULL for default.
#' @return The absolute error variance.
#'
#' @keywords internal
compute_absolute_error <- function(vc, universe_spec, error_spec = NULL, agg_facets = NULL, residual_is = NULL) {
if (!is.null(error_spec)) {
vc_filtered <- vc %>%
filter(component %in% error_spec)

sum(vc_filtered$var, na.rm = TRUE)
} else {
vc_filtered <- vc %>%
mutate(
is_agg_main = if (!is.null(agg_facets)) {
purrr::map_lgl(component, function(c) {
if (c == "Residual") return(FALSE)
facets <- parse_component_facets(c)
length(facets) == 1 && facets %in% agg_facets
})
} else {
FALSE
}
) %>%
filter(
!(component %in% universe_spec),
!is_agg_main
) %>%
select(-is_agg_main)

sum(vc_filtered$var, na.rm = TRUE)
}
}

#' Compute Generalizability Coefficient
#'
#' Calculates the generalizability coefficient (G coefficient), which represents
#' the proportion of observed variance attributable to universe score variance
#' for relative decisions.
#'
#' @param variance_components A tibble of variance components.
#' @param n A named list specifying the number of levels for each facet.
#' @param object Specification for object of measurement (character, vector, or formula).
#' @param error Specification for error components, or NULL for default.
#' @param aggregation Character vector of facets to aggregate over.
#' @param residual_is Character string specifying which facets make up the residual.
#' @param n_provided Logical indicating if n was explicitly provided.
#' @return The G coefficient as a numeric value.
#'
#' @keywords internal
compute_g_coefficient <- function(variance_components, n, object, error = NULL,
                                  aggregation = NULL, residual_is = NULL, n_provided = FALSE) {
  d_vc <- calculate_dstudy_variance(variance_components, n, object, aggregation, n_provided, residual_is)

  coefs <- calculate_coefficients(d_vc, n, object, error, aggregation, residual_is)

  coefs$g
}

#' Compute Dependability Coefficient
#'
#' Calculates the dependability coefficient (Phi coefficient), which represents
#' the proportion of observed variance attributable to universe score variance
#' for absolute decisions.
#'
#' @param variance_components A tibble of variance components.
#' @param n A named list specifying the number of levels for each facet.
#' @param object Specification for object of measurement (character, vector, or formula).
#' @param error Specification for error components, or NULL for default.
#' @param aggregation Character vector of facets to aggregate over.
#' @param residual_is Character string specifying which facets make up the residual.
#' @param n_provided Logical indicating if n was explicitly provided.
#' @return The Phi coefficient as a numeric value.
#'
#' @keywords internal
compute_phi_coefficient <- function(variance_components, n, object, error = NULL,
                                    aggregation = NULL, residual_is = NULL, n_provided = FALSE) {
  d_vc <- calculate_dstudy_variance(variance_components, n, object, aggregation, n_provided, residual_is)

  coefs <- calculate_coefficients(d_vc, n, object, error, aggregation, residual_is)

  coefs$phi
}

#' Compute Standard Error of Measurement
#'
#' Calculates the standard error of measurement (SEM) for a D-study design.
#'
#' @param error_variance The error variance (relative or absolute).
#' @return The standard error of measurement.
#'
#' @keywords internal
compute_sem <- function(error_variance) {
  sqrt(error_variance)
}

#' Calculate Coefficients Using Posterior Draws
#'
#' Computes generalizability and dependability coefficients using full posterior
#' distributions from a brms model. This provides proper uncertainty quantification
#' for D-study coefficients.
#'
#' @param gstudy_obj A gstudy object fitted with backend = "brms".
#' @param n Named list of sample sizes for each facet.
#' @param object Specification for object of measurement.
#' @param error Specification for error components.
#'   Note: If the same facet is specified in both `object` and `error`, an error is raised.
#' @param aggregation Character vector of facets to aggregate over.
#'   Note: If the same facet is specified in both `object` and `aggregation`, an error is raised.
#' @param residual_is Character string specifying residual composition.
#' @param is_sweep Logical indicating if this is a sweep over multiple sample sizes.
#' @param n_grid If is_sweep, a data frame with all combinations of sample sizes.
#'
#' @return A list with:
#'   \item{coefficients}{Data frame with coefficient means (and sweep combinations if applicable)}
#'   \item{posterior}{List of posterior distribution vectors (or list of lists for sweep)}
#'
#' @keywords internal
calculate_coefficients_posterior <- function(gstudy_obj, n, object = NULL, universe = NULL, error = NULL,
aggregation = NULL, residual_is = NULL,
is_sweep = FALSE, n_grid = NULL,
n_provided = FALSE, use_scaled = TRUE) {
if (!requireNamespace("brms", quietly = TRUE)) {
stop("Package 'brms' is required for posterior estimation.", call. = FALSE)
}

# Extract posterior draws
draws <- brms::as_draws_matrix(gstudy_obj$model)

# Get variance component names and their draws
vc_draws <- extract_variance_draws(gstudy_obj, draws)

# Get object of measurement (always from G-study)
if (is.null(object)) {
object <- gstudy_obj$object
}
object_spec <- parse_specification(object)

# Process universe specification
# Default: universe = object only
if (is.null(universe) || length(universe) == 0) {
universe_spec <- object_spec
} else {
universe_spec <- parse_specification(universe)
}

# Parse error specification
if (!is.null(error)) {
error_spec <- parse_specification(error)
} else {
error_spec <- NULL
}

# Parse aggregation specification
if (!is.null(aggregation)) {
agg_facets <- parse_specification(aggregation)
} else {
agg_facets <- NULL
}

# Validate universe vs error
if (!is.null(error_spec) && length(error_spec) > 0) {
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

  is_mv <- gstudy_obj$is_multivariate
  mv_dims <- if (is_mv) gstudy_obj$dimensions else character(0)

  if (is_sweep) {
    # Calculate for each sweep combination
    results <- vector("list", nrow(n_grid))
    posterior_list <- vector("list", nrow(n_grid))

    for (i in seq_len(nrow(n_grid))) {
      n_current <- as.list(n_grid[i, , drop = FALSE])
      n_current <- lapply(n_current, as.numeric)

      if (is_mv) {
        # Multivariate: iterate over dimensions
        dim_results <- vector("list", length(mv_dims))
        names(dim_results) <- mv_dims
        dim_posteriors <- vector("list", length(mv_dims))
        names(dim_posteriors) <- mv_dims

for (d in mv_dims) {
combo_result <- calculate_single_posterior(
vc_draws = vc_draws[[d]],
n = n_current,
object_spec = object_spec,
universe_spec = universe_spec,
error_spec = error_spec,
agg_facets = agg_facets,
residual_is = residual_is,
gstudy_obj = gstudy_obj,
n_provided = n_provided,
use_scaled = use_scaled
)
dim_results[[d]] <- combo_result$summary
dim_posteriors[[d]] <- combo_result$distributions
}

        # Combine results with dimension column
        results[[i]] <- cbind(
          data.frame(n_current, stringsAsFactors = FALSE),
          do.call(rbind, lapply(names(dim_results), function(d) {
            cbind(dim_results[[d]], dim = d)
          }))
        )
        posterior_list[[i]] <- dim_posteriors
} else {
# Univariate: single dimension
combo_result <- calculate_single_posterior(
vc_draws = vc_draws,
n = n_current,
object_spec = object_spec,
universe_spec = universe_spec,
error_spec = error_spec,
agg_facets = agg_facets,
residual_is = residual_is,
gstudy_obj = gstudy_obj,
n_provided = n_provided,
use_scaled = use_scaled
)

        results[[i]] <- cbind(
          data.frame(n_current, stringsAsFactors = FALSE),
          combo_result$summary
        )
        posterior_list[[i]] <- combo_result$distributions
      }
    }

    coefficients <- do.call(rbind, results)
    coefficients <- tibble::as_tibble(coefficients)

    return(list(
      coefficients = coefficients,
      posterior = posterior_list
    ))
  } else {
    # Single sample size
    if (is_mv) {
      dim_results <- vector("list", length(mv_dims))
      names(dim_results) <- mv_dims
      dim_posteriors <- vector("list", length(mv_dims))
      names(dim_posteriors) <- mv_dims

for (d in mv_dims) {
result <- calculate_single_posterior(
vc_draws = vc_draws[[d]],
n = n,
object_spec = object_spec,
universe_spec = universe_spec,
error_spec = error_spec,
agg_facets = agg_facets,
residual_is = residual_is,
gstudy_obj = gstudy_obj,
n_provided = n_provided,
use_scaled = use_scaled
)
dim_results[[d]] <- cbind(result$summary, dim = d)
dim_posteriors[[d]] <- result$distributions
}

      coefficients <- tibble::as_tibble(do.call(rbind, lapply(mv_dims, function(d) dim_results[[d]])))

      return(list(
        coefficients = coefficients,
        posterior = dim_posteriors
      ))
} else {
result <- calculate_single_posterior(
vc_draws = vc_draws,
n = n,
object_spec = object_spec,
universe_spec = universe_spec,
error_spec = error_spec,
agg_facets = agg_facets,
residual_is = residual_is,
gstudy_obj = gstudy_obj,
n_provided = n_provided,
use_scaled = use_scaled
)

      return(list(
        coefficients = tibble::as_tibble(result$summary),
        posterior = result$distributions
      ))
    }
  }
}

#' Extract Variance Component Draws from brms Model
#'
#' Extracts posterior draws of variance components (SD^2) from a brms model.
#' For multivariate models, returns a nested list structure by dimension.
#'
#' @param gstudy_obj A gstudy object with brms backend.
#' @param draws A posterior draws_matrix object (e.g., from brms::as_draws_matrix()).
#'
#' @return A named list. For univariate: list of variance draws per component.
#'   For multivariate: nested list with dimensions as names, each containing
#'   variance draws for that dimension.
#'
#' @keywords internal
extract_variance_draws <- function(gstudy_obj, draws) {
  model <- gstudy_obj$model
  vc <- gstudy_obj$variance_components
  is_mv <- gstudy_obj$is_multivariate
  
  if (is_mv) {
    dims <- unique(vc$dim)
    vc_draws <- vector("list", length(dims))
    names(vc_draws) <- dims
    
    for (d in dims) {
      vc_dim <- vc[vc$dim == d, ]
      vc_draws[[d]] <- list()
      
      for (comp in vc_dim$component) {
        if (comp == "Residual") {
          param_name <- paste0("sigma_", d)
        } else {
          clean_comp <- gsub(":", "__", comp)
          param_name <- paste0("sd_", clean_comp, "__", d, "_Intercept")
        }
        
        if (param_name %in% colnames(draws)) {
          sd_draws <- posterior::extract_variable(draws, param_name)
          vc_draws[[d]][[comp]] <- sd_draws^2
        } else {
          var_estimate <- vc_dim$var[vc_dim$component == comp]
          sd_draws <- rep(sqrt(var_estimate), nrow(draws))
          vc_draws[[d]][[comp]] <- sd_draws^2
        }
      }
    }
  } else {
    vc_draws <- list()
    
    for (comp in vc$component) {
      if (comp == "Residual") {
        param_name <- "sigma"
      } else {
        clean_comp <- gsub(":", "__", comp)
        param_name <- paste0("sd_", clean_comp, "__Intercept")
      }
      
      if (param_name %in% colnames(draws)) {
        sd_draws <- posterior::extract_variable(draws, param_name)
        vc_draws[[comp]] <- sd_draws^2
      } else {
        param_name_alt <- paste0("sd_", clean_comp)
        if (param_name_alt %in% colnames(draws)) {
          sd_draws <- posterior::extract_variable(draws, param_name_alt)
          vc_draws[[comp]] <- sd_draws^2
        } else {
          matching_params <- grep(paste0("^sd_", clean_comp), colnames(draws), value = TRUE)
          if (length(matching_params) > 0) {
            sd_draws <- posterior::extract_variable(draws, matching_params[1])
            vc_draws[[comp]] <- sd_draws^2
          } else {
            var_estimate <- vc$var[vc$component == comp]
            sd_draws <- rep(sqrt(var_estimate), nrow(draws))
            vc_draws[[comp]] <- sd_draws^2
          }
        }
      }
    }
  }

  vc_draws
}

#' Calculate Posterior Coefficients for a Single Sample Size
#'
#' @param vc_draws Named list of variance component draws.
#' @param n Named list of sample sizes.
#' @param object_spec Character vector of object components.
#' @param universe_spec Character vector of universe components.
#' @param error_spec Character vector of error components (or NULL).
#' @param agg_facets Character vector of aggregation facets (or NULL).
#' @param residual_is Character string for residual composition.
#' @param gstudy_obj Original gstudy object.
#' @param n_provided Logical indicating if n was explicitly provided.
#' @param use_scaled Logical indicating whether to use scaled variance draws.
#'   When FALSE, uses original variance draws and does not scale universe components.
#'   When TRUE, uses scaled variance draws and scales non-object universe components.
#'
#' @return List with 'summary' (data frame of means) and 'distributions' (list of vectors).
#'
#' @keywords internal
calculate_single_posterior <- function(vc_draws, n, object_spec, universe_spec,
error_spec, agg_facets, residual_is, gstudy_obj,
n_provided = FALSE, use_scaled = TRUE) {
n_draws <- length(vc_draws[[1]])

# Default universe to object if not specified
if (is.null(universe_spec) || length(universe_spec) == 0) {
universe_spec <- object_spec
}

# Scale variance components for D-study (only if use_scaled)
if (use_scaled) {
scaled_draws <- scale_variance_draws(vc_draws, n, object_spec, agg_facets,
residual_is, gstudy_obj, n_provided)
} else {
scaled_draws <- vc_draws
}

# Calculate universe score variance (uni)
# Object components are NEVER scaled, non-object universe components are scaled only if use_scaled
uni <- rep(0, n_draws)
for (comp in universe_spec) {
if (comp %in% names(scaled_draws)) {
if (comp %in% object_spec) {
# Object component: NEVER scaled
uni <- uni + scaled_draws[[comp]]
} else {
# Non-object universe component
if (use_scaled) {
# Scale by sample sizes
scale_factor <- get_universe_scale_factor_for_component(comp, n, agg_facets, residual_is, comp)
uni <- uni + scaled_draws[[comp]] / scale_factor
} else {
# No scaling for unscaled estimate
uni <- uni + scaled_draws[[comp]]
}
}
}
}

# Calculate relative error variance (sigma2_delta)
sigma2_delta <- compute_relative_error_draws(scaled_draws, universe_spec,
error_spec, agg_facets, residual_is)

# Calculate absolute error variance (sigma2_delta_abs)
sigma2_delta_abs <- compute_absolute_error_draws(scaled_draws, universe_spec,
error_spec, agg_facets, residual_is)

# Calculate coefficients
g <- uni / (uni + sigma2_delta)
phi <- uni / (uni + sigma2_delta_abs)
sem_rel <- sqrt(sigma2_delta)
sem_abs <- sqrt(sigma2_delta_abs)

# Handle infinite/NaN cases
g[is.nan(g) | is.infinite(g)] <- NA
phi[is.nan(phi) | is.infinite(phi)] <- NA

# Return summary (means) and distributions
summary_df <- data.frame(
uni = mean(uni, na.rm = TRUE),
sigma2_delta = mean(sigma2_delta, na.rm = TRUE),
sigma2_delta_abs = mean(sigma2_delta_abs, na.rm = TRUE),
g = mean(g, na.rm = TRUE),
phi = mean(phi, na.rm = TRUE),
sem_rel = mean(sem_rel, na.rm = TRUE),
sem_abs = mean(sem_abs, na.rm = TRUE)
)

distributions <- list(
uni = uni,
sigma2_delta = sigma2_delta,
sigma2_delta_abs = sigma2_delta_abs,
g = g,
phi = phi,
sem_rel = sem_rel,
sem_abs = sem_abs
)

list(
summary = summary_df,
distributions = distributions
)
}

#' Scale Variance Draws for D-Study
#'
#' Applies D-study scaling to posterior draws of variance components.
#'
#' @param vc_draws Named list of variance component draws.
#' @param n Named list of sample sizes.
#' @param object_spec Character vector of object components.
#' @param agg_facets Character vector of aggregation facets.
#' @param residual_is Character string for residual composition.
#' @param gstudy_obj Original gstudy object.
#' @param n_provided Logical indicating if n was explicitly provided.
#'
#' @return Named list of scaled variance draws.
#'
#' @keywords internal
scale_variance_draws <- function(vc_draws, n, object_spec, agg_facets,
  residual_is, gstudy_obj, n_provided = FALSE) {
  scaled_draws <- vc_draws

  do_rescaling <- !is.null(agg_facets) || (n_provided && length(n) > 0)

  if (!do_rescaling) {
    return(scaled_draws)
  }

  total_n <- if (length(n) > 0) prod(unlist(n)) else 1

  for (comp in names(scaled_draws)) {
    if (comp %in% object_spec) {
      next
    }

    if (comp == "Residual" || grepl("^Residual", comp)) {
      if (!is.null(agg_facets)) {
        agg_n <- 1
        for (f in agg_facets) {
          if (f %in% names(n)) agg_n <- agg_n * n[[f]]
        }
        scaled_draws[[comp]] <- scaled_draws[[comp]] / agg_n
      } else {
        scaled_draws[[comp]] <- scaled_draws[[comp]] / total_n
      }
      next
    }

    comp_facets <- parse_component_facets(comp)

    has_agg <- !is.null(agg_facets) && any(agg_facets %in% comp_facets)

    if (has_agg) {
      agg_in_comp <- comp_facets[comp_facets %in% agg_facets]
      agg_n <- 1
      for (f in agg_in_comp) {
        if (f %in% names(n)) agg_n <- agg_n * n[[f]]
      }
      if (agg_n > 1) {
        scaled_draws[[comp]] <- scaled_draws[[comp]] / agg_n
      }
    } else {
      scale_factor <- 1
      for (facet in comp_facets) {
        if (facet %in% names(n)) {
          scale_factor <- scale_factor * n[[facet]]
        }
      }
      if (scale_factor > 1) {
        scaled_draws[[comp]] <- scaled_draws[[comp]] / scale_factor
      }
    }
  }

  if (!is.null(agg_facets) && length(agg_facets) > 0 && !is.null(residual_is)) {
    for (comp in names(scaled_draws)) {
      if (comp == "Residual" || grepl("^Residual", comp)) {
        residual_facets <- parse_component_facets(residual_is)
        intersect_facets <- agg_facets[agg_facets %in% residual_facets]
        if (length(intersect_facets) > 0) {
          agg_factor <- 1
          for (facet in intersect_facets) {
            if (facet %in% names(n)) {
              agg_factor <- agg_factor * n[[facet]]
            }
          }
          if (agg_factor > 1) {
            scaled_draws[[comp]] <- scaled_draws[[comp]] * agg_factor
          }
        }
      }
    }
  }

  scaled_draws
}

#' Compute Relative Error Variance from Draws
#'
#' @param scaled_draws Named list of scaled variance draws.
#' @param universe_spec Character vector of universe components.
#' @param error_spec Character vector of error components (or NULL).
#' @param agg_facets Character vector of aggregation facets.
#'
#' @return Numeric vector of relative error variances.
#'
#' @keywords internal
compute_relative_error_draws <- function(scaled_draws, universe_spec, error_spec,
agg_facets, residual_is = NULL) {
n_draws <- length(scaled_draws[[1]])

if (!is.null(error_spec)) {
result <- rep(0, n_draws)
for (comp in error_spec) {
if (comp %in% names(scaled_draws)) {
result <- result + scaled_draws[[comp]]
}
}
return(result)
}

result <- rep(0, n_draws)

for (comp in names(scaled_draws)) {
if (comp %in% universe_spec) {
next
}

if (comp == "Residual" || grepl("^Residual", comp)) {
result <- result + scaled_draws[[comp]]
next
}

comp_facets <- parse_component_facets(comp)
is_agg_main <- length(comp_facets) == 1 && comp_facets %in% agg_facets
if (is_agg_main) {
next
}

if (any(universe_spec %in% comp_facets)) {
result <- result + scaled_draws[[comp]]
}
}

result
}

#' Compute Absolute Error Variance from Draws
#'
#' @param scaled_draws Named list of scaled variance draws.
#' @param universe_spec Character vector of universe components.
#' @param error_spec Character vector of error components (or NULL).
#' @param agg_facets Character vector of aggregation facets.
#'
#' @return Numeric vector of absolute error variances.
#'
#' @keywords internal
compute_absolute_error_draws <- function(scaled_draws, universe_spec, error_spec,
agg_facets, residual_is = NULL) {
n_draws <- length(scaled_draws[[1]])

if (!is.null(error_spec)) {
result <- rep(0, n_draws)
for (comp in error_spec) {
if (comp %in% names(scaled_draws)) {
result <- result + scaled_draws[[comp]]
}
}
return(result)
}

result <- rep(0, n_draws)

for (comp in names(scaled_draws)) {
if (comp %in% universe_spec) {
next
}

    if (comp == "Residual" || grepl("^Residual", comp)) {
      result <- result + scaled_draws[[comp]]
      next
    }

    comp_facets <- parse_component_facets(comp)
    is_agg_main <- length(comp_facets) == 1 && comp_facets %in% agg_facets
    if (is_agg_main) {
      next
    }

    result <- result + scaled_draws[[comp]]
  }

  result
}
