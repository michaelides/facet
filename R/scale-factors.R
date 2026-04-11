#' Scale Factor Functions
#'
#' Functions for computing scale factors in D-study variance component calculations.
#' These functions handle the scaling of variance components based on sample sizes.
#'
#' @name scale-factors
#' @keywords internal
NULL

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

  comp_facets <- parse_component_facets(component)

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

  component_for_scale <- if (!is.null(actual_comp)) actual_comp else comp

  if (component_for_scale == "Residual") {
    if (!is.null(residual_is) && residual_is != "") {
      comp_facets <- parse_component_facets(residual_is)
    } else {
      return(1)
    }
  } else {
    comp_facets <- parse_component_facets(component_for_scale)
  }

  scale_factor <- 1
  for (facet in comp_facets) {
    if (facet %in% names(n)) {
      scale_factor <- scale_factor * n[[facet]]
    }
  }

  scale_factor
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
compute_scale_factor_from_facets <- function(facets, n) {
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

#' Compute Scale Factor for a Component
#'
#' Computes the scale factor for a variance component.
#' Used in composite variance calculations.
#'
#' @param comp The variance component name.
#' @param n Named list of sample sizes.
#' @param object_spec Character vector of object components.
#' @param n_provided Logical indicating if n was explicitly provided.
#' @return The scale factor for this component.
#'
#' @keywords internal
compute_component_scale_factor <- function(comp, n, object_spec, n_provided) {
  if (comp %in% object_spec) {
    return(1)
  }

  total_n <- if (length(n) > 0) prod(unlist(n)) else 1

  if (n_provided && length(n) > 0) {
    if (comp == "Residual") {
      return(total_n)
    }

    comp_facets <- parse_component_facets(comp)
    scale_factor <- 1
    for (f in comp_facets) {
      if (f %in% names(n)) {
        scale_factor <- scale_factor * n[[f]]
      }
    }
    return(scale_factor)
  }

  return(1)
}

#' Compute Scale Factor for Composite Variance
#'
#' Computes the scale factor for scaling covariances in composite variance
#' calculation. The scale factor depends on the component type:
#' - For main effects (e.g., Item): scale_factor = n_facet
#' - For Residual: scale_factor = product of all error facet sample sizes
#'
#' @param component The variance component name.
#' @param n Named list of sample sizes.
#' @param object_spec Object of measurement specification.
#' @param universe_spec Universe components specification.
#'
#' @return The scale factor for this component.
#'
#' @keywords internal
compute_scale_factor <- function(component, n, object_spec, universe_spec) {
  if (is.null(n) || length(n) == 0) {
    return(1)
  }

  if (component == "Residual") {
    scale_factor <- 1
    for (facet_name in names(n)) {
      if (!(facet_name %in% object_spec)) {
        scale_factor <- scale_factor * n[[facet_name]]
      }
    }
    return(scale_factor)
  }

  facets <- parse_component_facets(component)

  if (length(facets) == 0) {
    return(1)
  }

  if (length(facets) == 1) {
    facet_name <- facets[1]
    if (facet_name %in% names(n)) {
      return(n[[facet_name]])
    }
    return(1)
  }

  scale_factor <- 1
  for (facet in facets) {
    if (facet %in% names(n) && !(facet %in% object_spec)) {
      scale_factor <- scale_factor * n[[facet]]
    }
  }
  scale_factor
}
