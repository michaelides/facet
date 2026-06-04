#' Coefficient Helper Functions
#'
#' Parsing and predicate functions used across coefficient calculations.
#' These are internal functions that support the main coefficient computation.
#'
#' @name coefficient-helpers
#' @keywords internal
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

#' Parse Specification for Object or Error Components
#'
#' Parses a specification for object of measurement or error components.
#' Accepts:
#' - A single character string: "p"
#' - A character vector: c("p", "p:d")
#' - A formula: obj ~ p + p:d (LHS is ignored for error spec, or can be used for object name)
#' - A one-sided formula: ~ p + p:d (extracts variables from RHS)
#'
#' @param x Specification (character, character vector, or formula).
#' @return A character vector of component names.
#'
#' @keywords internal
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
