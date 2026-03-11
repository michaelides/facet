#' S3 Class Definitions for mgt Package
#'
#' This file defines the S3 classes used throughout the mgt package.
#' The main classes are "gstudy" for G-study results and "dstudy" for D-study results.
#'
#' @name mgt-classes
#' @keywords internal
NULL

#' Create a gstudy Object
#'
#' Constructor function for the "gstudy" class. This is an internal low-level
#' constructor. Users should use [gstudy()] to create gstudy objects.
#'
#' @param model The fitted model object.
#' @param variance_components A tibble of variance components.
#' @param facets Character vector of facet names.
#' @param object Character string naming the object of measurement.
#' @param backend Character string indicating the backend used.
#' @param is_multivariate Logical indicating if the model is multivariate.
#' @param formula The formula used for fitting.
#' @param data The original data.
#' @return An object of class "gstudy".
#'
#' @keywords internal
new_gstudy <- function(model, variance_components, facets, object, backend,
                       is_multivariate, formula, data) {
  structure(
    list(
      model = model,
      variance_components = variance_components,
      facets = facets,
      object = object,
      backend = backend,
      is_multivariate = is_multivariate,
      formula = formula,
      data = data,
      n_obs = nrow(data)
    ),
    class = "gstudy"
  )
}

#' Create a dstudy Object
#'
#' Constructor function for the "dstudy" class. This is an internal low-level
#' constructor. Users should use [dstudy()] to create dstudy objects.
#'
#' @param gstudy The original G-study object.
#' @param variance_components A tibble of variance components for the D-study.
#' @param coefficients A tibble with G and D coefficients.
#' @param n A named list of facet levels.
#' @param object Character string naming the object of measurement.
#' @param is_sweep Logical indicating if this is a sweep over multiple sample sizes.
#' @return An object of class "dstudy".
#'
#' @keywords internal
new_dstudy <- function(gstudy, variance_components, coefficients, n, object,
                       is_sweep = FALSE) {
  structure(
    list(
      gstudy = gstudy,
      variance_components = variance_components,
      coefficients = coefficients,
      n = n,
      object = object,
      is_sweep = is_sweep
    ),
    class = "dstudy"
  )
}

#' Validate a gstudy Object
#'
#' Checks that a gstudy object has the required components and they are valid.
#'
#' @param x An object to validate.
#' @return TRUE if valid, otherwise raises an error.
#'
#' @keywords internal
validate_gstudy <- function(x) {
  if (!is.list(x)) {
    stop("gstudy object must be a list", call. = FALSE)
  }

  required <- c(
    "model", "variance_components", "facets", "backend",
    "is_multivariate", "formula", "data", "n_obs"
  )
  missing <- setdiff(required, names(x))

  if (length(missing) > 0) {
    stop("gstudy object is missing required components: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  # Validate variance_components is a tibble with required columns
  if (!tibble::is_tibble(x$variance_components)) {
    stop("variance_components must be a tibble", call. = FALSE)
  }

  required_vc_cols <- c("component", "var", "pct")
  missing_vc_cols <- setdiff(required_vc_cols, names(x$variance_components))
  if (length(missing_vc_cols) > 0) {
    stop("variance_components is missing required columns: ",
         paste(missing_vc_cols, collapse = ", "),
         call. = FALSE)
  }

  if (!"dim" %in% names(x$variance_components)) {
    warning(
      "variance_components is missing 'dim' column. ",
      "This is expected for older gstudy objects. ",
      "Consider recreating the gstudy object.",
      call. = FALSE
    )
  }

  required_vc_cols <- c("component", "var", "pct")
  missing_vc_cols <- setdiff(required_vc_cols, names(x$variance_components))
  if (length(missing_vc_cols) > 0) {
    stop("variance_components is missing required columns: ",
      paste(missing_vc_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Validate backend
  if (!x$backend %in% c("lme4", "brms", "mom")) {
    stop("backend must be 'lme4', 'brms', or 'mom'", call. = FALSE)
  }

  # Validate facet_n if present
  if (!is.null(x$facet_n)) {
    if (!is.numeric(x$facet_n)) {
      stop("facet_n must be numeric", call. = FALSE)
    }
  }

  TRUE
}

#' Validate a dstudy Object
#'
#' Checks that a dstudy object has the required components and they are valid.
#'
#' @param x An object to validate.
#' @return TRUE if valid, otherwise raises an error.
#'
#' @keywords internal
validate_dstudy <- function(x) {
  if (!is.list(x)) {
    stop("dstudy object must be a list", call. = FALSE)
  }

  required <- c("gstudy", "variance_components", "coefficients", "n", "object")
  missing <- setdiff(required, names(x))

  if (length(missing) > 0) {
    stop("dstudy object is missing required components: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  # Validate gstudy is a gstudy object
  if (!inherits(x$gstudy, "gstudy")) {
    stop("'gstudy' component must be a gstudy object", call. = FALSE)
  }

  # Validate coefficients is a data frame
  if (!is.data.frame(x$coefficients)) {
    stop("coefficients must be a data frame", call. = FALSE)
  }

  TRUE
}

#' Check if Object is a gstudy Object
#'
#' @param x An object to check.
#' @return TRUE if x is a gstudy object, FALSE otherwise.
#'
#' @export
is.gstudy <- function(x) {
  inherits(x, "gstudy")
}

#' Create an mgstudy Object
#'
#' Constructor function for the "mgstudy" class (multivariate G-study).
#' This is an internal low-level constructor.
#' Users should use [gstudy()] to create mgstudy objects.
#'
#' @param model The fitted model object.
#' @param variance_components A tibble of variance components with 'dim' column.
#' @param facets Character vector of facet names.
#' @param object Character string naming the object of measurement.
#' @param backend Character string indicating the backend used.
#' @param formula The formula used for fitting.
#' @param data The original data.
#' @param dimensions Character vector of dimension (response variable) names.
#' @return An object of class "mgstudy".
#'
#' @keywords internal
new_mgstudy <- function(model, variance_components, facets, object, backend,
                        formula, data, dimensions) {
  structure(
    list(
      model = model,
      variance_components = variance_components,
      facets = facets,
      object = object,
      backend = backend,
      is_multivariate = TRUE,
      formula = formula,
      data = data,
      n_obs = nrow(data),
      dimensions = dimensions
    ),
    class = "mgstudy"
  )
}

#' Validate an mgstudy Object
#'
#' Checks that an mgstudy object has the required components and they are valid.
#'
#' @param x An object to validate.
#' @return TRUE if valid, otherwise raises an error.
#'
#' @keywords internal
validate_mgstudy <- function(x) {
  if (!is.list(x)) {
    stop("mgstudy object must be a list", call. = FALSE)
  }

  required <- c(
    "model", "variance_components", "facets", "backend",
    "formula", "data", "n_obs", "dimensions"
  )
  missing <- setdiff(required, names(x))

  if (length(missing) > 0) {
    stop("mgstudy object is missing required components: ",
         paste(missing, collapse = ", "),
         call. = FALSE)
  }

  if (!tibble::is_tibble(x$variance_components)) {
    stop("variance_components must be a tibble", call. = FALSE)
  }

  required_vc_cols <- c("component", "dim", "var", "pct")
  missing_vc_cols <- setdiff(required_vc_cols, names(x$variance_components))
  if (length(missing_vc_cols) > 0) {
    stop("variance_components is missing required columns: ",
         paste(missing_vc_cols, collapse = ", "),
         call. = FALSE)
  }

  if (!x$backend %in% c("lme4", "brms", "mom")) {
    stop("backend must be 'lme4', 'brms', or 'mom'", call. = FALSE)
  }

  if (!is.character(x$dimensions) || length(x$dimensions) < 2) {
    stop("dimensions must be a character vector with at least 2 response names",
         call. = FALSE)
  }

  TRUE
}

#' Check if Object is an mgstudy Object
#'
#' @param x An object to check.
#' @return TRUE if x is an mgstudy object, FALSE otherwise.
#'
#' @export
is.mgstudy <- function(x) {
  inherits(x, "mgstudy")
}

#' Check if Object is a dstudy Object
#'
#' @param x An object to check.
#' @return TRUE if x is a dstudy object, FALSE otherwise.
#'
#' @export
is.dstudy <- function(x) {
  inherits(x, "dstudy")
}
