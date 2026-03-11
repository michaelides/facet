#' Wrapper Functions for brms Formula Components
#'
#' These functions are wrappers around brms formula components to allow
#' specifying complex models (like multivariate ones) via the brms backend in gstudy.
#'
#' @name brms-formula-helpers
NULL

#' Set up brms formulas
#'
#' This function is a wrapper around [brms::bf()] to allow specifying
#' formulas for Bayesian models fit via the brms backend in gstudy.
#' See [brms::bf()] for full documentation.
#'
#' @param formula A formula object.
#' @param ... Additional arguments passed to [brms::bf()].
#'
#' @return A brmsformula object.
#'
#' @rdname bf
#' @export
bf <- function(formula, ...) {
    brms::bf(formula, ...)
}

#' Bind response variables for multivariate models
#'
#' This function is a wrapper around [brms::mvbind()] to allow specifying
#' multivariate models via the brms backend in gstudy.
#' See [brms::mvbind()] for full documentation.
#'
#' @param ... Unquoted names of variables to bind.
#'
#' @return A matrix-like object for brms formulas.
#'
#' @rdname mvbind
#' @export
mvbind <- function(...) {
    brms::mvbind(...)
}

#' Set residual correlation between response variables
#'
#' This function is a wrapper around [brms::set_rescor()] to allow specifying
#' residual correlations in multivariate models fit via the brms backend.
#' See [brms::set_rescor()] for full documentation.
#'
#' @param rescor Logical; whether to estimate residual correlations.
#'
#' @return A brmsformula helper object.
#'
#' @rdname set_rescor
#' @export
set_rescor <- function(rescor) {
    brms::set_rescor(rescor)
}
