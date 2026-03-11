#' Prior Wrapper Functions for brms
#'
#' These functions are wrappers around brms prior functions to allow
#' specifying priors for Bayesian models fit via the brms backend in gstudy.
#'
#' @name prior-functions
NULL

#' Define Priors for brms Models
#'
#' This function is a wrapper around [brms::set_prior()] to allow
#' specifying priors for Bayesian models fit via the brms backend in gstudy.
#' See [brms::set_prior()] for full documentation on parameter specifications.
#'
#' @param prior A character string defining a distribution in Stan language.
#' @param class The parameter class. Defaults to `"b"` (i.e., population-level effects).
#'   See Details in [brms::set_prior()] for other valid parameter classes.
#' @param coef For class `"b"`, the coefficient name.
#' @param group For class `"sd"` or `"cor"`, the grouping factor.
#' @param resp For multivariate models, the response variable.
#' @param dpar For distributional parameters, the parameter name.
#' @param nlpar For non-linear parameters, the parameter name.
#' @param lb Lower bound for parameters with lower bounds.
#' @param ub Upper bound for parameters with upper bounds.
#' @param check Logical; check if prior is valid for the model.
#' @param ... Additional arguments passed to [brms::set_prior()].
#'
#' @return A brmsprior object; a data frame containing the prior specifications.
#'
#' @references
#' Paul Bürkner (2017). brms: An R Package for Bayesian Multilevel Models.
#' \emph{Journal of Statistical Software}, 80(1), 1-28.
#' \doi{10.18637/jss.v080.i01}
#'
#' @rdname set_prior
#' @export
#'
#' @examples
#' \dontrun{
#' # Define a prior for standard deviations of random effects
#' prior <- set_prior("normal(0, 1)", class = "sd", group = "person")
#'
#' # Use with gstudy
#' g <- gstudy(score ~ (1 | person) + (1 | item),
#'   data = mydata,
#'   prior = prior,
#'   backend = "brms"
#' )
#' }
set_prior <- function(prior, class = "b", coef = "", group = "", resp = "",
                      dpar = "", nlpar = "", lb = NA, ub = NA, check = TRUE, ...) {
  brms::set_prior(
    prior = prior,
    class = class,
    coef = coef,
    group = group,
    resp = resp,
    dpar = dpar,
    nlpar = nlpar,
    lb = lb,
    ub = ub,
    check = check,
    ...
  )
}

#' Get Default Priors for Bayesian Models
#'
#' This function is a wrapper around [brms::default_prior()] to allow
#' inspecting default priors that would be applied to a model.
#' See [brms::default_prior()] for full documentation.
#'
#' @param object An object whose class determines the method. Typically a formula,
#'   brmsformula, or mvbrmsformula.
#' @param ... Additional arguments passed to [brms::default_prior()] including
#'   `data` (a data frame containing the model data).
#'
#' @return A brmsprior object; a data frame containing information about all
#'   parameters for which priors can be specified.
#'
#' @references
#' Paul Bürkner (2017). brms: An R Package for Bayesian Multilevel Models.
#' \emph{Journal of Statistical Software}, 80(1), 1-28.
#' \doi{10.18637/jss.v080.i01}
#'
#' @rdname default_prior
#' @export
#'
#' @examples
#' \dontrun{
#' # Get default priors for a gstudy model
#' prior_info <- default_prior(score ~ (1 | person) + (1 | item), data = mydata)
#' print(prior_info)
#' }
default_prior <- function(object, ...) {
  brms::default_prior(object, ...)
}

#' Define Priors Using Tidy Evaluation
#'
#' This function is a wrapper around [brms::prior()] to allow
#' specifying priors using tidy evaluation syntax.
#' See [brms::prior()] for full documentation.
#'
#' @param prior A call or expression defining a distribution.
#' @param ... Additional arguments passed to [brms::prior()].
#'
#' @return A brmsprior object.
#'
#' @references
#' Paul Bürkner (2017). brms: An R Package for Bayesian Multilevel Models.
#' \emph{Journal of Statistical Software}, 80(1), 1-28.
#' \doi{10.18637/jss.v080.i01}
#'
#' @rdname prior
#' @export
#'
#' @examples
#' \dontrun{
#' prior <- prior(normal(0, 1), class = sd)
#' g <- gstudy(score ~ (1 | person), data = mydata, prior = prior, backend = "brms")
#' }
prior <- function(prior, ...) {
  brms::prior(prior, ...)
}

#' Define Priors Using Formula Syntax
#'
#' This function is a wrapper around [brms::prior_()] to allow
#' specifying priors using formula syntax (with `~`).
#' See [brms::prior_()] for full documentation.
#'
#' @param prior A one-sided formula defining the prior distribution.
#' @param ... Additional arguments passed to [brms::prior_()].
#'
#' @return A brmsprior object.
#'
#' @references
#' Paul Bürkner (2017). brms: An R Package for Bayesian Multilevel Models.
#' \emph{Journal of Statistical Software}, 80(1), 1-28.
#' \doi{10.18637/jss.v080.i01}
#'
#' @rdname prior_
#' @export
#'
#' @examples
#' \dontrun{
#' prior <- prior_(~ normal(0, 1), class = "sd")
#' g <- gstudy(score ~ (1 | person), data = mydata, prior = prior, backend = "brms")
#' }
prior_ <- function(prior, ...) {
  brms::prior_(prior, ...)
}

#' Define Priors Using Character Strings
#'
#' This function is a wrapper around [brms::prior_string()] to allow
#' specifying priors using character strings.
#' See [brms::prior_string()] for full documentation.
#'
#' @param prior A character string defining the prior distribution.
#' @param ... Additional arguments passed to [brms::prior_string()].
#'
#' @return A brmsprior object.
#'
#' @references
#' Paul Bürkner (2017). brms: An R Package for Bayesian Multilevel Models.
#' \emph{Journal of Statistical Software}, 80(1), 1-28.
#' \doi{10.18637/jss.v080.i01}
#'
#' @rdname prior_string
#' @export
#'
#' @examples
#' \dontrun{
#' prior <- prior_string("normal(0, 1)", class = "sd", group = "person")
#' g <- gstudy(score ~ (1 | person), data = mydata, prior = prior, backend = "brms")
#' }
prior_string <- function(prior, ...) {
  brms::prior_string(prior, ...)
}

#' Create an Empty Prior Object
#'
#' This function is a wrapper around [brms::empty_prior()] to create
#' an empty prior object to which priors can be added.
#' See [brms::empty_prior()] for full documentation.
#'
#' @return A brmsprior object with no rows.
#'
#' @references
#' Paul Bürkner (2017). brms: An R Package for Bayesian Multilevel Models.
#' \emph{Journal of Statistical Software}, 80(1), 1-28.
#' \doi{10.18637/jss.v080.i01}
#'
#' @rdname empty_prior
#' @export
#'
#' @examples
#' \dontrun{
#' # Start with empty prior and add custom priors
#' my_prior <- empty_prior()
#' my_prior <- rbind(my_prior, set_prior("normal(0, 1)", class = "sd", group = "person"))
#' g <- gstudy(score ~ (1 | person), data = mydata, prior = my_prior, backend = "brms")
#' }
empty_prior <- function() {
  brms::empty_prior()
}
