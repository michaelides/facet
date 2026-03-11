#' @keywords internal
#' @details
#' The facet package provides tools for conducting generalizability theory analyses
#' using variance component models. It supports multiple backends including lme4
#' for frequentist analyses and brms for Bayesian analyses.
#'
#' Key functions:
#' - [gstudy()]: Conduct a G-study to estimate variance components
#' - [dstudy()]: Conduct a D-study to compute reliability coefficients
#'
#' @seealso
#' Useful links:
#' - <https://github.com/yourusername/mgt>
#'
#' @importFrom lme4 lmer
#' @importFrom reformulas findbars
#' @importFrom brms brm bf mvbind set_rescor set_prior default_prior prior prior_ prior_string empty_prior
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr filter mutate select bind_rows case_when
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom rlang as_name as_label is_formula
#' @importFrom stats as.formula formula terms
#' @importFrom methods is
#' @importFrom ggplot2 ggplot aes geom_bar geom_col geom_errorbar labs theme_minimal coord_flip
#' @importFrom ggplot2 geom_pointrange geom_hline theme element_text facet_wrap
#' @importFrom broom tidy glance
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
