# Conduct a Generalizability (G) Study

Estimate variance components for a generalizability theory analysis
using mixed effects models. Supports frequentist (lme4), Bayesian
(brms), and method of moments (aov) estimators.

## Usage

``` r
gstudy(
  formula,
  data,
  estimator = c("auto", "lme4", "brms", "aov"),
  facets = NULL,
  nested = NULL,
  unbalanced = FALSE,
  ci_method = c("none", "profile", "boot"),
  nsim = 1000,
  boot.type = c("perc", "basic", "norm"),
  prior = NULL,
  ...
)
```

## Arguments

- formula:

  A formula specifying the G-study model. The left-hand side should be
  the outcome variable, and the right-hand side should specify the
  variance components using lme4-style syntax (e.g.,
  `y ~ (1|person) + (1|item)`). For multivariate models, use brms
  syntax: `mvbind(y1, y2) ~ (1|person)`.

- data:

  A data frame containing the variables in the formula.

- estimator:

  Character string specifying the estimator to use. One of "auto"
  (default, chooses based on formula type), "lme4", "brms", or "aov".
  "aov" uses the method of moments (ANOVA-based) for variance component
  estimation.

- facets:

  Character vector of facet names (optional, auto-detected from formula
  if NULL).

- nested:

  Optional named list specifying nesting relationships. Names are nested
  facets and values are nesting facets. For example,
  `list(task = "rater")` means task is nested within rater. If NULL
  (default), nesting is auto-detected from the data structure. For the
  aov estimator, this is used by
  [`adapt_aov_formula()`](https://github.com/yourorg/facet/reference/adapt_aov_formula.md)
  to rewrite the aov Error() decomposition so designs with nested facets
  (e.g. the brennan dataset, where Rater is nested within Task) fit
  without spurious "Error() model is singular" warnings from
  [`stats::aov()`](https://rdrr.io/r/stats/aov.html). The lme4 and brms
  estimators do not require this adaptation.

- unbalanced:

  Logical indicating whether to enable unbalanced multivariate
  estimation. Default is FALSE. When TRUE:

  - For **aov estimator**: Each dimension is analyzed with its available
    data using Henderson's Method III for variance component estimation.
    Correlations are computed using pairwise complete cases.

  - For **brms estimator**: Not implemented. Use long-format
    specification with
    `bf(Score ~ 0 + Dimension + (0+Dimension|Facet), sigma ~ 0 + Dimension)`
    for sparse multivariate data.

  - For **lme4 estimator**: Not applicable (lme4 does not support
    multivariate models).

- ci_method:

  Character string specifying the method for confidence intervals for
  the lme4 estimator. One of "none" (default, no CIs), "profile" (more
  accurate, slower), or "boot" (bootstrap, most accurate, slowest). Only
  applicable for lme4 estimator; brms and mom provide CIs automatically.

- nsim:

  Integer: number of bootstrap simulations (only for ci_method =
  "boot"). Default is 1000.

- boot.type:

  Character: bootstrap type, "perc" (percentile), "basic", or "norm"
  (normal-theory) (only for ci_method = "boot"). Default is "perc".

- prior:

  A brmsprior object or list of priors created by
  [`set_prior()`](https://github.com/yourorg/facet/reference/set_prior.md)
  or related functions. Only applicable when using brms estimator. Use
  [`default_prior()`](https://github.com/yourorg/facet/reference/default_prior.md)
  to see available parameters for priors.

- ...:

  Additional arguments passed to the estimator fitting function (e.g.,
  `lmer()` or `brm()`) and to `confint.merMod` for confidence intervals.

## Value

An object of class "gstudy" containing:

- model:

  The fitted model object from the estimator

- variance_components:

  A tibble of estimated variance components

- facets:

  Character vector of facet names

- facet_n:

  Named numeric vector of sample sizes for each main effect facet

- sample_size_info:

  Comprehensive sample size information including main effects,
  interactions, residual, and nested effects

- estimator:

  The estimator used for fitting

- is_multivariate:

  Logical indicating if the model is multivariate

- is_unbalanced:

  Logical indicating if the multivariate model was fit with
  `unbalanced = TRUE`. Only set for multivariate models. When TRUE, the
  per-dimension totals are stored in `n_per_dim` and
  [`dstudy()`](https://github.com/yourorg/facet/reference/dstudy.md)
  will accept a per-dimension `n` tibble.

- n_per_dim:

  Named list of per-dimension observation counts (only populated for
  multivariate models fit with `unbalanced = TRUE`).

- sample_size_info_per_dim:

  Named list of per-dimension sample size information (only populated
  for multivariate models fit with `unbalanced = TRUE`). Each element is
  itself a list keyed by facet/component, mirroring the structure of
  `sample_size_info`.

- formula:

  The formula used

- data:

  The original data

- n_obs:

  Number of observations

## See also

[`dstudy()`](https://github.com/yourorg/facet/reference/dstudy.md) for
conducting D-studies

## Examples

``` r
# Basic univariate G-study with lme4 (default).
# The brennan data is Person x Task x (Rater:Task) with one observation
# per Person x Rater:Task cell, so the canonical "all possible variance
# components" formula has the three main effects plus the only estimable
# interaction (Person:Task); the Person x Rater:Task interaction is
# confounded with the residual.
g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
data = brennan
)

# G-study with profile confidence intervals
if (FALSE) { # \dontrun{
g_prof <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
  data = brennan,
  ci_method = "profile"
)
} # }

# Method of moments G-study (ANOVA-based)
g_mom <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
data = brennan,
estimator = "aov"
)

# \donttest{
# Bayesian G-study with brms
g_bayes <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
  data = brennan,
  estimator = "brms",
  iter = 2000, cores = 4, refresh = 1000
)
#> Compiling Stan program...
#> Start sampling
#> Warning: There were 32 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems

# Bayesian G-study with custom priors
my_prior <- set_prior("normal(0, 1)", class = "sd", group = "Person")
g_bayes_prior <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
  data = brennan,
  prior = my_prior,
  estimator = "brms",
  iter = 2000, cores = 4, refresh = 1000
)
#> Compiling Stan program...
#> Start sampling
#>  
#> Warning: There were 13 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems

# Multivariate G-study (automatically uses brms)
# Formatting data for multivariate example
b_wide <- tidyr::pivot_wider(brennan, names_from = Task, values_from = Score,
  names_prefix = "Task")
g_multi <- gstudy(mvbind(Task1, Task2) ~ (1 | Person) + (1 | Rater),
  data = b_wide,
  iter = 2000, cores = 4, refresh = 1000
)
#> Missing values detected in response variables.
#> 
#>  - Original observations: 120
#>  - Complete cases used: 0 (120 rows removed)
#> 
#> Facet levels remain consistent across all dimensions.
#> Setting 'rescor' to TRUE by default for this model
#> Warning: In the future, 'rescor' will be set to FALSE by default for all models. It is thus recommended to explicitely set 'rescor' via 'set_rescor' instead of using the default.
#> Warning: Rows containing NAs were excluded from the model.
#> Error: Error fitting model with brms: All observations in the data were removed presumably because of NA values.
# }
```
