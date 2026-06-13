# Conduct a Generalizability (G) Study

Estimate variance components for a generalizability theory analysis
using mixed effects models. Supports frequentist (lme4), Bayesian
(brms), and method of moments (mom) backends.

## Usage

``` r
gstudy(
  formula,
  data,
  backend = c("auto", "lme4", "brms", "mom"),
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

- backend:

  Character string specifying the backend to use. One of "auto"
  (default, chooses based on formula type), "lme4", "brms", or "mom".
  "mom" uses the method of moments (ANOVA-based) for variance component
  estimation.

- facets:

  Character vector of facet names (optional, auto-detected from formula
  if NULL).

- nested:

  Optional named list specifying nesting relationships. Names are nested
  facets and values are nesting facets. For example,
  `list(task = "rater")` means task is nested within rater. If NULL
  (default), nesting is auto-detected from the data structure.

- unbalanced:

  Logical indicating whether to enable unbalanced multivariate
  estimation. Default is FALSE. When TRUE:

  - For **mom backend**: Each dimension is analyzed with its available
    data using Henderson's Method III for variance component estimation.
    Correlations are computed using pairwise complete cases.

  - For **brms backend**: Not implemented. Use long-format specification
    with
    `bf(Score ~ 0 + Dimension + (0+Dimension|Facet), sigma ~ 0 + Dimension)`
    for sparse multivariate data.

  - For **lme4 backend**: Not applicable (lme4 does not support
    multivariate models).

- ci_method:

  Character string specifying the method for confidence intervals for
  the lme4 backend. One of "none" (default, no CIs), "profile" (more
  accurate, slower), or "boot" (bootstrap, most accurate, slowest). Only
  applicable for lme4 backend; brms and mom provide CIs automatically.

- nsim:

  Integer: number of bootstrap simulations (only for ci_method =
  "boot"). Default is 1000.

- boot.type:

  Character: bootstrap type, "perc" (percentile), "basic", or "norm"
  (normal-theory) (only for ci_method = "boot"). Default is "perc".

- prior:

  A brmsprior object or list of priors created by
  [`set_prior()`](https://github.com/yourorg/facet/reference/set_prior.md)
  or related functions. Only applicable when using brms backend. Use
  [`default_prior()`](https://github.com/yourorg/facet/reference/default_prior.md)
  to see available parameters for priors.

- ...:

  Additional arguments passed to the backend fitting function (e.g.,
  `lmer()` or `brm()`) and to `confint.merMod` for confidence intervals.

## Value

An object of class "gstudy" containing:

- model:

  The fitted model object from the backend

- variance_components:

  A tibble of estimated variance components

- facets:

  Character vector of facet names

- facet_n:

  Named numeric vector of sample sizes for each main effect facet

- sample_size_info:

  Comprehensive sample size information including main effects,
  interactions, residual, and nested effects

- backend:

  The backend used for fitting

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
# Basic univariate G-study with lme4 (default)
g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task) + (1 | Person:Rater) + (1 | Task:Rater),
data = brennan
)
#> Error: Interaction term(s) have as many or more levels than observations, making them confounded with the residual:
#>  - Person:Rater: 120 levels (120 observations)
#> 
#> This occurs when an interaction has a unique combination for each observation, meaning the interaction variance cannot be distinguished from residual variance.
#> 
#> Solutions:
#>  1. Remove the problematic interaction term from your model
#>  2. Use method = "mom" for ANOVA-based estimation which handles this differently
#>  3. Use backend = "brms" for Bayesian estimation
#> 
#> Note: In a fully crossed design where all facets are crossed with each other, the highest-order interaction is typically confounded with the residual and should not be included.

# G-study with profile confidence intervals
if (FALSE) { # \dontrun{
g_prof <- gstudy(Score ~ (1 | Person) + (1 | Task),
  data = brennan,
  ci_method = "profile"
)
} # }

# Method of moments G-study (ANOVA-based)
g_mom <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task) + (1 | Person:Rater) + (1 | Task:Rater),
data = brennan,
backend = "mom"
)
#> Warning: Error() model is singular

# \donttest{
# Bayesian G-study with brms
g_bayes <- gstudy(Score ~ (1 | Person) + (1 | Task),
  data = brennan,
  backend = "brms"
)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.6e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.26 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.123 seconds (Warm-up)
#> Chain 1:                0.148 seconds (Sampling)
#> Chain 1:                0.271 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 6e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.131 seconds (Warm-up)
#> Chain 2:                0.095 seconds (Sampling)
#> Chain 2:                0.226 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 1.7e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.17 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.136 seconds (Warm-up)
#> Chain 3:                0.123 seconds (Sampling)
#> Chain 3:                0.259 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 7e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.07 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.125 seconds (Warm-up)
#> Chain 4:                0.109 seconds (Sampling)
#> Chain 4:                0.234 seconds (Total)
#> Chain 4: 
#> Warning: There were 12 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems

# Bayesian G-study with custom priors
my_prior <- set_prior("normal(0, 1)", class = "sd", group = "Person")
g_bayes_prior <- gstudy(Score ~ (1 | Person) + (1 | Task),
  data = brennan,
  prior = my_prior,
  backend = "brms"
)
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.27 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.129 seconds (Warm-up)
#> Chain 1:                0.132 seconds (Sampling)
#> Chain 1:                0.261 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 9e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.128 seconds (Warm-up)
#> Chain 2:                0.116 seconds (Sampling)
#> Chain 2:                0.244 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 6e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.114 seconds (Warm-up)
#> Chain 3:                0.107 seconds (Sampling)
#> Chain 3:                0.221 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 6e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.134 seconds (Warm-up)
#> Chain 4:                0.109 seconds (Sampling)
#> Chain 4:                0.243 seconds (Total)
#> Chain 4: 
#> Warning: There were 15 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

# Multivariate G-study (automatically uses brms)
# Formatting data for multivariate example
b_wide <- tidyr::pivot_wider(brennan, names_from = Task, values_from = Score,
  names_prefix = "Task")
g_multi <- gstudy(mvbind(Task1, Task2) ~ (1 | Person) + (1 | Rater),
  data = b_wide
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
