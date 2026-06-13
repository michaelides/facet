# Conduct a Decision (D) Study

Compute generalizability and dependability coefficients based on the
results of a G-study. Allows for exploring different measurement designs
by modifying the number of levels for each facet.

## Usage

``` r
dstudy(
  gstudy_obj,
  n = list(),
  universe = NULL,
  error = NULL,
  aggregation = NULL,
  residual_is = NULL,
  estimation = NULL,
  cut_score = NULL,
  ci = NULL,
  probs = c(0.025, 0.975),
  weights = NULL,
  ...
)
```

## Arguments

- gstudy_obj:

  An object of class "gstudy" from
  [`gstudy()`](https://github.com/yourorg/facet/reference/gstudy.md).

- n:

  A named list specifying the number of levels for each facet in the
  D-study design (e.g., `list(items = 10, raters = 3)`). Can also be a
  list of vectors to explore multiple sample sizes. For multivariate
  G-studies fit with `unbalanced = TRUE`, this may instead be a tibble
  with columns `dim`, `facet`, and `n`, giving a different sample size
  per (dimension, facet). Supplying multiple rows for the same
  (dimension, facet) triggers a sweep, in which case the sweep is
  performed independently within each dimension.

- universe:

  Specification for components that contribute to the universe score.
  The object of measurement (first component from the G-study formula)
  is always included in the universe. Can be:

  - NULL (default): universe includes only the object of measurement

  - A character string: "person:item" (object is automatically added)

  - A character vector: c("person", "person:item") (object should be
    included)

  - A formula: ~ person + person:item Non-object components in the
    universe are estimated from unscaled (G-study) variance components.
    For example, if the object is "p" and universe includes "p:o", the
    universe score variance is computed as: var(p) + var(p:o).

- error:

  Specification for error components. Can be:

  - NULL (default): all components not in universe (and not in
    aggregation)

  - A character string: "person:item"

  - A character vector: c("person:item", "person:rater")

  - A formula: ~ person:item + person:rater Note: If a component is
    specified in both `universe` and `error`, an error is raised. Note:
    If a component is specified in both `aggregation` and `error`, a
    warning is issued and the component is removed from error
    (interaction terms are preserved).

- aggregation:

  Character vector of facets to aggregate over. Components containing
  these facets will be divided by the sample size. For example,
  `aggregation = "item"` divides item and person:item by n_item. When
  aggregation is specified, the main effects of aggregation facets are
  automatically excluded from error variance (but interaction terms are
  kept). Default is NULL (no aggregation).

- residual_is:

  Character string specifying which facets make up the residual. For
  example, "person:item:rater" means the residual represents the
  three-way interaction. Required when aggregation is specified and you
  want the residual to be rescaled. Default is NULL (residual is not
  rescaled).

- estimation:

  Character string specifying how to calculate coefficients:

  - "simple": uses point estimates from variance components (for
    lme4/mom backends)

  - "posterior": uses full posterior distributions (for brms backend)
    For brms backend, posterior estimation is always used to ensure
    consistency between variance component estimates and coefficient
    calculations. This avoids Jensen's inequality bias that would occur
    if variance estimates were computed as mean(SD)^2 rather than
    mean(SD^2). When estimation = "posterior" is requested with non-brms
    backend, a warning is issued and the gstudy model is refit with
    backend = "brms".

- cut_score:

  Optional numeric value specifying a cutoff score for
  criterion-referenced decisions. When provided, calculates phi-cut
  coefficient (phi_cut) in addition to standard phi coefficient. For
  multivariate models, can be a single value applied to all dimensions.
  Default is NULL (no phi-cut calculation).

- ci:

  Character vector specifying which coefficients to compute credible
  intervals for. Options: "g", "phi", "phi-cut". Can specify multiple:
  `ci = c("g", "phi")`. Credible intervals are only available when using
  the brms backend. Default is NULL (no credible intervals computed).

- probs:

  Numeric vector of length 2 specifying the quantiles for credible
  interval calculation. Default is `c(0.025, 0.975)` for a 95% credible
  interval. Works like the
  [`quantile()`](https://rdrr.io/r/stats/quantile.html) function. Only
  used when `ci` is specified.

- weights:

  Numeric vector of weights for computing composite coefficients. Length
  must match the number of dimensions in the multivariate design.
  Default is NULL, which uses equal weights (1 for each dimension). Only
  applicable for multivariate G-studies (mgstudy objects).

- ...:

  Additional arguments (currently unused).

## Value

An object of class "dstudy" containing:

- gstudy:

  The original G-study object

- variance_components:

  A tibble of variance components for the D-study with columns:

  - component: Name of variance component

  - var_unscaled: Unscaled variance estimate (from G-study)

  - pct_unscaled: Percentage of total unscaled variance

  - var_scaled: Scaled variance (divided by D-study sample sizes)

  - pct_scaled: Percentage of total scaled variance

  - dim: Dimension/response variable (for multivariate models)

  For multivariate models with posterior estimation, variance_components
  includes additional rows with dim = "Composite" showing the weighted
  composite variance for each component type, computed with full
  posterior propagation.

- coefficients:

  A tibble with G and D coefficients. For multivariate designs, the
  coefficients tibble includes an additional row with dim = "Composite"
  showing the weighted composite G and Phi coefficients.

- n:

  The number of levels for each facet

- object:

  The object of measurement (first component from G-study)

- universe:

  The universe components specification

- error:

  The error specification (if provided)

- aggregation:

  The aggregation specification (if provided)

- residual_is:

  The residual specification (if provided)

- residual_composition:

  The facets that make up the residual (from parse_residual_facets)

- estimation:

  The estimation method used ("simple" or "posterior")

- posterior:

  List of posterior distributions (only when estimation = "posterior")

- composite_posterior:

  Posterior draws of composite variance components for downstream
  analysis (only for multivariate models with posterior estimation)

- cut_score:

  The cutoff score used for phi-cut calculation (if provided)

- mu_y:

  The grand mean(s) used for phi-cut calculation (if cut_score provided)

- ci:

  The credible interval specification (if provided)

- probs:

  The probability levels used for credible interval calculation (if ci
  provided)

## Details

### Universe, Error, and Aggregation

By default, the *universe score* variance contains only the object of
measurement (e.g., `Person`), and all remaining components contribute to
*error* variance.

- Use `universe` to include additional components in the universe score
  (e.g., a facet that is fixed in your applied setting).

- Use `error` to restrict which components count as error (e.g., when
  some facets are considered fixed and not a source of measurement error
  for your inference).

- Use `aggregation` when the measurement procedure involves averaging
  over a facet (e.g., averaging across multiple raters), so that the
  variance components for that facet are divided by their respective
  sample sizes before computing coefficients.

### Composite Coefficients for Multivariate Designs

For multivariate G-studies (mgstudy objects), composite coefficients
provide an overall reliability estimate for a weighted sum of dimension
scores. The composite score is:

`Y_composite = sum(w_d * Y_d)`

The variance of this composite for any variance component is computed
as:

`sigma^2_composite = w' * Sigma * w`

Where:

- `w` is the vector of weights (one per dimension)

- `Sigma` is the variance-covariance matrix for that component

Expanded, this equals:

`sigma^2_composite = sum(w_d^2 * sigma^2_d) + 2 * sum(w_d * w_d' * sigma_dd')`

The covariances (off-diagonal elements) contribute when dimensions are
correlated. This means composite reliability often exceeds the average
of dimension-specific reliabilities, as the composite "borrows strength"
from correlations.

The `weights` parameter controls the contribution of each dimension:

- Default: `NULL` uses equal weights (1 for each dimension)

- Custom: Provide a numeric vector with length matching the number of
  dimensions

For posterior estimation with brms backend, composite coefficients are
computed for each posterior draw, properly propagating uncertainty to
the final estimates.

## See also

[`gstudy()`](https://github.com/yourorg/facet/reference/gstudy.md) for
conducting G-studies

Other decision studies:
[`coefficients`](https://github.com/yourorg/facet/reference/coefficients.md)

## Examples

``` r
# First conduct a G-study using the brennan dataset
# (Person crossed with Task and Rater)
g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task) + (1 | Person:Rater) + (1 | Task:Rater),
data = brennan)
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

# D-study with specific sample sizes
d <- dstudy(g, n = list(Task = 3, Rater = 4))
#> Error: object 'g' not found
print(d)
#> Error: object 'd' not found

# D-study exploring multiple sample sizes (sweep)
d_sweep <- dstudy(g, n = list(Task = c(3, 5, 10), Rater = c(2, 4, 8)))
#> Error: object 'g' not found
print(d_sweep)
#> Error: object 'd_sweep' not found

# D-study with aggregation (averaging over Raters)
d_agg <- dstudy(g, n = list(Task = 3, Rater = 4),
  aggregation = "Rater",
  residual_is = "Person:Task:Rater")
#> Error: object 'g' not found

# \donttest{
# D-study with posterior estimation (requires brms backend)
g_brms <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater),
  data = brennan, backend = "brms")
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.25 seconds.
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
#> Chain 1:  Elapsed Time: 0.285 seconds (Warm-up)
#> Chain 1:                0.232 seconds (Sampling)
#> Chain 1:                0.517 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.1 seconds.
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
#> Chain 2:  Elapsed Time: 0.222 seconds (Warm-up)
#> Chain 2:                0.201 seconds (Sampling)
#> Chain 2:                0.423 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 9e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
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
#> Chain 3:  Elapsed Time: 0.275 seconds (Warm-up)
#> Chain 3:                0.228 seconds (Sampling)
#> Chain 3:                0.503 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 9e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
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
#> Chain 4:  Elapsed Time: 0.255 seconds (Warm-up)
#> Chain 4:                0.223 seconds (Sampling)
#> Chain 4:                0.478 seconds (Total)
#> Chain 4: 
#> Warning: There were 20 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
d_post <- dstudy(g_brms, n = list(Task = 3, Rater = 4))

# D-study with credible intervals (requires brms backend)
d_ci <- dstudy(g_brms, n = list(Task = 3, Rater = 4), ci = c("g", "phi"))
print(d_ci)
#> Decision Study (D-Study)
#> ========================
#> 
#> Based on G-Study with brms backend
#> Object of measurement: Person 
#> Universe components: Person 
#> Error components for relative error (sigma2_delta): Person:Rater (Residual) 
#> Error components for absolute error (sigma2_delta_abs): Task, Rater, Person:Rater (Residual) 
#> 
#> Sample Sizes:
#>  Task: 3
#>  Rater: 4
#> 
#> Variance Components:
#> # A tibble: 4 × 6
#>   component dim   var_unscaled pct_unscaled var_scaled pct_scaled
#>   <chr>     <chr>        <dbl>        <dbl>      <dbl>      <dbl>
#> 1 Person    Score        0.949         14.1      0.949      37.1 
#> 2 Task      Score        1.96          29.0      0.652      25.5 
#> 3 Rater     Score        0.933         13.8      0.233       9.12
#> 4 Residual  Score        2.90          43.0      0.725      28.3 
#> 
#> Coefficients:
#> # A tibble: 1 × 9
#>     uni sigma2_delta sigma2_delta_abs     g   phi  g_LL  g_UL phi_LL phi_UL
#>   <dbl>        <dbl>            <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1 0.949        0.242             1.13 0.733 0.466 0.360 0.927 0.0888  0.817

# Custom probability levels (90% credible interval)
d_ci_90 <- dstudy(g_brms, n = list(Task = 3, Rater = 4),
  ci = "g", probs = c(0.05, 0.95))
# }
```
