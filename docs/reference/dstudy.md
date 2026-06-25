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
# (Person crossed with Task; Rater nested in Task).
# Canonical "all possible variance components" formula:
g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
data = brennan)

# D-study with specific sample sizes
d <- dstudy(g, n = list(Task = 3, Rater = 4))
print(d)
#> Decision Study (D-Study)
#> ========================
#> 
#> Based on G Study with lme4 backend
#> Object of measurement: Person 
#> Universe components: Person 
#> Error components for relative error (sigma2_delta): Person:Task, Person:Rater (Residual) 
#> Error components for absolute error (sigma2_delta_abs): Task, Rater, Person:Task, Person:Rater (Residual) 
#> 
#> Sample Sizes:
#>  Task: 3
#>  Rater: 4
#> 
#> Variance Components:
#> # A tibble: 5 × 6
#>   component   dim   var_unscaled pct_unscaled var_scaled pct_scaled
#>   <chr>       <chr>        <dbl>        <dbl>      <dbl>      <dbl>
#> 1 Person      Score       0.473       10.7851     0.473     31.0181
#> 2 Task        Score       0.3253       7.4173     0.1084     7.1108
#> 3 Rater       Score       0.6475      14.7639     0.1619    10.6153
#> 4 Person:Task Score       0.5596      12.7597     0.1865    12.2324
#> 5 Residual    Score       2.3803      54.2741     0.5951    39.0234
#> 
#> Coefficients:
#>    dim   uni sigma2_delta sigma2_delta_abs     g    phi
#>  Score 0.473       0.7816            1.052 0.377 0.3102

# D-study exploring multiple sample sizes (sweep)
d_sweep <- dstudy(g, n = list(Task = c(3, 5, 10), Rater = c(2, 4, 8)))
print(d_sweep)
#> Decision Study (D-Study)
#> ========================
#> 
#> Based on G Study with lme4 backend
#> Object of measurement: Person 
#> Universe components: Person 
#> Error components for relative error (sigma2_delta): Person:Task, Person:Rater (Residual) 
#> Error components for absolute error (sigma2_delta_abs): Task, Rater, Person:Task, Person:Rater (Residual) 
#> 
#> Sample Size Sweep (by dimension):
#> 
#> Dimension: Score
#> ---------------------------------------- 
#> # A tibble: 9 × 8
#>    Task Rater dim     uni sigma2_delta sigma2_delta_abs     g   phi
#>   <dbl> <dbl> <chr> <dbl>        <dbl>            <dbl> <dbl> <dbl>
#> 1     3     2 Score 0.473        1.38             1.81  0.256 0.207
#> 2     5     2 Score 0.473        1.30             1.69  0.266 0.219
#> 3    10     2 Score 0.473        1.25             1.60  0.275 0.228
#> 4     3     4 Score 0.473        0.782            1.05  0.377 0.310
#> 5     5     4 Score 0.473        0.707            0.934 0.401 0.336
#> 6    10     4 Score 0.473        0.651            0.845 0.421 0.359
#> 7     3     8 Score 0.473        0.484            0.673 0.494 0.413
#> 8     5     8 Score 0.473        0.409            0.555 0.536 0.460
#> 9    10     8 Score 0.473        0.353            0.467 0.572 0.503
#> 

# D-study with aggregation (averaging over Raters)
d_agg <- dstudy(g, n = list(Task = 3, Rater = 4),
  aggregation = "Rater",
  residual_is = "Person:Task:Rater"
)

# \donttest{
# D-study with posterior estimation (requires brms backend)
g_brms <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
  data = brennan, backend = "brms",
  iter = 2000, cores = 4, refresh = 1000)
#> Compiling Stan program...
#> Trying to compile a simple C file
#> Running /Library/Frameworks/R.framework/Resources/bin/R CMD SHLIB foo.c
#> using C compiler: ‘Apple clang version 21.0.0 (clang-2100.1.1.101)’
#> using SDK: ‘MacOSX26.5.sdk’
#> clang -arch arm64 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I"/Users/afp18axu/Library/R/arm64/4.6/library/Rcpp/include/"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/RcppEigen/include/"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/RcppEigen/include/unsupported"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/BH/include" -I"/Users/afp18axu/Library/R/arm64/4.6/library/StanHeaders/include/src/"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/StanHeaders/include/"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/RcppParallel/include/"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DUSE_STANC3 -DSTRICT_R_HEADERS  -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION  -D_HAS_AUTO_PTR_ETC=0  -include '/Users/afp18axu/Library/R/arm64/4.6/library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c foo.c -o foo.o
#> In file included from <built-in>:1:
#> In file included from /Users/afp18axu/Library/R/arm64/4.6/library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp:22:
#> In file included from /Users/afp18axu/Library/R/arm64/4.6/library/RcppEigen/include/Eigen/Dense:1:
#> In file included from /Users/afp18axu/Library/R/arm64/4.6/library/RcppEigen/include/Eigen/Core:19:
#> /Users/afp18axu/Library/R/arm64/4.6/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:679:10: fatal error: 'cmath' file not found
#>   679 | #include <cmath>
#>       |          ^~~~~~~
#> 1 error generated.
#> make: *** [foo.o] Error 1
#> Start sampling
#> Warning: There were 39 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
d_post <- dstudy(g_brms, n = list(Task = 3, Rater = 4))
#> Warning: Low tail effective sample size (ESS < 400) detected for: Task. The g and phi coefficients cannot be trusted.

# D-study with credible intervals (requires brms backend)
d_ci <- dstudy(g_brms, n = list(Task = 3, Rater = 4), ci = c("g", "phi"))
#> Warning: Low tail effective sample size (ESS < 400) detected for: Task. The g and phi coefficients cannot be trusted.
print(d_ci)
#> Decision Study (D-Study)
#> ========================
#> 
#> Based on G Study with brms backend
#> Object of measurement: Person 
#> Universe components: Person 
#> Error components for relative error (sigma2_delta): Person:Task, Person:Rater (Residual) 
#> Error components for absolute error (sigma2_delta_abs): Task, Rater, Person:Task, Person:Rater (Residual) 
#> 
#> Sample Sizes:
#>  Task: 3
#>  Rater: 4
#> 
#> Variance Components:
#> # A tibble: 5 × 6
#>   component   dim   var_unscaled pct_unscaled var_scaled pct_scaled
#>   <chr>       <chr>        <dbl>        <dbl>      <dbl>      <dbl>
#> 1 Person      Score       0.6891       9.6062     0.6891    26.9261
#> 2 Task        Score       2.2899      31.9217     0.7633    29.8254
#> 3 Rater       Score       0.9688      13.5053     0.2422     9.4638
#> 4 Person:Task Score       0.6984       9.7358     0.2328     9.0965
#> 5 Residual    Score       2.5273      35.2311     0.6318    24.6881
#> 
#> Coefficients:
#> # A tibble: 1 × 9
#>     uni sigma2_delta sigma2_delta_abs     g   phi    g_LL  g_UL  phi_LL phi_UL
#>   <dbl>        <dbl>            <dbl> <dbl> <dbl>   <dbl> <dbl>   <dbl>  <dbl>
#> 1 0.689        0.865             1.87 0.358 0.254 0.00827 0.763 0.00391  0.655

# Custom probability levels (90% credible interval)
d_ci_90 <- dstudy(g_brms, n = list(Task = 3, Rater = 4),
  ci = "g", probs = c(0.05, 0.95))
#> Warning: Low tail effective sample size (ESS < 400) detected for: Task. The g and phi coefficients cannot be trusted.
# }
```
