# Extract Posterior Draws from gstudy or dstudy Objects

Extracts pre-calculated posterior draws from gstudy or dstudy objects
and returns them as a wide-format tibble where each variance component
or coefficient is a separate column.

## Usage

``` r
gt_draws(object, ...)

# S3 method for class 'gstudy'
gt_draws(object, components = NULL, dims = NULL, n = NULL, ...)

# S3 method for class 'mgstudy'
gt_draws(object, components = NULL, dims = NULL, n = NULL, ...)

# S3 method for class 'dstudy'
gt_draws(
  object,
  what = c("all", "coefficients", "composite", "var", "variance"),
  coefficients = NULL,
  dims = NULL,
  ...
)
```

## Arguments

- object:

  A gstudy or dstudy object

- ...:

  Additional arguments passed to methods

- components:

  Character vector of variance components to extract. If NULL (default),
  all components are extracted.

- dims:

  Character vector of dimensions to extract (for multivariate models).
  If NULL (default), all dimensions are extracted.

- n:

  Named list of sample sizes for D-study calculation. If NULL (default),
  sample sizes are extracted from the G-study data.

- what:

  Character string specifying what to extract:

  "all"

  :   Extract all available draws (default)

  "coefficients"

  :   Extract coefficient draws (uni, sigma2_delta, g, phi, etc.)

  "composite"

  :   Extract composite posterior draws (multivariate only)

  "var"

  :   Extract VAR and PRMSE draws (multivariate only)

  "variance"

  :   Extract variance component draws from underlying gstudy

- coefficients:

  Character vector of coefficients to extract. If NULL (default), all
  coefficients are extracted.

## Value

A tibble or named list of tibbles with posterior draws in wide format

For univariate gstudy objects, a tibble with columns:

- `draw`: Integer draw index

- Variance component columns (e.g., `person`, `item`, `Residual`)

- Coefficient columns: `uni`, `sigma2_delta`, `sigma2_delta_abs`, `g`,
  `phi`, `sem_rel`, `sem_abs`

For multivariate gstudy objects (mgstudy), a named list of tibbles, one
per dimension, each with the same column structure as univariate.

For univariate dstudy objects, a tibble with columns:

- `draw`: Integer draw index

- Coefficient columns: `uni`, `sigma2_delta`, `sigma2_delta_abs`, `g`,
  `phi`, `sem_rel`, `sem_abs`

For multivariate dstudy objects, a named list of tibbles:

- One tibble per dimension (with coefficients and VAR columns merged)

- `$composite`: Composite posterior draws (if available)

## Details

For gstudy objects, the function extracts:

- Variance component draws for each random effect and residual

- Derived coefficients (uni, sigma2_delta, g, phi, etc.) calculated
  using the sample sizes from the original G-study or provided via `n`

For dstudy objects, the function extracts pre-calculated posterior draws
that were computed when `estimation = "posterior"` was specified in the
[`dstudy()`](https://github.com/yourorg/facet/reference/dstudy.md) call.

For multivariate dstudy objects, VAR/PRMSE columns are merged with
coefficient columns in each dimension's tibble.

## Examples

``` r
if (FALSE) { # \dontrun{
# Univariate gstudy
g <- gstudy(score ~ (1 | person) + (1 | item), data = my_data, backend = "brms",
  iter = 2000, cores = 4, refresh = 1000)
draws <- gt_draws(g)

# Filter specific components
draws <- gt_draws(g, components = c("person", "Residual"))

# Multivariate gstudy - returns named list
mg <- gstudy(cbind(score1, score2) ~ (1 | person) + (1 | item),
  data = my_data, backend = "brms", iter = 2000, cores = 4, refresh = 1000)
draws <- gt_draws(mg)
draws$score1  # Access dimension-specific draws
draws$score2

# Filter to specific dimensions
draws <- gt_draws(mg, dims = "score1")
} # }
if (FALSE) { # \dontrun{
# D-study with posterior estimation
g <- gstudy(score ~ (1 | person) + (1 | item), data = my_data, backend = "brms",
  iter = 2000, cores = 4, refresh = 1000)
d <- dstudy(g, n = list(item = 10), estimation = "posterior")

# Extract all draws
draws <- gt_draws(d)

# Extract only coefficient draws
draws <- gt_draws(d, what = "coefficients")

# Multivariate - returns named list
mg <- gstudy(cbind(score1, score2) ~ (1 | person) + (1 | item),
  data = my_data, backend = "brms", iter = 2000, cores = 4, refresh = 1000)
md <- dstudy(mg, n = list(item = 10), estimation = "posterior")
draws <- gt_draws(md)
draws$score1  # Per-dimension draws with VAR columns
draws$composite  # Composite draws
} # }
```
