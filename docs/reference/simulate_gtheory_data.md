# Simulate Data for Generalizability Theory Analysis

Generate simulated data for univariate and multivariate generalizability
theory (G-theory) analyses. Supports both wide-format and long-format
multivariate designs, with optional correlation structures for residual
and random effects.

## Usage

``` r
simulate_gtheory_data(
  facets = NULL,
  formula = NULL,
  vc = NULL,
  sd_residual = 1,
  multivariate = FALSE,
  n_dims = 2,
  dim_names = NULL,
  format = c("wide", "long"),
  dim_var = "dimension",
  response_var = NULL,
  residual_cor = NULL,
  re_cor = NULL,
  nested = NULL,
  seed = NULL,
  prefix = ""
)
```

## Arguments

- facets:

  A named list specifying the number of levels for each facet. For
  example, `list(person = 50, item = 10)` creates a crossed design with
  50 persons and 10 items. Names should be valid R variable names.

- formula:

  An optional lme4-style formula from which to derive the facet
  structure. If provided, overrides `facets`. The formula should use
  random effects syntax, e.g., `y ~ (1|person) + (1|item)`.

- vc:

  A named list of variance components (as variances, not SDs). Names
  should match facet names or interactions (use `:` for interactions).
  For example, `list(person = 1, item = 0.5, "person:item" = 0.3)`.

- sd_residual:

  Standard deviation of the residual error term. Default is 1.

- multivariate:

  Logical. If TRUE, generates multivariate data. Default is FALSE.

- n_dims:

  Number of dimensions (outcome variables) for multivariate data.
  Default is 2. Ignored if `multivariate = FALSE`.

- dim_names:

  Character vector of dimension names. If NULL, dimensions are named
  `y1`, `y2`, etc. for wide format.

- format:

  Character string specifying the output format: "wide" or "long". Wide
  format has multiple outcome columns; long format has a dimension
  indicator variable. Default is "wide".

- dim_var:

  Name of the dimension indicator variable for long format. Default is
  "dimension".

- response_var:

  Name of the response variable. Default is "y" for univariate and
  "score" for multivariate long format.

- residual_cor:

  A correlation matrix for residual errors across dimensions. Must be a
  symmetric positive-definite matrix with dimensions matching `n_dims`.
  NULL (default) implies independent residuals (identity matrix).

- re_cor:

  A named list of correlation matrices for random effects. Names should
  match facet names. Each matrix should be symmetric positive-definite
  with dimensions matching `n_dims`. NULL (default) implies independent
  random effects.

- nested:

  A named list specifying nesting relationships. Names are nested
  facets, values are the facets they are nested within. For example,
  `list(item = "person")` means items are nested within persons.

- seed:

  Optional random seed for reproducibility.

- prefix:

  Prefix for factor level names. Default is empty string, producing
  levels like "1", "2", etc.

## Value

A data frame with:

- Columns for each facet (as factors)

- One or more outcome columns (univariate: `y`, wide multivariate: `y1`,
  `y2`, ..., long multivariate: `score` with `dimension` column)

- For long format: a dimension indicator column

## Details

### Variance Component Specification

Variance components can be specified using either the `vc` parameter or
individual `sd_*` parameters. The `vc` parameter takes a named list
where names correspond to facet names or interactions (use `:` to
separate).

### Univariate Simulation

For univariate designs, the function generates data according to:
\$\$y\_{ijk...} = \mu + \epsilon_p + \epsilon_i + \epsilon\_{pi} +
\epsilon\_{residual}\$\$ where each \\\epsilon\\ term is drawn from a
normal distribution with the specified variance component.

### Multivariate Simulation

For multivariate designs, correlations between dimensions can be
specified:

- **Residual correlations**: Correlations between residual errors across
  dimensions

- **Random effect correlations**: Correlations between random effects
  for a specific facet

Correlations are implemented using Cholesky decomposition to ensure
proper multivariate normal distributions.

### Format Options

- **Wide format**: Multiple outcome columns (`y1`, `y2`, ...) - suitable
  for [`mvbind()`](https://github.com/yourorg/facet/reference/mvbind.md)
  syntax in
  [`gstudy()`](https://github.com/yourorg/facet/reference/gstudy.md)

- **Long format**: Single outcome column with dimension indicator -
  suitable for brms long-format models like
  `score ~ 0 + dimension + (0+dimension|facet)`

## See also

[`gstudy()`](https://github.com/yourorg/facet/reference/gstudy.md) for
conducting G-studies on the simulated data

## Examples

``` r
# Univariate p x i design
data_uni <- simulate_gtheory_data(
  facets = list(person = 20, item = 10),
  vc = list(person = 1, item = 0.5, "person:item" = 0.3),
  sd_residual = 1
)

# Wide-format multivariate (2 dimensions)
data_wide <- simulate_gtheory_data(
  facets = list(person = 50, item = 10),
  vc = list(person = 1, item = 0.5),
  multivariate = TRUE,
  n_dims = 2,
  format = "wide"
)

# Long-format multivariate with correlations
cor_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
data_long <- simulate_gtheory_data(
  facets = list(person = 50, item = 10),
  vc = list(person = 1, item = 0.5),
  multivariate = TRUE,
  n_dims = 2,
  format = "long",
  residual_cor = cor_matrix
)
```
