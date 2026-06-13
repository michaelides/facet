# Extract Variance Components from brms Long-Format Model

Extracts variance components for long-format multivariate models where
dimensions are defined by a dimension variable.

## Usage

``` r
extract_vc_brms_long_format(
  model,
  conf_level = 0.95,
  formula = NULL,
  dimension_var = NULL,
  data = NULL
)
```

## Arguments

- model:

  A brmsfit object.

- conf_level:

  Credible interval level (default 0.95).

- formula:

  The original formula.

- dimension_var:

  Name of the dimension variable.

- data:

  The original data.

## Value

A tibble with variance components per dimension.
