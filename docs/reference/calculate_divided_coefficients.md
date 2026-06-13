# Calculate Coefficients Using Divided Variance Components

Computes coefficients using the "divided" variance estimates where each
component is divided only by the sample sizes of non-object facets.

## Usage

``` r
calculate_divided_coefficients(
  vc_divided,
  object,
  error = NULL,
  aggregation = NULL,
  residual_is = NULL,
  universe = NULL,
  cut_score = NULL,
  mu_y = NULL
)
```

## Arguments

- vc_divided:

  Variance components tibble with 'var_divided' column

- object:

  Specification for object of measurement

- error:

  Specification for error components (optional)

- aggregation:

  Aggregation facets (optional)

- residual_is:

  Residual composition specification

- universe:

  Universe specification (optional)

## Value

Data frame with coefficient estimates
