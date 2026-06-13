# Append Composite Coefficients if Multivariate

Conditionally appends composite coefficient rows to the coefficients
tibble when the model is multivariate with more than one dimension.

## Usage

``` r
maybe_append_composite(
  coefficients,
  vc,
  n,
  weights,
  object,
  error,
  aggregation,
  residual_is_effective,
  universe_spec,
  gstudy_obj,
  cut_score = NULL,
  mu_y = NULL,
  estimate_label = NULL
)
```

## Arguments

- coefficients:

  Existing coefficients tibble

- vc:

  Variance components tibble

- n:

  Sample sizes

- weights:

  Dimension weights

- object:

  Object of measurement name

- error:

  Error specification

- aggregation:

  Aggregation specification

- residual_is_effective:

  Effective residual specification

- universe_spec:

  Universe specification

- gstudy_obj:

  The G-study object (for correlations, data, dimension_var)

- cut_score:

  Optional cut score

- mu_y:

  Grand mean

- estimate_label:

  Optional label for the estimate column

## Value

A list with coefficients (possibly with composite row) and var_results
