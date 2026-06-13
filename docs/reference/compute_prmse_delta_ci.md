# Compute PRMSE Confidence Intervals via Delta Method

Computes confidence intervals for PRMSE metrics using the delta method,
propagating uncertainty from variance component estimates.

## Usage

``` r
compute_prmse_delta_ci(
  dstudy_obj,
  gstudy_obj,
  metrics,
  probs,
  n,
  weights,
  object,
  universe,
  error,
  aggregation,
  residual_is
)
```

## Arguments

- dstudy_obj:

  A dstudy object

- gstudy_obj:

  The associated gstudy object

- metrics:

  Character vector: "prmse", "var", or both

- probs:

  Numeric vector of length 2 for CI bounds

- n:

  Named list of D-study sample sizes

- weights:

  Named numeric vector of dimension weights

- object:

  Object of measurement specification

- universe:

  Universe components specification

- error:

  Error components specification

- aggregation:

  Aggregation specification

- residual_is:

  Residual specification

## Value

A list with CI bounds for each requested metric

## Details

The delta method approximates the variance of a function using:
Var(f(X)) = sum((df/dx_i)^2 \* Var(x_i))

For PRMSE metrics, this involves computing partial derivatives with
respect to each variance component and propagating their SEs.
