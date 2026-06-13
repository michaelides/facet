# Compute PRMSE Confidence Intervals via Bootstrap

Computes confidence intervals for PRMSE metrics using parametric
bootstrap, resampling variance components from their estimated
distributions.

## Usage

``` r
compute_prmse_bootstrap_ci(
  dstudy_obj,
  gstudy_obj,
  metrics,
  probs,
  n_bootstrap,
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

- n_bootstrap:

  Number of bootstrap samples

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

For each bootstrap iteration:

1.  Resample variance components from N(VC, SE^2)

2.  Recompute PRMSE metrics

3.  Use percentiles for CI bounds
