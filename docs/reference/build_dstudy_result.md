# Build D-Study Result Object

Constructs the final dstudy result list from computed components.

## Usage

``` r
build_dstudy_result(
  gstudy_obj,
  d_vc,
  coefficients,
  n,
  n_tibble,
  n_per_dim,
  object,
  universe_spec,
  error,
  aggregation,
  residual_is,
  residual_composition,
  is_sweep,
  estimation,
  posterior,
  composite_post,
  var_results,
  is_multivariate,
  cut_score,
  mu_y,
  ci,
  probs,
  weights
)
```

## Arguments

- gstudy_obj:

  The G-study object

- d_vc:

  D-study variance components

- coefficients:

  Coefficients tibble

- n:

  Sample sizes

- n_tibble:

  Per-dimension sample sizes tibble

- n_per_dim:

  Per-dimension sample sizes list

- object:

  Object of measurement

- universe_spec:

  Universe specification

- error:

  Error specification

- aggregation:

  Aggregation specification

- residual_is:

  Original residual_is parameter

- residual_composition:

  Residual composition

- is_sweep:

  Whether this is a sweep

- estimation:

  Estimation method used

- posterior:

  Posterior draws (NULL if simple)

- composite_post:

  Composite posterior draws

- var_results:

  VAR results

- is_multivariate:

  Whether multivariate

- cut_score:

  Cut score (NULL if not provided)

- mu_y:

  Grand mean (NULL if no cut_score)

- ci:

  CI specification

- probs:

  Probability levels

- weights:

  Dimension weights

## Value

A dstudy object
