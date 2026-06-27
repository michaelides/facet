# Calculate Composite Coefficients from Posterior Draws

Computes composite reliability coefficients from posterior draws of
variance components and covariances. This ensures proper uncertainty
propagation for multivariate D-studies with brms estimator.

## Usage

``` r
calculate_composite_coefficients_from_draws(
  vc_draws,
  cov_draws,
  weights,
  scale_factors,
  object_spec,
  universe_spec,
  error_spec,
  agg_facets,
  residual_is,
  cut_score = NULL,
  mu_y = NULL,
  ci = NULL,
  probs = c(0.025, 0.975),
  gstudy_data = NULL,
  dimension_var = NULL
)
```

## Arguments

- vc_draws:

  Named list (by dimension) of named lists (by component) of variance
  draws

- cov_draws:

  Named list from extract_covariance_draws()

- weights:

  Named numeric vector of weights (names = dimensions)

- scale_factors:

  Named list of scale factors per component

- object_spec:

  Parsed object specification

- universe_spec:

  Parsed universe specification

- error_spec:

  Parsed error specification (or NULL)

- agg_facets:

  Parsed aggregation specification (or NULL)

- residual_is:

  Residual composition specification

- cut_score:

  Optional cut score for phi-cut

- mu_y:

  Named list of grand means (or single value)

- ci:

  Optional CI specification (character vector)

- probs:

  Probability levels for credible intervals

## Value

List with:

- summary:

  Data frame with coefficient point estimates

- distributions:

  Named list of posterior draws for each coefficient
