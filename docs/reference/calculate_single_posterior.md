# Calculate Posterior Coefficients for a Single Sample Size

Calculate Posterior Coefficients for a Single Sample Size

## Usage

``` r
calculate_single_posterior(
  vc_draws,
  n,
  object_spec,
  universe_spec,
  error_spec,
  agg_facets,
  residual_is,
  gstudy_obj,
  n_provided = FALSE,
  use_scaled = TRUE,
  cut_score = NULL,
  mu_y = NULL,
  ci = NULL,
  probs = c(0.025, 0.975)
)
```

## Arguments

- vc_draws:

  Named list of variance component draws.

- n:

  Named list of sample sizes.

- object_spec:

  Character vector of object components.

- universe_spec:

  Character vector of universe components.

- error_spec:

  Character vector of error components (or NULL).

- agg_facets:

  Character vector of aggregation facets (or NULL).

- residual_is:

  Character string for residual composition.

- gstudy_obj:

  Original gstudy object.

- n_provided:

  Logical indicating if n was explicitly provided.

- use_scaled:

  Logical indicating whether to use scaled variance draws. When FALSE,
  uses original variance draws for both universe and error. When TRUE,
  uses scaled variance draws for error components only; universe
  components always use unscaled (G-study) variance draws.

- cut_score:

  Optional cut score for phi-cut calculation.

- mu_y:

  Optional grand mean for phi-cut calculation.

## Value

List with 'summary' (data frame of means) and 'distributions' (list of
vectors).
