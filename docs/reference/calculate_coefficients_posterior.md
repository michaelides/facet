# Calculate Coefficients Using Posterior Draws

Computes generalizability and dependability coefficients using full
posterior distributions from a brms model. This provides proper
uncertainty quantification for D-study coefficients.

## Usage

``` r
calculate_coefficients_posterior(
  gstudy_obj,
  n,
  object = NULL,
  universe = NULL,
  error = NULL,
  aggregation = NULL,
  residual_is = NULL,
  is_sweep = FALSE,
  n_grid = NULL,
  n_provided = FALSE,
  use_scaled = TRUE,
  cut_score = NULL,
  mu_y = NULL,
  ci = NULL,
  probs = c(0.025, 0.975),
  weights = NULL,
  n_per_dim = NULL
)
```

## Arguments

- gstudy_obj:

  A gstudy object fitted with estimator = "brms".

- n:

  Named list of sample sizes for each facet.

- object:

  Specification for object of measurement.

- error:

  Specification for error components. Note: If the same facet is
  specified in both `object` and `error`, an error is raised.

- aggregation:

  Character vector of facets to aggregate over. Note: If the same facet
  is specified in both `object` and `aggregation`, an error is raised.

- residual_is:

  Character string specifying residual composition.

- is_sweep:

  Logical indicating if this is a sweep over multiple sample sizes.

- n_grid:

  If is_sweep, a data frame with all combinations of sample sizes.

- cut_score:

  Optional cut score for phi-cut calculation.

- mu_y:

  Optional grand mean for phi-cut calculation.

## Value

A list with:

- coefficients:

  Data frame with coefficient means (and sweep combinations if
  applicable)

- posterior:

  List of posterior distribution vectors (or list of lists for sweep)
