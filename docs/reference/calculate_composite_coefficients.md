# Calculate Composite G and Phi Coefficients

Computes composite reliability coefficients for multivariate designs by
combining variance-covariance information across dimensions using
weights.

## Usage

``` r
calculate_composite_coefficients(
  vc,
  n,
  weights,
  object,
  error = NULL,
  aggregation = NULL,
  residual_is = NULL,
  universe = NULL,
  correlations = NULL,
  cut_score = NULL,
  mu_y = NULL,
  gstudy_data = NULL,
  dimension_var = NULL
)
```

## Arguments

- vc:

  A variance components tibble (scaled for D-study)

- n:

  Named list of sample sizes for each facet

- weights:

  Named numeric vector of weights (names = dimensions)

- object:

  The object of measurement specification

- error:

  Optional error component specification

- aggregation:

  Optional aggregation specification

- residual_is:

  Optional residual specification

- universe:

  Optional universe specification

- correlations:

  List of covariance information from the G-study

- cut_score:

  Optional cut-score for phi-cut calculation

- mu_y:

  Optional grand mean(s) for phi-cut calculation

## Value

A data.frame with one row containing:

- dim: "Composite"

- uni: Universe score variance (composite)

- sigma2_delta: Relative error variance (composite)

- sigma2_delta_abs: Absolute error variance (composite)

- g: G coefficient (composite)

- phi: Phi coefficient (composite)

- sem_rel: Standard error of measurement (relative)

- sem_abs: Standard error of measurement (absolute)

- phi_cut: Phi-cut coefficient (if cut_score provided)
