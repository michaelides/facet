# Calculate G and D Coefficients

Computes generalizability (G) and dependability (D/Phi) coefficients
from D-study variance components.

## Usage

``` r
calculate_coefficients(
  vc,
  n = NULL,
  object = NULL,
  error = NULL,
  aggregation = NULL,
  residual_is = NULL,
  universe = NULL,
  cut_score = NULL,
  mu_y = NULL
)
```

## Arguments

- vc:

  D-study variance components tibble.

- n:

  Named list of sample sizes for each facet (required for aggregation).

- object:

  Specification for object of measurement. Can be:

  - NULL (default): uses first component in vc

  - A character string: "p"

  - A character vector: c("p", "p:d")

  - A formula: obj ~ p + p:d

- error:

  Specification for error components. Can be:

  - NULL (default): all components not in object specification

  - A character string: "p:i"

  - A character vector: c("p:i", "p:r")

  - A formula: ~ p:i + p:r Note: If the same facet is specified in both
    `object` and `error`, an error is raised.

- aggregation:

  Character vector of facets to aggregate over. Components containing
  these facets will be divided by the sample size. Default is NULL (no
  aggregation). Note: If the same facet is specified in both `object`
  and `aggregation`, an error is raised.

- residual_is:

  Character string specifying which facets make up the residual. For
  example, "p:i:r" means the residual represents the p:i:r interaction.
  Default is NULL (residual is not rescaled).

## Value

A data frame with columns:

- uni:

  Universe score variance. Represents the true individual differences in
  the trait being measured. It reflects the systematic variance
  attributable to the object of measurement (e.g., persons) after
  accounting for measurement error. In classical test theory terms, this
  corresponds to the "true score" variance. A larger universe score
  variance indicates greater heterogeneity among objects of measurement,
  which generally improves the precision of relative comparisons.

- sigma2_delta:

  Relative error variance. Error variance affecting relative decisions
  (rankings). Includes only variance components that interact with the
  object of measurement (i.e., components that contain the object
  facet). For example, in a person x rater design with person as object,
  this includes person:rater and residual components.

- sigma2_delta_abs:

  Absolute error variance. Error variance affecting absolute decisions.
  Includes ALL variance components except the object of measurement. For
  example, in a person x rater design with person as object, this
  includes rater, person:rater, and residual components.

- g:

  Generalizability coefficient. Quantifies the reliability of
  measurements for relative decisions—situations where the focus is on
  rank-ordering or comparing objects of measurement (e.g., identifying
  high vs. low performers). Calculated as: g = uni / (uni +
  sigma2_delta). The G coefficient ranges from 0 to 1, with higher
  values indicating better reliability. Values \>= 0.70 are typically
  considered acceptable for group-level decisions, while \>= 0.80 is
  preferred for individual-level decisions.

- phi:

  Dependability coefficient (Phi coefficient). Quantifies the
  reliability of measurements for absolute decisions—situations where
  the focus is on the absolute level of a measurement (e.g., determining
  if a score exceeds a cutoff). Calculated as: phi = uni / (uni +
  sigma2_delta_abs). The Phi coefficient is always \<= G coefficient
  because absolute error includes all error sources, not just those
  affecting relative rankings.

- sem_rel:

  Standard error of measurement for relative decisions. Calculated as
  sqrt(sigma2_delta). Indicates the typical error when interpreting
  relative standing or rank-order among objects of measurement.

- sem_abs:

  Standard error of measurement for absolute decisions. Calculated as
  sqrt(sigma2_delta_abs). Indicates the typical error when interpreting
  absolute score levels against criterion-referenced standards.
