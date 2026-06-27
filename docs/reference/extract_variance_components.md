# Extract Variance Components from a Fitted Model

Generic function that dispatches to the appropriate estimator-specific
extractor.

## Usage

``` r
extract_variance_components(
  model,
  estimator,
  ci_method = "none",
  nsim = 1000,
  boot.type = "perc",
  formula = NULL,
  ...
)
```

## Arguments

- model:

  A fitted model object.

- estimator:

  Character string indicating the estimator used.

- ci_method:

  Character string specifying the CI method for lme4 estimator. One of
  "none", "profile", or "boot". Default is "none".

- nsim:

  Integer: number of bootstrap simulations (only for ci_method =
  "boot"). Default is 1000.

- boot.type:

  Character: bootstrap type, "perc" or "basic" (only for ci_method =
  "boot"). Default is "perc".

- ...:

  Additional arguments passed to confint.merMod.

## Value

A tibble of variance components with columns:

- component:

  Name of the variance component

- facet:

  Associated facet name

- type:

  Type: "main", "interaction", or "residual"

- var:

  Point estimate of variance

- pct:

  Percentage of total variance

For lme4 with ci_method != "none" or brms estimator, also includes:

- lower:

  Lower confidence interval bound

- upper:

  Upper confidence interval bound

For brms estimator, also includes:

- error:

  Standard error of the estimate

- sd:

  Standard deviation

- Rhat:

  Convergence diagnostic

- ESS:

  Effective sample size
