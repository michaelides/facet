# Extract Variance Components from brms Model

Extracts variance components from a brms model using posterior samples.
Computes variance estimates as mean(SD^2) from posterior draws, ensuring
proper uncertainty quantification and avoiding Jensen's inequality bias
that would occur with mean(SD)^2.

## Usage

``` r
extract_vc_brms(model, conf_level = 0.95, formula = NULL)
```

## Arguments

- model:

  A brmsfit object from brms::brm().

- conf_level:

  Credible interval level (default 0.95).

- formula:

  The original formula (used to extract response names).

## Value

A tibble with columns:

- component: Name of the variance component

- dim: Response variable name (for multivariate models)

- type: "main", "interaction", or "residual"

- var: Mean of variance draws (mean(SD^2))

- error: Standard error of variance (sd(variance_draws))

- lower: 2.5th percentile of variance draws

- upper: 97.5th percentile of variance draws

- sd: Mean of SD draws (for display purposes)

- Rhat: Convergence diagnostic

- Bulk_ESS: Effective sample size (bulk)

- Tail_ESS: Effective sample size (tail)

- pct: Percentage of total variance
