# Extract Variance Components from lme4 Model

Extracts variance components from an lme4 model using VarCorr.
Calculates confidence intervals using confint.merMod.

## Usage

``` r
extract_vc_lme4(
  model,
  ci_method = "none",
  conf_level = 0.95,
  nsim = 1000,
  boot.type = "perc",
  formula = NULL,
  ...
)
```

## Arguments

- model:

  An lmerMod object from lme4::lmer().

- ci_method:

  Method for confidence intervals: "none", "profile", or "boot".

  - "none": No confidence intervals (default).

  - "profile": Profile likelihood (more accurate, slower).

  - "boot": Parametric bootstrap (most accurate, slowest).

- conf_level:

  Confidence level for intervals (default 0.95).

- nsim:

  Number of bootstrap simulations (if ci_method = "boot", default 1000).

- boot.type:

  Bootstrap type: "perc" (percentile), "basic", or "norm"
  (normal-theory) (default "perc").

- formula:

  The original formula (used to extract response name).

- ...:

  Additional arguments passed to confint.merMod.

## Value

A tibble with columns: component, dim, type, var, pct. If ci_method !=
"none", also includes: lower, upper.
