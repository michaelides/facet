# Fit a Model Using lme4

Fits a mixed effects model using lme4::lmer(). Converts brms-style
formulas to lme4 format if necessary.

## Usage

``` r
fit_lme4(formula, data, ...)
```

## Arguments

- formula:

  A formula for the model (lme4 style).

- data:

  A data frame.

- ...:

  Additional arguments passed to
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).

## Value

A fitted lmerMod object.
