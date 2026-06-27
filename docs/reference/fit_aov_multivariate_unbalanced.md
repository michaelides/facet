# Fit Unbalanced Multivariate Model Using ANOVA-based Estimation (Henderson's Method III)

Implements Henderson's Method III for variance component estimation in
unbalanced multivariate designs. Each dimension is analyzed separately
using its available data, and correlations are computed using pairwise
complete cases.

## Usage

``` r
fit_aov_multivariate_unbalanced(formula, data, nested = NULL, ...)
```

## Arguments

- formula:

  A formula for the model

- data:

  A data frame containing the variables

- ...:

  Additional arguments

## Value

A aovfit object with unbalanced estimation
