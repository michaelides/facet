# Fit Unbalanced Multivariate Model Using Method of Moments (Henderson's Method III)

Implements Henderson's Method III for variance component estimation in
unbalanced multivariate designs. Each dimension is analyzed separately
using its available data, and correlations are computed using pairwise
complete cases.

## Usage

``` r
fit_mom_multivariate_unbalanced(formula, data, ...)
```

## Arguments

- formula:

  A formula for the model

- data:

  A data frame containing the variables

- ...:

  Additional arguments

## Value

A momfit object with unbalanced estimation
