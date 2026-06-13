# Detect Method of Moments Estimation Issues

Checks for truncated variance estimates (when MS \< residual MS,
indicating negative variance estimates that were truncated to zero).

## Usage

``` r
detect_mom_issues(model)
```

## Arguments

- model:

  A momfit object from fit_mom().

## Value

A list with components:

- truncated_variance:

  Logical indicating if variance was truncated

- truncated_components:

  Character vector of components with truncated variance
