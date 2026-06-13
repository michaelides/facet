# Detect if Formula is Multivariate

Checks whether a formula represents a multivariate model by looking for
brms-specific multivariate syntax (mvbind or set_rescor).

## Usage

``` r
is_multivariate(formula)
```

## Arguments

- formula:

  A formula or brmsformula object.

## Value

Logical indicating whether the formula is multivariate.
