# Select the Appropriate Backend

Determines which backend to use based on formula characteristics and
user preference. If backend is "auto", selects brms for multivariate
formulas and lme4 for univariate.

## Usage

``` r
select_backend(formula, backend = c("auto", "lme4", "brms", "mom"))
```

## Arguments

- formula:

  The model formula.

- backend:

  Character string specifying the desired backend: "auto", "lme4", or
  "brms".

## Value

Character string indicating the backend to use.
