# Extract Response Variable Names from Formula

Extracts the names of response (dependent) variables from a formula.
Works for both univariate and multivariate (brms mvbind) formulas.

## Usage

``` r
extract_response_names(formula)
```

## Arguments

- formula:

  A formula or brmsformula object.

## Value

Character vector of response variable names. For univariate formulas,
returns a single name. For multivariate formulas, returns all response
names.
