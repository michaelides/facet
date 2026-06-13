# Extract Random Effect Correlation Specifications from Formula

Parses the formula to detect random effects with correlation structures.
In brms, this is specified as (1\|cor_name\|facet) where cor_name is the
name of a correlation matrix.

## Usage

``` r
extract_random_effect_cor(formula)
```

## Arguments

- formula:

  A formula object.

## Value

Named list where names are facet names and values are correlation matrix
names. For example: list(person = "cor_p", item = "cor_i")
