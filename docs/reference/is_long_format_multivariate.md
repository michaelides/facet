# Detect Long-Format Multivariate Models

Checks if a brms formula represents a long-format multivariate model,
characterized by an auxiliary sigma formula and dimension variable in
random effects.

## Usage

``` r
is_long_format_multivariate(formula)
```

## Arguments

- formula:

  A formula or brmsformula object.

## Value

A list with:

- is_long: Logical indicating if this is a long-format multivariate
  model

- dimension_var: Character name of the dimension variable, or NULL
