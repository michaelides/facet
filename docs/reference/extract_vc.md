# Extract Variance Components from a Model

Generic function for extracting variance components from different model
types.

## Usage

``` r
extract_vc(model, ...)

# S3 method for class 'merMod'
extract_vc(model, ...)

# S3 method for class 'brmsfit'
extract_vc(model, ...)
```

## Arguments

- model:

  A fitted model object.

- ...:

  Additional arguments passed to methods.

## Value

A tibble with columns: component, variance, percent.
