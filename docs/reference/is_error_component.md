# Check if Component is Part of Error Specification

Determines whether a variance component is part of the error
specification.

## Usage

``` r
is_error_component(component, error_spec, object_spec, vc_all)
```

## Arguments

- component:

  Character string naming the variance component.

- error_spec:

  Character vector of error component names, or NULL for default.

- object_spec:

  Character vector of object component names.

- vc_all:

  Variance components tibble (for default error calculation).

## Value

TRUE if the component is part of the error specification.
