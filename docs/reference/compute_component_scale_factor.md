# Compute Scale Factor for a Component

Computes the scale factor for a variance component. Used in composite
variance calculations.

## Usage

``` r
compute_component_scale_factor(comp, n, object_spec, n_provided)
```

## Arguments

- comp:

  The variance component name.

- n:

  Named list of sample sizes.

- object_spec:

  Character vector of object components.

- n_provided:

  Logical indicating if n was explicitly provided.

## Value

The scale factor for this component.
