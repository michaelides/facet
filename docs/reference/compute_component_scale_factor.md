# Compute Scale Factor for a Component

Computes the scale factor for a variance component. Used in composite
variance calculations. The object of measurement is always excluded from
the divisor (per generalizability theory). For the residual,
`residual_is` is used to decompose the residual into its constituent
facets before computing the divisor.

## Usage

``` r
compute_component_scale_factor(
  comp,
  n,
  object_spec,
  n_provided,
  residual_is = NULL
)
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

- residual_is:

  Character string specifying residual composition (e.g.,
  "Person:Rater"). Used only when `comp == "Residual"`.

## Value

The scale factor for this component.
