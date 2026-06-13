# Compute Scale Factors for Viable Recalculation

Computes scale factors for multiple variance components. Universe
components always receive a scale factor of 1 (unscaled). Non-universe
components are scaled by the product of non-object facet sample sizes in
that component.

## Usage

``` r
compute_scale_factors_for_viable(components, n, universe_spec, object_spec)
```

## Arguments

- components:

  Character vector of variance component names

- n:

  Named list of sample sizes

- universe_spec:

  Universe components specification

- object_spec:

  Object of measurement specification
