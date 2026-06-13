# Create a dstudy Object

Constructor function for the "dstudy" class. This is an internal
low-level constructor. Users should use
[`dstudy()`](https://github.com/yourorg/facet/reference/dstudy.md) to
create dstudy objects.

## Usage

``` r
new_dstudy(
  gstudy,
  variance_components,
  coefficients,
  n,
  object,
  is_sweep = FALSE
)
```

## Arguments

- gstudy:

  The original G-study object.

- variance_components:

  A tibble of variance components for the D-study.

- coefficients:

  A tibble with G and D coefficients.

- n:

  A named list of facet levels.

- object:

  Character string naming the object of measurement.

- is_sweep:

  Logical indicating if this is a sweep over multiple sample sizes.

## Value

An object of class "dstudy".
