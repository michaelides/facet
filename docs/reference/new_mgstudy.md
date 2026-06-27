# Create an mgstudy Object

Constructor function for the "mgstudy" class (multivariate G-study).
This is an internal low-level constructor. Users should use
[`gstudy()`](https://github.com/yourorg/facet/reference/gstudy.md) to
create mgstudy objects.

## Usage

``` r
new_mgstudy(
  model,
  variance_components,
  facets,
  object,
  estimator,
  formula,
  data,
  dimensions
)
```

## Arguments

- model:

  The fitted model object.

- variance_components:

  A tibble of variance components with 'dim' column.

- facets:

  Character vector of facet names.

- object:

  Character string naming the object of measurement.

- estimator:

  Character string indicating the estimator used.

- formula:

  The formula used for fitting.

- data:

  The original data.

- dimensions:

  Character vector of dimension (response variable) names.

## Value

An object of class "mgstudy".
