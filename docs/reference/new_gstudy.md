# Create a gstudy Object

Constructor function for the "gstudy" class. This is an internal
low-level constructor. Users should use
[`gstudy()`](https://github.com/yourorg/facet/reference/gstudy.md) to
create gstudy objects.

## Usage

``` r
new_gstudy(
  model,
  variance_components,
  facets,
  object,
  estimator,
  is_multivariate,
  formula,
  data
)
```

## Arguments

- model:

  The fitted model object.

- variance_components:

  A tibble of variance components.

- facets:

  Character vector of facet names.

- object:

  Character string naming the object of measurement.

- estimator:

  Character string indicating the estimator used.

- is_multivariate:

  Logical indicating if the model is multivariate.

- formula:

  The formula used for fitting.

- data:

  The original data.

## Value

An object of class "gstudy".
