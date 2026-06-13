# Scale Variance Draws for D-Study

Applies D-study scaling to posterior draws of variance components.

## Usage

``` r
scale_variance_draws(
  vc_draws,
  n,
  object_spec,
  agg_facets,
  residual_is,
  gstudy_obj,
  n_provided = FALSE
)
```

## Arguments

- vc_draws:

  Named list of variance component draws.

- n:

  Named list of sample sizes.

- object_spec:

  Character vector of object components.

- agg_facets:

  Character vector of aggregation facets.

- residual_is:

  Character string for residual composition.

- gstudy_obj:

  Original gstudy object.

- n_provided:

  Logical indicating if n was explicitly provided.

## Value

Named list of scaled variance draws.
