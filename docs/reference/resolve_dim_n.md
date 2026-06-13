# Resolve Per-Dimension Sample Size for One Dimension

Extracts a per-dim named list of facet -\> n from the `n_per_dim`
structure built by
[`resolve_dstudy_sample_sizes()`](https://github.com/yourorg/facet/reference/resolve_dstudy_sample_sizes.md).
Falls back to the global `n` for facets not specified per-dim. Returns a
list keyed by dim when iterating in a loop. Mirrors the per-dim override
pattern in
[`calculate_coefficients_posterior()`](https://github.com/yourorg/facet/reference/calculate_coefficients_posterior.md)
so non-posterior and posterior paths agree on what "per-dim n" means.

## Usage

``` r
resolve_dim_n(dim_name, n_per_dim, n_global)
```

## Arguments

- dim_name:

  Character: the dimension to extract.

- n_per_dim:

  The `n_per_dim` slot of
  [`resolve_dstudy_sample_sizes()`](https://github.com/yourorg/facet/reference/resolve_dstudy_sample_sizes.md)
  output.

- n_global:

  The flat `n` list (per-facet, possibly scalar) used as fallback.

## Value

A named list of facet -\> n. Returns `n_global` unchanged when no
per-dim data is available for `dim_name`.
