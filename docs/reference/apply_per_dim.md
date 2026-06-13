# Apply a Variance Function Per Dimension

When
[`needs_per_dim_vc_loop()`](https://github.com/yourorg/facet/reference/needs_per_dim_vc_loop.md)
is TRUE, splits the variance components tibble by `dim` and calls
`fn(vc_dim, n_dim, ...)` for each dimension, then row-binds the results.
Otherwise calls `fn(vc, n, ...)` once with the full vc and global n.
Used to give the non-posterior D-study paths per-dim awareness without
changing their signatures.

## Usage

``` r
apply_per_dim(vc, n, n_per_dim, n_tibble, fn, ...)
```

## Arguments

- vc:

  Variance components tibble with a `dim` column.

- n:

  Named list of facet -\> n (global n; per-dim n overrides apply).

- n_per_dim:

  The `n_per_dim` slot from
  [`resolve_dstudy_sample_sizes()`](https://github.com/yourorg/facet/reference/resolve_dstudy_sample_sizes.md).

- n_tibble:

  The `n_tibble` slot from
  [`resolve_dstudy_sample_sizes()`](https://github.com/yourorg/facet/reference/resolve_dstudy_sample_sizes.md).

- fn:

  A variance function such as
  [`calculate_divided_variance()`](https://github.com/yourorg/facet/reference/calculate_divided_variance.md)
  or
  [`calculate_dstudy_variance()`](https://github.com/yourorg/facet/reference/calculate_dstudy_variance.md).

- ...:

  Additional arguments passed to `fn`.

## Value

Tibble of variance components (per-dim, row-bound).
