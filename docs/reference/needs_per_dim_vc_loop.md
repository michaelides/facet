# Should the D-Study Path Iterate Per Dimension?

Returns TRUE when the mgstudy object has per-dim sample-size information
and the user-provided `n` is a per-dim tibble (i.e., we should call the
variance functions separately for each dim rather than once with scalar
n).

## Usage

``` r
needs_per_dim_vc_loop(is_multivariate, n_per_dim, n_tibble)
```

## Arguments

- is_multivariate:

  Logical: whether the model is multivariate.

- n_per_dim:

  The `n_per_dim` slot from
  [`resolve_dstudy_sample_sizes()`](https://github.com/yourorg/facet/reference/resolve_dstudy_sample_sizes.md).

- n_tibble:

  The `n_tibble` slot from
  [`resolve_dstudy_sample_sizes()`](https://github.com/yourorg/facet/reference/resolve_dstudy_sample_sizes.md).

## Value

Logical scalar.
