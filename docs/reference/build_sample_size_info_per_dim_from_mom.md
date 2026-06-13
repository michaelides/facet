# Build Per-Dimension Sample Size Info from a mom Unbalanced Fit

For wide-format mom fits with `unbalanced = TRUE`, the per-response N is
stored as a named list on the momfit object. This adapter reconstructs a
`sample_size_info` object for each dimension by re-running
[`calculate_sample_size_info()`](https://github.com/yourorg/facet/reference/calculate_sample_size_info.md)
on the non-NA subset of the original data for that response. The result
has the same shape as the long-format `sample_size_info_per_dim`, so
downstream code can consume it uniformly.

## Usage

``` r
build_sample_size_info_per_dim_from_mom(mom_model)
```

## Arguments

- mom_model:

  A `momfit` object produced by
  [`fit_mom_multivariate_unbalanced()`](https://github.com/yourorg/facet/reference/fit_mom_multivariate_unbalanced.md).
  Must have `$n_per_dim`, `$data`, `$formula`, and `$responses`.

## Value

A named list of `sample_size_info` objects, one per dimension. Returns
`NULL` if `mom_model` is not an unbalanced momfit.
