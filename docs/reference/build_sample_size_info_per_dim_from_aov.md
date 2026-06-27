# Build Per-Dimension Sample Size Info from a mom Unbalanced Fit

For wide-format aov fits with `unbalanced = TRUE`, the per-response N is
stored as a named list on the aovfit object. This adapter reconstructs a
`sample_size_info` object for each dimension by re-running
[`calculate_sample_size_info()`](https://github.com/yourorg/facet/reference/calculate_sample_size_info.md)
on the non-NA subset of the original data for that response. The result
has the same shape as the long-format `sample_size_info_per_dim`, so
downstream code can consume it uniformly.

## Usage

``` r
build_sample_size_info_per_dim_from_aov(aov_model)
```

## Arguments

- aov_model:

  An `aovfit` object produced by
  [`fit_aov_multivariate_unbalanced()`](https://github.com/yourorg/facet/reference/fit_aov_multivariate_unbalanced.md).
  Must have `$n_per_dim`, `$data`, `$formula`, and `$responses`.

## Value

A named list of `sample_size_info` objects, one per dimension. Returns
`NULL` if `aov_model` is not an unbalanced aovfit.
