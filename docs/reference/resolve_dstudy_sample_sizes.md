# Resolve D-Study Sample Sizes

Determines sample sizes from multiple possible sources: user-provided n,
per-dimension sample sizes, or extracted from the G-study. Also
determines whether this is a sweep (multiple sample sizes per facet).

## Usage

``` r
resolve_dstudy_sample_sizes(gstudy_obj, n = list(), is_multivariate = FALSE)
```

## Arguments

- gstudy_obj:

  A gstudy or mgstudy object

- n:

  User-provided sample sizes

- is_multivariate:

  Whether this is a multivariate model

## Value

A list with n, n_provided, n_per_dim, n_tibble, is_sweep
