# Generate Wide-Format Multivariate Data

Generate Wide-Format Multivariate Data

## Usage

``` r
generate_multivariate_wide(
  facets,
  vc,
  n_dims,
  dim_names,
  sd_residual,
  residual_cor,
  re_cor,
  nested,
  prefix
)
```

## Arguments

- facets:

  Named list of facet levels

- vc:

  Named list of variance components

- n_dims:

  Number of dimensions

- dim_names:

  Names of dimensions

- sd_residual:

  Residual standard deviation (or vector for heterogeneous)

- residual_cor:

  Residual correlation matrix

- re_cor:

  Named list of random effect correlation matrices

- nested:

  Named list of nesting relationships

- prefix:

  Prefix for level names

## Value

Data frame with simulated multivariate data
