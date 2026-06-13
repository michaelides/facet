# Simulate Data with Individual SD Parameters (Legacy Style)

A convenience wrapper that uses individual `sd_*` parameters similar to
the `makedata()` function style. Useful for complex designs with many
variance components.

## Usage

``` r
simulate_gtheory_legacy(
  n_p = 200,
  n_i = 25,
  n_d = 5,
  n_l = 20,
  sd_res = 1,
  sd_p = 1,
  sd_i = 0.5,
  sd_d = 0.3,
  sd_l = 0.5,
  sd_pd = 0.3,
  sd_li = 0.3,
  sd_ld = 0.3,
  multivariate = FALSE,
  n_dims = 2,
  seed = NULL
)
```

## Arguments

- n_p:

  Number of persons

- n_i:

  Number of items

- n_d:

  Number of domains (optional)

- n_l:

  Number of levels (optional)

- sd_res:

  Residual standard deviation

- sd_p:

  Person standard deviation

- sd_i:

  Item standard deviation

- sd_d:

  Domain standard deviation

- sd_l:

  Level standard deviation

- sd_pd:

  Person-by-domain interaction SD

- sd_li:

  Level-by-item interaction SD

- sd_ld:

  Level-by-domain interaction SD

- multivariate:

  Logical. If TRUE, generates multivariate data

- n_dims:

  Number of dimensions for multivariate data

- seed:

  Random seed

## Value

A data frame with simulated data

## Examples

``` r
# Simulate data similar to the original makedata() example
data <- simulate_gtheory_legacy(
  n_p = 200, n_i = 25, n_d = 5, n_l = 20,
  sd_res = 1, sd_p = 1, sd_i = .5, sd_d = .3, sd_l = .5,
  sd_pd = .3, sd_li = .3, sd_ld = .3
)
```
