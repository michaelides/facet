# Calculate Divided Variance Components

Divides each variance component by the sample sizes of non-object facets
that appear in that component. This provides an alternative scaling
where:

- Object component is NOT scaled (same as standard D-study)

- For non-object components, divides only by n for non-object facets in
  component

- For Residual, divides by product of all non-object facet sample sizes

## Usage

``` r
calculate_divided_variance(vc, n, object, residual_is = NULL)
```

## Arguments

- vc:

  Variance components tibble (unscaled, from G-study)

- n:

  Named list of sample sizes for each facet

- object:

  Specification for object of measurement

- residual_is:

  Character string specifying which facets make up the residual

## Value

Tibble with additional 'var_divided' column

## Details

This differs from standard D-study variance scaling in how interactions
are handled. Standard D-study scales interactions by all facets in the
component, while this method scales only by the non-object facets.
