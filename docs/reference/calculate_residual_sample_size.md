# Calculate Sample Size for Residual

Calculates the number of unique observations for the residual component.

## Usage

``` r
calculate_residual_sample_size(data, residual_facets)
```

## Arguments

- data:

  A data frame.

- residual_facets:

  Character string of facets making up the residual, separated by colons
  (e.g., "person:rater:item").

## Value

Integer number of unique combinations.
