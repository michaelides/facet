# Calculate Sample Sizes for Nested Effects

For each nesting relationship, calculates the number of nested units per
nesting unit, including mean and harmonic mean for unbalanced designs.

## Usage

``` r
calculate_nested_sample_sizes(data, nesting_info)
```

## Arguments

- data:

  A data frame.

- nesting_info:

  List of nesting relationships from detect_nesting_patterns().

## Value

A list where each element contains:

- nested_facet:

  Name of the nested facet

- nesting_facet:

  Name of the nesting facet

- n_groups:

  Number of groups (levels of nesting facet)

- mean_per_group:

  Mean number of nested units per group

- harmonic_mean_per_group:

  Harmonic mean of nested units per group

- counts_per_group:

  Named vector of counts for each group
