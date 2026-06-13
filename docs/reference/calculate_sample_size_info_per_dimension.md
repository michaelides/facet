# Calculate Sample Size Information Per Dimension

Calculates sample sizes separately for each level of the dimension
variable in long-format multivariate models.

## Usage

``` r
calculate_sample_size_info_per_dimension(formula, data, dimension_var)
```

## Arguments

- formula:

  A formula object.

- data:

  A data frame containing the variables.

- dimension_var:

  Character string naming the dimension variable.

## Value

A named list where each element is a sample_size_info for that
dimension.
