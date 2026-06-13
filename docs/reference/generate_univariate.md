# Generate Univariate Data

Generate Univariate Data

## Usage

``` r
generate_univariate(
  facets,
  vc,
  sd_residual,
  response_var,
  nested = NULL,
  prefix = ""
)
```

## Arguments

- facets:

  Named list of facet levels

- vc:

  Named list of variance components

- sd_residual:

  Residual standard deviation

- response_var:

  Name of response variable

- nested:

  Named list of nesting relationships

- prefix:

  Prefix for level names

## Value

Data frame with simulated data
