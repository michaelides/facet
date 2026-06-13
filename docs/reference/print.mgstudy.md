# Print Method for mgstudy Objects

Print Method for mgstudy Objects

## Usage

``` r
# S3 method for class 'mgstudy'
print(
  x,
  digits = 3,
  scale = c("variance", "sd"),
  cor_format = c("long", "matrix"),
  vc_format = c("dimension", "facet"),
  type = c("correlation", "covariance"),
  ...
)
```

## Arguments

- x:

  An mgstudy object.

- digits:

  Number of digits to display.

- scale:

  Scale for displaying results: "variance" (default) or "sd".

- cor_format:

  Format for displaying correlations/covariances: "long" (default) or
  "matrix".

- vc_format:

  Format for displaying variance components: "dimension" (default) or
  "facet". When "facet", variance components are grouped by facet with
  inline correlations.

- type:

  Type of association to display: "correlation" (default) or
  "covariance".

- ...:

  Additional arguments (ignored).

## Value

Invisibly returns x.
