# Summary Method for mgstudy Objects

Summary Method for mgstudy Objects

## Usage

``` r
# S3 method for class 'mgstudy'
summary(
  object,
  scale = c("variance", "sd"),
  digits = 3,
  cor_format = c("long", "matrix"),
  vc_format = c("dimension", "facet"),
  type = c("correlation", "covariance"),
  ...
)
```

## Arguments

- object:

  An mgstudy object.

- scale:

  Scale for displaying results: "variance" (default) or "sd".

- digits:

  Number of digits to display (default 3).

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

Invisibly returns object.
