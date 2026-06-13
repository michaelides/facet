# Plot Method for dstudy Objects

Plot Method for dstudy Objects

## Usage

``` r
# S3 method for class 'dstudy'
plot(
  x,
  type = c("coefficients", "sweep"),
  coefficient = c("both", "g", "phi"),
  ...
)
```

## Arguments

- x:

  A dstudy object.

- type:

  Type of plot: "coefficients" or "sweep".

- coefficient:

  Which coefficient to plot: "g", "phi", or "both".

- ...:

  Additional arguments passed to plotting functions.

## Value

A ggplot object (invisibly).
