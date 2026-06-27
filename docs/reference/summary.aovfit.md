# Summary Method for aovfit Objects

Prints a summary of a aovfit object including variance components.

## Usage

``` r
# S3 method for class 'aovfit'
summary(object, ...)
```

## Arguments

- object:

  A aovfit object from
  [`fit_aov()`](https://github.com/yourorg/facet/reference/fit_aov.md).

- ...:

  Additional arguments passed to other methods.

## Value

Side effect: prints summary to console. Returns the object invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, each = 5)),
  item = factor(rep(1:5, times = 20))
)
model <- fit_aov(score ~ (1 | person) + (1 | item), data)
summary(model)
} # }
```
