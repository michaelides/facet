# Variance-Covariance Matrix Extractor Generic

Generic function for extracting variance-covariance matrices from fitted
models. Dispatches to the appropriate method based on the class of the
object.

## Usage

``` r
VarCorr(x, ...)
```

## Arguments

- x:

  An object from which to extract variance components.

- ...:

  Additional arguments passed to methods.

## Value

Variance-covariance matrices or posterior summaries, depending on the
method.
