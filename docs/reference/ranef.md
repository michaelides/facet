# Random Effects Extractor Generic

Generic function for extracting random effects (BLUPs/conditional modes)
from fitted models. Dispatches to the appropriate method based on the
class of the object.

## Usage

``` r
ranef(object, ...)
```

## Arguments

- object:

  An object from which to extract random effects.

- ...:

  Additional arguments passed to methods.

## Value

Random effects from the fitted model.
