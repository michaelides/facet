# Extract set_rescor Setting from Formula

Parses the formula to detect if set_rescor() is specified and returns
its value. Returns FALSE by default (no residual correlations) if not
specified.

## Usage

``` r
extract_rescor_setting(formula)
```

## Arguments

- formula:

  A formula object.

## Value

Logical indicating whether to estimate residual correlations (default
FALSE).
