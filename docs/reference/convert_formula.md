# Convert Formula to Backend-Specific Format

Converts a G-study formula to the format required by a specific backend.
For lme4, this removes any brms-specific syntax.

## Usage

``` r
convert_formula(formula, backend)
```

## Arguments

- formula:

  A formula object.

- backend:

  Character string indicating the target backend.

## Value

A formula object in the backend-specific format.
