# Compute Residual Divisor Excluding the Object of Measurement

Calculates the divisor for the residual variance component when scaling
for a D-study, excluding the object of measurement from the divisor.
This implements the generalizability-theory rule that the object of
measurement is NOT averaged over when rescaling the residual.

## Usage

``` r
compute_residual_divisor(residual_is, n, object_spec)
```

## Arguments

- residual_is:

  Character string specifying residual composition (e.g.,
  "Person:Rater"). If NULL or empty, all facets in `n` are used
  (excluding object).

- n:

  Named list of sample sizes for each facet.

- object_spec:

  Character vector of object component names.

## Value

Numeric divisor (product of `n` for non-object facets in the residual).

## Details

The residual is decomposed into its constituent facets (e.g.,
`"Person:Rater"` → `c("Person", "Rater")`). The divisor is the product
of `n[[f]]` for facets that are NOT the object of measurement. For
example, with `residual_is = "Person:Rater"`, `object = "Person"`, and
`n = list(Rater = 4)`, the divisor is `4` (Rater only; Person excluded).
