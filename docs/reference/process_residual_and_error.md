# Process Residual and Error Specification Overlap

Determines residual composition from formula/data, handles overlap
between error and aggregation specifications, and removes residual from
error spec if double-counting would occur.

## Usage

``` r
process_residual_and_error(
  gstudy_obj,
  error = NULL,
  aggregation = NULL,
  residual_is = NULL
)
```

## Arguments

- gstudy_obj:

  A gstudy or mgstudy object

- error:

  Error specification

- aggregation:

  Aggregation specification

- residual_is:

  User-specified residual composition

## Value

A list with residual_composition, residual_is_effective, error (cleaned)
