# Parse D-Study Specifications

Extracts variance components, object of measurement, and parses
universe/error/aggregation specifications. Validates that universe and
error don't overlap and warns about non-interacting universe components.

## Usage

``` r
parse_dstudy_specifications(
  gstudy_obj,
  universe = NULL,
  error = NULL,
  aggregation = NULL
)
```

## Arguments

- gstudy_obj:

  A gstudy or mgstudy object

- universe:

  Universe specification

- error:

  Error specification

- aggregation:

  Aggregation specification

## Value

A list with vc, object, object_spec, universe_spec, error_spec, agg_spec
