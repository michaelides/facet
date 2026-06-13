# Check Multivariate Design Balance

Detects and reports unbalanced designs in wide-format multivariate
models. An unbalanced design occurs when facet levels differ across
dimensions after missing value handling.

## Usage

``` r
check_multivariate_balance(
  formula,
  data,
  is_long_format = FALSE,
  unbalanced = FALSE
)
```

## Arguments

- formula:

  A formula object (potentially brmsformula or mvbrmsformula).

- data:

  Data frame containing the variables.

- is_long_format:

  Logical indicating if this is a long-format model.

- unbalanced:

  Logical indicating whether unbalanced estimation is enabled. When
  TRUE, warnings are suppressed and an informational message is
  provided.

## Value

A list with:

- is_balanced: TRUE if no missing values OR missing values with
  consistent levels

- has_missing: TRUE if any missing values in response variables

- has_different_levels: TRUE if facet levels differ across dimensions

- is_genuine_missing: TRUE if missing values present but levels
  consistent

- n_original: Original number of observations

- n_complete: Number of complete cases

- n_removed: Number of rows removed

- dimensions: Dimension/response names

- facet_levels_per_dim: Named list of facet levels per dimension

- warning_message: Warning message for unbalanced design (or NULL)

- info_message: Informational message for genuine missing (or NULL)
