# Calculate Comprehensive Sample Size Information

Calculates sample sizes for all variance components in a G-study design,
including main effects, interactions, residual, and nested effects.

## Usage

``` r
calculate_sample_size_info(formula, data, nested = NULL)
```

## Arguments

- formula:

  A formula object with lme4-style random effects.

- data:

  A data frame containing the variables in the formula.

- nested:

  Optional list specifying nesting relationships. Named list where names
  are nested facets and values are nesting facets. For example,
  `list(task = "rater")` means task is nested within rater. If NULL,
  nesting is auto-detected from the data.

## Value

A list with components:

- main:

  Named vector of sample sizes for main effects

- interactions:

  Named vector of sample sizes for interaction terms specified in
  formula

- residual:

  List with facets (character) and n (numeric) for the residual

- nested:

  List of nested effect information, or NULL if none detected
