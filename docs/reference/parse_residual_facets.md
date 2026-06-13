# Parse Residual Facets from Formula

Determines which facets make up the residual (lowest-level) variance
component based on the formula and data structure. The residual is
typically the interaction of all non-nested facets.

## Usage

``` r
parse_residual_facets(formula, data = NULL)
```

## Arguments

- formula:

  A formula object with lme4-style random effects (e.g.,
  `score ~ (1|p) + (1|i:d)`).

- data:

  Optional data frame to detect actual nesting patterns from the data.
  If NULL, assumes all facets are crossed (most conservative).

## Value

A character string representing the facets that make up the residual,
with facets separated by colons (e.g., `"p:i:d"`).

## Details

The function determines the residual composition by:

1.  Parsing the formula to extract all unique facets from random effects

2.  If `data` is provided, detecting which facets are nested within
    others

3.  Returning the interaction of all non-nested (crossed) facets

### Nesting Detection

When `data` is provided, the function analyzes the actual data structure
to determine if facets appearing in interaction terms (like `i:d`) are
nested or crossed:

- **Nested**: If each level of `i` appears with only one level of `d`,
  then `i` is nested in `d`. The nested facet is absorbed and the
  residual excludes the nesting facet (e.g., residual = `"p:i"` instead
  of `"p:i:d"`).

- **Crossed**: If levels of `i` appear with multiple levels of `d`, the
  facets are crossed. The residual includes all facets (e.g.,
  `"p:i:d"`).

## See also

[`detect_crossing_patterns_from_data`](https://github.com/yourorg/facet/reference/detect_crossing_patterns_from_data.md)
for the underlying nesting detection
[`parse_g_formula`](https://github.com/yourorg/facet/reference/parse_g_formula.md)
for parsing formula structure

## Examples

``` r
# Fully crossed design: p, i, d all crossed with each other
f1 <- score ~ (1 | p) + (1 | i) + (1 | d)
parse_residual_facets(f1) # Returns "p:i:d"
#> Error in parse_residual_facets(f1): could not find function "parse_residual_facets"

# With nested data: i nested in d
nested_data <- data.frame(
  score = rnorm(100),
  p = factor(rep(1:10, 10)),
  i = factor(rep(1:5, each = 20)), # i levels repeat within d
  d = factor(rep(1:2, each = 50))
)
f2 <- score ~ (1 | p) + (1 | i:d) + (1 | d)
parse_residual_facets(f2, nested_data) # May return "p:i" if i nested in d
#> Error in parse_residual_facets(f2, nested_data): could not find function "parse_residual_facets"
```
