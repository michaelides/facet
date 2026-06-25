# Adapt aov Formula for the Method-of-Moments Backend

Inspects the formula and data to detect facets that are nested within
other facets in the data, and rewrites the aov Error() strata so the
canonical G-study specification works without spurious rank-deficiency
warnings. This is only used by the mom backend; lme4 and brms are
unaffected.

## Usage

``` r
adapt_mom_aov_formula(formula, data, nested = NULL)
```

## Arguments

- formula:

  A formula object.

- data:

  A data frame.

- nested:

  Optional named list overriding auto-detected nesting (e.g.,
  `list(Rater = "Task")`).

## Value

A list with components:

- aov_formula:

  The rewritten aov() formula string.

- name_mapping:

  Named list mapping the rewritten ANOVA source name to the
  user-supplied component name.

- nested_info:

  The (possibly user-overridden) nesting relationships.

- random_facet_specs:

  The (possibly rewritten) random-effect specs.

## Details

Behaviour:

- Parses the formula with
  [`parse_g_formula()`](https://github.com/yourorg/facet/reference/parse_g_formula.md)
  to obtain the user-supplied random-effect specifications.

- Calls
  [`detect_nesting_patterns()`](https://github.com/yourorg/facet/reference/detect_nesting_patterns.md)
  (in R/sample-sizes.R) to find facets that are nested in other facets
  in the data. The user may override auto-detection by passing `nested`
  (the same argument accepted by
  [`gstudy()`](https://github.com/yourorg/facet/reference/gstudy.md)).

- For each detected nested pair (A nested in B):

  - If the user wrote the explicit nested form `A:B` in the formula, the
    helper leaves the spec alone and records the identity mapping
    `A:B -> A:B` so downstream code can label the variance component.

  - If the user wrote the un-nested main effect `A` (the canonical
    G-study shorthand for the A-within-B effect when labels are already
    unique within B), the helper substitutes `A` with `A:B` in the
    Error() term and records the mapping `A:B -> A` so the variance
    component is labelled with the user-supplied name even though the
    Error() decomposition uses the nested form.

- Returns the rewritten aov formula string and a name mapping that
  [`compute_mom_vc_from_results()`](https://github.com/yourorg/facet/reference/compute_mom_vc_from_results.md)
  applies to the resulting ANOVA table.

Any `Error() model is singular` warning from
[`aov()`](https://rdrr.io/r/stats/aov.html) is suppressed by the caller;
for these designs the warning is spurious — the variance components are
still computed correctly, matching lme4 to several decimals.
