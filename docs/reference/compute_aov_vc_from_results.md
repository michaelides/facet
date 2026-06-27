# Compute Variance Components from ANOVA Results

Compute Variance Components from ANOVA Results

## Usage

``` r
compute_aov_vc_from_results(
  anova_results,
  data,
  response,
  random_facets,
  random_facet_specs = NULL,
  name_mapping = NULL
)
```

## Arguments

- anova_results:

  List of ANOVA results per stratum.

- data:

  The data frame used to fit the model.

- response:

  Name of the response variable.

- random_facets:

  Character vector of individual facet names from random effects.

- random_facet_specs:

  Character vector of facet specifications as user specified (e.g.,
  "Rater:Task").

- name_mapping:

  Optional named list mapping rewritten ANOVA source names (e.g.
  "Rater:Task") to user-supplied component names (e.g. "Rater"), as
  produced by
  [`adapt_aov_formula()`](https://github.com/yourorg/facet/reference/adapt_aov_formula.md).
