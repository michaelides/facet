# Detect Crossing Patterns Between Facets in Data

Analyzes the data to determine if pairs of facets are nested, crossed,
or partially crossed.

- Nested: each level of A appears with exactly 1 level of B

- Crossed: each level of A appears with all levels of B

- Partial: some levels appear with 1, others with multiple

## Usage

``` r
detect_crossing_patterns_from_data(data, facet_specs)
```

## Arguments

- data:

  The data frame

- facet_specs:

  Character vector of facet specifications (e.g., c("Person", "Task",
  "Rater:Task"))

## Value

List containing:

- pairwise: list of crossing patterns for each facet pair

- facets: list with crossing info for each facet

- pair_counts: number of unique observed pairs for each combination
