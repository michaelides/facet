# Extract Facet Specifications Preserving User's Order

Extracts facet specifications as the user specified them (e.g.,
"Rater:Task" vs "Task:Rater"), preserving the original order.

## Usage

``` r
extract_facet_specs(formula = NULL)
```

## Arguments

- formula:

  A formula object.

## Value

Character vector of facet specifications in user-specified order.
