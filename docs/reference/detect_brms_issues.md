# Detect brms Estimation Issues

Checks for low ESS, high R-hat, and treedepth issues in brms models. ESS
threshold is 100 per chain.

## Usage

``` r
detect_brms_issues(model, vc = NULL)
```

## Arguments

- model:

  A brmsfit object from brms::brm().

- vc:

  Optional variance components tibble (to avoid re-extraction).

## Value

A list with components:

- high_rhat:

  Logical indicating if R-hat \> 1.05 was detected

- high_rhat_components:

  Character vector of components with high R-hat

- low_bulk_ess:

  Logical indicating if low Bulk ESS was detected

- low_bulk_ess_components:

  Character vector of components with low Bulk ESS

- low_tail_ess:

  Logical indicating if low Tail ESS was detected

- low_tail_ess_components:

  Character vector of components with low Tail ESS

- treedepth_exceeded:

  Logical indicating if max treedepth was exceeded

- treedepth_count:

  Number of iterations that hit max treedepth
