# Extract Correlation Matrices from Multivariate brms Models

Extracts residual correlations and correlated random effects from
multivariate brms models.

## Usage

``` r
extract_correlations_brms(model)
```

## Arguments

- model:

  A brmsfit object from brms::brm().

## Value

A list containing:

- residual_cor:

  Tibble with columns: dim1, dim2, estimate, se, lower, upper, Rhat,
  Bulk_ESS, Tail_ESS

- random_effect_cor:

  Named list of tibbles for random effects (same columns as
  residual_cor)

- residual_cor_matrix:

  Residual correlation matrix (for backward compatibility)

- random_effect_cor_matrix:

  Named list of correlation matrices (for backward compatibility)
