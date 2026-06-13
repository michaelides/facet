# Generate Correlated Effects Using Cholesky Decomposition

Generate Correlated Effects Using Cholesky Decomposition

## Usage

``` r
generate_correlated_effects(n, cor_matrix, sds)
```

## Arguments

- n:

  Number of observations

- cor_matrix:

  Correlation matrix

- sds:

  Standard deviations (vector of length ncol(cor_matrix))

## Value

Matrix of correlated effects (n x ncol(cor_matrix))
