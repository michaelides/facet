# Solve Generalized Eigenvalue Problem for Composite Optimization

Finds optimal weights that maximize composite reliability by solving the
generalized eigenvalue problem: Sigma_tau \* w = lambda \* Sigma_total
\* w

## Usage

``` r
solve_generalized_eigenvalue(Sigma_tau, Sigma_total, dimensions)
```

## Arguments

- Sigma_tau:

  Universe score variance-covariance matrix

- Sigma_total:

  Total variance-covariance matrix (Sigma_tau + Sigma_delta)

- dimensions:

  Character vector of dimension names

## Value

List with optimal weights, eigenvalue, and explained variance
