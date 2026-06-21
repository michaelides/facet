# Compute VAR Using Haberman Formula (Per Draw)

Computes the Value Added Ratio (VAR) for each subscale using the
Haberman (2008) formula, which evaluates whether a subscale's own scores
better estimate its universe score than the composite does.

## Usage

``` r
compute_var_haberman_draws(
  uni_cov_draws,
  total_rel_cov_draws,
  total_abs_cov_draws,
  dimensions,
  weights
)
```

## Arguments

- uni_cov_draws:

  3D array (n_draws x n_dims x n_dims) of universe score covariances.

- total_rel_cov_draws:

  3D array (n_draws x n_dims x n_dims) of total observed covariances for
  relative error.

- total_abs_cov_draws:

  3D array (n_draws x n_dims x n_dims) of total observed covariances for
  absolute error.

- dimensions:

  Character vector of dimension names.

- weights:

  Named numeric vector of weights for the composite scores.

## Value

A list containing:

- prmse_s_rel, prmse_s_abs:

  Univariate subscale reliabilities.

- prmse_c_rel, prmse_c_abs:

  PRMSE when projecting from weighted composite.

- prmse_p_rel, prmse_p_abs:

  PRMSE when projecting from entire profile (diagonal slice).

- var_rel, var_abs:

  Value Added Ratios (PRMSE_S / PRMSE_C).

- prmse_mv_rel, prmse_mv_abs:

  Global multivariate PRMSE (trace formula).

## Details

The formula is: \$\$VAR(S_i) = \frac{PRMSE(S_i)}{PRMSE(C \rightarrow
S_i)}\$\$

Where PRMSE(C -\> S_i) is the Haberman (2008) correct form: \$\$PRMSE(C
\rightarrow S_i) = \frac{\[\mathrm{Cov}(\tau_i,
C)\]^2}{\mathrm{Var}(\tau_i) \cdot \mathrm{Rel}(C) \cdot
\mathrm{Var}(C)}\$\$ with \\\mathrm{Cov}(\tau_i, C) = \sum_j w_j \\
\mathrm{Cov}(\tau_i, \tau_j)\\ computed from the universe-score
covariance matrix. This is the corrected form of what Vispoel et al.
(2023) labeled as "Equation 38" (the printed equation contains a
typesetting error in the denominator; see
[`prmse`](https://github.com/yourorg/facet/reference/prmse.md) for
details).

## References

Haberman, S. J. (2008). When can subscores have value? *Journal of
Educational and Behavioral Statistics*, 33(2), 204-229.

Vispoel, W. P., Lee, H., Hong, H., & Chen, T. (2023). Applying
multivariate generalizability theory to psychological assessments.
*Psychological Methods*, 28(4), 847-870.
