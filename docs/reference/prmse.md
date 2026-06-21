# Tidy PRMSE and VAR Results

**Experimental.** `prmse()` is experimental: testing has been limited and
results should be verified before being cited in published work. The
function and its outputs may change in future releases.

Extract PRMSE(C→S_i) and VAR results as a data frame from a dstudy
object. For univariate models, returns the G and Phi coefficients as
PRMSE metrics.

## Usage

``` r
prmse(
  dstudy_obj,
  include_composite = FALSE,
  include_profile = TRUE,
  ci = NULL,
  ci_method = c("delta", "bootstrap"),
  n_bootstrap = 1000,
  probs = c(0.025, 0.975),
  weights = NULL,
  optimize = NULL,
  optimize_target = "rel",
  grid_resolution = 0.1,
  subscale = NULL
)
```

## Arguments

- dstudy_obj:

  A dstudy object from
  [`dstudy()`](https://github.com/yourorg/facet/reference/dstudy.md)

- include_composite:

  Logical. If TRUE, adds a Composite row with summary metrics. Default
  FALSE. The Composite row shows composite reliability (G/Phi) in
  prmse_s_rel/abs columns, with prmse_c_rel/abs and var_rel/abs set to
  1.0 (trivially, the composite perfectly predicts itself).

- include_profile:

  Logical. If TRUE (default), includes prmse_p_rel and prmse_p_abs
  columns showing how well the entire profile predicts each subscale. If
  FALSE, these columns are omitted for a more compact output.

- ci:

  Character vector specifying which metrics to compute CIs for. Options:
  "prmse", "var", or both. Default NULL (no CIs). **Note:** Credible
  intervals are only available for brms backend with posterior
  estimation. For mom backend, this parameter is ignored with a warning.

- ci_method:

  Unused. Kept for backward compatibility.

- n_bootstrap:

  Unused. Kept for backward compatibility.

- probs:

  Numeric vector of length 2 specifying the quantile probabilities for
  credible intervals (brms backend only). Default is `c(0.025, 0.975)`
  for 95% CI.

- weights:

  Numeric vector of weights for computing composite coefficients. Length
  must match the number of dimensions. Default NULL uses weights from
  dstudy. When alternative weights are provided, all calculations are
  performed using posterior draws to properly propagate uncertainty.

- optimize:

  Character string specifying optimization method. One of:

  "composite"

  :   Maximize composite reliability via generalized eigenvalue

  "subscale"

  :   Maximize VAR for subscales (per-subscale and minimax)

  "tuning"

  :   Grid search guided by VAR diagnostics

  Default NULL (no optimization).

- optimize_target:

  Character string specifying which VAR to optimize. Either "rel"
  (G-based) or "abs" (Phi-based). Default "rel".

- grid_resolution:

  Numeric step size for grid search in "tuning" method. Must be between
  0.01 and 0.5. Default 0.1.

- subscale:

  Character string specifying target subscale for "subscale" method. If
  NULL (default), returns all per-subscale results plus minimax.

## Value

For univariate models, a tibble with:

- dim:

  The response variable name

- prmse_rel:

  The G coefficient (relative reliability)

- prmse_abs:

  The Phi coefficient (absolute reliability)

- prmse_rel_LL, prmse_rel_UL, prmse_abs_LL, prmse_abs_UL:

  Confidence intervals (if ci specified)

For multivariate models, a tibble with one row per dimension and:

- dim:

  Dimension/subscale name

- prmse_s_rel, prmse_s_abs:

  Univariate subscale reliabilities (G and Phi coefficients). These are
  baseline measures that ignore information from other dimensions.

- prmse_c_rel, prmse_c_abs:

  PRMSE when predicting subscale true scores from the weighted
  composite. Used for subscale viability analysis: if PRMSE_S \>
  PRMSE_C, the subscale provides unique information beyond the
  composite.

- prmse_p_rel, prmse_p_abs:

  PRMSE when projecting from the entire profile (all dimensions). Can
  exceed PRMSE_S when dimensions are correlated by borrowing predictive
  information from related dimensions. Only included when
  `include_profile = TRUE` (default).

- var_rel, var_abs:

  Value Added Ratio (PRMSE_S / PRMSE_C).

- ...\_LL, ...\_UL:

  Confidence interval limits for requested metrics.

When `include_composite = TRUE`, a Composite row is added with:

- prmse_s_rel/abs: Composite reliability (G/Phi coefficients)

- prmse_c_rel/abs: 1.0 (composite trivially predicts itself)

- prmse_p_rel/abs: Multivariate reliability (if include_profile = TRUE)

- var_rel/abs: 1.0

When `optimize` is specified, returns a list with optimization results
(see `viable()` for details).

## Details

### Univariate Models

For univariate models, PRMSE_rel and PRMSE_abs correspond directly to
the G and Phi coefficients, respectively:

- `prmse_rel` = G coefficient (relative decision reliability)

- `prmse_abs` = Phi coefficient (absolute decision reliability)

These are the same metrics but labeled differently to maintain
consistency with the multivariate output structure.

### Multivariate Models: Three Types of PRMSE

For multivariate models, three types of PRMSE metrics are computed, each
answering a different question about prediction accuracy:

#### PRMSE_S (Subscale Reliability)

The univariate G coefficient for each dimension. This measures how well
the subscale's own observed scores predict its true scores, ignoring
information from other dimensions.

\$\$\text{PRMSE}\_S = G =
\frac{\sigma^2\_{\text{universe}}}{\sigma^2\_{\text{universe}} +
\sigma^2\_{\text{relative error}}}\$\$

This is the baseline reliability for each dimension when considered in
isolation.

#### PRMSE_C (Composite Prediction)

How well the weighted composite score predicts each subscale's true
scores (Equation 38 from Brennan, 2001):

\$\$\text{PRMSE}\_{(C)} =
\frac{\left(\hat{\sigma}^2\_{(\text{Subscale}\_i)} \cdot
\text{Rel}\_{(\text{Subscale}\_i)} + \sum\_{j \neq i}
\hat{\sigma}\_{(\text{Subscale}\_i,\\
\text{Subscale}\_j)}\right)^2}{\hat{\sigma}^2\_{(\text{Subscale}\_i)}
\cdot \text{Rel}\_{(C)} \cdot \hat{\sigma}^2\_{(C)}}\$\$

where:

- \\\hat{\sigma}^2\_{(\text{Subscale}\_i)}\\: Observed score variance of
  subscale i

- \\\text{Rel}\_{(\text{Subscale}\_i)}\\: Reliability of subscale i (G
  coefficient)

- \\\hat{\sigma}\_{(\text{Subscale}\_i, \text{Subscale}\_j)}\\: Observed
  covariance between subscales i and j

- \\\text{Rel}\_{(C)}\\: Reliability of the composite score

- \\\hat{\sigma}^2\_{(C)}\\: Observed score variance of the composite

**Viability Analysis:** Subscale viability is supported when
\\\text{PRMSE}\_S \> \text{PRMSE}\_C\\, meaning the subscale predicts
its own true scores better than the composite does. The Value Added
Ratio (VAR) formalizes this: \\\text{VAR} = \text{PRMSE}\_S /
\text{PRMSE}\_C\\.

#### PRMSE_P (Profile Projection)

How well the entire profile (all dimensions simultaneously) predicts
each subscale's true scores. This uses the multivariate structure to
potentially improve prediction accuracy by borrowing information from
correlated dimensions:

\$\$\text{PRMSE}\_P\[d\] = \frac{\[\Sigma\_\tau
\Sigma\_{\text{obs}}^{-1} \Sigma\_\tau\]\_{dd}}{\sigma^2\_{\tau,d}}\$\$

where \\\Sigma\_\tau\\ is the true-score covariance matrix and
\\\Sigma\_{\text{obs}}\\ is the observed covariance matrix.

**Relationship to PRMSE_S:** When dimensions are uncorrelated (all
off-diagonal covariances = 0), the matrices become diagonal and PRMSE_P
simplifies to PRMSE_S. In this case, `prmse_p_rel` will approximately
equal `prmse_s_rel`.

These metrics will only differ meaningfully when:

- Dimensions are correlated with each other

- The correlations are substantial enough to provide predictive value

A large difference between PRMSE_P and PRMSE_S indicates that the
multivariate structure provides additional predictive information beyond
what each dimension offers alone. This can occur when one dimension has
high reliability and is correlated with a less reliable dimension,
allowing the profile to "borrow strength" for improved prediction.

#### Summary Comparison

|         |                                                         |
|---------|---------------------------------------------------------|
| Metric  | Question Answered                                       |
| PRMSE_S | "How reliable is this subscale alone?"                  |
| PRMSE_C | "How well does the composite predict this subscale?"    |
| PRMSE_P | "How well does the full profile predict this subscale?" |

### Composite Row

When `include_composite = TRUE`, a Composite row is added to the output.
This row summarizes the composite score's reliability:

- `prmse_s_rel/abs`: The composite's G and Phi coefficients

- `prmse_c_rel/abs`: Always 1.0 (the composite perfectly predicts
  itself)

- `prmse_p_rel/abs`: Multivariate reliability (when include_profile =
  TRUE)

- `var_rel/abs`: Always 1.0

### Weight Optimization

The three optimization methods address a fundamental tension in
multivariate generalizability theory:

- **"composite"**: Finds weights that maximize composite reliability
  (PRMSE). Uses analytical solution via generalized eigenvalue
  decomposition. Best for selection/norm-referencing decisions.

- **"subscale"**: Finds weights that maximize VAR for subscales. Returns
  both per-subscale optimizations and minimax optimization. Best for
  profile interpretation where all subscales matter.

- **"tuning"**: Exhaustive grid search over weight combinations. Useful
  for diagnostic purposes and understanding weight-VAR relationships.

### Confidence Intervals

Credible intervals are only available when using the brms backend with
posterior estimation. These are computed directly from the posterior
draws of the variance components and PRMSE metrics.

For the mom backend (method of moments), credible intervals are not
available because the method produces only point estimates without a
posterior distribution. If credible intervals are needed, consider using
the brms backend:

    g_mv <- gstudy(formula, data, backend = "brms")
    d_mv <- dstudy(g_mv, n = ...)
    prmse(d_mv, ci = "prmse") # CIs from posterior draws

## See also

[`dstudy()`](https://github.com/yourorg/facet/reference/dstudy.md) for
computing D-study results.

## Examples

``` r
if (FALSE) { # \dontrun{
# Univariate model
g_uni <- gstudy(Score ~ (1 | Person) + (1 | Item), data = mydata)
d_uni <- dstudy(g_uni, n = list(Item = 5))
prmse(d_uni) # Returns G and Phi as prmse_rel and prmse_abs

# Multivariate model with mom backend (no CIs available)
g_mv <- gstudy(Score ~ 0 + Subtest + (0 + Subtest | Person),
  data = data, dimension_var = "Subtest", backend = "mom")
d_mv <- dstudy(g_mv, n = list(Person = 5))
prmse(d_mv)

# Include Composite row
prmse(d_mv, include_composite = TRUE)

# Exclude profile projection metrics for compact output
prmse(d_mv, include_profile = FALSE)

# Understanding PRMSE_S vs PRMSE_P:
# When dimensions are uncorrelated, prmse_p_rel ≈ prmse_s_rel
# When correlated, prmse_p_rel may exceed prmse_s_rel
result <- prmse(d_mv)
diff <- result$prmse_s_rel - result$prmse_p_rel
if (all(abs(diff) < 0.01)) {
  message("Dimensions appear uncorrelated: PRMSE_P ≈ PRMSE_S")
} else {
  message("Dimensions are correlated: PRMSE_P differs from PRMSE_S")
}

# Weight optimization
prmse(d_mv, optimize = "composite")
prmse(d_mv, optimize = "subscale")

# Multivariate model with brms backend (credible intervals available)
library(brms)
g_mv <- gstudy(
  bf(Score ~ 0 + Subtest + (0 + Subtest | r | Person)),
  data = data, backend = "brms"
)
d_mv <- dstudy(g_mv, n = list(Person = 5))
prmse(d_mv, ci = "prmse") # Returns point estimates with credible intervals
} # }
```
