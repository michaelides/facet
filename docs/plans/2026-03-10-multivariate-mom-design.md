# Multivariate G-Study with Method of Moments Backend

## Overview

Implement full multivariate gstudy support for the method of moments (mom) backend, matching the brms backend's functionality and output format.

## Model Specification

Same syntax as brms:
```r
g <- gstudy(mvbind(y1, y2) ~ (1 | person) + (1 | rater), 
            data = my_data, 
            backend = "mom")
```

## Implementation Details

### 1. Formula Parsing
- Extract response variables from `mvbind(y1, y2)` syntax
- Extract facet specifications from random effects
- Detect nesting patterns from data

### 2. MANOVA Fitting
- Use base R's `manova()` with all responses
- Use Error() terms for facet strata
- Extract multivariate summary tables (SS matrices, MS matrices)

### 3. Variance Component Estimation
For each facet × response:
- Extract Expected Mean Squares (EMS) from MANOVA
- Solve for variance components: σ² = (MS_between - MS_within) / n_effective
- Output format: tibble with columns `component`, `dim`, `type`, `var`, `pct`

### 4. Random Effect Correlations
For each facet, compute correlations between random effects across responses:
- Extract random effect estimates for each response
- Compute Pearson correlation with Fisher's z-transformation for CIs

### 5. Residual Correlations
From MANOVA residual (Within) stratum:
- Compute residual covariance matrix
- Extract response residual correlations

### 6. Confidence Intervals
- Variance components: asymptotic approximation (SE = √(2 × MS² / df))
- Correlations: Fisher's z-transformation

## D-Study Integration

Ensure dstudy works with multivariate mom:
- Extract variance components per dimension (using `dim` column)
- Compute G and D coefficients for each response separately
- Output coefficients with `dim` column for identification

## Output Format

Matches brms exactly:
- Variance components: tibble with `dim` column identifying response
- Correlations: stored in `$correlations` list with `residual_cor` tibble and `random_effect_cor` list

## Files to Modify

1. `R/backend-mom.R` - Implement multivariate fitting
2. `R/variance-components.R` - Handle multivariate mom output
3. `R/gstudy.R` - Remove limited support warning
4. `R/backends.R` - Allow mom for multivariate selection
5. `tests/testthat/` - Add comprehensive tests

## Approval

- Design approved by user
- Implementation to follow
