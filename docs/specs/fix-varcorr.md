# Spec: fix-varcorr

Scope: feature

# Feature: Fix VarCorr.gstudy for Multivariate brms Models

## Overview
Fix VarCorr.gstudy method to properly handle multivariate brms models by detecting multivariate structure and processing brms::VarCorr output appropriately.

## Problem Statement
- VarCorr.gstudy works with univariate brms and lme4 backends
- VarCorr.gstudy fails when gstudy is fitted with multivariate brms (mvbrmsformula)
- Root cause: brms::VarCorr returns different output structure for multivariate models

## Technical Details

### Current Behavior (Broken)
```r
VarCorr.gstudy <- function(x, ...) {
  model <- x$model
  if (x$backend == "brms") {
    brms::VarCorr(model, ...)  # Fails for multivariate models
  }
}
```

### Expected Behavior
- VarCorr.gstudy should detect multivariate brms models using `x$is_multivariate` field
- Should process brms::VarCorr output appropriately for multivariate structure
- Should return variance components organized by response variable

### Implementation Approach
1. Check `x$is_multivariate` field in VarCorr.gstudy
2. For multivariate models, iterate over response variables using `brms::VarCorr(model, resp = ...)`
3. Restructure output to match expected format for multivariate models

## Acceptance Criteria
1. VarCorr.gstudy works with univariate lme4 backend (existing behavior preserved)
2. VarCorr.gstudy works with univariate brms backend (existing behavior preserved)  
3. VarCorr.gstudy works with multivariate brms backend (new behavior)
4. Returns properly structured output for multivariate models (variance-covariance matrices per response)