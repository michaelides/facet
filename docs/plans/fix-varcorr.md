---
plan name: fix-varcorr
plan description: Fix VarCorr for multivariate brms
plan status: active
---

## Idea
Fix VarCorr.gstudy to properly handle multivariate brms models by adapting how it calls and processes brms::VarCorr output

## Implementation
- 1. Modify VarCorr.gstudy in R/methods.R to detect multivariate brms models using x$is_multivariate field
- 2. For multivariate brms models, pass appropriate arguments to brms::VarCorr to get variance components per response variable
- 3. Restructure the output to match expected format for multivariate models (similar to extract_vc_brms logic)
- 4. Add tests for VarCorr.gstudy with multivariate brms models
- 5. Verify fix works with both univariate and multivariate brms models

## Required Specs
<!-- SPECS_START -->
- fix-varcorr
<!-- SPECS_END -->