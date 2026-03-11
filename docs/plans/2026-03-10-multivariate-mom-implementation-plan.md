# Multivariate G-Study with Mom Backend Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement full multivariate gstudy support for the method of moments (mom) backend, matching brms backend's functionality and output format, including dstudy integration.

**Architecture:** Use MANOVA-based estimation for variance components and correlations. Fit multivariate ANOVA with all responses, extract variance-covariance matrices from MANOVA strata, compute correlations from random effect estimates and residuals.

**Tech Stack:** R base functions (manova, aov), facet package structure

---

## Task 0: Preparation

**Files:**
- Verify: `R/gstudy.R`
- Verify: `R/backend-mom.R`
- Verify: `R/backends.R`
- Verify: `R/variance-components.R`
- Verify: `R/dstudy.R`

**Step 1: Check current state of multivariate support in mom backend**

Run in R:
```r
library(facet)
data <- data.frame(
  y1 = rnorm(100),
  y2 = rnorm(100),
  person = factor(rep(1:20, 5)),
  rater = factor(rep(1:5, each = 20))
)
# This should currently show warning about limited support
result <- tryCatch(
  gstudy(mvbind(y1, y2) ~ (1 | person) + (1 | rater), data = data, backend = "mom"),
  error = function(e) e$message
)
print(result)
```

---

## Task 1: Remove Limited Multivariate Support Warning

**Files:**
- Modify: `R/gstudy.R:135-142` - Remove warning about limited multivariate support
- Modify: `R/formula-helpers.R:173-179` - Remove warning in validate_formula

**Step 1: Remove warning in gstudy.R**

Find lines 135-142 in gstudy.R and remove:
```r
# 7c. Check multivariate compatibility for mom backend
if (is_mv && selected_backend == "mom") {
  warning(
    "Method of moments backend has limited multivariate support. ",
    "For full multivariate functionality, consider using brms backend.",
    call. = FALSE
  )
}
```

**Step 2: Remove warning in formula-helpers.R**

Find lines 173-179 and remove:
```r
if (is_mv && backend == "mom") {
  warning(
    "Method of moments backend has limited multivariate support.\n",
    "For full multivariate functionality, consider using brms backend.",
    call. = FALSE
  )
}
```

**Step 3: Test that warning is removed**

Run: `Rscript -e "library(facet); data <- data.frame(y1=rnorm(100), y2=rnorm(100), person=factor(rep(1:20,5)), rater=factor(rep(1:5,each=20))); gstudy(mvbind(y1,y2)~(1|person)+(1|rater), data=data, backend='mom')"`
Expected: No warning about limited support

---

## Task 2: Update Backend Selection for Mom Multivariate

**Files:**
- Modify: `R/backends.R:20-74` - Update select_backend function

**Step 1: Modify select_backend to allow mom for multivariate**

In `select_backend()` function, after checking backend availability, modify the multivariate check:

Current code at lines 65-71:
```r
if (is_mv && backend == "lme4") {
  stop(
    "Multivariate formulas (containing mvbind or set_rescor) require brms backend.\n",
    "Use: backend = 'brms'",
    call. = FALSE
  )
}
```

Change to:
```r
if (is_mv && backend == "lme4") {
  stop(
    "Multivariate formulas (containing mvbind or set_rescor) require brms or mom backend.\n",
    "Use: backend = 'brms' or backend = 'mom'",
    call. = FALSE
  )
}
```

**Step 2: Also update the auto selection logic**

In the auto selection block (lines 26-50), ensure mom is considered for multivariate:
```r
if (backend == "auto") {
  if (is_mv) {
    # Prefer brms for multivariate but allow mom
    if (requireNamespace("brms", quietly = TRUE)) {
      return("brms")
    } else if (requireNamespace("lme4", quietly = TRUE)) {
      return("lme4")
    } else {
      return("mom")  # mom is always available
    }
  } else {
    # Default to lme4 for univariate
    if (requireNamespace("lme4", quietly = TRUE)) {
      return("lme4")
    } else if (requireNamespace("brms", quietly = TRUE)) {
      return("brms")
    } else {
      return("mom")
    }
  }
}
```

**Step 3: Run test to verify backend selection**

```r
library(facet)
# Test that backend selection works
data <- data.frame(y1=rnorm(100), y2=rnorm(100), person=factor(rep(1:20,5)), rater=factor(rep(1:5,each=20)))
# Should work without error now
```

---

## Task 3: Implement Multivariate MANOVA Fitting

**Files:**
- Modify: `R/backend-mom.R:655-840` - Rewrite fit_mom_multivariate function

**Step 1: Write failing test first**

Create test file `tests/testthat/test-mom-multivariate.R`:
```r
test_that("mom backend fits multivariate model", {
  skip_on_cran()
  
  set.seed(123)
  data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  result <- gstudy(mvbind(y1, y2) ~ (1 | person) + (1 | rater), 
                   data = data, 
                   backend = "mom")
  
  # Check basic structure
  expect_s3_class(result, "mgstudy")
  expect_true("variance_components" %in% names(result))
  
  # Check that we have variance components for each dimension
  vc <- result$variance_components
  expect_true("dim" %in% names(vc))
  expect_equal(sort(unique(vc$dim)), c("y1", "y2"))
  
  # Check that we have variance components for each facet per dimension
  expect_true(all(c("y1", "y2") %in% vc$dim))
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e "devtools::test(filter='test-mom-multivariate')"`
Expected: FAIL - function not fully implemented

**Step 3: Implement MANOVA-based fit_mom_multivariate**

Replace `fit_mom_multivariate` function (lines 655-724 in backend-mom.R) with:

```r
#' Fit Multivariate Model Using Method of Moments
#'
#' @keywords internal
fit_mom_multivariate <- function(formula, data, ...) {
  # Extract response variables from mvbind
  formula_char <- deparse(formula)
  
  resp_match <- regmatches(formula_char, regexpr("mvbind\\s*\\([^)]+\\)", formula_char))
  if (length(resp_match) == 0 || resp_match == "") {
    stop("Could not parse multivariate formula for method of moments",
      call. = FALSE
    )
  }
  
  # Extract variable names from mvbind(...)
  resp_str <- sub("mvbind\\s*\\(", "", resp_match)
  resp_str <- sub("\\)\\s*$", "", resp_str)
  responses <- unlist(strsplit(resp_str, "\\s*,\\s*"))
  responses <- trimws(responses)
  
  # Parse the formula to get facets
  parsed <- parse_g_formula(formula)
  random_facets <- parsed$random_facets
  random_facet_specs <- parsed$random_facet_specs
  
  # Build multivariate formula for manova
  # Response side: cbind(y1, y2)
  rhs_formula <- sub("mvbind\\s*\\([^)]+\\)\\s*~", "", formula_char)
  rhs_formula <- sub("^\\s*~", "", rhs_formula)
  
  # Build manova formula: cbind(y1, y2) ~ facet1 + facet2 + ...
  facets_in_formula <- paste(random_facets, collapse = " + ")
  manova_formula <- paste0("cbind(", paste(responses, collapse = ", "), ") ~ ", facets_in_formula)
  
  # Fit MANOVA model with Error terms for each facet
  # Build aovlist-style formula with Error strata
  aov_formula <- build_aov_formula_mv(formula, data, responses)
  
  # Fit the MANOVA model
  manova_model <- tryCatch(
    {
      manova(as.formula(aov_formula), data = data)
    },
    error = function(e) {
      stop("Error fitting multivariate model with method of moments: ",
        e$message,
        call. = FALSE
      )
    }
  )
  
  # Extract variance components for each response
  all_vc <- list()
  for (resp in responses) {
    # Build univariate formula for this response
    univ_formula <- paste(resp, "~", rhs_formula)
    
    # Fit univariate aov for this response
    aov_univ <- tryCatch(
      {
        aov(as.formula(aov_formula), data =      error = function data)
      },
(e) NULL
    )
    
    if (!is.null(aov_univ)) {
      # Extract ANOVA results
      anova_results <- extract_anova_from_aovlist(aov_univ, random_facets)
      
      # Compute variance components
      vc <- compute_mom_vc_from_results(
        anova_results, data, resp, random_facets, random_facet_specs
      )
      
      # Add dimension column
      vc$dim <- resp
      all_vc[[resp]] <- vc
    }
  }
  
  # Combine variance components
  combined_vc <- do.call(rbind, all_vc)
  
  # Compute correlations (both residual and random effect)
  correlations <- compute_mom_correlations(formula, data, responses, random_facets)
  
  # Create result object
  result <- structure(
    list(
      formula = formula,
      data = data,
      manova_model = manova_model,
      variance_components = combined_vc,
      responses = responses,
      random_facets = random_facets,
      random_facet_specs = random_facet_specs,
      correlations = correlations,
      is_multivariate = TRUE
    ),
    class = "momfit"
  )
  
  result
}
```

**Step 4: Add helper function build_aov_formula_mv**

Add new function after `build_aov_formula` (around line 172):

```r
#' Build aov Formula for Multivariate Model
#'
#' @keywords internal
build_aov_formula_mv <- function(formula, data, responses) {
  parsed <- parse_g_formula(formula)
  random_facets <- parsed$random_facets
  
  if (length(random_facets) == 0) {
    stop("No random facets found in multivariate formula", call. = FALSE)
  }
  
  # Use Error() terms similar to univariate
  if (length(random_facets) == 1) {
    aov_formula <- paste(
      paste(responses, collapse = " + "), 
      "~", random_facets[1],
      "+ Error(", random_facets[1], ")"
    )
  } else {
    error_strata <- paste(random_facets, collapse = " + ")
    aov_formula <- paste(
      paste(responses, collapse = " + "),
      "~", paste(random_facets, collapse = " + "),
      "+ Error(", error_strata, ")"
    )
  }
  
  aov_formula
}
```

**Step 5: Add correlation computation function**

Add new function after `compute_mom_variance_components`:

```r
#' Compute Correlations for Multivariate Mom Model
#'
#' @keywords internal
compute_mom_correlations <- function(formula, data, responses, random_facets) {
  if (length(responses) < 2) {
    return(NULL)
  }
  
  result <- list(
    residual_cor = NULL,
    random_effect_cor = list(),
    residual_cor_matrix = NULL,
    random_effect_cor_matrix = list()
  )
  
  # Compute residual correlations from data
  # For each pair of responses, compute correlation of residuals
  resid_cor_list <- list()
  
  for (i in seq_along(responses)) {
    for (j in seq_along(responses)) {
      if (i < j) {
        resp_i <- responses[i]
        resp_j <- responses[j]
        
        # Fit main effects model and get residuals
        rhs <- paste(random_facets, collapse = " + ")
        
        # Simple approach: compute correlation of responses directly
        # (approximation for balanced designs)
        cor_ij <- cor(data[[resp_i]], data[[resp_j]])
        n <- nrow(data)
        se_ij <- sqrt((1 - cor_ij^2) / (n - 2))
        
        # Fisher's z transformation for CI
        z_ij <- atanh(cor_ij)
        se_z <- 1 / sqrt(n - 3)
        z_lower <- z_ij - 1.96 * se_z
        z_upper <- z_ij + 1.96 * se_z
        lower <- tanh(z_lower)
        upper <- tanh(z_upper)
        
        resid_cor_list[[length(resid_cor_list) + 1]] <- data.frame(
          dim1 = resp_i,
          dim2 = resp_j,
          estimate = cor_ij,
          se = se_ij,
          lower = lower,
          upper = upper,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(resid_cor_list) > 0) {
    result$residual_cor <- do.call(rbind, resid_cor_list)
    
    # Build matrix format
    n_resp <- length(responses)
    cor_matrix <- matrix(NA, n_resp, n_resp)
    rownames(cor_matrix) <- responses
    colnames(cor_matrix) <- responses
    for (i in seq_len(nrow(result$residual_cor))) {
      r <- result$residual_cor[i, ]
      cor_matrix[r$dim1, r$dim2] <- r$estimate
      cor_matrix[r$dim2, r$dim1] <- r$estimate
    }
    diag(cor_matrix) <- 1
    result$residual_cor_matrix <- cor_matrix
  }
  
  # For random effect correlations, compute from facet means
  # For each facet, compute correlation of means across responses
  for (facet in random_facets) {
    facet_means <- aggregate(
      data[, responses, drop = FALSE],
      by = list(facet = data[[facet]]),
      FUN = mean
    )
    
    if (nrow(facet_means) > 2) {
      facet_cor_list <- list()
      
      for (i in seq_along(responses)) {
        for (j in seq_along(responses)) {
          if (i < j) {
            resp_i <- responses[i]
            resp_j <- responses[j]
            
            cor_ij <- cor(facet_means[[resp_i]], facet_means[[resp_j]])
            n <- nrow(facet_means)
            se_ij <- sqrt((1 - cor_ij^2) / (n - 2))
            
            z_ij <- atanh(cor_ij)
            se_z <- 1 / sqrt(n - 3)
            z_lower <- z_ij - 1.96 * se_z
            z_upper <- z_ij + 1.96 * se_z
            lower <- tanh(z_lower)
            upper <- tanh(z_upper)
            
            facet_cor_list[[length(facet_cor_list) + 1]] <- data.frame(
              dim1 = resp_i,
              dim2 = resp_j,
              estimate = cor_ij,
              se = se_ij,
              lower = lower,
              upper = upper,
              stringsAsFactors = FALSE
            )
          }
        }
      }
      
      if (length(facet_cor_list) > 0) {
        result$random_effect_cor[[facet]] <- do.call(rbind, facet_cor_list)
        
        # Build matrix
        n_resp <- length(responses)
        cor_matrix <- matrix(NA, n_resp, n_resp)
        rownames(cor_matrix) <- responses
        colnames(cor_matrix) <- responses
        for (i in seq_len(nrow(result$random_effect_cor[[facet]]))) {
          r <- result$random_effect_cor[[facet]][i, ]
          cor_matrix[r$dim1, r$dim2] <- r$estimate
          cor_matrix[r$dim2, r$dim1] <- r$estimate
        }
        diag(cor_matrix) <- 1
        result$random_effect_cor_matrix[[facet]] <- cor_matrix
      }
    }
  }
  
  result
}
```

**Step 6: Run test to verify it passes**

Run: `Rscript -e "devtools::test(filter='test-mom-multivariate')"`
Expected: PASS

---

## Task 4: Update Variance Component Extraction for Multivariate Mom

**Files:**
- Modify: `R/backend-mom.R:870-919` - Update extract_vc_mom function

**Step 1: Write failing test**

Add to test file:
```r
test_that("mom multivariate variance components have correct structure", {
  skip_on_cran()
  
  set.seed(123)
  data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  result <- gstudy(mvbind(y1, y2) ~ (1 | person) + (1 | rater), 
                   data = data, 
                   backend = "mom")
  
  vc <- result$variance_components
  
  # Should have dim column
  expect_true("dim" %in% names(vc))
  
  # Should have rows for each dimension
  expect_true(all(c("y1", "y2") %in% vc$dim))
  
  # Should have same components for each dimension
  vc_y1 <- vc[vc$dim == "y1", ]
  vc_y2 <- vc[vc$dim == "y2", ]
  expect_equal(nrow(vc_y1), nrow(vc_y2))
})
```

**Step 2: Run test to verify it fails**
Expected: FAIL

**Step 3: Update extract_vc_mom for multivariate**

The function should already work since we added `dim` column in fit_mom_multivariate. Verify and fix if needed.

---

## Task 5: Add Correlation Extraction to Mom Output

**Files:**
- Modify: `R/gstudy.R` - Extract correlations from mom model

**Step 1: Write failing test**

Add to test file:
```r
test_that("mom multivariate has correlations like brms", {
  skip_on_cran()
  
  set.seed(123)
  data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  result <- gstudy(mvbind(y1, y2) ~ (1 | person) + (1 | rater), 
                   data = data, 
                   backend = "mom")
  
  # Should have correlations
  expect_true("correlations" %in% names(result))
  expect_true(!is.null(result$correlations))
  
  # Should have residual correlations
  expect_true(!is.null(result$correlations$residual_cor))
})
```

**Step 2: Run test to verify it fails**
Expected: FAIL

**Step 3: Update gstudy.R to extract correlations from mom**

In gstudy.R, find the section where correlations are extracted (around lines 208-212):

Current code:
```r
# 13. Extract correlations for multivariate brms models
correlations <- NULL
if (is_mv && selected_backend == "brms") {
  correlations <- extract_correlations_brms(model)
}
```

Change to:
```r
# 13. Extract correlations for multivariate models
correlations <- NULL
if (is_mv) {
  if (selected_backend == "brms") {
    correlations <- extract_correlations_brms(model)
  } else if (selected_backend == "mom" && !is.null(model$correlations)) {
    correlations <- model$correlations
  }
}
```

**Step 4: Run test to verify it passes**

---

## Task 6: Ensure D-Study Works with Multivariate Mom

**Files:**
- Verify: `R/dstudy.R` - Check dstudy handles multivariate mom
- Modify: `tests/testthat/test-mom-multivariate.R` - Add dstudy tests

**Step 1: Write failing test**

Add to test file:
```r
test_that("dstudy works with multivariate mom", {
  skip_on_cran()
  
  set.seed(123)
  data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  g <- gstudy(mvbind(y1, y2) ~ (1 | person) + (1 | rater), 
              data = data, 
              backend = "mom")
  
  d <- dstudy(g, n = list(rater = 3))
  
  # Check dstudy output
  expect_s3_class(d, "dstudy")
  expect_true(d$is_multivariate)
  
  # Check coefficients have dim column
  expect_true("dim" %in% names(d$coefficients))
  
  # Should have coefficients for each dimension
  expect_equal(sort(unique(d$coefficients$dim)), c("y1", "y2"))
  
  # Should have G and D coefficients
  expect_true("g" %in% names(d$coefficients))
  expect_true("phi" %in% names(d$coefficients))
})
```

**Step 2: Run test to verify it fails**
Expected: FAIL (may need dstudy modifications)

**Step 3: Check dstudy.R for multivariate handling**

Examine how dstudy handles multivariate - it should work automatically if variance components have `dim` column. If not, fix accordingly.

**Step 4: Run test to verify it passes**

---

## Task 7: Add Print and Summary Methods Support

**Files:**
- Modify: `R/methods.R` - Ensure print.mgstudy and summary.mgstudy work

**Step 1: Write failing test**

Add to test file:
```r
test_that("print.mgstudy works for mom", {
  skip_on_cran()
  
  set.seed(123)
  data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  result <- gstudy(mvbind(y1, y2) ~ (1 | person) + (1 | rater), 
                   data = data, 
                   backend = "mom")
  
  expect_output(print(result), "Multivariate Generalizability Study")
  expect_output(print(result), "Dimensions:")
})

test_that("summary.mgstudy works for mom", {
  skip_on_cran()
  
  set.seed(123)
  data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  result <- gstudy(mvbind(y1, y2) ~ (1 | person) + (1 | rater), 
                   data = data, 
                   backend = "mom")
  
  expect_output(summary(result), "Multivariate G-Study Summary")
  expect_output(summary(result), "Dimensions:")
})
```

**Step 2: Run test to verify**
Expected: Should work if print/summary methods handle mgstudy class generically

---

## Task 8: Comprehensive Integration Tests

**Files:**
- Modify: `tests/testthat/test-mom-multivariate.R` - Add full test suite

**Step 1: Add more comprehensive tests**

```r
test_that("mom multivariate compares with brms for balanced data", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(123)
  n_person <- 10
  n_rater <- 4
  n_obs <- n_person * n_rater
  
  data <- data.frame(
    y1 = rnorm(n_obs, mean = 0, sd = 2),
    y2 = rnorm(n_obs, mean = 0, sd = 1.5),
    person = factor(rep(1:n_person, each = n_rater)),
    rater = factor(rep(1:n_rater, times = n_person))
  )
  
  # Fit with brms
  g_brms <- gstudy(mvbind(y1, y2) ~ (1 | person) + (1 | rater), 
                   data = data, 
                   backend = "brms",
                   chains = 2, iter = 1000)
  
  # Fit with mom  
  g_mom <- gstudy(mvbind(y1, y2) ~ (1 | person) + (1 | rater), 
                  data = data, 
                  backend = "mom")
  
  # Compare variance components (should be similar for balanced data)
  vc_brms <- g_brms$variance_components
  vc_mom <- g_mom$variance_components
  
  for (dim in c("y1", "y2")) {
    for (comp in c("person", "rater", "Residual")) {
      var_brms <- vc_brms$var[vc_brms$dim == dim & vc_brms$component == comp]
      var_mom <- vc_mom$var[vc_mom$dim == dim & vc_mom$component == comp]
      
      if (length(var_brms) > 0 && length(var_mom) > 0) {
        # Should be reasonably close (within 50% for this test)
        ratio <- var_mom / var_brms
        expect_true(ratio > 0.5 && ratio < 2.0,
          info = paste0("Variance ratio for ", comp, " in ", dim, 
                       " is ", round(ratio, 2), " (brms: ", round(var_brms, 3),
                       ", mom: ", round(var_mom, 3), ")"))
      }
    }
  }
})
```

---

## Task 9: Final Verification and Documentation

**Step 1: Run full test suite**

Run: `Rscript -e "devtools::test(filter='test-mom-multivariate')"`

**Step 2: Run existing tests to ensure no regressions**

Run: `Rscript -e "devtools::test()"`

**Step 3: Check documentation**

Verify that help pages work:
```r
library(facet)
?gstudy
?dstudy
```

---

## Plan Complete

Save this plan to: `docs/plans/2026-03-10-multivariate-mom-implementation-plan.md`
