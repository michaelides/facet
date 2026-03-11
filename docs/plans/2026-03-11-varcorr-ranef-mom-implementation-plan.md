# VarCorr and ranef for mom Backend Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement VarCorr and ranef methods for gstudy/mgstudy objects with mom backend, returning list of variance matrices and means-based BLUPs respectively.

**Architecture:** Add mom backend handling to existing VarCorr.gstudy and ranef.gstudy functions in R/methods.R. Create helper functions to extract variance-covariance matrices from momfit objects and compute empirical BLUPs.

**Tech Stack:** R, lme4 (for reference API), method of moments backend

---

### Task 1: Write failing tests for VarCorr.mom

**Files:**
- Modify: `tests/testthat/test-gstudy.R`

**Step 1: Add test for VarCorr with univariate mom backend**

Add this test after existing VarCorr tests (around line 204):

```r
test_that("VarCorr.gstudy works with univariate mom backend", {
  result <- gstudy(score ~ (1 | person) + (1 | item), 
                   data = test_data, backend = "mom")
  vc <- VarCorr.gstudy(result)
  expect_type(vc, "list")
  expect_true("person" %in% names(vc))
  expect_true("item" %in% names(vc))
  # Check matrix structure
  expect_type(vc[["person"]], "double")
  expect_equal(dim(vc[["person"]]), c(1, 1))
})
```

**Step 2: Run test to verify it fails**

Run: `devtools::test(filter = "mom")` or manually run test
Expected: FAIL with "Unknown backend: mom"

**Step 3: Commit**
```bash
git add tests/testthat/test-gstudy.R
git commit -m "test: add failing test for VarCorr with mom backend"
```

---

### Task 2: Implement VarCorr for univariate mom

**Files:**
- Modify: `R/methods.R:1000-1036`

**Step 1: Add mom case to VarCorr.gstudy**

Replace the else clause in VarCorr.gstudy (line 1033-1035):

```r
} else if (x$backend == "mom") {
  # VarCorr for method of moments backend
  model <- x$model
  
  # Extract variance components
  vc <- model$variance_components
  
  # Filter to non-residual components
  random_effects <- vc[vc$component != "Residual", ]
  
  # Build list of variance matrices
  result <- list()
  for (i in seq_len(nrow(random_effects))) {
    comp_name <- random_effects$component[i]
    var_est <- random_effects$var[i]
    # Ensure non-negative
    var_est <- max(0, var_est)
    result[[comp_name]] <- matrix(var_est, nrow = 1, ncol = 1,
                                   dimnames = list(NULL, comp_name))
  }
  
  # Add residual variance
  resid_vc <- vc[vc$component == "Residual", ]
  if (nrow(resid_vc) > 0) {
    resid_var <- max(0, resid_vc$var[1])
    result[["Residual"]] <- matrix(resid_var, nrow = 1, ncol = 1,
                                    dimnames = list(NULL, "Residual"))
  }
  
  return(result)
}
```

**Step 2: Run test to verify it passes**

Run test again
Expected: PASS

**Step 3: Commit**
```bash
git add R/methods.R
git commit -m "feat: add VarCorr support for univariate mom backend"
```

---

### Task 3: Write failing test for VarCorr multivariate mom

**Files:**
- Modify: `tests/testthat/test-gstudy.R`

**Step 1: Add test for multivariate mom**

```r
test_that("VarCorr.gstudy works with multivariate mom backend", {
  skip_if_not_installed("brms")  # For mvbind syntax
  test_data_mv <- data.frame(
    score1 = rnorm(100),
    score2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:10, each = 10))
  )
  result <- gstudy(mvbind(score1, score2) ~ (1 | person) + (1 | item), 
                   data = test_data_mv, backend = "mom")
  vc <- VarCorr.gstudy(result)
  expect_type(vc, "list")
  expect_true("score1" %in% names(vc))
  expect_true("score2" %in% names(vc))
  expect_true("person" %in% names(vc[["score1"]]))
})
```

**Step 2: Run test to verify it fails**
Expected: FAIL (current code doesn't handle multivariate mom)

**Step 3: Commit**
```bash
git add tests/testthat/test-gstudy.R
git commit -m "test: add failing test for VarCorr with multivariate mom"
```

---

### Task 4: Implement VarCorr for multivariate mom

**Files:**
- Modify: `R/methods.R`

**Step 1: Modify VarCorr.gstudy to handle multivariate mom**

In VarCorr.gstudy, after the univariate mom case, add:

```r
} else if (x$backend == "mom" && x$is_multivariate) {
  # Multivariate method of moments
  model <- x$model
  vc <- model$variance_components
  dimensions <- x$dimensions
  random_facets <- model$random_facets
  
  result <- list()
  
  for (dim in dimensions) {
    # Filter VC for this dimension
    vc_dim <- vc[vc$dim == dim, ]
    random_effects <- vc_dim[vc_dim$component != "Residual", ]
    
    dim_result <- list()
    
    for (i in seq_len(nrow(random_effects))) {
      comp_name <- random_effects$component[i]
      var_est <- max(0, random_effects$var[i])
      dim_result[[comp_name]] <- matrix(var_est, nrow = 1, ncol = 1,
                                         dimnames = list(NULL, comp_name))
    }
    
    # Add residual
    resid_vc <- vc_dim[vc_dim$component == "Residual", ]
    if (nrow(resid_vc) > 0) {
      resid_var <- max(0, resid_vc$var[1])
      dim_result[["Residual"]] <- matrix(resid_var, nrow = 1, ncol = 1,
                                          dimnames = list(NULL, "Residual"))
    }
    
    result[[dim]] <- dim_result
  }
  
  # Add correlation info if available
  if (!is.null(model$correlations)) {
    result$random_effect_cor <- model$correlations$random_effect_cor
    result$residual_cor <- model$correlations$residual_cor
  }
  
  return(result)
}
```

**Step 2: Run test to verify it passes**
Expected: PASS

**Step 3: Commit**
```bash
git add R/methods.R
git commit -m "feat: add VarCorr support for multivariate mom backend"
```

---

### Task 5: Write failing test for ranef mom

**Files:**
- Modify: `tests/testthat/test-gstudy.R`

**Step 1: Add test for ranef with mom backend**

```r
test_that("ranef.gstudy works with univariate mom backend", {
  result <- gstudy(score ~ (1 | person) + (1 | item), 
                   data = test_data, backend = "mom")
  re <- ranef.gstudy(result)
  expect_type(re, "list")
  expect_true("person" %in% names(re))
  expect_true("item" %in% names(re))
  # Check structure - should be named vectors
  expect_type(re[["person"]], "double")
})
```

**Step 2: Run test to verify it fails**
Expected: FAIL (current code returns NULL with message)

**Step 3: Commit**
```bash
git add tests/testthat/test-gstudy.R
git commit -m "test: add failing test for ranef with mom backend"
```

---

### Task 6: Implement ranef for univariate mom

**Files:**
- Modify: `R/methods.R:1070-1087`

**Step 1: Replace the mom case in ranef.gstudy**

Replace lines 1081-1083:
```r
} else if (object$backend == "mom") {
  # Empirical BLUPs for method of moments
  model <- object$model
  data <- object$data
  response <- model$response
  random_facets <- model$random_facets
  
  # Compute grand mean
  grand_mean <- mean(data[[response]], na.rm = TRUE)
  
  result <- list()
  
  for (facet in random_facets) {
    if (!facet %in% names(data)) next
    
    # Compute means by facet level
    facet_means <- aggregate(
      data[[response]],
      by = list(facet = data[[facet]]),
      FUN = mean, na.rm = TRUE
    )
    
    # Empirical BLUPs = level mean - grand mean
    eblup <- setNames(
      facet_means$x - grand_mean,
      facet_means$facet
    )
    
    result[[facet]] <- eblup
  }
  
  return(result)
}
```

**Step 2: Run test to verify it passes**
Expected: PASS

**Step 3: Commit**
```bash
git add R/methods.R
git commit -m "feat: add ranef support for univariate mom backend"
```

---

### Task 7: Write test for ranef multivariate mom

**Files:**
- Modify: `tests/testthat/test-gstudy.R`

**Step 1: Add test**

```r
test_that("ranef.gstudy works with multivariate mom backend", {
  test_data_mv <- data.frame(
    score1 = rnorm(100),
    score2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:10, each = 10))
  )
  result <- gstudy(mvbind(score1, score2) ~ (1 | person) + (1 | item), 
                   data = test_data_mv, backend = "mom")
  re <- ranef.gstudy(result)
  expect_type(re, "list")
  expect_true("score1" %in% names(re))
  expect_true("score2" %in% names(re))
  expect_true("person" %in% names(re[["score1"]]))
})
```

**Step 2: Run test to verify it fails**
Expected: FAIL

**Step 3: Commit**
```bash
git add tests/testthat/test-gstudy.R
git commit -m "test: add failing test for ranef multivariate mom"
```

---

### Task 8: Implement ranef for multivariate mom

**Files:**
- Modify: `R/methods.R`

**Step 1: Add multivariate mom case in ranef.gstudy**

After the univariate mom case, add:

```r
} else if (object$backend == "mom" && object$is_multivariate) {
  # Multivariate method of moments - empirical BLUPs
  model <- object$model
  data <- object$data
  responses <- model$responses
  random_facets <- model$random_facets
  
  result <- list()
  
  for (resp in responses) {
    # Compute grand mean for this response
    grand_mean <- mean(data[[resp]], na.rm = TRUE)
    
    resp_result <- list()
    
    for (facet in random_facets) {
      if (!facet %in% names(data)) next
      
      # Compute means by facet level
      facet_means <- aggregate(
        data[[resp]],
        by = list(facet = data[[facet]]),
        FUN = mean, na.rm = TRUE
      )
      
      # Empirical BLUPs
      eblup <- setNames(
        facet_means$x - grand_mean,
        facet_means$facet
      )
      
      resp_result[[facet]] <- eblup
    }
    
    result[[resp]] <- resp_result
  }
  
  return(result)
}
```

**Step 2: Run test to verify it passes**
Expected: PASS

**Step 3: Commit**
```bash
git add R/methods.R
git commit -m "feat: add ranef support for multivariate mom backend"
```

---

### Task 9: Run full test suite

**Step 1: Run all gstudy tests**

```bash
devtools::test(filter = "gstudy")
```

Expected: All tests pass

**Step 2: Commit**
```bash
git commit -m "test: run full test suite for gstudy"
```

---

**Plan complete.**
