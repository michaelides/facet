# Tests for viable function

# =============================================================================
# Input validation tests
# =============================================================================

test_that("viable requires a dstudy object", {
  expect_error(viable(list()), "must be a dstudy object")
  expect_error(viable(NULL), "must be a dstudy object")
  expect_error(viable("not a dstudy"), "must be a dstudy object")
})

test_that("viable requires multivariate dstudy object", {
  skip_if_not_installed("lme4")
  data(brennan)
  g <- gstudy(Score ~ (1 | Person) + (1 | Task), data = brennan)
  d <- dstudy(g, n = list(Task = 3))
  expect_error(viable(d), "VAR is only computed for multivariate models")
})

test_that("viable requires posterior estimation", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("brms")
  
  set.seed(123)
  n_persons <- 20
  data <- data.frame(
    person = 1:n_persons,
    A = rnorm(n_persons),
    B = rnorm(n_persons)
  )
  
  g <- gstudy(mvbind(A, B) ~ (1 | person), data = data, backend = "mom")
  d <- dstudy(g, n = list(person = 10))
  
  expect_error(viable(d), "requires posterior estimation")
})

test_that("viable validates ci parameter", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(123)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  expect_error(viable(d, ci = "invalid"), "'arg' should be one of")
})

test_that("viable validates probs parameter", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(123)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  expect_error(viable(d, ci = "prmse", probs = c(0.1)), "must have exactly 2 elements")
  expect_error(viable(d, ci = "prmse", probs = c(0.9, 0.1)), "must be in increasing order")
  expect_error(viable(d, ci = "prmse", probs = c(-0.1, 0.9)), "must be between 0 and 1")
})

test_that("viable validates weights length", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(123)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  expect_error(viable(d, weights = c(1)), "must have length")
})

# =============================================================================
# Output structure tests
# =============================================================================

test_that("viable returns tibble with correct columns (no ci)", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(123)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  result <- viable(d)
  expect_s3_class(result, "tbl_df")
  expected_cols <- c("dim", "prmse_c_rel", "prmse_c_abs", "var_rel", "var_abs")
  expect_true(all(expected_cols %in% names(result)))
  expect_false("prmse_c_rel_LL" %in% names(result))
  expect_false("var_rel_LL" %in% names(result))
})

test_that("viable returns CIs for prmse when ci='prmse'", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(123)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  result <- viable(d, ci = "prmse")
  expect_true("prmse_c_rel_LL" %in% names(result))
  expect_true("prmse_c_rel_UL" %in% names(result))
  expect_true("prmse_c_abs_LL" %in% names(result))
  expect_true("prmse_c_abs_UL" %in% names(result))
  expect_false("var_rel_LL" %in% names(result))
})

test_that("viable returns CIs for var when ci='var'", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(123)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  result <- viable(d, ci = "var")
  expect_true("var_rel_LL" %in% names(result))
  expect_true("var_rel_UL" %in% names(result))
  expect_true("var_abs_LL" %in% names(result))
  expect_true("var_abs_UL" %in% names(result))
  expect_false("prmse_c_rel_LL" %in% names(result))
})

test_that("viable returns both metric CIs when both specified", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(123)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  result <- viable(d, ci = c("prmse", "var"))
  expect_true("prmse_c_rel_LL" %in% names(result))
  expect_true("var_rel_LL" %in% names(result))
})

test_that("viable CI values are ordered correctly", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(123)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  result <- viable(d, ci = c("prmse", "var"))
  expect_true(all(result$var_rel_LL < result$var_rel_UL))
  expect_true(all(result$prmse_c_rel_LL < result$prmse_c_rel_UL))
})

# =============================================================================
# Alternative weights tests
# =============================================================================

test_that("viable recalculates with alternative weights", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(456)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  result_default <- viable(d)
  result_custom <- viable(d, weights = c(A = 2, B = 1))
  
  expect_false(isTRUE(all.equal(result_default$var_rel, result_custom$var_rel)))
})

test_that("viable with custom weights uses posterior draws", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(456)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  result <- viable(d, weights = c(A = 2, B = 1), ci = "var")
  expect_s3_class(result, "tbl_df")
  expect_true("var_rel_LL" %in% names(result))
  expect_true(all(result$var_rel > 0))
})

# =============================================================================
# Integration tests
# =============================================================================

test_that("viable results match prmse() when weights match", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(789)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  viable_result <- viable(d, probs = c(0.025, 0.975), ci = c("prmse", "var"))
  prmse_result <- prmse(d, probs = c(0.025, 0.975))
  
  expect_equal(viable_result$var_rel, prmse_result$var_rel, tolerance = 1e-10)
  expect_equal(viable_result$prmse_c_rel, prmse_result$prmse_c_rel, tolerance = 1e-10)
})

# =============================================================================
# Sweep mode tests
# =============================================================================

test_that("viable works with sweep dstudy objects", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(123)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  
  d_sweep <- dstudy(gu, n = list(person = c(5, 10, 15)), weights = c(A = 1, B = 1))
  expect_true(d_sweep$is_sweep)
  
  result <- viable(d_sweep)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_true(all(c("dim", "prmse_c_rel", "prmse_c_abs", "var_rel", "var_abs") %in% names(result)))
})

test_that("viable sweep uses actual sample sizes from gstudy", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(456)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  
  actual_n <- gu$facet_n
  
  d_sweep <- dstudy(gu, n = list(person = c(5, 10, 15)), weights = c(A = 1, B = 1))
  
  d_actual <- dstudy(gu, n = as.list(actual_n), weights = c(A = 1, B = 1))
  
  result_sweep <- viable(d_sweep)
  result_actual <- viable(d_actual)
  
  expect_equal(result_sweep$var_rel, result_actual$var_rel, tolerance = 0.1)
  expect_equal(result_sweep$prmse_c_rel, result_actual$prmse_c_rel, tolerance = 0.1)
})

test_that("viable sweep with CIs produces CI columns", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(789)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  
  d_sweep <- dstudy(gu, n = list(person = c(5, 10, 15)), weights = c(A = 1, B = 1))
  
  result <- viable(d_sweep, ci = c("prmse", "var"))
  expect_true("prmse_c_rel_LL" %in% names(result))
  expect_true("var_rel_LL" %in% names(result))
  expect_true(all(result$var_rel_LL < result$var_rel_UL))
})

test_that("viable works with custom probs", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(789)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  result_95 <- viable(d, ci = "var", probs = c(0.025, 0.975))
  result_90 <- viable(d, ci = "var", probs = c(0.05, 0.95))
  
  width_95 <- result_95$var_rel_UL - result_95$var_rel_LL
  width_90 <- result_90$var_rel_UL - result_90$var_rel_LL
  
  expect_true(all(width_90 < width_95))
})

test_that("viable returns positive VAR values", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  set.seed(123)
  n_persons <- 30
  person_effect <- rnorm(n_persons, 0, 2)
  data <- data.frame(
    person = 1:n_persons,
    A = person_effect + rnorm(n_persons, 0, 1),
    B = 0.7 * person_effect + rnorm(n_persons, 0, 1)
  )
  
  gu <- gstudy(
    mvbind(A, B) ~ (1 | person),
    data = data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))
  
  result <- viable(d)
  expect_true(all(result$var_rel > 0))
  expect_true(all(result$var_abs > 0))
})
