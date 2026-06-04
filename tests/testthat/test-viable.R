# Tests for prmse viability analysis functionality (include_composite = TRUE)

# =============================================================================
# Output structure tests
# =============================================================================

test_that("prmse with include_composite returns tibble with correct columns (no ci)", {
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

  result <- prmse(d, include_composite = TRUE)
  expect_s3_class(result, "tbl_df")
  expected_cols <- c("dim", "prmse_c_rel", "prmse_c_abs", "var_rel", "var_abs")
  expect_true(all(expected_cols %in% names(result)))
  expect_false("prmse_c_rel_LL" %in% names(result))
  expect_false("var_rel_LL" %in% names(result))
})

test_that("prmse with include_composite returns CIs for prmse when ci='prmse'", {
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

  result <- prmse(d, include_composite = TRUE, ci = "prmse")
  expect_true("prmse_c_rel_LL" %in% names(result))
  expect_true("prmse_c_rel_UL" %in% names(result))
  expect_true("prmse_c_abs_LL" %in% names(result))
  expect_true("prmse_c_abs_UL" %in% names(result))
  expect_false("var_rel_LL" %in% names(result))
})

test_that("prmse with include_composite returns CIs for var when ci='var'", {
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

  result <- prmse(d, include_composite = TRUE, ci = "var")
  expect_true("var_rel_LL" %in% names(result))
  expect_true("var_rel_UL" %in% names(result))
  expect_true("var_abs_LL" %in% names(result))
  expect_true("var_abs_UL" %in% names(result))
  expect_false("prmse_c_rel_LL" %in% names(result))
})

test_that("prmse with include_composite returns both metric CIs when both specified", {
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

  result <- prmse(d, include_composite = TRUE, ci = c("prmse", "var"))
  expect_true("prmse_c_rel_LL" %in% names(result))
  expect_true("var_rel_LL" %in% names(result))
})

test_that("prmse with include_composite CI values are ordered correctly", {
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

  result <- prmse(d, include_composite = TRUE, ci = c("prmse", "var"))
  expect_true(all(result$var_rel_LL < result$var_rel_UL))
  expect_true(all(result$prmse_c_rel_LL < result$prmse_c_rel_UL))
})

# =============================================================================
# Alternative weights tests
# =============================================================================

test_that("prmse recalculates with alternative weights", {
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

  result_default <- prmse(d, include_composite = TRUE)
  result_custom <- prmse(d, include_composite = TRUE, weights = c(A = 2, B = 1))

  expect_false(isTRUE(all.equal(result_default$var_rel, result_custom$var_rel)))
})

test_that("prmse with custom weights uses posterior draws", {
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

  result <- prmse(d, include_composite = TRUE, weights = c(A = 2, B = 1), ci = "var")
  expect_s3_class(result, "tbl_df")
  expect_true("var_rel_LL" %in% names(result))
  expect_true(all(result$var_rel > 0))
})

# =============================================================================
# Sweep mode tests
# =============================================================================

test_that("prmse with include_composite works with sweep dstudy objects", {
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

  result <- prmse(d_sweep, include_composite = TRUE)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_true(all(c("dim", "prmse_c_rel", "prmse_c_abs", "var_rel", "var_abs") %in% names(result)))
})

test_that("prmse sweep uses actual sample sizes from gstudy", {
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

  result_sweep <- prmse(d_sweep, include_composite = TRUE)
  result_actual <- prmse(d_actual, include_composite = TRUE)

  expect_equal(result_sweep$var_rel, result_actual$var_rel, tolerance = 0.1)
  expect_equal(result_sweep$prmse_c_rel, result_actual$prmse_c_rel, tolerance = 0.1)
})

test_that("prmse sweep with CIs produces CI columns", {
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

  result <- prmse(d_sweep, include_composite = TRUE, ci = c("prmse", "var"))
  expect_true("prmse_c_rel_LL" %in% names(result))
  expect_true("var_rel_LL" %in% names(result))
  expect_true(all(result$var_rel_LL < result$var_rel_UL))
})

test_that("prmse works with custom probs", {
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

  result_95 <- prmse(d, include_composite = TRUE, ci = "var", probs = c(0.025, 0.975))
  result_90 <- prmse(d, include_composite = TRUE, ci = "var", probs = c(0.05, 0.95))

  width_95 <- result_95$var_rel_UL - result_95$var_rel_LL
  width_90 <- result_90$var_rel_UL - result_90$var_rel_LL

  expect_true(all(width_90 < width_95))
})

test_that("prmse returns positive VAR values", {
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

  result <- prmse(d, include_composite = TRUE)
  expect_true(all(result$var_rel > 0))
  expect_true(all(result$var_abs > 0))
})

# =============================================================================
# Optimization tests
# =============================================================================

test_that("prmse validates optimize parameter", {
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

  expect_error(prmse(d, optimize = "invalid"), "'arg' should be one of")
})

test_that("prmse validates optimize_target parameter", {
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

  expect_error(prmse(d, optimize = "composite", optimize_target = "invalid"),
               "'arg' should be one of")
})

test_that("prmse validates grid_resolution parameter", {
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

  expect_error(prmse(d, optimize = "tuning", grid_resolution = 0),
               "must be between")
  expect_error(prmse(d, optimize = "tuning", grid_resolution = 1),
               "must be between")
})

test_that("prmse composite optimization returns valid structure", {
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

  result <- prmse(d, optimize = "composite")

  expect_true(is.list(result))
  expect_equal(result$method, "composite")
  expect_true(!is.null(result$weights))
  expect_equal(sum(result$weights), 1, tolerance = 1e-10)
  expect_true(all(result$weights >= 0))
  expect_true(!is.null(result$metrics))
  expect_true(!is.null(result$composite_reliability))
})

test_that("prmse composite weights improve reliability", {
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

  result_opt <- prmse(d, optimize = "composite")

  k <- length(d$gstudy$dimensions)
  result_equal <- recalculate_var_with_weights(
    d, weights = rep(1/k, k), dims = d$gstudy$dimensions,
    ci = NULL, probs = c(0.025, 0.975), n = NULL
  )

  expect_gte(result_opt$composite_reliability, 0)
})

test_that("prmse subscale optimization returns valid structure", {
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

  result <- prmse(d, optimize = "subscale")

  expect_true(is.list(result))
  expect_equal(result$method, "subscale")
  expect_true(!is.null(result$minimax))
  expect_true(!is.null(result$per_subscale))
  expect_true(!is.null(result$comparison))
})

test_that("prmse subscale with target returns single optimization", {
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

  result <- prmse(d, optimize = "subscale", subscale = "A")

  expect_true(is.list(result))
  expect_equal(result$method, "subscale")
  expect_equal(result$target, "A")
  expect_true(!is.null(result$weights))
  expect_equal(sum(result$weights), 1, tolerance = 1e-10)
})

test_that("prmse tuning optimization returns valid structure", {
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

  result <- prmse(d, optimize = "tuning", grid_resolution = 0.2)

  expect_true(is.list(result))
  expect_equal(result$method, "tuning")
  expect_true(!is.null(result$best))
  expect_true(!is.null(result$grid_results))
  expect_true(!is.null(result$best$weights))
  expect_equal(sum(result$best$weights), 1, tolerance = 1e-10)
  expect_true(result$resolution == 0.2)
})

test_that("prmse tuning grid_resolution affects grid size", {
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

  result_01 <- suppressWarnings(prmse(d, optimize = "tuning", grid_resolution = 0.1))
  result_02 <- suppressWarnings(prmse(d, optimize = "tuning", grid_resolution = 0.2))

  expect_gt(nrow(result_01$grid_results), nrow(result_02$grid_results))
})

test_that("prmse optimize_target affects subscale results", {
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

  result_rel <- suppressWarnings(prmse(d, optimize = "composite", optimize_target = "rel"))
  result_abs <- suppressWarnings(prmse(d, optimize = "composite", optimize_target = "abs"))

  expect_equal(sum(result_rel$weights), 1, tolerance = 1e-10)
  expect_equal(sum(result_abs$weights), 1, tolerance = 1e-10)
})
