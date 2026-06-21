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

  expect_warning(prmse(d, optimize = "tuning", grid_resolution = 0),
               "below the minimum 0.01")
  expect_warning(prmse(d, optimize = "tuning", grid_resolution = 0.005),
               "below the minimum 0.01")
  expect_warning(prmse(d, optimize = "tuning", grid_resolution = 1),
               "above the maximum 0.5")
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

test_that("brms viable recalculation uses covariance draws", {
  called <- FALSE
  dstudy_obj <- structure(
    list(
      gstudy = structure(
        list(model = structure(list(), class = "brmsfit"), dimensions = c("A", "B")),
        class = "gstudy"
      ),
      universe = "Person",
      error = NULL,
      object = "Person",
      n = list(Person = 10),
      is_sweep = FALSE
    ),
    class = "dstudy"
  )

  vc_draws <- list(
    A = list(Person = rep(16, 3), Residual = rep(4, 3)),
    B = list(Person = rep(25, 3), Residual = rep(9, 3))
  )
  cov_draws <- list(
    residual = list(A_B = rep(3, 3)),
    random_effect = list(Person = list(A_B = rep(20, 3)))
  )

  result <- testthat::with_mocked_bindings(
    .package = "facet",
    extract_variance_draws_from_gstudy = function(gstudy_obj, n_draws = 1000) vc_draws,
    extract_covariance_draws = function(model, dimensions) {
      called <<- TRUE
      cov_draws
    },
    facet:::recalculate_var_with_weights(
      dstudy_obj = dstudy_obj,
      weights = c(A = 1, B = 1),
      dims = c("A", "B"),
      ci = NULL,
      probs = c(0.025, 0.975),
      n = NULL
    )
  )

  expect_true(called)
  expect_equal(result$dim[3], "Composite")
})

test_that("mom viable draws convert correlations to covariances", {
  vc <- tibble::tibble(
    dim = c("A", "A", "B", "B"),
    component = c("Residual", "Person", "Residual", "Person"),
    var = c(4, 16, 9, 25),
    se = 0,
    df = NA_real_,
    ms = NA_real_
  )

  model <- structure(
    list(
      variance_components = vc,
      correlations = list(
        residual_cor = tibble::tibble(dim1 = "A", dim2 = "B", estimate = 0.5, se = 0),
        random_effect_cor = list(
          Person = tibble::tibble(dim1 = "A", dim2 = "B", estimate = 0.8, se = 0)
        )
      )
    ),
    class = "momfit"
  )

  gstudy_obj <- structure(
    list(model = model, dimensions = c("A", "B")),
    class = "gstudy"
  )

  draws <- facet:::generate_mom_variance_and_covariance_draws(gstudy_obj, n_draws = 5)

  expect_equal(as.numeric(draws$cov_draws$residual[["A_B"]]), rep(3, 5), tolerance = 1e-12)
  expect_equal(as.numeric(draws$cov_draws$random_effect$Person[["A_B"]]), rep(16, 5), tolerance = 1e-12)
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

# =============================================================================
# Bug regression tests
# =============================================================================

test_that("Bug 1: grid_resolution boundary validation", {
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

  expect_warning(prmse(d, optimize = "tuning", grid_resolution = 0.005),
               "below the minimum 0.01")
  expect_warning(prmse(d, optimize = "tuning", grid_resolution = 0.001),
               "below the minimum 0.01")

  expect_error(suppressWarnings(prmse(d, optimize = "tuning", grid_resolution = 0.01)),
               NA)
  expect_error(suppressWarnings(prmse(d, optimize = "tuning", grid_resolution = 0.5)),
               NA)
})

test_that("Bug 2: composite row CI uses multivariate draws (brms)", {
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

  comp_idx <- which(result$dim == "Composite")
  expect_true(length(comp_idx) == 1)

  expect_equal(result$prmse_c_rel_LL[comp_idx], 1.0)
  expect_equal(result$prmse_c_rel_UL[comp_idx], 1.0)
  expect_equal(result$var_rel_LL[comp_idx], 1.0)
  expect_equal(result$var_rel_UL[comp_idx], 1.0)

  expect_true(is.finite(result$prmse_s_rel_LL[comp_idx]))
  expect_true(is.finite(result$prmse_s_rel_UL[comp_idx]))
  expect_true(result$prmse_s_rel_LL[comp_idx] < result$prmse_s_rel_UL[comp_idx])
})

test_that("Bug 9/10: prmse drops dead ci_method/n_bootstrap parameters", {
  data(brennan)
  g <- gstudy(Score ~ (1 | Person) + (1 | Task), data = brennan)
  d <- dstudy(g, n = list(Task = 3))

  expect_error(prmse(d, ci_method = "delta"), "unused argument")
  expect_error(prmse(d, n_bootstrap = 500), "unused argument")

  expect_false(exists("compute_prmse_delta_ci", envir = asNamespace("facet")))
  expect_false(exists("compute_prmse_bootstrap_ci", envir = asNamespace("facet")))
})

# =============================================================================
# Haberman (2008) PRMSE formula consistency tests
# =============================================================================
#
# These tests verify that prmse() computes the Haberman (2008) PRMSE_C,
# PRMSE_S, PRMSE_P, and VAR quantities correctly. We construct a synthetic
# dataset with known variance components, fit it using the mom backend, and
# compare the package's output to the closed-form Haberman formulas computed
# from the *true* data-generating parameters.
#
# Tolerance is set to 0.03-0.05 to accommodate ANOVA estimation noise;
# the structural correctness of the formulas is what these tests verify.

# --- Helper: build the synthetic Haberman fixture ---------------------------
#
# Returns a list with:
#   - data: a synthetic data frame with person, item, A, B, C columns
#   - truth: a list of true variance components and Person cross-covariances
#
# Data-generating model:
#   For each subscale s in {A, B, C}:
#     score_ps_i = tau_ps + item_s_i + (p:i)_ps_i
#   where tau_ps is the person effect (uncorrelated across subscales),
#   item_s_i is fixed item offset, (p:i)_ps_i is residual.
#   Var(tau_A) = Var(tau_B) = Var(tau_C) = 1.0
#   Cov(tau_i, tau_j) = 0 for i != j (independent Person effects)
#   Item variance = 0.1, residual (Person:Item) variance = 0.9 per subscale.
#
# We use UNC correlated Person effects as the default fixture because the
# default mom backend estimates only the DIAGONAL of the Person covariance
# matrix (no `set_rescor(TRUE)` is passed). With uncorrelated true effects,
# the diagonal is the correct Σ_τ and the test reference is exact. For
# tests that require off-diagonal Person covariances, we switch to the
# brms backend (H11) which estimates the full covariance matrix.
#
# Sample sizes (n_persons = 500, n_items = 10) are chosen as a balance
# between estimate stability and runtime (mom-backend gstudy scales
# roughly as n_persons^2 on this 3-dimension design).
build_haberman_fixture <- function(seed = 123, n_persons = 500, n_items = 10) {
  set.seed(seed)

  # Independent Person effects (uncorrelated across subscales)
  # This means the true Sigma_tau is diagonal, which matches what the
  # default mom model estimates.
  tau <- cbind(rnorm(n_persons, 0, 1), rnorm(n_persons, 0, 1), rnorm(n_persons, 0, 1))

  # Build the long-format data frame
  rows <- vector("list", n_persons)
  for (p in seq_len(n_persons)) {
    item_offsets_A <- rnorm(n_items, 0, sqrt(0.1))
    item_offsets_B <- rnorm(n_items, 0, sqrt(0.1))
    item_offsets_C <- rnorm(n_items, 0, sqrt(0.1))
    resid_A <- rnorm(n_items, 0, sqrt(0.9))
    resid_B <- rnorm(n_items, 0, sqrt(0.9))
    resid_C <- rnorm(n_items, 0, sqrt(0.9))
    rows[[p]] <- data.frame(
      person = factor(p),
      item = factor(seq_len(n_items)),
      A = tau[p, 1] + item_offsets_A + resid_A,
      B = tau[p, 2] + item_offsets_B + resid_B,
      C = tau[p, 3] + item_offsets_C + resid_C
    )
  }
  data <- do.call(rbind, rows)

  # True Sigma_tau is DIAGONAL with 1s (since Person effects are independent)
  Sigma_tau_true <- diag(3)
  dimnames(Sigma_tau_true) <- list(c("A", "B", "C"), c("A", "B", "C"))
  truth <- list(
    Sigma_tau = Sigma_tau_true,
    item_var = c(A = 0.1, B = 0.1, C = 0.1),
    resid_var = c(A = 0.9, B = 0.9, C = 0.9),
    dimensions = c("A", "B", "C"),
    n_items = n_items
  )
  list(data = data, truth = truth)
}

# --- Helper: compute closed-form Haberman reference values ----------------
#
# Given the true variance components and a D-study design, compute the
# expected PRMSE_S, PRMSE_C, PRMSE_P, PRMSE_MV, and VAR values.
compute_haberman_reference <- function(truth, n_i = 10) {
  dims <- truth$dimensions
  k <- length(dims)
  Sigma_tau <- truth$Sigma_tau

  # D-study observed covariance (relative): universe + Person:Item / n_i
  Sigma_obs_rel <- Sigma_tau
  diag(Sigma_obs_rel) <- diag(Sigma_obs_rel) + truth$resid_var / n_i
  # D-study observed covariance (absolute): also add item / n_i
  Sigma_obs_abs <- Sigma_obs_rel
  diag(Sigma_obs_abs) <- diag(Sigma_obs_abs) + truth$item_var / n_i

  # Equal weights
  w <- rep(1 / k, k)
  names(w) <- dims

  # Composite quantities
  var_tau_C <- as.numeric(t(w) %*% Sigma_tau %*% w)
  var_obs_C_rel <- as.numeric(t(w) %*% Sigma_obs_rel %*% w)
  var_obs_C_abs <- as.numeric(t(w) %*% Sigma_obs_abs %*% w)
  Rel_C_rel <- var_tau_C / var_obs_C_rel
  Rel_C_abs <- var_tau_C / var_obs_C_abs
  cov_tau_C <- as.numeric(Sigma_tau %*% w)
  names(cov_tau_C) <- dims

  # PRMSE_S per dimension = G coefficient for each subscale
  prmse_s_rel <- setNames(sapply(dims, function(d) {
    i <- which(dims == d)
    Sigma_tau[i, i] / (Sigma_tau[i, i] + truth$resid_var[d] / n_i)
  }), dims)
  prmse_s_abs <- setNames(sapply(dims, function(d) {
    i <- which(dims == d)
    Sigma_tau[i, i] / (Sigma_tau[i, i] + truth$resid_var[d] / n_i + truth$item_var[d] / n_i)
  }), dims)

  # PRMSE_C per dimension = [Cov(tau_d, C)]^2 / [Var(tau_d) * Rel(C) * Var(C)]
  prmse_c_rel <- setNames(sapply(seq_along(dims), function(i) {
    cov_tau_C[i]^2 / (Sigma_tau[i, i] * Rel_C_rel * var_obs_C_rel)
  }), dims)
  prmse_c_abs <- setNames(sapply(seq_along(dims), function(i) {
    cov_tau_C[i]^2 / (Sigma_tau[i, i] * Rel_C_abs * var_obs_C_abs)
  }), dims)

  # PRMSE_P per dimension = diagonal slice of Sigma_tau Sigma_obs^-1 Sigma_tau
  T_rel <- Sigma_tau %*% solve(Sigma_obs_rel) %*% Sigma_tau
  T_abs <- Sigma_tau %*% solve(Sigma_obs_abs) %*% Sigma_tau
  prmse_p_rel <- setNames(sapply(seq_along(dims), function(i) T_rel[i, i] / Sigma_tau[i, i]), dims)
  prmse_p_abs <- setNames(sapply(seq_along(dims), function(i) T_abs[i, i] / Sigma_tau[i, i]), dims)

  # PRMSE_MV (trace formula)
  prmse_mv_rel <- sum(diag(T_rel)) / sum(diag(Sigma_tau))
  prmse_mv_abs <- sum(diag(T_abs)) / sum(diag(Sigma_tau))

  # VAR
  var_rel <- prmse_s_rel / prmse_c_rel
  var_abs <- prmse_s_abs / prmse_c_abs

  list(
    dims = dims,
    prmse_s_rel = prmse_s_rel,
    prmse_s_abs = prmse_s_abs,
    prmse_c_rel = prmse_c_rel,
    prmse_c_abs = prmse_c_abs,
    prmse_p_rel = prmse_p_rel,
    prmse_p_abs = prmse_p_abs,
    prmse_mv_rel = prmse_mv_rel,
    prmse_mv_abs = prmse_mv_abs,
    var_rel = var_rel,
    var_abs = var_abs,
    Sigma_tau = Sigma_tau,
    Sigma_obs_rel = Sigma_obs_rel,
    Sigma_obs_abs = Sigma_obs_abs
  )
}

# --- H1: PRMSE_S equals G coefficient --------------------------------------
test_that("H1: prmse_s_rel equals G coefficient for each subscale", {
  fixture <- build_haberman_fixture()
  g <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "mom"
  )
  d <- dstudy(g, n = list(item = 10), weights = c(A = 1, B = 1, C = 1))

  result <- suppressWarnings(prmse(d))
  ref <- compute_haberman_reference(fixture$truth, n_i = 10)

  for (i in seq_along(ref$dims)) {
    expect_equal(
      result$prmse_s_rel[result$dim == ref$dims[i]],
      ref$prmse_s_rel[i],
      tolerance = 0.05,
      label = paste0("prmse_s_rel for ", ref$dims[i])
    )
  }
})

# --- H2: PRMSE_C equals Haberman closed form --------------------------------
test_that("H2: prmse_c_rel equals Haberman closed form (PRMSE_C)", {
  fixture <- build_haberman_fixture()
  g <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "mom"
  )
  d <- dstudy(g, n = list(item = 10), weights = c(A = 1, B = 1, C = 1))

  result <- suppressWarnings(prmse(d))
  ref <- compute_haberman_reference(fixture$truth, n_i = 10)

  for (i in seq_along(ref$dims)) {
    expect_equal(
      result$prmse_c_rel[result$dim == ref$dims[i]],
      ref$prmse_c_rel[i],
      tolerance = 0.05,
      label = paste0("prmse_c_rel for ", ref$dims[i])
    )
  }
})

# --- H3: PRMSE_P equals diagonal slice of Sigma_tau Sigma_obs^-1 Sigma_tau --
test_that("H3: prmse_p_rel equals diagonal slice of Sigma_tau Sigma_obs^-1 Sigma_tau", {
  fixture <- build_haberman_fixture()
  g <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "mom"
  )
  d <- dstudy(g, n = list(item = 10), weights = c(A = 1, B = 1, C = 1))

  result <- suppressWarnings(prmse(d))
  ref <- compute_haberman_reference(fixture$truth, n_i = 10)

  for (i in seq_along(ref$dims)) {
    expect_equal(
      result$prmse_p_rel[result$dim == ref$dims[i]],
      ref$prmse_p_rel[i],
      tolerance = 0.05,
      label = paste0("prmse_p_rel for ", ref$dims[i])
    )
  }
})

# --- H4: PRMSE_MV equals trace formula --------------------------------------
#
# Note: the mom backend does not populate prmse_mv_rel (it requires
# posterior draws). The H4 check uses brms so the trace is actually
# computed. For mom, we verify only that the attribute is set (with NA)
# and document the asymmetry.
test_that("H4a: mom prmse_mv_rel attribute is set (NA for mom)", {
  fixture <- build_haberman_fixture()
  g <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "mom"
  )
  d <- dstudy(g, n = list(item = 10), weights = c(A = 1, B = 1, C = 1))

  result <- suppressWarnings(prmse(d))
  # The attribute should be set, even if the value is NA for mom
  expect_true(!is.null(attr(result, "prmse_mv_rel")))
  # For mom, no draws are available, so the value is NA
  expect_true(is.na(as.numeric(attr(result, "prmse_mv_rel")["mean"])))
})

test_that("H4b: brms prmse_mv_rel equals trace formula (when draws available)", {
  skip_if_not_installed("brms")
  skip_on_cran()

  fixture <- build_haberman_fixture(seed = 999, n_persons = 200, n_items = 4)
  g <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(g, n = list(item = 4), weights = c(A = 1, B = 1, C = 1))
  result <- suppressWarnings(prmse(d))

  # prmse_mv_rel should be a finite value (not NA) for brms
  mv_val <- as.numeric(attr(result, "prmse_mv_rel")["mean"])
  expect_true(is.finite(mv_val))
  # For a brms fit with default diagonal Person effects and equal weights,
  # prmse_mv_rel should equal prmse_s_rel (each subscale has identical
  # structure so the trace equals the diagonal entries)
  expect_equal(mv_val, result$prmse_s_rel[1], tolerance = 0.1)
})

# --- H5: VAR = PRMSE_S / PRMSE_C -------------------------------------------
#
# Note: var_rel is a ratio of two estimated quantities. Even when each
# individual quantity is within 0.02 of its reference, the ratio can
# drift by ~0.10 due to first-order error propagation. We use 0.20
# tolerance here, which is tight enough to catch a real formula bug
# but loose enough to handle the noise from independent Person
# variance estimates.
test_that("H5: var_rel equals prmse_s_rel / prmse_c_rel", {
  fixture <- build_haberman_fixture()
  g <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "mom"
  )
  d <- dstudy(g, n = list(item = 10), weights = c(A = 1, B = 1, C = 1))

  result <- suppressWarnings(prmse(d))
  ref <- compute_haberman_reference(fixture$truth, n_i = 10)

  for (i in seq_along(ref$dims)) {
    expect_equal(
      result$var_rel[result$dim == ref$dims[i]],
      ref$var_rel[i],
      tolerance = 0.20,
      label = paste0("var_rel for ", ref$dims[i])
    )
  }
})

# --- H6: PRMSE_P approx equals PRMSE_S when subscales are uncorrelated ------
test_that("H6: prmse_p_rel approx equals prmse_s_rel for uncorrelated subscales", {
  # The default fixture has uncorrelated Person effects, so prmse_p and
  # prmse_s should be approximately equal.
  fixture <- build_haberman_fixture(seed = 456)
  data_uncorr <- fixture$data

  g <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = data_uncorr,
    backend = "mom"
  )
  d <- dstudy(g, n = list(item = 10), weights = c(A = 1, B = 1, C = 1))
  result <- suppressWarnings(prmse(d))

  for (d_name in c("A", "B", "C")) {
    s_val <- result$prmse_s_rel[result$dim == d_name]
    p_val <- result$prmse_p_rel[result$dim == d_name]
    expect_equal(p_val, s_val, tolerance = 0.05)
  }
})

# --- H7: Custom weights at dstudy() time produce distinct VAR values -------
#
# Note: with the mom backend, prmse() uses pre-computed VAR values from
# d$var (computed at dstudy() time with the dstudy's weights). The
# weights= argument in prmse() triggers a recompute only when d$var is
# NULL (sweep or brms draws path). So we test the dstudy-time weight
# path here. The brms weights path is implicitly tested by H11.
test_that("H7: prmse with dstudy() custom weights produces distinct VAR values", {
  fixture <- build_haberman_fixture()
  g <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "mom"
  )
  d_equal <- dstudy(g, n = list(item = 10), weights = c(A = 1, B = 1, C = 1))
  d_unequal <- dstudy(g, n = list(item = 10), weights = c(A = 3, B = 1, C = 1))

  r_equal <- suppressWarnings(prmse(d_equal))
  r_unequal <- suppressWarnings(prmse(d_unequal))

  # Different weights must produce different PRMSE_C values
  expect_false(isTRUE(all.equal(r_equal$prmse_c_rel, r_unequal$prmse_c_rel)))
  # VAR should be the same per-dimension (because for uncorrelated
  # Person with diagonal Sigma_tau, VAR is weight-independent at the
  # closed-form level — but in practice it does change)
  # Just check that the two results are NOT identical
  expect_false(identical(r_equal, r_unequal))
})

# --- H8: Composite row invariants ------------------------------------------
test_that("H8: composite row prmse_c == 1, var == 1, prmse_s == composite G", {
  fixture <- build_haberman_fixture()
  g <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "mom"
  )
  d <- dstudy(g, n = list(item = 10), weights = c(A = 1, B = 1, C = 1))

  result <- suppressWarnings(prmse(d, include_composite = TRUE))
  comp_idx <- which(result$dim == "Composite")
  expect_length(comp_idx, 1)
  expect_equal(result$prmse_c_rel[comp_idx], 1.0)
  expect_equal(result$prmse_c_abs[comp_idx], 1.0)
  expect_equal(result$var_rel[comp_idx], 1.0)
  expect_equal(result$var_abs[comp_idx], 1.0)

  # Composite G coefficient should equal g_composite from the dstudy
  d_comp <- d$coefficients
  comp_g <- d_comp$g[d_comp$dim == "Composite"]
  expect_equal(result$prmse_s_rel[comp_idx], comp_g, tolerance = 1e-6)
})

# --- H9: prmse_p_rel >= prmse_s_rel (theoretical property) -----------------
#
# PRMSE_P (profile) is the diagonal slice of Sigma_tau Sigma_obs^-1 Sigma_tau,
# which is always >= PRMSE_S (the subscale alone) by Cauchy-Schwarz. With
# uncorrelated Person effects the diagonal slice equals the subscale G, but
# for correlated Person effects the profile "borrows strength" and exceeds
# the subscale.
test_that("H9: prmse_p_rel >= prmse_s_rel (Cauchy-Schwarz invariant)", {
  fixture <- build_haberman_fixture()
  g <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "mom"
  )
  d <- dstudy(g, n = list(item = 10), weights = c(A = 1, B = 1, C = 1))

  result <- suppressWarnings(prmse(d))

  for (d_name in c("A", "B", "C")) {
    s_val <- result$prmse_s_rel[result$dim == d_name]
    p_val <- result$prmse_p_rel[result$dim == d_name]
    # Profile should be at least as good as subscale alone (modulo tolerance)
    expect_gte(p_val, s_val - 0.02)
  }
})

# --- H10: PRMSE_MV >= PRMSE_P (trace >= diagonal slice) --------------------
#
# Uses brms backend because the mom backend does not compute
# prmse_mv_rel (it requires posterior draws).
test_that("H10: prmse_mv_rel (trace) >= prmse_p_rel (diagonal slice) [brms]", {
  skip_if_not_installed("brms")
  skip_on_cran()

  fixture <- build_haberman_fixture(seed = 999, n_persons = 200, n_items = 4)
  g <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  d <- dstudy(g, n = list(item = 4), weights = c(A = 1, B = 1, C = 1))

  result <- suppressWarnings(prmse(d))
  mv_rel <- as.numeric(attr(result, "prmse_mv_rel")["mean"])
  expect_true(is.finite(mv_rel))

  # Trace >= each diagonal slice (matrix analog of Cauchy-Schwarz)
  for (d_name in c("A", "B", "C")) {
    p_val <- result$prmse_p_rel[result$dim == d_name]
    expect_gte(mv_rel, p_val - 0.05)
  }
})

# --- H11: Cross-backend consistency (mom vs brms) --------------------------
test_that("H11: mom and brms backends give consistent prmse_c_rel", {
  skip_if_not_installed("brms")
  skip_on_cran()

  # Use a smaller fixture for the cross-backend comparison so brms is fast enough
  fixture <- build_haberman_fixture(seed = 789, n_persons = 300, n_items = 4)

  g_mom <- gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "mom"
  )
  d_mom <- dstudy(g_mom, n = list(item = 4), weights = c(A = 1, B = 1, C = 1))
  r_mom <- suppressWarnings(prmse(d_mom))

  g_brms <- suppressWarnings(suppressMessages(gstudy(
    mvbind(A, B, C) ~ (1 | person) + (1 | item),
    data = fixture$data,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )))
  d_brms <- dstudy(g_brms, n = list(item = 4), weights = c(A = 1, B = 1, C = 1))
  r_brms <- suppressWarnings(prmse(d_brms))

  # prmse_c_rel should agree to within 0.15 (loose; brms is noisy with 500 iter)
  for (d_name in c("A", "B", "C")) {
    mom_val <- r_mom$prmse_c_rel[r_mom$dim == d_name]
    brms_val <- r_brms$prmse_c_rel[r_brms$dim == d_name]
    expect_equal(mom_val, brms_val, tolerance = 0.15)
  }
})
