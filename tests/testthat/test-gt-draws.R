# Tests for gt_draws function

# Helper function to create test data
create_test_data <- function(n_persons = 20, n_items = 10, seed = 123) {
  set.seed(seed)
  data.frame(
    person = factor(rep(1:n_persons, each = n_items)),
    item = factor(rep(1:n_items, n_persons)),
    score = rnorm(n_persons * n_items, mean = 50, sd = 10)
  )
}

# Test: gt_draws.gstudy basic structure (univariate)
test_that("gt_draws.gstudy returns correct structure for univariate model", {
  skip_if_not_installed("brms")

  test_data <- create_test_data()
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "brms", refresh = 0)

  result <- gt_draws(g)

  expect_s3_class(result, "tbl_df")
  expect_true("draw" %in% names(result))
  expect_true("person" %in% names(result))
  expect_true("item" %in% names(result))
  expect_true("Residual" %in% names(result))
  expect_true("uni" %in% names(result))
  expect_true("sigma2_delta" %in% names(result))
  expect_true("sigma2_delta_abs" %in% names(result))
  expect_true("g" %in% names(result))
  expect_true("phi" %in% names(result))
  expect_true("sem_rel" %in% names(result))
  expect_true("sem_abs" %in% names(result))

  expect_true(nrow(result) > 0)
})

# Test: gt_draws.gstudy filtering components
test_that("gt_draws.gstudy filters by components", {
  skip_if_not_installed("brms")

  test_data <- create_test_data()
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "brms", refresh = 0)

  result <- gt_draws(g, components = c("person", "Residual"))

  expect_s3_class(result, "tbl_df")
  expect_true("person" %in% names(result))
  expect_true("Residual" %in% names(result))
  expect_false("item" %in% names(result))
})

# Test: gt_draws.dstudy coefficients (univariate)
test_that("gt_draws.dstudy extracts coefficient draws for univariate", {
  skip_if_not_installed("brms")

  test_data <- create_test_data()
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "brms", refresh = 0)
  d <- dstudy(g, n = list(item = 5), estimation = "posterior")

  result <- gt_draws(d, what = "coefficients")

  expect_s3_class(result, "tbl_df")
  expect_true("draw" %in% names(result))
  expect_true("uni" %in% names(result))
  expect_true("sigma2_delta" %in% names(result))
  expect_true("g" %in% names(result))
  expect_true("phi" %in% names(result))

  expect_true(nrow(result) > 0)
})

# Test: gt_draws.dstudy coefficients match means (univariate)
test_that("gt_draws.dstudy coefficient means match summary", {
  skip_if_not_installed("brms")

  test_data <- create_test_data()
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "brms", refresh = 0)
  d <- dstudy(g, n = list(item = 5), estimation = "posterior")

  result <- gt_draws(d, what = "coefficients")

  expect_equal(mean(result$g), d$coefficients$g, tolerance = 0.001)
  expect_equal(mean(result$phi), d$coefficients$phi, tolerance = 0.001)
})

# Test: gt_draws.dstudy error for non-posterior
test_that("gt_draws.dstudy errors for non-posterior estimation", {
  test_data <- create_test_data()
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "lme4")
  d <- dstudy(g, n = list(item = 5), estimation = "simple")

  expect_error(gt_draws(d), "estimation = 'posterior'")
})

# Test: gt_draws.gstudy errors for non-brms
test_that("gt_draws.gstudy errors for non-brms estimator", {
  test_data <- create_test_data()
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "lme4")

  expect_error(gt_draws(g), "brms estimator")
})

# Test: gt_draws.gstudy variance values are positive
test_that("gt_draws.gstudy variance draws are positive", {
  skip_if_not_installed("brms")

  test_data <- create_test_data()
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "brms", refresh = 0)

  result <- gt_draws(g)

  expect_true(all(result$person > 0))
  expect_true(all(result$item > 0))
  expect_true(all(result$Residual > 0))
})

# Test: gt_draws.gstudy g and phi are between 0 and 1
test_that("gt_draws.gstudy g and phi are between 0 and 1", {
  skip_if_not_installed("brms")

  test_data <- create_test_data()
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "brms", refresh = 0)

  result <- gt_draws(g)

  expect_true(all(result$g >= 0 & result$g <= 1, na.rm = TRUE))
  expect_true(all(result$phi >= 0 & result$phi <= 1, na.rm = TRUE))
})

# Test: gt_draws.dstudy filtering coefficients
test_that("gt_draws.dstudy filters by coefficients", {
  skip_if_not_installed("brms")

  test_data <- create_test_data()
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "brms", refresh = 0)
  d <- dstudy(g, n = list(item = 5), estimation = "posterior")

  result <- gt_draws(d, what = "coefficients", coefficients = c("g", "phi"))

  expect_s3_class(result, "tbl_df")
  expect_true("g" %in% names(result))
  expect_true("phi" %in% names(result))
  expect_false("uni" %in% names(result))
})

# Test: gt_draws.gstudy multivariate returns named list
test_that("gt_draws.gstudy returns named list for multivariate model", {
  skip_if_not_installed("brms")
  skip_if_not_installed("tidyr")

  library(tidyr)
  data(rajaratnam)
  raj <- rajaratnam
  raj$Subtest <- paste0("test", raj$Subtest)
  rajl <- raj |> pivot_wider(
    names_from = Subtest,
    names_sep = ".",
    values_from = c(Score)
  )

  g <- gstudy(bf(mvbind(test1, test2, test3) ~ (1|pr|Person) + (1|ItemId)) + set_rescor(FALSE),
              data = rajl, estimator = "brms", cores = 4, refresh = 0, silent = 2)

  result <- gt_draws(g)

  expect_type(result, "list")
  expect_true("test1" %in% names(result))
  expect_true("test2" %in% names(result))
  expect_true("test3" %in% names(result))

  expect_s3_class(result$test1, "tbl_df")
  expect_true("draw" %in% names(result$test1))
  expect_true("Person" %in% names(result$test1))
  expect_true("ItemId" %in% names(result$test1))
  expect_true("Residual" %in% names(result$test1))
  expect_true("uni" %in% names(result$test1))
  expect_true("g" %in% names(result$test1))
  expect_true("phi" %in% names(result$test1))
})

# Test: gt_draws.gstudy multivariate dims filter returns named list
test_that("gt_draws.gstudy dims filter returns named list", {
  skip_if_not_installed("brms")
  skip_if_not_installed("tidyr")

  library(tidyr)
  data(rajaratnam)
  raj <- rajaratnam
  raj$Subtest <- paste0("test", raj$Subtest)
  rajl <- raj |> pivot_wider(
    names_from = Subtest,
    names_sep = ".",
    values_from = c(Score)
  )

  g <- gstudy(bf(mvbind(test1, test2, test3) ~ (1|pr|Person) + (1|ItemId)) + set_rescor(FALSE),
              data = rajl, estimator = "brms", cores = 4, refresh = 0, silent = 2)

  result <- gt_draws(g, dims = "test1")

  expect_type(result, "list")
  expect_true("test1" %in% names(result))
  expect_false("test2" %in% names(result))
})

# Test: gt_draws.dstudy multivariate returns named list
test_that("gt_draws.dstudy returns named list for multivariate model", {
  skip_if_not_installed("brms")
  skip_if_not_installed("tidyr")

  library(tidyr)
  data(rajaratnam)
  raj <- rajaratnam
  raj$Subtest <- paste0("test", raj$Subtest)
  rajl <- raj |> pivot_wider(
    names_from = Subtest,
    names_sep = ".",
    values_from = c(Score)
  )

  g <- gstudy(bf(mvbind(test1, test2, test3) ~ (1|pr|Person) + (1|ItemId)) + set_rescor(FALSE),
              data = rajl, estimator = "brms", cores = 4, refresh = 0, silent = 2)

  d <- dstudy(g, silent = TRUE)

  result <- gt_draws(d, what = "coefficients")

  expect_type(result, "list")
  expect_true("test1" %in% names(result))
  expect_true("test2" %in% names(result))
  expect_true("test3" %in% names(result))

  expect_s3_class(result$test1, "tbl_df")
  expect_true("draw" %in% names(result$test1))
  expect_true("uni" %in% names(result$test1))
  expect_true("g" %in% names(result$test1))
  expect_true("phi" %in% names(result$test1))

  expect_false("dim" %in% names(result$test1))
  expect_false("type" %in% names(result$test1))
})

# Test: gt_draws.dstudy multivariate what = "all" includes VAR columns
test_that("gt_draws.dstudy what = 'all' includes VAR columns merged", {
  skip_if_not_installed("brms")
  skip_if_not_installed("tidyr")

  library(tidyr)
  data(rajaratnam)
  raj <- rajaratnam
  raj$Subtest <- paste0("test", raj$Subtest)
  rajl <- raj |> pivot_wider(
    names_from = Subtest,
    names_sep = ".",
    values_from = c(Score)
  )

  g <- gstudy(bf(mvbind(test1, test2, test3) ~ (1|pr|Person) + (1|ItemId)) + set_rescor(FALSE),
              data = rajl, estimator = "brms", cores = 4, refresh = 0, silent = 2)

  d <- dstudy(g, silent = TRUE)

  result <- gt_draws(d, what = "all")

  expect_type(result, "list")
  expect_true("test1" %in% names(result))
  expect_true("test2" %in% names(result))
  expect_true("test3" %in% names(result))
  expect_true("composite" %in% names(result))

  expect_true("var_rel" %in% names(result$test1))
  expect_true("var_abs" %in% names(result$test1))
  expect_true("prmse_c_rel" %in% names(result$test1))
  expect_true("prmse_c_abs" %in% names(result$test1))

  expect_false("var_rel" %in% names(result$composite))
})

# Test: gt_draws.dstudy multivariate composite structure
test_that("gt_draws.dstudy composite has correct structure", {
  skip_if_not_installed("brms")
  skip_if_not_installed("tidyr")

  library(tidyr)
  data(rajaratnam)
  raj <- rajaratnam
  raj$Subtest <- paste0("test", raj$Subtest)
  rajl <- raj |> pivot_wider(
    names_from = Subtest,
    names_sep = ".",
    values_from = c(Score)
  )

  g <- gstudy(bf(mvbind(test1, test2, test3) ~ (1|pr|Person) + (1|ItemId)) + set_rescor(FALSE),
              data = rajl, estimator = "brms", cores = 4, refresh = 0, silent = 2)

  d <- dstudy(g, silent = TRUE)

  result <- gt_draws(d, what = "composite")

  expect_s3_class(result, "tbl_df")
  expect_true("draw" %in% names(result))
  expect_true("uni" %in% names(result))
  expect_true("g" %in% names(result))
  expect_true("phi" %in% names(result))

  expect_false("var_rel" %in% names(result))
})

# Test: gt_draws.dstudy multivariate what = "variance"
test_that("gt_draws.dstudy what = 'variance' returns gstudy draws", {
  skip_if_not_installed("brms")
  skip_if_not_installed("tidyr")

  library(tidyr)
  data(rajaratnam)
  raj <- rajaratnam
  raj$Subtest <- paste0("test", raj$Subtest)
  rajl <- raj |> pivot_wider(
    names_from = Subtest,
    names_sep = ".",
    values_from = c(Score)
  )

  g <- gstudy(bf(mvbind(test1, test2, test3) ~ (1|pr|Person) + (1|ItemId)) + set_rescor(FALSE),
              data = rajl, estimator = "brms", cores = 4, refresh = 0, silent = 2)

  d <- dstudy(g, silent = TRUE)

  result <- gt_draws(d, what = "variance")

  expect_type(result, "list")
  expect_true("test1" %in% names(result))
  expect_true("test2" %in% names(result))
  expect_true("test3" %in% names(result))

  expect_true("Person" %in% names(result$test1))
  expect_true("ItemId" %in% names(result$test1))
  expect_true("Residual" %in% names(result$test1))
})

# Test: gt_draws.dstudy multivariate no duplicate names
test_that("gt_draws.dstudy what = 'all' has no duplicate names", {
  skip_if_not_installed("brms")
  skip_if_not_installed("tidyr")

  library(tidyr)
  data(rajaratnam)
  raj <- rajaratnam
  raj$Subtest <- paste0("test", raj$Subtest)
  rajl <- raj |> pivot_wider(
    names_from = Subtest,
    names_sep = ".",
    values_from = c(Score)
  )

  g <- gstudy(bf(mvbind(test1, test2, test3) ~ (1|pr|Person) + (1|ItemId)) + set_rescor(FALSE),
              data = rajl, estimator = "brms", cores = 4, refresh = 0, silent = 2)

  d <- dstudy(g, silent = TRUE)

  result <- gt_draws(d, what = "all")

  nms <- names(result)
  expect_equal(length(nms), length(unique(nms)))
})
