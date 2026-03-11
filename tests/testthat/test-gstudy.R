# Tests for gstudy function

# Test data for reuse
test_data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, 5)),
  rater = factor(rep(1:5, each = 20)),
  item = factor(rep(1:10, each = 10))
)

# =============================================================================
# Basic gstudy tests
# =============================================================================

test_that("gstudy returns a gstudy object", {
  skip_if_not_installed("lme4")
  result <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  expect_s3_class(result, "gstudy")
})

test_that("gstudy validates backend argument", {
  expect_error(
    gstudy(score ~ (1 | person), data = test_data, backend = "invalid"),
    "should be one of"
  )
})

test_that("gstudy works with lme4 backend", {
  skip_if_not_installed("lme4")
  result <- gstudy(
    score ~ (1 | person) + (1 | rater), 
    data = test_data, 
    backend = "lme4"
  )
  expect_s3_class(result, "gstudy")
  expect_equal(result$backend, "lme4")
})

test_that("gstudy works with auto backend selection", {
  skip_if_not_installed("lme4")
  result <- gstudy(
    score ~ (1 | person) + (1 | rater), 
    data = test_data, 
    backend = "auto"
  )
  expect_s3_class(result, "gstudy")
  expect_true(result$backend %in% c("lme4", "brms"))
})

# =============================================================================
# Input validation tests
# =============================================================================

test_that("gstudy requires data frame", {
  expect_error(
    gstudy(score ~ (1 | person), data = "not a data frame"),
    "must be a data frame"
  )
})

test_that("gstudy validates formula has random effects", {
  expect_error(
    gstudy(score ~ item, data = test_data),
    "random effect"
  )
})

# =============================================================================
# Output structure tests
# =============================================================================

test_that("gstudy stores variance components", {
  skip_if_not_installed("lme4")
  result <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  expect_true(!is.null(result$variance_components))
  expect_s3_class(result$variance_components, "tbl_df")
})

test_that("gstudy stores facets", {
  skip_if_not_installed("lme4")
  result <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  expect_true(!is.null(result$facets))
  expect_true("person" %in% result$facets)
  expect_true("rater" %in% result$facets)
})

test_that("gstudy stores formula", {
  skip_if_not_installed("lme4")
  f <- score ~ (1 | person) + (1 | rater)
  result <- gstudy(f, data = test_data)
  expect_equal(result$formula, f)
})

test_that("gstudy stores number of observations", {
  skip_if_not_installed("lme4")
  result <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  expect_equal(result$n_obs, nrow(test_data))
})

test_that("gstudy stores backend used", {
  skip_if_not_installed("lme4")
  result <- gstudy(score ~ (1 | person), data = test_data, backend = "lme4")
  expect_equal(result$backend, "lme4")
})

test_that("gstudy stores is_multivariate flag", {
  skip_if_not_installed("lme4")
  result <- gstudy(score ~ (1 | person), data = test_data)
  expect_false(result$is_multivariate)
})

# =============================================================================
# Object of measurement tests
# =============================================================================

test_that("gstudy auto-detects object of measurement as first facet", {
skip_if_not_installed("lme4")
result <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
expect_equal(result$object, "person")
})

test_that("gstudy correctly identifies object of measurement with double-bar syntax", {
skip_if_not_installed("brms")
result <- gstudy(score ~ (1 | cor_p | person) + (1 | rater), data = test_data, backend = "brms")
expect_equal(result$object, "person")
})

test_that("gstudy extracts correct facets with double-bar syntax", {
skip_if_not_installed("brms")
result <- gstudy(score ~ (1 | cor_p | person) + (1 | rater), data = test_data, backend = "brms")
expect_true("person" %in% result$facets)
expect_true("rater" %in% result$facets)
expect_false("cor_p" %in% result$facets)
})

# =============================================================================
# Custom facets tests
# =============================================================================

test_that("gstudy accepts custom facets", {
  skip_if_not_installed("lme4")
  result <- gstudy(
    score ~ (1 | person) + (1 | rater), 
    data = test_data, 
    facets = c("person", "rater")
  )
  expect_equal(result$facets, c("person", "rater"))
})

# =============================================================================
# is.gstudy test
# =============================================================================

test_that("is.gstudy returns TRUE for gstudy objects", {
  skip_if_not_installed("lme4")
  result <- gstudy(score ~ (1 | person), data = test_data)
  expect_true(is.gstudy(result))
})

test_that("is.gstudy returns FALSE for non-gstudy objects", {
  expect_false(is.gstudy(list()))
  expect_false(is.gstudy(NULL))
  expect_false(is.gstudy("not a gstudy"))
})

# =============================================================================
# VarCorr.gstudy tests
# =============================================================================

test_that("VarCorr.gstudy works with lme4 backend", {
  skip_if_not_installed("lme4")
  result <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, backend = "lme4")
  vc <- VarCorr.gstudy(result)
  expect_type(vc, "list")
  expect_true("person" %in% vc$Name[vc$Name != ""])
})

test_that("VarCorr.gstudy works with univariate brms backend", {
  skip_if_not_installed("brms")
  result <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, backend = "brms", chains = 1, iter = 500)
  vc <- VarCorr.gstudy(result)
  expect_type(vc, "list")
})

test_that("VarCorr.gstudy works with multivariate brms backend", {
  skip_if_not_installed("brms")
  test_data_mv <- data.frame(
    score1 = rnorm(100),
    score2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:10, each = 10))
  )
  result <- gstudy(mvbind(score1, score2) ~ (1 | person) + (1 | item), 
                   data = test_data_mv, backend = "brms", chains = 1, iter = 500)
  vc <- VarCorr.gstudy(result)
  expect_type(vc, "list")
  expect_true(length(vc) == 2)
  expect_true("score1" %in% names(vc))
  expect_true("score2" %in% names(vc))
})

test_that("VarCorr.gstudy fails with non-gstudy object", {
  expect_error(VarCorr.gstudy("not a gstudy"), "must be a gstudy object")
})

test_that("VarCorr.gstudy works with univariate mom backend", {
  result <- gstudy(score ~ (1 | person) + (1 | item), 
                   data = test_data, backend = "mom")
  vc <- VarCorr.gstudy(result)
  expect_s3_class(vc, "data.frame")
  expect_true("person" %in% vc$Group)
  expect_true("item" %in% vc$Group)
  expect_true("Residual" %in% vc$Group)
  expect_true("Std.Dev." %in% names(vc))
  expect_true("Var." %in% names(vc))
})

test_that("VarCorr.gstudy works with multivariate mom backend", {
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
  # Each dimension should be a data.frame
  expect_s3_class(vc[["score1"]], "data.frame")
  expect_true("person" %in% vc[["score1"]]$Group)
})

# =============================================================================
# ranef.gstudy tests
# =============================================================================

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
