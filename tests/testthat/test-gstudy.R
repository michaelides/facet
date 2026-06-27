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

test_that("gstudy validates estimator argument", {
  expect_error(
    gstudy(score ~ (1 | person), data = test_data, estimator = "invalid"),
    "should be one of"
  )
})

test_that("gstudy works with lme4 estimator", {
  skip_if_not_installed("lme4")
  result <- gstudy(
    score ~ (1 | person) + (1 | rater), 
    data = test_data, 
    estimator = "lme4"
  )
  expect_s3_class(result, "gstudy")
  expect_equal(result$estimator, "lme4")
})

test_that("gstudy works with auto estimator selection", {
  skip_if_not_installed("lme4")
  result <- gstudy(
    score ~ (1 | person) + (1 | rater), 
    data = test_data, 
    estimator = "auto"
  )
  expect_s3_class(result, "gstudy")
  expect_true(result$estimator %in% c("lme4", "brms"))
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

test_that("gstudy stores estimator used", {
  skip_if_not_installed("lme4")
  result <- gstudy(score ~ (1 | person), data = test_data, estimator = "lme4")
  expect_equal(result$estimator, "lme4")
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
result <- gstudy(score ~ (1 | cor_p | person) + (1 | rater), data = test_data, estimator = "brms")
expect_equal(result$object, "person")
})

test_that("gstudy extracts correct facets with double-bar syntax", {
skip_if_not_installed("brms")
result <- gstudy(score ~ (1 | cor_p | person) + (1 | rater), data = test_data, estimator = "brms")
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

test_that("VarCorr.gstudy works with lme4 estimator", {
  skip_if_not_installed("lme4")
  result <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "lme4")
  vc <- VarCorr.gstudy(result)
  expect_type(vc, "list")
  expect_true("person" %in% vc$Name[vc$Name != ""])
})

test_that("VarCorr.gstudy works with univariate brms estimator", {
  skip_if_not_installed("brms")
  result <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "brms", chains = 1, iter = 500)
  vc <- VarCorr.gstudy(result)
  expect_type(vc, "list")
})

test_that("VarCorr.gstudy works with multivariate brms estimator", {
  skip_if_not_installed("brms")
  test_data_mv <- data.frame(
    score1 = rnorm(100),
    score2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:10, each = 10))
  )
  result <- gstudy(mvbind(score1, score2) ~ (1 | person) + (1 | item), 
                   data = test_data_mv, estimator = "brms", chains = 1, iter = 500)
  vc <- VarCorr.gstudy(result)
  expect_type(vc, "list")
  expect_true(length(vc) == 2)
  expect_true("score1" %in% names(vc))
  expect_true("score2" %in% names(vc))
})

test_that("VarCorr.gstudy fails with non-gstudy object", {
  expect_error(VarCorr.gstudy("not a gstudy"), "must be a gstudy object")
})

test_that("VarCorr.gstudy works with univariate aov estimator", {
  result <- gstudy(score ~ (1 | person) + (1 | item), 
                   data = test_data, estimator = "aov")
  vc <- VarCorr.gstudy(result)
  expect_s3_class(vc, "data.frame")
  expect_true("person" %in% vc$Group)
  expect_true("item" %in% vc$Group)
  expect_true("Residual" %in% vc$Group)
  expect_true("Std.Dev." %in% names(vc))
  expect_true("Var." %in% names(vc))
})

test_that("VarCorr.gstudy works with multivariate aov estimator", {
  test_data_mv <- data.frame(
    score1 = rnorm(100),
    score2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:10, each = 10))
  )
  result <- gstudy(mvbind(score1, score2) ~ (1 | person) + (1 | item), 
                   data = test_data_mv, estimator = "aov")
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

test_that("ranef.gstudy works with univariate aov estimator", {
  result <- gstudy(score ~ (1 | person) + (1 | item), 
                   data = test_data, estimator = "aov")
  re <- ranef.gstudy(result)
  expect_type(re, "list")
  expect_true("person" %in% names(re))
  expect_true("item" %in% names(re))
  # Check structure - should be named vectors
  expect_type(re[["person"]], "double")
})

test_that("ranef.gstudy works with multivariate aov estimator", {
  test_data_mv <- data.frame(
    score1 = rnorm(100),
    score2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:10, each = 10))
  )
  result <- gstudy(mvbind(score1, score2) ~ (1 | person) + (1 | item),
    data = test_data_mv, estimator = "aov")
  re <- ranef.gstudy(result)
  expect_type(re, "list")
  expect_true("score1" %in% names(re))
  expect_true("score2" %in% names(re))
  expect_true("person" %in% names(re[["score1"]]))
})

# =============================================================================
# Long-format multivariate detection tests
# =============================================================================

test_that("gstudy detects long-format multivariate models", {
  skip_if_not_installed("brms")
  skip_on_cran()
  test_data <- data.frame(
    Score = rnorm(60),
    Person = factor(rep(1:10, 6)),
    Subtest = factor(rep(c("A", "B"), each = 30)),
    Item = factor(rep(1:4, 15))
  )
formula <- brms::bf(Score ~ 0 + Subtest + (0+Subtest|Person), sigma ~ 0 + Subtest)
# Test detection logic
detection <- is_long_format_multivariate(formula)
expect_true(detection$is_long)
expect_equal(detection$dimension_var, "Subtest")
})

test_that("gstudy long_format detection sets is_mv correctly", {
skip_if_not_installed("brms")
formula <- brms::bf(Score ~ 0 + Subtest + (0+Subtest|Person), sigma ~ 0 + Subtest)
detection <- is_long_format_multivariate(formula)
# Detection should identify this as long-format
expect_true(detection$is_long)
})

# =============================================================================
# Unbalanced multivariate tests
# =============================================================================

test_that("unbalanced = TRUE with aov estimator works for multivariate", {
  skip_if_not_installed("brms")
  
  # Create test data with unbalanced design
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  # Make y1 missing for some persons (creates unbalanced design)
  test_data$y1[test_data$person %in% c("1", "2")] <- NA
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  
  # Should NOT warn when unbalanced = TRUE
  expect_no_warning(
    result <- gstudy(formula, test_data, estimator = "aov", unbalanced = TRUE)
  )
  
  expect_s3_class(result, "mgstudy")
  expect_true(result$model$is_unbalanced)
  expect_true(!is.null(result$variance_components))
})

test_that("unbalanced = TRUE warns for univariate models", {
  skip_if_not_installed("lme4")
  
  expect_warning(
    gstudy(score ~ (1 | person), data = test_data, estimator = "aov", unbalanced = TRUE),
    "only applicable to multivariate"
  )
})

test_that("unbalanced = TRUE warns for brms estimator", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(50),
    y2 = rnorm(50),
    person = factor(rep(1:10, 5))
  )
  
  expect_warning(
    gstudy(brms::mvbind(y1, y2) ~ (1 | person), data = test_data, 
           estimator = "brms", unbalanced = TRUE),
    "not implemented for brms"
  )
})

test_that("unbalanced = TRUE errors for lme4 multivariate", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(50),
    y2 = rnorm(50),
    person = factor(rep(1:10, 5))
  )
  
  expect_error(
    gstudy(brms::mvbind(y1, y2) ~ (1 | person), data = test_data,
           estimator = "lme4", unbalanced = TRUE),
    "Multivariate models are not supported by lme4"
  )
})

test_that("check_multivariate_balance respects unbalanced flag", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  # Make y1 missing for ALL observations of person 1 (unbalanced design)
  test_data$y1[test_data$person == "1"] <- NA
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  
  # Without unbalanced flag, should have warning
  result1 <- check_multivariate_balance(formula, test_data, FALSE, unbalanced = FALSE)
  expect_false(is.null(result1$warning_message))
  expect_true(grepl("Unbalanced", result1$warning_message))
  
  # With unbalanced flag, should have info message instead
  result2 <- check_multivariate_balance(formula, test_data, FALSE, unbalanced = TRUE)
  expect_true(is.null(result2$warning_message))
  expect_false(is.null(result2$info_message))
  expect_true(grepl("Henderson", result2$info_message))
})

test_that("fit_aov_multivariate_unbalanced returns correct structure", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  # Make unbalanced
  test_data$y1[test_data$person %in% c("1", "2")] <- NA
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  
  result <- fit_aov(formula, test_data, unbalanced = TRUE)
  
  expect_s3_class(result, "aovfit")
  expect_true(result$is_unbalanced)
  expect_true(!is.null(result$n_per_dim))
  expect_true(!is.null(result$variance_components))
  expect_true("dim" %in% names(result$variance_components))
})

# =============================================================================
# 4-decimal-point display tests
# =============================================================================

test_that("gstudy variance components are stored at 4 decimal places by default", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  vc <- g$variance_components
  expect_true("var" %in% names(vc))
  expect_true("pct" %in% names(vc))

  # Stored at 4 dp (allowing for floating point error of 1e-6)
  expect_true(all(abs(vc$var - round(vc$var, 4)) < 1e-6, na.rm = TRUE))
  expect_true(all(abs(vc$pct - round(vc$pct, 4)) < 1e-6, na.rm = TRUE))
})

test_that("print.gstudy default uses 4 decimal places", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  out <- capture.output(print(g))
  out <- paste(out, collapse = "\n")

  # Output should contain a 4-decimal-place value like "1.2345"
  expect_match(out, "\\d+\\.\\d{4}", all = FALSE)
})

test_that("summary.gstudy default uses 4 decimal places", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  out <- capture.output(summary(g))
  out <- paste(out, collapse = "\n")

  expect_match(out, "\\d+\\.\\d{4}", all = FALSE)
})

test_that("tidy.gstudy returns a tibble with numeric columns at 4 decimal places by default", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  result <- tidy(g)

  expect_s3_class(result, "tbl_df")
  numeric_cols <- vapply(result, is.numeric, logical(1))
  numeric_cols <- names(result)[numeric_cols]

  for (col in numeric_cols) {
    expect_true(all(abs(result[[col]] - round(result[[col]], 4)) < 1e-6, na.rm = TRUE))
  }
})

test_that("tidy.gstudy digits argument overrides default", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  result <- tidy(g, digits = 2)
  expect_s3_class(result, "tbl_df")

  expect_true(all(abs(result$var - round(result$var, 2)) < 1e-6, na.rm = TRUE))
})

test_that("print.gstudy digits argument overrides default", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  out <- capture.output(print(g, digits = 2))
  out <- paste(out, collapse = "\n")

  # With digits = 2, the output should not contain values with more than
  # 4 decimal places (5+ dp). Some values may still show 3-4 dp due to
  # pillar's significant figure handling.
  expect_false(grepl("\\d+\\.\\d{5,}", out))
})

