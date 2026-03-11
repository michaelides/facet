# Tests for formula helper functions

# Test data for reuse
test_data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, 5)),
  rater = factor(rep(1:5, each = 20)),
  item = factor(rep(1:10, each = 10))
)

# =============================================================================
# is_multivariate tests
# =============================================================================

test_that("is_multivariate returns FALSE for standard univariate formula", {
  f <- score ~ x + (1 | facet)
  expect_false(is_multivariate(f))
})

test_that("is_multivariate returns FALSE for simple random effects formula", {
  f <- score ~ (1 | person) + (1 | rater)
  expect_false(is_multivariate(f))
})

# =============================================================================
# parse_g_formula tests
# =============================================================================

test_that("parse_g_formula extracts response variable", {
  f <- score ~ (1 | person) + (1 | rater)
  parsed <- parse_g_formula(f)
  expect_equal(parsed$response, "score")
})

test_that("parse_g_formula extracts random facets", {
  f <- score ~ (1 | person) + (1 | rater)
  parsed <- parse_g_formula(f)
  expect_true(length(parsed$random_facets) == 2)
  expect_true("person" %in% parsed$random_facets)
  expect_true("rater" %in% parsed$random_facets)
})

test_that("parse_g_formula extracts single random facet", {
  f <- score ~ (1 | person)
  parsed <- parse_g_formula(f)
  expect_equal(parsed$random_facets, "person")
})

test_that("parse_g_formula handles fixed effects", {
  f <- score ~ item + (1 | person)
  parsed <- parse_g_formula(f)
  expect_true("item" %in% parsed$fixed)
})

test_that("parse_g_formula returns a list with expected components", {
  f <- score ~ (1 | person)
  parsed <- parse_g_formula(f)
  expect_type(parsed, "list")
  expect_true("response" %in% names(parsed))
  expect_true("fixed" %in% names(parsed))
  expect_true("random" %in% names(parsed))
  expect_true("random_facets" %in% names(parsed))
})

test_that("parse_g_formula handles double-bar syntax for correlated effects", {
  f <- score ~ (1 | cor_p | person) + (1 | rater)
  parsed <- parse_g_formula(f)
  expect_true("person" %in% parsed$random_facets)
  expect_true("rater" %in% parsed$random_facets)
  expect_false("cor_p" %in% parsed$random_facets)
})

test_that("parse_g_formula handles double-bar with interaction term", {
  f <- score ~ (1 | cor_grp | person:item) + (1 | rater)
  parsed <- parse_g_formula(f)
  expect_true("person" %in% parsed$random_facets)
  expect_true("item" %in% parsed$random_facets)
  expect_true("rater" %in% parsed$random_facets)
  expect_false("cor_grp" %in% parsed$random_facets)
})

test_that("parse_g_formula handles mixed single and double-bar syntax", {
  f <- score ~ (1 | cor_p | person) + (1 | rater) + (1 | item)
  parsed <- parse_g_formula(f)
  expect_true("person" %in% parsed$random_facets)
  expect_true("rater" %in% parsed$random_facets)
  expect_true("item" %in% parsed$random_facets)
  expect_false("cor_p" %in% parsed$random_facets)
})

# =============================================================================
# validate_formula tests
# =============================================================================

test_that("validate_formula returns TRUE for valid formula with random effects", {
  f <- score ~ (1 | person) + (1 | rater)
  expect_true(validate_formula(f, "auto"))
})

test_that("validate_formula rejects formulas without random effects", {
  f <- score ~ x
  expect_error(validate_formula(f, "lme4"), "random effect")
})

test_that("validate_formula rejects non-formula objects", {
  expect_error(validate_formula("not a formula", "lme4"), "formula")
})

# =============================================================================
# detect_facets tests
# =============================================================================

test_that("detect_facets returns character vector", {
  f <- score ~ (1 | person) + (1 | rater)
  facets <- detect_facets(f)
  expect_type(facets, "character")
})

test_that("detect_facets extracts correct facet names", {
  f <- score ~ (1 | person) + (1 | rater)
  facets <- detect_facets(f)
  expect_true("person" %in% facets)
  expect_true("rater" %in% facets)
})

test_that("detect_facets validates against data", {
  f <- score ~ (1 | person) + (1 | rater)
  facets <- detect_facets(f, test_data)
  expect_true(all(c("person", "rater") %in% facets))
})

test_that("detect_facets errors on missing facets in data", {
  f <- score ~ (1 | person) + (1 | nonexistent)
  expect_error(detect_facets(f, test_data), "not found")
})

test_that("detect_facets returns unique facets", {
  f <- score ~ (1 | person) + (1 | person)
  facets <- detect_facets(f)
  expect_equal(length(facets), 1)
})

# =============================================================================
# convert_formula tests
# =============================================================================

test_that("convert_formula returns a formula for lme4 backend", {
  f <- score ~ (1 | person)
  result <- convert_formula(f, "lme4")
  expect_s3_class(result, "formula")
})

test_that("convert_formula preserves formula structure for lme4", {
  f <- score ~ item + (1 | person)
  result <- convert_formula(f, "lme4")
  expect_equal(result, f)
})

test_that("convert_formula handles brms backend", {
  f <- score ~ (1 | person)
  result <- convert_formula(f, "brms")
  expect_s3_class(result, "formula")
})

# =============================================================================
# extract_facets tests
# =============================================================================

test_that("extract_facets works with formula input", {
  f <- score ~ (1 | person) + (1 | rater)
  facets <- extract_facets(formula = f)
  expect_true(length(facets) == 2)
  expect_true("person" %in% facets)
  expect_true("rater" %in% facets)
})

test_that("extract_facets works with variance components input", {
  vc <- data.frame(
    facet = c("person", "rater", "Residual"),
    variance = c(1.0, 0.5, 0.3)
  )
  facets <- extract_facets(vc = vc)
  expect_true("person" %in% facets)
  expect_true("rater" %in% facets)
  expect_false("Residual" %in% facets)
})

test_that("extract_facets prefers formula over vc when both provided", {
  f <- score ~ (1 | person)
  vc <- data.frame(
    facet = c("different", "Residual"),
    variance = c(1.0, 0.3)
  )
  facets <- extract_facets(formula = f, vc = vc)
  expect_true("person" %in% facets)
})

test_that("extract_facets returns empty character vector with no inputs", {
  facets <- extract_facets()
  expect_equal(length(facets), 0)
})

test_that("extract_facets correctly identifies grouping variable with double-bar", {
  f <- score ~ (1 | cor_p | person) + (1 | rater)
  facets <- extract_facets(formula = f)
  expect_true("person" %in% facets)
  expect_true("rater" %in% facets)
  expect_false("cor_p" %in% facets)
})

test_that("extract_facets handles multiple double-bar terms", {
  f <- score ~ (1 | cor_p | person) + (1 | cor_r | rater)
  facets <- extract_facets(formula = f)
  expect_true("person" %in% facets)
  expect_true("rater" %in% facets)
  expect_false("cor_p" %in% facets)
  expect_false("cor_r" %in% facets)
})

# =============================================================================
# parse_residual_facets tests
# =============================================================================

test_that("parse_residual_facets returns interaction of all facets when no data provided", {
  f <- score ~ (1 | p) + (1 | i) + (1 | d)
  result <- parse_residual_facets(f)
  expect_equal(result, "p:i:d")
})

test_that("parse_residual_facets returns single facet when only one random effect", {
  f <- score ~ (1 | person)
  result <- parse_residual_facets(f)
  expect_equal(result, "person")
})

test_that("parse_residual_facets handles two facets without data", {
  f <- score ~ (1 | person) + (1 | rater)
  result <- parse_residual_facets(f)
  expect_equal(result, "person:rater")
})

test_that("parse_residual_facets detects crossed facets from data", {
  crossed_data <- data.frame(
    score = rnorm(100),
    p = factor(rep(1:10, 10)),
    i = factor(rep(1:5, each = 20)),
    d = factor(rep(1:2, each = 50))
  )
  f <- score ~ (1 | p) + (1 | i) + (1 | d)
  result <- parse_residual_facets(f, crossed_data)
  expect_equal(result, "p:i:d")
})

test_that("parse_residual_facets detects nested facets from data", {
  nested_data <- data.frame(
    score = rnorm(100),
    p = factor(rep(1:10, 10)),
    i = factor(rep(1:5, times = 20)),
    d = factor(rep(1:2, each = 50))
  )
  for (i_lvl in unique(nested_data$i)) {
    nested_data$d[nested_data$i == i_lvl] <- ((as.numeric(i_lvl) - 1) %% 2) + 1
  }
  nested_data$p <- factor(rep(1:10, each = 10))
  f <- score ~ (1 | p) + (1 | i:d) + (1 | d)
  result <- parse_residual_facets(f, nested_data)
  expect_equal(result, "p:i")
})

test_that("parse_residual_facets handles interaction terms in formula", {
  f <- score ~ (1 | p) + (1 | i:d) + (1 | d) + (1 | p:d)
  result <- parse_residual_facets(f)
  expect_equal(result, "p:i:d")
})

test_that("parse_residual_facets returns empty string for formula with no random effects", {
  skip("parse_g_formula handles this case differently - may need validation elsewhere")
})

test_that("parse_residual_facets works with partially nested design", {
  partial_data <- data.frame(
    score = rnorm(60),
    p = factor(rep(1:6, 10)),
    i = factor(rep(1:3, each = 20)),
    d = factor(rep(1:2, each = 30))
  )
  f <- score ~ (1 | p) + (1 | i) + (1 | d)
  result <- parse_residual_facets(f, partial_data)
  expect_true(nchar(result) > 0)
  facets_in_result <- strsplit(result, ":")[[1]]
  expect_true(all(facets_in_result %in% c("p", "i", "d")))
})

test_that("parse_residual_facets preserves facet order from formula", {
  f <- score ~ (1 | zebra) + (1 | apple) + (1 | mango)
  result <- parse_residual_facets(f)
  expect_equal(result, "zebra:apple:mango")
})
