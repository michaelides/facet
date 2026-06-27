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

test_that("convert_formula returns a formula for lme4 estimator", {
  f <- score ~ (1 | person)
  result <- convert_formula(f, "lme4")
  expect_s3_class(result, "formula")
})

test_that("convert_formula preserves formula structure for lme4", {
  f <- score ~ item + (1 | person)
  result <- convert_formula(f, "lme4")
  expect_equal(result, f)
})

test_that("convert_formula handles brms estimator", {
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

# =============================================================================
# is_long_format_multivariate tests
# =============================================================================

test_that("is_long_format_multivariate returns FALSE for non-brms formulas", {
  result <- is_long_format_multivariate(Score ~ (1|Person))
  expect_false(result$is_long)
})

test_that("is_long_format_multivariate returns FALSE for wide-format multivariate", {
  skip_if_not_installed("brms")
  formula <- brms::mvbind(score1, score2) ~ (1|Person)
  result <- is_long_format_multivariate(formula)
  expect_false(result$is_long)
})

test_that("is_long_format_multivariate detects sigma auxiliary formula", {
  skip_if_not_installed("brms")
  formula <- brms::bf(Score ~ 0 + Subtest + (0+Subtest|Person), sigma ~ 0 + Subtest)
  result <- is_long_format_multivariate(formula)
  expect_true(result$is_long)
  expect_equal(result$dimension_var, "Subtest")
})

test_that("is_long_format_multivariate requires dimension var in random effects", {
  skip_if_not_installed("brms")
  formula <- brms::bf(Score ~ (1|Person), sigma ~ 1)
  result <- is_long_format_multivariate(formula)
  expect_false(result$is_long)
})

# =============================================================================
# validate_interaction_levels tests
# =============================================================================

test_that("validate_interaction_levels passes for valid interaction terms", {
  data <- data.frame(
    score = rnorm(120),
    person = factor(rep(1:10, 12)),
    item = factor(rep(1:6, each = 20))
  )
  # person:item has 10*6 = 60 levels, which is < 120 observations
  f <- score ~ (1 | person) + (1 | item) + (1 | person:item)
  expect_silent(validate_interaction_levels(f, data, "lme4"))
})

test_that("validate_interaction_levels errors when interaction has same levels as observations", {
  data <- data.frame(
    score = rnorm(20),
    person = factor(rep(1:10, 2)),
    item = factor(rep(1:2, each = 10))
  )
  f <- score ~ (1 | person) + (1 | item) + (1 | person:item)
  expect_error(
    validate_interaction_levels(f, data, "lme4"),
    "as many or more levels than observations"
  )
})

test_that("validate_interaction_levels does not error for non-lme4 estimators", {
  data <- data.frame(
    score = rnorm(20),
    person = factor(rep(1:10, 2)),
    item = factor(rep(1:2, each = 10))
  )
  f <- score ~ (1 | person) + (1 | item) + (1 | person:item)
  expect_silent(validate_interaction_levels(f, data, "brms"))
  expect_silent(validate_interaction_levels(f, data, "mom"))
})

test_that("validate_interaction_levels provides helpful error message", {
  data <- data.frame(
    score = rnorm(20),
    person = factor(rep(1:10, 2)),
    item = factor(rep(1:2, each = 10))
  )
  f <- score ~ (1 | person) + (1 | item) + (1 | person:item)
  expect_error(
    validate_interaction_levels(f, data, "lme4"),
    "Solutions"
  )
})

test_that("validate_interaction_levels handles multiple problematic terms", {
  data <- data.frame(
    score = rnorm(12),
    a = factor(rep(1:3, 4)),
    b = factor(rep(1:4, each = 3)),
    c = factor(rep(1:2, each = 6))
  )
  # Both a:b and a:c would have 12 levels
  f <- score ~ (1 | a) + (1 | b) + (1 | c) + (1 | a:b) + (1 | a:c)
  expect_error(
    validate_interaction_levels(f, data, "lme4"),
    "Interaction term"
  )
})

test_that("validate_interaction_levels passes for nested designs", {
  data <- data.frame(
    score = rnorm(24),
    person = factor(rep(1:6, 4)),
    task = factor(rep(1:2, each = 12)),
    rater = factor(paste0(rep(1:2, each = 12), ".", rep(1:2, each = 12)))
  )
  # Rater labels are unique within each task (e.g. "1.1", "1.2", "2.1", "2.2"),
  # so rater:task has 4 levels, not 24. The formula uses the nested form here
  # for the unit test of validate_interaction_levels itself; in canonical G-study
  # formulas the un-nested form (1 | rater) is preferred when the rater labels
  # are already unique within the nesting facet.
  f <- score ~ (1 | person) + (1 | task) + (1 | rater:task)
  expect_silent(validate_interaction_levels(f, data, "lme4"))
})

# =============================================================================
# check_multivariate_balance tests
# =============================================================================

test_that("check_multivariate_balance returns balanced for no missing values", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  result <- check_multivariate_balance(formula, test_data, FALSE)
  
  expect_true(result$is_balanced)
  expect_false(result$has_missing)
  expect_false(result$has_different_levels)
  expect_null(result$warning_message)
  expect_null(result$info_message)
})

test_that("check_multivariate_balance detects genuine missing values with consistent levels", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  # Add missing values that don't affect facet levels
  test_data$y1[c(1, 50)] <- NA
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  result <- check_multivariate_balance(formula, test_data, FALSE)
  
  expect_true(result$is_balanced)
  expect_true(result$has_missing)
  expect_false(result$has_different_levels)
  expect_true(result$is_genuine_missing)
  expect_null(result$warning_message)
  expect_false(is.null(result$info_message))
  expect_true(grepl("Missing values detected", result$info_message))
})

test_that("check_multivariate_balance detects unbalanced design with different levels", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  # Make y1 missing for ALL observations of person 1 (5 rows)
  test_data$y1[test_data$person == "1"] <- NA
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  result <- check_multivariate_balance(formula, test_data, FALSE)
  
  expect_false(result$is_balanced)
  expect_true(result$has_missing)
  expect_true(result$has_different_levels)
  expect_false(result$is_genuine_missing)
  expect_false(is.null(result$warning_message))
  expect_true(grepl("Unbalanced multivariate design", result$warning_message))
})

test_that("check_multivariate_balance returns balanced for long-format", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  result <- check_multivariate_balance(formula, test_data, is_long_format = TRUE)
  
  expect_true(result$is_balanced)
  expect_false(result$has_missing)
  expect_null(result$warning_message)
  expect_null(result$info_message)
})

test_that("check_multivariate_balance returns balanced for univariate formula", {
  test_data <- data.frame(
    y = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  formula <- y ~ (1 | person) + (1 | rater)
  result <- check_multivariate_balance(formula, test_data, FALSE)
  
  expect_true(result$is_balanced)
  expect_false(result$has_missing)
})

test_that("check_multivariate_balance counts observations correctly", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  # Add 10 missing values
  test_data$y1[1:10] <- NA
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  result <- check_multivariate_balance(formula, test_data, FALSE)
  
  expect_equal(result$n_original, 100)
  expect_equal(result$n_removed, 10)
  expect_equal(result$n_complete, 90)
})

test_that("check_multivariate_balance extracts dimensions correctly", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    score1 = rnorm(50),
    score2 = rnorm(50),
    score3 = rnorm(50),
    person = factor(rep(1:10, 5))
  )
  
  formula <- brms::mvbind(score1, score2, score3) ~ (1 | person)
  result <- check_multivariate_balance(formula, test_data, FALSE)
  
  expect_equal(sort(result$dimensions), c("score1", "score2", "score3"))
})

test_that("check_multivariate_balance facet_levels_per_dim has correct structure", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  result <- check_multivariate_balance(formula, test_data, FALSE)
  
  expect_true(is.list(result$facet_levels_per_dim))
  expect_true("y1" %in% names(result$facet_levels_per_dim))
  expect_true("y2" %in% names(result$facet_levels_per_dim))
  expect_true("person" %in% names(result$facet_levels_per_dim$y1))
  expect_true("rater" %in% names(result$facet_levels_per_dim$y1))
})

test_that("check_multivariate_balance handles multiple facets with different levels", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  # Make y1 missing for all of person 1 and rater 1
  test_data$y1[test_data$person == "1" | test_data$rater == "1"] <- NA
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  result <- check_multivariate_balance(formula, test_data, FALSE)
  
  expect_false(result$is_balanced)
  expect_true(result$has_different_levels)
  # Should report both person and rater as having different levels
  expect_true(grepl("person", result$warning_message))
  expect_true(grepl("rater", result$warning_message))
})

test_that("check_multivariate_balance info message has correct format", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  test_data$y1[1:5] <- NA
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  result <- check_multivariate_balance(formula, test_data, FALSE)
  
  expect_true(grepl("Original observations: 100", result$info_message))
  expect_true(grepl("Complete cases used: 95", result$info_message))
  expect_true(grepl("5 rows removed", result$info_message))
  expect_true(grepl("consistent", result$info_message))
})

test_that("check_multivariate_balance warning message has correct format", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  
  test_data$y1[test_data$person == "1"] <- NA
  
  formula <- brms::mvbind(y1, y2) ~ (1 | person) + (1 | rater)
  result <- check_multivariate_balance(formula, test_data, FALSE)
  
  expect_true(grepl("Unbalanced multivariate design", result$warning_message))
  expect_true(grepl("person:", result$warning_message))
  expect_true(grepl("long-format", result$warning_message))
  expect_true(grepl("bf\\(", result$warning_message))
})

test_that("check_multivariate_balance works with mvbrmsformula and set_rescor", {
  skip_if_not_installed("brms")
  
  test_data <- data.frame(
    test1 = rnorm(100),
    test2 = rnorm(100),
    test3 = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:4, each = 25))
  )
  
  # Make test1 missing for items 3 and 4
  test_data$test1[test_data$item %in% c("3", "4")] <- NA
  
  formula <- brms::bf(brms::mvbind(test1, test2, test3) ~ (1 | person) + (1 | item)) + 
    brms::set_rescor(FALSE)
  result <- check_multivariate_balance(formula, test_data, FALSE)
  
  expect_equal(sort(result$dimensions), c("test1", "test2", "test3"))
  expect_false(result$is_balanced)
  expect_true(result$has_different_levels)
  expect_true(grepl("item:", result$warning_message))
})

test_that("extract_response_names works with mvbrmsformula objects", {
  skip_if_not_installed("brms")
  
  formula <- brms::bf(brms::mvbind(y1, y2, y3) ~ (1 | person)) + brms::set_rescor(TRUE)
  result <- extract_response_names(formula)
  
  expect_equal(sort(result), c("y1", "y2", "y3"))
})
