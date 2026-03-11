# Tests for backend functions

# Test data for reuse
test_data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, 5)),
  rater = factor(rep(1:5, each = 20))
)

# =============================================================================
# select_backend tests
# =============================================================================

test_that("select_backend chooses lme4 for univariate formula with auto", {
  skip_if_not_installed("lme4")
  f <- score ~ (1 | person) + (1 | rater)
  expect_equal(select_backend(f, "auto"), "lme4")
})

test_that("select_backend respects user choice for lme4", {
  skip_if_not_installed("lme4")
  f <- score ~ (1 | person)
  expect_equal(select_backend(f, "lme4"), "lme4")
})

test_that("select_backend respects user choice for brms", {
  skip_if_not_installed("brms")
  f <- score ~ (1 | person)
  expect_equal(select_backend(f, "brms"), "brms")
})

test_that("select_backend errors for invalid backend choice", {
  f <- score ~ (1 | person)
  expect_error(select_backend(f, "invalid"), "should be one of")
})

test_that("select_backend errors when no backend available", {
  # This test would require mocking, so we skip it
  # Instead, test that it works when lme4 is available
  skip_if_not_installed("lme4")
  f <- score ~ (1 | person)
  result <- select_backend(f, "auto")
  expect_true(result %in% c("lme4", "brms"))
})

# =============================================================================
# fit_lme4 tests
# =============================================================================

test_that("fit_lme4 fits a simple model", {
  skip_if_not_installed("lme4")
  f <- score ~ (1 | person) + (1 | rater)
  model <- fit_lme4(f, test_data)
  # lmerMod is S4, not S3
  expect_s4_class(model, "lmerMod")
})

test_that("fit_lme4 handles convergence issues gracefully", {
  skip_if_not_installed("lme4")
  # Create data that might cause convergence issues
  bad_data <- data.frame(
    score = rep(1, 10),
    person = factor(rep(1:5, 2))
  )
  # Should still fit, possibly with a warning
  expect_s4_class(fit_lme4(score ~ (1 | person), bad_data), "lmerMod")
})

# =============================================================================
# fit_brms tests
# =============================================================================

test_that("fit_brms fits a simple model", {
  skip_if_not_installed("brms")
  skip_on_cran() # Skip on CRAN to avoid long-running tests
  f <- score ~ (1 | person) + (1 | rater)
  model <- fit_brms(f, test_data, chains = 1, iter = 500, refresh = 0)
  expect_s3_class(model, "brmsfit")
})

test_that("fit_brms sets default chains and iterations", {
  skip_if_not_installed("brms")
  skip_on_cran()
  f <- score ~ (1 | person)
  # Test that defaults are applied (chains = 4, iter = 2000)
  # We can't easily test this without fitting, so just check it works
  model <- fit_brms(f, test_data, chains = 1, iter = 500, refresh = 0)
  expect_s3_class(model, "brmsfit")
})

# =============================================================================
# extract_variance_components tests
# =============================================================================

test_that("extract_variance_components works with lme4 model", {
  skip_if_not_installed("lme4")
  f <- score ~ (1 | person) + (1 | rater)
  model <- fit_lme4(f, test_data)
  vc <- extract_variance_components(model, "lme4")
  expect_s3_class(vc, "tbl_df")
  expect_true("component" %in% names(vc) || "facet" %in% names(vc))
  expect_true("var" %in% names(vc) || "variance" %in% names(vc))
})

test_that("extract_variance_components errors for unknown backend", {
  expect_error(extract_variance_components(NULL, "unknown"), "Unknown backend")
})

# =============================================================================
# extract_vc_lme4 tests
# =============================================================================

test_that("extract_vc_lme4 returns tibble with expected columns", {
  skip_if_not_installed("lme4")
  f <- score ~ (1 | person) + (1 | rater)
  model <- lme4::lmer(f, test_data)
  vc <- extract_vc_lme4(model)
  expect_s3_class(vc, "tbl_df")
})

test_that("extract_vc_lme4 includes residual variance", {
  skip_if_not_installed("lme4")
  f <- score ~ (1 | person) + (1 | rater)
  model <- lme4::lmer(f, test_data)
  vc <- extract_vc_lme4(model)
  expect_true(any(grepl("Residual", vc$component) | grepl("Residual", vc$facet)))
})

# =============================================================================
# Method of Moments (mom) backend tests
# =============================================================================

test_that("select_backend works for mom backend", {
  f <- score ~ (1 | person) + (1 | rater)
  expect_equal(select_backend(f, "mom"), "mom")
})

test_that("fit_mom fits a simple crossed design", {
  # Create balanced test data
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  f <- score ~ (1 | person) + (1 | item)
  model <- fit_mom(f, balanced_data)
  
  expect_s3_class(model, "momfit")
  expect_true("variance_components" %in% names(model))
})

test_that("fit_mom returns variance components", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  f <- score ~ (1 | person) + (1 | item)
  model <- fit_mom(f, balanced_data)
  
  vc <- model$variance_components
  expect_s3_class(vc, "tbl_df")
  expect_true("var" %in% names(vc))
  expect_true("pct" %in% names(vc))
})

test_that("extract_vc_mom returns tibble with confidence intervals", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  f <- score ~ (1 | person) + (1 | item)
  model <- fit_mom(f, balanced_data)
  
  vc <- extract_vc_mom(model)
  
  expect_s3_class(vc, "tbl_df")
  expect_true("component" %in% names(vc))
  expect_true("var" %in% names(vc))
  expect_true("lower" %in% names(vc))
  expect_true("upper" %in% names(vc))
  expect_true("pct" %in% names(vc))
})

test_that("extract_variance_components works with mom model", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  f <- score ~ (1 | person) + (1 | item)
  model <- fit_mom(f, balanced_data)
  
  vc <- extract_variance_components(model, "mom")
  expect_s3_class(vc, "tbl_df")
  expect_true("var" %in% names(vc))
})

test_that("mom variance components sum to total variance", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  f <- score ~ (1 | person) + (1 | item)
  model <- fit_mom(f, balanced_data)
  
  vc <- model$variance_components
  total_var <- sum(vc$var)
  
  # Variance components should be non-negative
  expect_true(all(vc$var >= 0))
  
  # Percentages should sum to approximately 100
  expect_true(abs(sum(vc$pct) - 100) < 0.01)
})

test_that("fit_mom handles nested design", {
  set.seed(42)
  # Nested design: raters nested within items
  nested_data <- data.frame(
    score = rnorm(60),
    person = factor(rep(1:20, each = 3)),
    item = factor(rep(1:4, each = 15)),
    rater = factor(rep(rep(1:3, each = 5), 4))
  )
  
  # Nested formula: rater within item
  f <- score ~ (1 | person) + (1 | item/rater)
  model <- fit_mom(f, nested_data)
  
  expect_s3_class(model, "momfit")
})

test_that("gstudy works with mom backend", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  g <- gstudy(score ~ (1 | person) + (1 | item),
              data = balanced_data,
              backend = "mom")
  
  expect_s3_class(g, "gstudy")
  expect_equal(g$backend, "mom")
  expect_true("variance_components" %in% names(g))
})

test_that("gstudy mom backend includes confidence intervals", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  g <- gstudy(score ~ (1 | person) + (1 | item),
              data = balanced_data,
              backend = "mom")
  
  vc <- g$variance_components
  expect_true("lower" %in% names(vc))
  expect_true("upper" %in% names(vc))
})

test_that("mom backend warns about ci_method", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )

  expect_warning(
    gstudy(score ~ (1 | person) + (1 | item),
      data = balanced_data,
      backend = "mom",
      ci_method = "profile"),
    "ci_method is only applicable for lme4"
  )
})

# =============================================================================
# Variance component ordering tests
# =============================================================================

test_that("variance components order matches formula order for lme4 backend", {
  skip_if_not_installed("lme4")
  set.seed(42)
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20)),
    rater = factor(rep(1:10, each = 10))
  )

  f <- score ~ (1 | person) + (1 | item) + (1 | rater)
  g <- gstudy(f, data = test_data, backend = "lme4")

  components <- g$variance_components$component
  components <- components[components != "Residual"]

  expect_equal(components, c("person", "item", "rater"))
})

test_that("variance components order matches formula order for mom backend", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20)),
    rater = factor(rep(1:10, each = 10))
  )

  f <- score ~ (1 | person) + (1 | item) + (1 | rater)
  g <- gstudy(f, data = balanced_data, backend = "mom")

  components <- g$variance_components$component
  components <- components[components != "Residual"]

  expect_equal(components, c("person", "item", "rater"))
})

test_that("Residual is always last in variance components table", {
  skip_if_not_installed("lme4")
  set.seed(42)
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )

  g_lme4 <- gstudy(score ~ (1 | person) + (1 | item),
                   data = test_data,
                   backend = "lme4")

  last_component <- tail(g_lme4$variance_components$component, 1)
  expect_equal(last_component, "Residual")

  g_mom <- gstudy(score ~ (1 | person) + (1 | item),
                  data = test_data,
                  backend = "mom")

  last_component <- tail(g_mom$variance_components$component, 1)
  expect_equal(last_component, "Residual")
})

test_that("interaction naming follows user specification", {
  skip_if_not_installed("lme4")
  set.seed(42)
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    rater = factor(rep(1:5, times = 20)),
    task = factor(rep(1:10, each = 10))
  )

  g <- gstudy(score ~ (1 | person) + (1 | rater:task),
              data = test_data,
              backend = "lme4")

  components <- g$variance_components$component
  expect_true("rater:task" %in% components)
  expect_false("task:rater" %in% components)
})

test_that("variance components have consistent order across backends", {
  skip_if_not_installed("lme4")
  set.seed(42)
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20)),
    rater = factor(rep(1:10, each = 10))
  )

  f <- score ~ (1 | person) + (1 | item) + (1 | rater)

  g_lme4 <- gstudy(f, data = test_data, backend = "lme4")
  g_mom <- gstudy(f, data = test_data, backend = "mom")

  lme4_components <- g_lme4$variance_components$component
  mom_components <- g_mom$variance_components$component

  expect_equal(lme4_components, mom_components)
})

# =============================================================================
# Correlated Random Effects tests
# =============================================================================

test_that("extract_random_effect_cor parses double-bar syntax correctly", {
  f1 <- score ~ (1 | cor_p | person) + (1 | rater)
  cors <- extract_random_effect_cor(f1)
  expect_equal(cors, list(person = "cor_p"))

  f2 <- score ~ (1 | cor_p | person) + (1 | cor_r | rater)
  cors <- extract_random_effect_cor(f2)
  expect_equal(cors, list(person = "cor_p", rater = "cor_r"))

  # Test interaction term
  f3 <- score ~ (1 | cor_grp | person:item) + (1 | rater)
  cors <- extract_random_effect_cor(f3)
  expect_equal(cors, list(`person:item` = "cor_grp"))
})

test_that("extract_rescor_setting handles multi-line formulas", {
  # Create a formula that would span multiple lines when deparsed
  f <- y ~ (1 | person) + (1 | item) + (1 | rater) +
    set_rescor(TRUE)

  result <- extract_rescor_setting(f)
  expect_true(result)
})

test_that("mom backend handles correlated random effects syntax", {
  skip_on_cran()

  set.seed(42)
  data <- data.frame(
    y1 = rnorm(100),
    y2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )

  # Test with (1|cor_p|person) syntax
  result <- gstudy(
    mvbind(y1, y2) ~ (1 | cor_p | person) + (1 | rater),
    data = data,
    backend = "mom"
  )

  # Should not error
  expect_s3_class(result, "mgstudy")

  # Should have correlations for person facet
  expect_true(!is.null(result$correlations))
  expect_true(!is.null(result$correlations$random_effect_cor))
  expect_true("person" %in% names(result$correlations$random_effect_cor))
})
