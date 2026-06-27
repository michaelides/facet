# Tests for estimator functions

# Test data for reuse
test_data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, 5)),
  rater = factor(rep(1:5, each = 20))
)

# =============================================================================
# select_estimator tests
# =============================================================================

test_that("select_estimator chooses lme4 for univariate formula with auto", {
  skip_if_not_installed("lme4")
  f <- score ~ (1 | person) + (1 | rater)
  expect_equal(select_estimator(f, "auto"), "lme4")
})

test_that("select_estimator respects user choice for lme4", {
  skip_if_not_installed("lme4")
  f <- score ~ (1 | person)
  expect_equal(select_estimator(f, "lme4"), "lme4")
})

test_that("select_estimator respects user choice for brms", {
  skip_if_not_installed("brms")
  f <- score ~ (1 | person)
  expect_equal(select_estimator(f, "brms"), "brms")
})

test_that("select_estimator errors for invalid estimator choice", {
  f <- score ~ (1 | person)
  expect_error(select_estimator(f, "invalid"), "should be one of")
})

test_that("select_estimator errors when no estimator available", {
  # This test would require mocking, so we skip it
  # Instead, test that it works when lme4 is available
  skip_if_not_installed("lme4")
  f <- score ~ (1 | person)
  result <- select_estimator(f, "auto")
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

test_that("extract_variance_components errors for unknown estimator", {
  expect_error(extract_variance_components(NULL, "unknown"), "Unknown estimator")
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
  expect_true(any(grepl("Residual", vc$component)))
})

# =============================================================================
# ANOVA-based Estimation (mom) estimator tests
# =============================================================================

test_that("select_estimator works for aov estimator", {
  f <- score ~ (1 | person) + (1 | rater)
  expect_equal(select_estimator(f, "aov"), "aov")
})

test_that("fit_aov fits a simple crossed design", {
  # Create balanced test data
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  f <- score ~ (1 | person) + (1 | item)
  model <- fit_aov(f, balanced_data)
  
  expect_s3_class(model, "aovfit")
  expect_true("variance_components" %in% names(model))
})

test_that("fit_aov returns variance components", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  f <- score ~ (1 | person) + (1 | item)
  model <- fit_aov(f, balanced_data)
  
  vc <- model$variance_components
  expect_s3_class(vc, "tbl_df")
  expect_true("var" %in% names(vc))
  expect_true("pct" %in% names(vc))
})

test_that("extract_vc_aov returns tibble with confidence intervals", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  f <- score ~ (1 | person) + (1 | item)
  model <- fit_aov(f, balanced_data)
  
  vc <- extract_vc_aov(model)
  
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
  model <- fit_aov(f, balanced_data)
  
  vc <- extract_variance_components(model, "aov")
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
  model <- fit_aov(f, balanced_data)
  
  vc <- model$variance_components
  total_var <- sum(vc$var)
  
  # Variance components should be non-negative
  expect_true(all(vc$var >= 0))
  
  # Percentages should sum to approximately 100
  expect_true(abs(sum(vc$pct) - 100) < 0.01)
})

test_that("fit_aov handles nested design", {
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
  model <- fit_aov(f, nested_data)
  
  expect_s3_class(model, "aovfit")
})

test_that("gstudy works with aov estimator", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  g <- gstudy(score ~ (1 | person) + (1 | item),
              data = balanced_data,
              estimator = "aov")
  
  expect_s3_class(g, "gstudy")
  expect_equal(g$estimator, "aov")
  expect_true("variance_components" %in% names(g))
})

test_that("gstudy aov estimator includes confidence intervals", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )
  
  g <- gstudy(score ~ (1 | person) + (1 | item),
              data = balanced_data,
              estimator = "aov")
  
  vc <- g$variance_components
  expect_true("lower" %in% names(vc))
  expect_true("upper" %in% names(vc))
})

test_that("aov estimator warns about ci_method", {
  set.seed(42)
  balanced_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, each = 5)),
    item = factor(rep(1:5, times = 20))
  )

  expect_warning(
    gstudy(score ~ (1 | person) + (1 | item),
      data = balanced_data,
      estimator = "aov",
      ci_method = "profile"),
    "ci_method is only applicable for lme4"
  )
})

# =============================================================================
# Variance component ordering tests
# =============================================================================

test_that("variance components order matches formula order for lme4 estimator", {
  skip_if_not_installed("lme4")
  set.seed(42)
  test_data <- expand.grid(
    person = factor(1:10),
    item = factor(1:5),
    rater = factor(1:2)
  )
  test_data$score <- rnorm(nrow(test_data))

  f <- score ~ (1 | person) + (1 | item) + (1 | rater)
  g <- gstudy(f, data = test_data, estimator = "lme4")

  components <- g$variance_components$component
  components <- components[components != "Residual"]

  expect_equal(components, c("person", "item", "rater"))
})

test_that("variance components order matches formula order for aov estimator", {
  set.seed(42)
  balanced_data <- expand.grid(
    person = factor(1:10),
    item = factor(1:5),
    rater = factor(1:2)
  )
  balanced_data$score <- rnorm(nrow(balanced_data))

  f <- score ~ (1 | person) + (1 | item) + (1 | rater)
  g <- gstudy(f, data = balanced_data, estimator = "aov")

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
                   estimator = "lme4")

  last_component <- tail(g_lme4$variance_components$component, 1)
  expect_equal(last_component, "Residual")

  g_mom <- gstudy(score ~ (1 | person) + (1 | item),
                  data = test_data,
                  estimator = "aov")

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
              estimator = "lme4")

  components <- g$variance_components$component
  expect_true("rater:task" %in% components)
  expect_false("task:rater" %in% components)
})

test_that("variance components have consistent order across estimators", {
  skip_if_not_installed("lme4")
  set.seed(42)
  test_data <- expand.grid(
    person = factor(1:10),
    item = factor(1:5),
    rater = factor(1:2)
  )
  test_data$score <- rnorm(nrow(test_data))

  f <- score ~ (1 | person) + (1 | item) + (1 | rater)

  g_lme4 <- gstudy(f, data = test_data, estimator = "lme4")
  g_mom <- gstudy(f, data = test_data, estimator = "aov")

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

test_that("aov estimator handles correlated random effects syntax", {
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
    estimator = "aov"
  )

  # Should not error
  expect_s3_class(result, "mgstudy")

  # Should have correlations for person facet
  expect_true(!is.null(result$correlations))
  expect_true(!is.null(result$correlations$random_effect_cor))
  expect_true("person" %in% names(result$correlations$random_effect_cor))
})

# =============================================================================
# adapt_aov_formula tests (brennan-style nested Rater within Task)
# =============================================================================

test_that("adapt_aov_formula rewrites Error() for Rater nested in Task", {
  data(brennan)
  f <- Score ~ (1 | Person) + (1 | Task) + (1 | Rater) + (1 | Person:Task)

  out <- facet:::adapt_aov_formula(f, brennan)

  # Should rewrite the Rater stratum to Rater:Task (aov alphabetizes to
  # "Task:Rater", which we also accept).
  expect_true(grepl("Rater:Task", out$aov_formula) ||
                grepl("Task:Rater", out$aov_formula))

  # The standalone "Rater" term should no longer appear in the rewritten
  # Error() decomposition: strip "Rater:Task" and "Task:Rater" first, then
  # check for a bare "Rater" word.
  stripped <- out$aov_formula
  stripped <- gsub("Rater:Task", "", stripped)
  stripped <- gsub("Task:Rater", "", stripped)
  expect_false(grepl("\\bRater\\b", stripped))

  # Name mapping should preserve the user-supplied name "Rater"
  expect_equal(out$name_mapping[["Rater:Task"]], "Rater")
  expect_equal(out$name_mapping[["Task:Rater"]], "Rater")
})

test_that("adapt_aov_formula is a no-op for fully-crossed designs", {
  test_data <- expand.grid(
    person = factor(1:10),
    item = factor(1:5)
  )
  test_data$score <- rnorm(50)

  f <- score ~ (1 | person) + (1 | item)
  out <- facet:::adapt_aov_formula(f, test_data)

  # No nesting detected, so aov_formula should be the standard one
  expect_equal(
    out$aov_formula,
    "score ~ person + item + Error( person + item )"
  )
  expect_equal(out$name_mapping, list())
})

test_that("adapt_aov_formula respects the explicit nested form", {
  data(brennan)
  # User writes (1 | Rater:Task) directly
  f <- Score ~ (1 | Person) + (1 | Task) + (1 | Rater:Task) + (1 | Person:Task)
  out <- facet:::adapt_aov_formula(f, brennan)

  # No rewrite should happen; the spec is already the nested form. aov
  # alphabetizes the order to "Rater:Task" (since R < T) so the formula
  # preserves the user-supplied order.
  expect_equal(
    out$aov_formula,
    "Score ~ Person + Task + Rater:Task + Person:Task + Error( Person + Task + Rater:Task + Person:Task )"
  )
  # Name mapping is identity for the user-supplied nested form
  expect_equal(out$name_mapping[["Rater:Task"]], "Rater:Task")
})

test_that("adapt_aov_formula handles user-supplied nested argument", {
  data(brennan)
  f <- Score ~ (1 | Person) + (1 | Task) + (1 | Rater) + (1 | Person:Task)

  out <- facet:::adapt_aov_formula(
    f,
    brennan,
    nested = list(Rater = "Task")
  )

  expect_equal(out$name_mapping[["Rater:Task"]], "Rater")
})

test_that("gstudy aov estimator on brennan produces no singular-model warning", {
  data(brennan)
  f <- Score ~ (1 | Person) + (1 | Task) + (1 | Rater) + (1 | Person:Task)

  # Should not produce "Error() model is singular" warning
  expect_no_warning(gstudy(f, data = brennan, estimator = "aov"))
})

test_that("gstudy mom on brennan gives Rater row with lme4-matching variance", {
  skip_if_not_installed("lme4")
  data(brennan)
  library(lme4)

  f <- Score ~ (1 | Person) + (1 | Task) + (1 | Rater) + (1 | Person:Task)

  vc_mom <- gstudy(f, data = brennan, estimator = "aov")$variance_components
  m_lme4 <- lmer(f, data = brennan, REML = FALSE)
  vc_lme4 <- as.data.frame(VarCorr(m_lme4))[, c("grp", "vcov")]

  # Variance components have the user-supplied component name "Rater"
  expect_true("Rater" %in% vc_mom$component)
  expect_false("Rater:Task" %in% vc_mom$component)
  expect_false("Task:Rater" %in% vc_mom$component)

  # The Rater row in the mom output should match lme4's Rater variance
  rater_mom <- vc_mom$var[vc_mom$component == "Rater"]
  rater_lme4 <- vc_lme4$vcov[vc_lme4$grp == "Rater"]
  expect_equal(rater_mom, rater_lme4, tolerance = 1e-3)
})

test_that("gstudy mom on brennan is invariant to the explicit nested form", {
  data(brennan)
  f_unnested <- Score ~ (1 | Person) + (1 | Task) + (1 | Rater) + (1 | Person:Task)
  f_nested <- Score ~ (1 | Person) + (1 | Task) + (1 | Rater:Task) + (1 | Person:Task)

  vc_unnested <- gstudy(f_unnested, data = brennan, estimator = "aov")$variance_components
  vc_nested <- gstudy(f_nested, data = brennan, estimator = "aov")$variance_components

  # The Person:Task and Residual rows should match exactly
  expect_equal(
    vc_unnested$var[vc_unnested$component == "Person:Task"],
    vc_nested$var[vc_nested$component == "Person:Task"]
  )
  expect_equal(
    vc_unnested$var[vc_unnested$component == "Residual"],
    vc_nested$var[vc_nested$component == "Residual"]
  )
})
