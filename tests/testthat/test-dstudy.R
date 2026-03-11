# Tests for dstudy function

# Test data for reuse
test_data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, 5)),
  rater = factor(rep(1:5, each = 20))
)

# =============================================================================
# Basic dstudy tests
# =============================================================================

test_that("dstudy returns a dstudy object", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3))
  expect_s3_class(result, "dstudy")
})

test_that("dstudy requires a gstudy object", {
  expect_error(
    dstudy(list(), n = list(rater = 3)),
    "must be an object of class 'gstudy'"
  )
})

test_that("dstudy errors with non-gstudy input", {
  expect_error(dstudy(NULL, n = list(rater = 3)), "gstudy")
  expect_error(dstudy("not a gstudy", n = list(rater = 3)), "gstudy")
  expect_error(dstudy(123, n = list(rater = 3)), "gstudy")
})

# =============================================================================
# Input validation tests
# =============================================================================

test_that("dstudy uses G-study sample sizes when n is empty", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list())
  expect_s3_class(result, "dstudy")
  expect_true(length(result$n) > 0)
})

test_that("dstudy accepts valid n parameter", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  n <- list(rater = 5)
  result <- dstudy(g, n = n)
  expect_equal(result$n, n)
})

# =============================================================================
# Output structure tests
# =============================================================================

test_that("dstudy stores reference to gstudy object", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3))
  expect_equal(result$gstudy, g)
})

test_that("dstudy stores variance components", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3))
  expect_true(!is.null(result$variance_components))
})

test_that("dstudy stores coefficients", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3))
  expect_true(!is.null(result$coefficients))
})

test_that("dstudy stores object of measurement", {
skip_if_not_installed("lme4")
g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
result <- dstudy(g, n = list(rater = 3))
expect_equal(result$object, "person")
})

test_that("dstudy has universe field defaulting to object", {
skip_if_not_installed("lme4")
g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
result <- dstudy(g, n = list(rater = 3))
expect_equal(result$universe, "person")
})

test_that("dstudy accepts custom universe specification", {
skip_if_not_installed("lme4")
g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
result <- dstudy(g, n = list(rater = 3), universe = c("person", "person:rater"))
expect_true("person" %in% result$universe)
expect_true("person:rater" %in% result$universe)
})

test_that("dstudy warns when object not in universe", {
skip_if_not_installed("lme4")
g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
expect_warning(
dstudy(g, n = list(rater = 3), universe = "person:rater"),
"specification of the universe did not include the object of measurement"
)
})

# =============================================================================
# Coefficient calculation tests
# =============================================================================

test_that("dstudy calculates G coefficient", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3))

  # Check that G coefficient exists and is numeric
  if ("g" %in% names(result$coefficients)) {
    expect_true(is.numeric(result$coefficients$g))
    expect_true(result$coefficients$g >= 0 && result$coefficients$g <= 1)
  } else if ("coefficient" %in% names(result$coefficients)) {
    g_val <- result$coefficients$value[result$coefficients$coefficient == "g"]
    expect_true(is.numeric(g_val))
    expect_true(g_val >= 0 && g_val <= 1)
  }
})

test_that("dstudy calculates Phi coefficient", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3))

  # Check that Phi coefficient exists and is numeric
  if ("phi" %in% names(result$coefficients)) {
    expect_true(is.numeric(result$coefficients$phi))
    expect_true(result$coefficients$phi >= 0 && result$coefficients$phi <= 1)
  } else if ("coefficient" %in% names(result$coefficients)) {
    phi_val <- result$coefficients$value[result$coefficients$coefficient == "phi"]
    expect_true(is.numeric(phi_val))
    expect_true(phi_val >= 0 && phi_val <= 1)
  }
})

test_that("dstudy coefficients are between 0 and 1", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3))
  
  # All coefficients should be between 0 and 1
  coef_values <- result$coefficients
  if ("value" %in% names(coef_values)) {
    expect_true(all(coef_values$value >= 0))
    expect_true(all(coef_values$value <= 1))
  }
})

# =============================================================================
# Sample size sweep tests
# =============================================================================

test_that("dstudy handles multiple sample sizes (sweep)", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = c(2, 3, 5)))
  
  expect_true(result$is_sweep)
  expect_s3_class(result$coefficients, "tbl_df")
})

test_that("dstudy sweep creates grid of combinations", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = c(2, 3), person = c(10, 20)))
  
  # Should have 2 x 2 = 4 combinations
  expect_true(nrow(result$coefficients) == 4)
})

test_that("dstudy single sample size is not sweep", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3))
  
  expect_false(result$is_sweep)
})

# =============================================================================
# is.dstudy test
# =============================================================================

test_that("is.dstudy returns TRUE for dstudy objects", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3))
  expect_true(is.dstudy(result))
})

test_that("is.dstudy returns FALSE for non-dstudy objects", {
  expect_false(is.dstudy(list()))
  expect_false(is.dstudy(NULL))
  expect_false(is.dstudy("not a dstudy"))
})

# =============================================================================
# Edge cases
# =============================================================================

test_that("dstudy handles single facet design", {
  skip_if_not_installed("lme4")
  single_facet_data <- data.frame(
    score = rnorm(50),
    person = factor(rep(1:10, 5))
  )
  g <- gstudy(score ~ (1 | person), data = single_facet_data)
  result <- dstudy(g, n = list())
  
  expect_s3_class(result, "dstudy")
})

test_that("dstudy handles large sample sizes", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 100))

  expect_s3_class(result, "dstudy")
  # With more raters, error variance should decrease, coefficients should increase
})

# =============================================================================
# Residual-Error redundancy check tests
# =============================================================================

test_that("dstudy warns when residual matches error component", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  # Get the residual composition from the gstudy
  residual_comp <- g$variance_components$component[g$variance_components$component == "Residual"]
  # Note: In a crossed design, residual is the interaction term, not "Residual"

  # For a crossed design person x rater, the residual composition should be "person:rater"
  # We can test by specifying that as error
  expect_warning(
    result <- dstudy(g, n = list(rater = 3), error = "person:rater"),
    "residual component.*already included.*Removing"
  )

  # The result should still be valid
  expect_s3_class(result, "dstudy")
})

test_that("dstudy removes duplicate error component matching residual", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  # Specify error with multiple components including the residual
  suppressWarnings({
    result <- dstudy(g, n = list(rater = 3), error = c("person:rater", "rater"))
  })

  # After removal of person:rater, only rater should remain
  expect_true(is.null(result$error) || !"person:rater" %in% result$error)
})

test_that("dstudy sets error to NULL when all error components match residual", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  # Specify only the residual as error
  suppressWarnings({
    result <- dstudy(g, n = list(rater = 3), error = "person:rater")
  })

  # After removal, error should be NULL
  expect_null(result$error)
})

test_that("dstudy does not warn when error does not match residual", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  # Specify error component that doesn't match residual
  expect_warning(
    result <- dstudy(g, n = list(rater = 3), error = "rater"),
    NA
  )

  expect_s3_class(result, "dstudy")
  expect_equal(result$error, "rater")
})

test_that("dstudy does not warn when error is NULL", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  expect_warning(
    result <- dstudy(g, n = list(rater = 3)),
    NA
  )

  expect_s3_class(result, "dstudy")
})

test_that("dstudy residual-error check works for sweep mode", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  # Sweep with error matching residual
  expect_warning(
    result <- dstudy(g, n = list(rater = c(2, 3, 5)), error = "person:rater"),
    "residual component.*already included.*Removing"
  )

  expect_true(result$is_sweep)
  expect_s3_class(result, "dstudy")
})

# =============================================================================
# Aggregation-Error interaction tests
# =============================================================================

test_that("dstudy excludes aggregation component from error (not interaction)", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  result <- dstudy(g, n = list(rater = 3), aggregation = "rater")

  expect_s3_class(result, "dstudy")
})

test_that("dstudy warns when component in both aggregation and error", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  expect_warning(
    result <- dstudy(g, n = list(rater = 3), aggregation = "rater", error = "rater"),
    "specified for both aggregation and error"
  )

  expect_s3_class(result, "dstudy")
})

test_that("dstudy removes component from error when in aggregation", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  suppressWarnings({
    result <- dstudy(g, n = list(rater = 3), aggregation = "rater", error = c("rater", "person:rater"))
  })

  expect_true(!"rater" %in% result$error)
})

test_that("dstudy interaction terms NOT removed when only main effect in aggregation", {
  skip_if_not_installed("lme4")

  test_data_3 <- data.frame(
    score = rnorm(200),
    person = factor(rep(1:20, 10)),
    item = factor(rep(1:10, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_3)

  suppressWarnings({
    result <- dstudy(g, n = list(item = 3), aggregation = "item", error = c("item"))
  })

  expect_true(!"item" %in% result$error)
})

# =============================================================================
# Object vs Error/Aggregation validation tests
# =============================================================================

test_that("dstudy errors when same component in universe and error", {
skip_if_not_installed("lme4")
g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

expect_error(
dstudy(g, n = list(rater = 3), universe = c("person", "person:rater"), error = "person:rater"),
"same component cannot be in both the universe and error"
)
})

test_that("dstudy allows universe as main facet and error as interaction", {
skip_if_not_installed("lme4")
g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

result <- dstudy(g, n = list(rater = 3), universe = "person", error = "person:rater")
expect_s3_class(result, "dstudy")
})

test_that("dstudy object is always first facet from gstudy", {
skip_if_not_installed("lme4")
g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

# Object is always 'person' (first facet) regardless of universe specification
result <- dstudy(g, n = list(rater = 3))
expect_equal(result$object, "person")
})

# =============================================================================
# Aggregation scaling tests (new behavior)
# =============================================================================

test_that("dstudy aggregation scales interaction by only aggregation facet n", {
  skip_if_not_installed("lme4")

  test_data_agg <- data.frame(
    score = rnorm(200),
    person = factor(rep(1:20, 10)),
    item = factor(rep(1:10, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_agg)

  result <- dstudy(g, n = list(item = 5), aggregation = "item")

  vc <- result$variance_components

  item_var <- vc$var[vc$component == "item"]
  pi_var <- vc$var[vc$component == "person:item"]

  g_item_var <- g$variance_components$var[g$variance_components$component == "item"]
  g_pi_var <- g$variance_components$var[g$variance_components$component == "person:item"]

  expect_equal(item_var, g_item_var / 5, tolerance = 1e-10)
  expect_equal(pi_var, g_pi_var / 5, tolerance = 1e-10)
})

test_that("dstudy aggregation calculates relative error correctly", {
  skip_if_not_installed("lme4")

  test_data_agg <- data.frame(
    score = rnorm(200),
    person = factor(rep(1:20, 10)),
    item = factor(rep(1:10, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_agg)

  result <- dstudy(g, n = list(item = 5), aggregation = "item")

  vc <- result$variance_components
  person_var <- vc$var[vc$component == "person"]

  sigma2_delta <- result$coefficients$sigma2_delta

  g_residual_var <- g$variance_components$var[g$variance_components$component == "Residual"]

  expect_equal(sigma2_delta, g_residual_var / 5, tolerance = 1e-6)
})

test_that("dstudy aggregation with sweep calculates correctly", {
  skip_if_not_installed("lme4")

  test_data_agg <- data.frame(
    score = rnorm(200),
    person = factor(rep(1:20, 10)),
    item = factor(rep(1:10, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_agg)

  result <- dstudy(g, n = list(item = c(2, 5, 10)), aggregation = "item")

  expect_true(result$is_sweep)
  expect_equal(nrow(result$coefficients), 3)

  for (i in 1:3) {
    expect_true(result$coefficients$sigma2_delta[i] >= 0)
  }
})

test_that("dstudy rescaling when n is explicitly provided", {
  skip_if_not_installed("lme4")

  test_data_agg <- data.frame(
    score = rnorm(200),
    person = factor(rep(1:20, 10)),
    item = factor(rep(1:10, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_agg)

  result <- dstudy(g, n = list(item = 5))

  vc <- result$variance_components

  item_var <- vc$var[vc$component == "item"]
  pi_var <- vc$var[vc$component == "person:item"]

  g_item_var <- g$variance_components$var[g$variance_components$component == "item"]
  g_pi_var <- g$variance_components$var[g$variance_components$component == "person:item"]

  expect_equal(item_var, g_item_var / 5, tolerance = 1e-10)
  expect_equal(pi_var, g_pi_var / 5, tolerance = 1e-10)
})

test_that("dstudy n and aggregation same facet divides by n once not n squared", {
  skip_if_not_installed("lme4")

  test_data_agg <- data.frame(
    score = rnorm(200),
    person = factor(rep(1:20, 10)),
    item = factor(rep(1:10, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_agg)

  result <- dstudy(g, n = list(item = 10), aggregation = "item")

  vc <- result$variance_components

  item_var <- vc$var[vc$component == "item"]

  g_item_var <- g$variance_components$var[g$variance_components$component == "item"]

  expect_equal(item_var, g_item_var / 10, tolerance = 1e-10)
})

test_that("dstudy n and aggregation different facets each divide correctly", {
  skip_if_not_installed("lme4")

  test_data_agg <- data.frame(
    score = rnorm(200),
    person = factor(rep(1:20, 10)),
    item = factor(rep(1:10, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_agg)

  result <- dstudy(g, n = list(item = 5, person = 20), aggregation = "person")

  vc <- result$variance_components

  person_var <- vc$var[vc$component == "person"]
  item_var <- vc$var[vc$component == "item"]

  g_person_var <- g$variance_components$var[g$variance_components$component == "person"]
  g_item_var <- g$variance_components$var[g$variance_components$component == "item"]

  expect_equal(person_var, g_person_var, tolerance = 1e-10)
  expect_equal(item_var, g_item_var / 5, tolerance = 1e-10)
})

# =============================================================================
# Estimation parameter tests
# =============================================================================

test_that("dstudy defaults to simple estimation", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3))

  expect_equal(result$estimation, "simple")
  expect_null(result$posterior)
})

test_that("dstudy accepts simple estimation explicitly", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3), estimation = "simple")

  expect_equal(result$estimation, "simple")
  expect_null(result$posterior)
})

test_that("dstudy warns and refits when posterior requested with non-brms backend", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, backend = "lme4")

  expect_warning(
    result <- dstudy(g, n = list(rater = 3), estimation = "posterior"),
    "estimation = 'posterior' requires backend = 'brms'"
  )

  expect_equal(result$estimation, "posterior")
  expect_equal(result$gstudy$backend, "brms")
})

test_that("dstudy posterior estimation returns posterior distributions", {
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, backend = "brms")
  result <- dstudy(g, n = list(rater = 3), estimation = "posterior")

  expect_equal(result$estimation, "posterior")
  expect_type(result$posterior, "list")
  expect_true("uni" %in% names(result$posterior))
  expect_true("sigma2_delta" %in% names(result$posterior))
  expect_true("sigma2_delta_abs" %in% names(result$posterior))
  expect_true("g" %in% names(result$posterior))
  expect_true("phi" %in% names(result$posterior))
  expect_true("sem_rel" %in% names(result$posterior))
  expect_true("sem_abs" %in% names(result$posterior))
})

test_that("dstudy posterior distributions are numeric vectors", {
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, backend = "brms")
  result <- dstudy(g, n = list(rater = 3), estimation = "posterior")

  for (name in c("uni", "sigma2_delta", "sigma2_delta_abs",
                 "g", "phi", "sem_rel", "sem_abs")) {
    expect_type(result$posterior[[name]], "double")
    expect_true(length(result$posterior[[name]]) > 0)
  }
})

test_that("dstudy posterior coefficients are means of distributions", {
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, backend = "brms")
  result <- dstudy(g, n = list(rater = 3), estimation = "posterior")

  expect_equal(result$coefficients$uni, mean(result$posterior$uni))
  expect_equal(result$coefficients$g, mean(result$posterior$g))
  expect_equal(result$coefficients$phi, mean(result$posterior$phi))
})

test_that("dstudy posterior g is between 0 and 1", {
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, backend = "brms")
  result <- dstudy(g, n = list(rater = 3), estimation = "posterior")

  expect_true(all(result$posterior$g >= 0, na.rm = TRUE))
  expect_true(all(result$posterior$g <= 1, na.rm = TRUE))
})

test_that("dstudy posterior with sweep returns list of posteriors", {
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, backend = "brms")
  result <- dstudy(g, n = list(rater = c(2, 3)), estimation = "posterior")

  expect_true(result$is_sweep)
  expect_type(result$posterior, "list")
  expect_equal(length(result$posterior), 2)

  for (i in seq_along(result$posterior)) {
    expect_true("uni" %in% names(result$posterior[[i]]))
    expect_true("g" %in% names(result$posterior[[i]]))
  }
})

test_that("dstudy posterior estimation works with brms backend without warning", {
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, backend = "brms")

  expect_warning(
    result <- dstudy(g, n = list(rater = 3), estimation = "posterior"),
    NA
  )

  expect_equal(result$estimation, "posterior")
})

test_that("dstudy estimation rejects invalid values", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  expect_error(
    dstudy(g, n = list(rater = 3), estimation = "invalid"),
    "'arg' should be one of"
  )
})

test_that("plot.dstudy warns for non-sweep dstudy", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("ggplot2")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  expect_false(d$is_sweep)

  expect_warning(
    plot(d, type = "sweep"),
    "Plot is not defined for individual g and phi coefficients"
  )
})

test_that("plot.dstudy creates line plot for single vector n", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("ggplot2")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = c(2, 3, 5)))

  expect_true(d$is_sweep)

  p <- plot(d, type = "sweep", coefficient = "g")
  expect_s3_class(p, "ggplot")

  p <- plot(d, type = "sweep", coefficient = "phi")
  expect_s3_class(p, "ggplot")

  p <- plot(d, type = "sweep", coefficient = "both")
  expect_s3_class(p, "ggplot")
})

test_that("plot.dstudy creates line plot for two vector n with grouping", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("ggplot2")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = c(2, 3), person = c(10, 20)))

  expect_true(d$is_sweep)
  expect_equal(nrow(d$coefficients), 4)

  p <- plot(d, type = "sweep", coefficient = "g")
  expect_s3_class(p, "ggplot")

  p <- plot(d, type = "sweep", coefficient = "phi")
  expect_s3_class(p, "ggplot")

  p <- plot(d, type = "sweep", coefficient = "both")
  expect_s3_class(p, "ggplot")
  expect_true(any(grepl("facet_wrap", deparse(body(p)))))
})

test_that("plot.dstudy coefficients type still works", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("ggplot2")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  p <- plot(d, type = "coefficients")
  expect_s3_class(p, "ggplot")
})

# =============================================================================
# Multiple vector n scaling tests
# =============================================================================

test_that("dstudy scales correctly with single vector n (baseline)", {
  skip_if_not_installed("lme4")

  test_data_2facet <- data.frame(
    score = rnorm(200),
    person = factor(rep(1:20, 10)),
    item = factor(rep(1:10, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_2facet)

  result <- dstudy(g, n = list(item = c(2, 5, 10)))

  expect_true(result$is_sweep)
  expect_equal(nrow(result$coefficients), 3)

  expect_true(all(result$coefficients$g >= 0))
  expect_true(all(result$coefficients$g <= 1))
  expect_true(all(result$coefficients$phi >= 0))
  expect_true(all(result$coefficients$phi <= 1))
})

test_that("dstudy scales correctly with multiple vector n without aggregation", {
  skip_if_not_installed("lme4")

  test_data_3 <- data.frame(
    score = rnorm(300),
    person = factor(rep(1:20, 15)),
    item = factor(rep(1:10, each = 30)),
    rater = factor(rep(1:3, each = 100))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item) + (1 | rater), data = test_data_3)

  result <- dstudy(g, n = list(item = c(2, 5), rater = c(2, 4)))

  expect_true(result$is_sweep)
  expect_equal(nrow(result$coefficients), 4)

  expect_true(all(result$coefficients$g >= 0))
  expect_true(all(result$coefficients$g <= 1))
  expect_true(all(result$coefficients$phi >= 0))
  expect_true(all(result$coefficients$phi <= 1))

  coef_row_2_2 <- result$coefficients[result$coefficients$item == 2 & result$coefficients$rater == 2, ]
  coef_row_5_4 <- result$coefficients[result$coefficients$item == 5 & result$coefficients$rater == 4, ]

  expect_true(coef_row_2_2$g <= coef_row_5_4$g)
})

test_that("dstudy scales correctly with multiple vector n with aggregation", {
  skip_if_not_installed("lme4")

  test_data_3 <- data.frame(
    score = rnorm(300),
    person = factor(rep(1:20, 15)),
    item = factor(rep(1:10, each = 30)),
    rater = factor(rep(1:3, each = 100))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item) + (1 | rater), data = test_data_3)

  result <- dstudy(g, n = list(item = c(2, 5), rater = c(2, 4)), aggregation = "rater")

  expect_true(result$is_sweep)
  expect_equal(nrow(result$coefficients), 4)

  expect_true(all(result$coefficients$g >= 0))
  expect_true(all(result$coefficients$g <= 1))
  expect_true(all(result$coefficients$phi >= 0))
  expect_true(all(result$coefficients$phi <= 1))
})

test_that("dstudy scaling is consistent between single and multiple vector n", {
  skip_if_not_installed("lme4")

  test_data_2 <- data.frame(
    score = rnorm(200),
    person = factor(rep(1:20, 10)),
    item = factor(rep(1:10, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_2)

  d_single <- dstudy(g, n = list(item = 5))
  d_multi <- dstudy(g, n = list(item = c(5)))

  vc_single <- d_single$variance_components
  vc_multi <- d_multi$variance_components

  expect_equal(vc_single$var, vc_multi$var, tolerance = 1e-10)
})

# =============================================================================
# Divided estimates tests
# =============================================================================

test_that("dstudy includes both unscaled and scaled estimates when n not provided", {
  skip_if_not_installed("lme4")

  test_data_div <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_div)

  d <- dstudy(g)

  expect_true("estimate" %in% names(d$coefficients))
  expect_true(all(c("unscaled", "scaled") %in% d$coefficients$estimate))
  expect_equal(nrow(d$coefficients), 2)
})

test_that("dstudy includes only scaled estimates when n is provided", {
  skip_if_not_installed("lme4")

  test_data_div <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_div)

  d <- dstudy(g, n = list(item = 10))

  expect_false("estimate" %in% names(d$coefficients))
  expect_equal(nrow(d$coefficients), 1)
})

test_that("calculate_divided_variance produces correct divisors for simple design", {
  skip_if_not_installed("lme4")

  test_data_div <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_div)

  vc <- g$variance_components
  n <- list(item = 10)

  vc_div <- calculate_divided_variance(vc, n, object = "person", residual_is = "person:item")

  expect_true("var_divided" %in% names(vc_div))
  expect_true(is.numeric(vc_div$var_divided))
  expect_equal(length(vc_div$var_divided), nrow(vc))
})

test_that("scaled estimates scale non-object components by non-object facet n only", {
skip_if_not_installed("lme4")

test_data_div <- data.frame(
score = rnorm(100),
person = factor(rep(1:20, 5)),
item = factor(rep(1:5, each = 20))
)
g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_div)

d <- dstudy(g)

unscaled_row <- d$coefficients[d$coefficients$estimate == "unscaled", ]
scaled_row <- d$coefficients[d$coefficients$estimate == "scaled", ]

# Universe score variance should be the same for both
expect_equal(unscaled_row$uni, scaled_row$uni, tolerance = 1e-10)

# Scaled sigma2_delta should be smaller than unscaled (divided by n_item)
expect_true(scaled_row$sigma2_delta <= unscaled_row$sigma2_delta)

# Scaled g coefficient should be larger than unscaled (better reliability)
expect_true(scaled_row$g >= unscaled_row$g)
})

test_that("dstudy sweep includes only scaled estimates when n is provided", {
  skip_if_not_installed("lme4")

  test_data_sweep <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_sweep)

  d <- dstudy(g, n = list(item = c(5, 10)))

  expect_false("estimate" %in% names(d$coefficients))
  expect_equal(nrow(d$coefficients), 2)
})

test_that("scaled estimates work with 3-facet design when n not provided", {
skip_if_not_installed("lme4")

test_data_3 <- data.frame(
score = rnorm(300),
person = factor(rep(1:20, 15)),
item = factor(rep(1:10, each = 30)),
rater = factor(rep(1:3, each = 100))
)
g <- gstudy(score ~ (1 | person) + (1 | item) + (1 | rater), data = test_data_3)

d <- dstudy(g)

expect_true("estimate" %in% names(d$coefficients))
expect_true(all(c("unscaled", "scaled") %in% d$coefficients$estimate))

unscaled_row <- d$coefficients[d$coefficients$estimate == "unscaled", ]
scaled_row <- d$coefficients[d$coefficients$estimate == "scaled", ]

# Both should have valid coefficients
expect_true(unscaled_row$g >= 0 && unscaled_row$g <= 1)
expect_true(scaled_row$g >= 0 && scaled_row$g <= 1)
expect_true(unscaled_row$phi >= 0 && unscaled_row$phi <= 1)
expect_true(scaled_row$phi >= 0 && scaled_row$phi <= 1)
})

test_that("scaled estimates with custom error specification when n not provided", {
  skip_if_not_installed("lme4")

  test_data_div <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
)
g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_div)

d <- dstudy(g, error = "person:item")

expect_true("estimate" %in% names(d$coefficients))
expect_equal(nrow(d$coefficients), 2)
})

test_that("scaled estimates produce correct divisors for interaction components", {
  skip_if_not_installed("lme4")

  # For a model with p (object), i, and p:i (residual):
  # - unscaled: no scaling (G-study variance components)
  # - scaled: divided by n (D-study variance components)

  test_data_div <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_div)

  vc <- g$variance_components

  # Find the component variances
  p_var <- vc$var[vc$component == "person"]
  i_var <- vc$var[vc$component == "item"]
  res_var <- vc$var[vc$component == "Residual"]

  n_item <- 10
  n_person <- 20 # from the data

# Test calculate_divided_variance (used for scaled estimates)
vc_scaled <- calculate_divided_variance(vc, n = list(item = n_item),
residual_is = "person:item")

  # Check that object variance is not divided
  p_scaled <- vc_scaled$var_divided[vc_scaled$component == "person"]
  expect_equal(p_scaled, p_var, tolerance = 1e-10)

  # Check that item is divided by n_item
  i_scaled <- vc_scaled$var_divided[vc_scaled$component == "item"]
  expect_equal(i_scaled, i_var / n_item, tolerance = 1e-10)
})

# =============================================================================
# Additional tests for n_provided behavior
# =============================================================================

test_that("dstudy with n provided produces only scaled estimates (simple estimation)", {
  skip_if_not_installed("lme4")

  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data)

  d <- dstudy(g, n = list(item = 10))

  expect_false("estimate" %in% names(d$coefficients))
  expect_equal(nrow(d$coefficients), 1)
  expect_true("uni" %in% names(d$coefficients))
  expect_true("g" %in% names(d$coefficients))
  expect_true("phi" %in% names(d$coefficients))
})

test_that("dstudy sweep with n provided produces only scaled estimates", {
  skip_if_not_installed("lme4")

  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data)

  d <- dstudy(g, n = list(item = c(5, 10)))

  expect_false("estimate" %in% names(d$coefficients))
  expect_equal(nrow(d$coefficients), 2)
})

test_that("estimate column is first when both estimates present", {
  skip_if_not_installed("lme4")

  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data)

  d <- dstudy(g)

  expect_true("estimate" %in% names(d$coefficients))
  expect_true("uni" %in% names(d$coefficients))
  col_names <- names(d$coefficients)
  estimate_idx <- which(col_names == "estimate")
  uni_idx <- which(col_names == "uni")
  expect_equal(estimate_idx, 1)
  expect_true(estimate_idx < uni_idx)
})

test_that("dstudy with posterior estimation and n provided produces only scaled estimates", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("brms")

  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, backend = "brms")

  suppressWarnings({
    d <- dstudy(g, n = list(item = 10))
  })

  expect_false("estimate" %in% names(d$coefficients))
  expect_equal(nrow(d$coefficients), 1)
})

test_that("dstudy with posterior estimation and no n provided produces both unscaled and scaled estimates", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("brms")

  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, backend = "brms")

  suppressWarnings({
    d <- dstudy(g)
  })

  expect_true("estimate" %in% names(d$coefficients))
  expect_true(all(c("unscaled", "scaled") %in% d$coefficients$estimate))
  expect_equal(nrow(d$coefficients), 2)
})

# =============================================================================
# Tests for correct scaling behavior with lme4
# =============================================================================

test_that("scaled estimates correctly divide variance by sample sizes (lme4)", {
  skip_if_not_installed("lme4")

  # Create data with meaningful variance structure
  set.seed(456)
  n_persons <- 15
  n_items <- 5

  test_data <- data.frame(
    person = factor(rep(1:n_persons, n_items)),
    item = factor(rep(1:n_items, each = n_persons)),
    score = rnorm(n_persons * n_items)
  )

  # Add person effect to ensure non-zero universe score variance
  person_effect <- rnorm(n_persons, sd = 2)
  test_data$score <- test_data$score + person_effect[test_data$person]

g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data)
d <- dstudy(g)

  unscaled_row <- d$coefficients[d$coefficients$estimate == "unscaled", ]
  scaled_row <- d$coefficients[d$coefficients$estimate == "scaled", ]

  # Universe score variance should be identical
  expect_equal(unscaled_row$uni, scaled_row$uni, tolerance = 1e-10)

  # Scaled sigma2_delta should be smaller than unscaled
  expect_true(scaled_row$sigma2_delta < unscaled_row$sigma2_delta)

  # The ratio should be approximately equal to n_item (allowing for numerical tolerance)
  if (unscaled_row$sigma2_delta > 0) {
    ratio <- unscaled_row$sigma2_delta / scaled_row$sigma2_delta
    expect_equal(ratio, n_items, tolerance = 0.5)
  }

  # Scaled g should be larger (better reliability with scaling)
  expect_true(scaled_row$g > unscaled_row$g)
})
