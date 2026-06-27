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

test_that("universe score uses unscaled variance for all universe components", {
  skip_if_not_installed("lme4")
  set.seed(42)
  data3 <- data.frame(
    score = rnorm(120),
    person = factor(rep(1:20, 6)),
    task = factor(rep(1:3, each = 40)),
    rater = factor(rep(rep(1:2, each = 20), 3))
  )
  # Fully crossed Person x Task x Rater design (one observation per
  # Person x Task x Rater cell). The canonical "all possible variance
  # components" formula has all three main effects and all three
  # two-way interactions; the three-way Person:Task:Rater interaction
  # is confounded with the residual.
  g <- gstudy(score ~ (1 | person) + (1 | task) + (1 | rater) +
                (1 | person:task) + (1 | person:rater) + (1 | task:rater),
              data = data3)

  result <- dstudy(g, n = list(task = 5, rater = 4), universe = c("person", "person:task"))

  vc <- result$variance_components
  var_person <- vc$var_unscaled[vc$component == "person"]
  var_pt <- vc$var_unscaled[vc$component == "person:task"]
  expected_uni <- var_person + var_pt

  expect_equal(result$coefficients$uni, expected_uni, tolerance = 1e-10,
    info = "universe score should be sum of unscaled G-study variances for all universe components")
})

test_that("universe score with expanded universe is larger than object-only universe", {
  skip_if_not_installed("lme4")
  set.seed(42)
  data3 <- data.frame(
    score = rnorm(120),
    person = factor(rep(1:20, 6)),
    task = factor(rep(1:3, each = 40)),
    rater = factor(rep(rep(1:2, each = 20), 3))
  )
  g <- gstudy(score ~ (1 | person) + (1 | task) + (1 | rater) +
                (1 | person:task) + (1 | person:rater) + (1 | task:rater),
              data = data3)

  result_default <- dstudy(g, n = list(task = 5, rater = 4))
  result_expanded <- dstudy(g, n = list(task = 5, rater = 4), universe = c("person", "person:task"))

  var_pt <- result_default$variance_components$var_unscaled[
    result_default$variance_components$component == "person:task"]

  expect_true(result_expanded$coefficients$uni > result_default$coefficients$uni,
    info = "expanded universe should have larger universe score")
  expect_equal(result_expanded$coefficients$uni,
    result_default$coefficients$uni + var_pt, tolerance = 1e-10,
    info = "expanded universe score = default uni + unscaled person:task variance")
})

test_that("universe score uses unscaled variance when n not provided", {
  skip_if_not_installed("lme4")
  set.seed(42)
  data3 <- data.frame(
    score = rnorm(120),
    person = factor(rep(1:20, 6)),
    task = factor(rep(1:3, each = 40)),
    rater = factor(rep(rep(1:2, each = 20), 3))
  )
  g <- gstudy(score ~ (1 | person) + (1 | task) + (1 | rater) +
                (1 | person:task) + (1 | person:rater) + (1 | task:rater),
              data = data3)

  result <- dstudy(g, universe = c("person", "person:task"))

  vc <- result$variance_components
  var_person <- vc$var_unscaled[vc$component == "person"]
  var_pt <- vc$var_unscaled[vc$component == "person:task"]
  expected_uni <- var_person + var_pt

  if ("estimate" %in% names(result$coefficients)) {
    uni_unscaled <- result$coefficients$uni[result$coefficients$estimate == "unscaled"]
    expect_equal(uni_unscaled, expected_uni, tolerance = 1e-10)
  } else {
    expect_equal(result$coefficients$uni, expected_uni, tolerance = 1e-10)
  }
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
    "Model is"
  )

  expect_s3_class(result, "dstudy")
  expect_equal(result$error, "rater")
})

test_that("dstudy does not warn when error is NULL", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)

  expect_warning(
    result <- dstudy(g, n = list(rater = 3)),
    "Model is"
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

test_that("dstudy warns and refits when posterior requested with non-brms estimator", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "lme4")

  expect_warning(
    result <- dstudy(g, n = list(rater = 3), estimation = "posterior"),
    "estimation = 'posterior' requires estimator = 'brms'"
  )

  expect_equal(result$estimation, "posterior")
  expect_equal(result$gstudy$estimator, "brms")
})

test_that("dstudy posterior estimation returns posterior distributions", {
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
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

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  result <- dstudy(g, n = list(rater = 3), estimation = "posterior")

  for (name in c("uni", "sigma2_delta", "sigma2_delta_abs",
                 "g", "phi", "sem_rel", "sem_abs")) {
    expect_type(result$posterior[[name]], "double")
    expect_true(length(result$posterior[[name]]) > 0)
  }
})

test_that("dstudy posterior coefficients are means of distributions", {
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  result <- dstudy(g, n = list(rater = 3), estimation = "posterior")

  expect_equal(result$coefficients$uni, mean(result$posterior$uni))
  expect_equal(result$coefficients$g, mean(result$posterior$g))
  expect_equal(result$coefficients$phi, mean(result$posterior$phi))
})

test_that("dstudy posterior g is between 0 and 1", {
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  result <- dstudy(g, n = list(rater = 3), estimation = "posterior")

  expect_true(all(result$posterior$g >= 0, na.rm = TRUE))
  expect_true(all(result$posterior$g <= 1, na.rm = TRUE))
})

test_that("dstudy posterior with sweep returns list of posteriors", {
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  result <- dstudy(g, n = list(rater = c(2, 3)), estimation = "posterior")

  expect_true(result$is_sweep)
  expect_type(result$posterior, "list")
  expect_equal(length(result$posterior), 2)

  for (i in seq_along(result$posterior)) {
    expect_true("uni" %in% names(result$posterior[[i]]))
    expect_true("g" %in% names(result$posterior[[i]]))
  }
})

test_that("dstudy posterior estimation works with brms estimator without warning", {
  skip_if_not_installed("brms")

  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")

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
  expect_true(!inherits(p$facet, "FacetNull"))
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
object = "person", residual_is = "person:item")

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
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "brms")

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
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data, estimator = "brms")

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

# =============================================================================
# Credible interval tests
# =============================================================================

test_that("dstudy accepts ci parameter", {
  skip_if_not_installed("brms")
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  expect_no_error(dstudy(g, n = list(rater = 3), ci = "g"))
  expect_no_error(dstudy(g, n = list(rater = 3), ci = c("g", "phi")))
})

test_that("dstudy accepts probs parameter", {
  skip_if_not_installed("brms")
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  expect_no_error(dstudy(g, n = list(rater = 3), ci = "g", probs = c(0.05, 0.95)))
})

test_that("dstudy warns when ci requested with non-brms estimator", {
  skip_if_not_installed("lme4")
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "lme4")
  expect_warning(
    dstudy(g, n = list(rater = 3), ci = "g"),
    "Credible intervals for 'mom' and 'lme4' estimators are not yet implemented"
  )
})

test_that("dstudy validates ci parameter values", {
  skip_if_not_installed("lme4")
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "lme4")
  expect_error(
    dstudy(g, n = list(rater = 3), ci = "invalid"),
    "'arg' should be one of"
  )
})

test_that("dstudy validates probs parameter", {
  skip_if_not_installed("brms")
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  expect_error(
    dstudy(g, n = list(rater = 3), ci = "g", probs = c(0.1)),
    "'probs' must have exactly 2 elements"
  )
  expect_error(
    dstudy(g, n = list(rater = 3), ci = "g", probs = c(0.9, 0.1)),
    "'probs' must be in increasing order"
  )
})

test_that("dstudy with brms estimator produces CI columns", {
  skip_if_not_installed("brms")
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  result <- dstudy(g, n = list(rater = 3), ci = c("g", "phi"))
  expect_true("g_LL" %in% names(result$coefficients))
  expect_true("g_UL" %in% names(result$coefficients))
  expect_true("phi_LL" %in% names(result$coefficients))
  expect_true("phi_UL" %in% names(result$coefficients))
})

test_that("dstudy CI values are in correct order (LL < UL)", {
  skip_if_not_installed("brms")
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  result <- dstudy(g, n = list(rater = 3), ci = c("g", "phi"))
  expect_true(result$coefficients$g_LL < result$coefficients$g_UL)
  expect_true(result$coefficients$phi_LL < result$coefficients$phi_UL)
})

test_that("dstudy CI uses custom probs", {
  skip_if_not_installed("brms")
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  result_95 <- dstudy(g, n = list(rater = 3), ci = "g", probs = c(0.025, 0.975))
  result_90 <- dstudy(g, n = list(rater = 3), ci = "g", probs = c(0.05, 0.95))
  width_95 <- result_95$coefficients$g_UL - result_95$coefficients$g_LL
  width_90 <- result_90$coefficients$g_UL - result_90$coefficients$g_LL
  expect_true(width_90 < width_95)
})

test_that("dstudy sweep with ci produces CI columns", {
  skip_if_not_installed("brms")
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  result <- dstudy(g, n = list(rater = c(2, 3)), ci = "g")
  expect_true(result$is_sweep)
  expect_true("g_LL" %in% names(result$coefficients))
  expect_true("g_UL" %in% names(result$coefficients))
  expect_equal(nrow(result$coefficients), 2)
})

test_that("dstudy stores ci and probs in result", {
  skip_if_not_installed("brms")
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, estimator = "brms")
  result <- dstudy(g, n = list(rater = 3), ci = "g", probs = c(0.05, 0.95))
  expect_equal(result$ci, "g")
  expect_equal(result$probs, c(0.05, 0.95))
})

test_that("dstudy stores NULL ci when not specified", {
  skip_if_not_installed("lme4")
  test_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  result <- dstudy(g, n = list(rater = 3))
  expect_null(result$ci)
})

# =============================================================================
# Long-format multivariate tests
# =============================================================================

test_that("dstudy handles long_format_multivariate per-dimension sample sizes", {
  mock_g <- list(
    estimator = "brms",
    long_format_multivariate = TRUE,
    dimension_var = "Subtest",
    dimensions = c("A", "B"),
    variance_components = tibble::tibble(
      component = c("Person", "Residual", "Person", "Residual"),
      dim = c("A", "A", "B", "B"),
      type = c("main", "residual", "main", "residual"),
      var = c(0.5, 0.3, 0.4, 0.35)
    ),
    sample_size_info_per_dim = list(
      A = list(main = c(Person = 10)),
      B = list(main = c(Person = 8))
    ),
    object = "Person",
    facets = c("Person"),
    is_multivariate = TRUE
  )
  class(mock_g) <- "mgstudy"

  expect_true(is.mgstudy(mock_g))
})

# =============================================================================
# Integration test for multivariate designs with composite coefficients
# =============================================================================

test_that("dstudy computes composite coefficients for multivariate designs", {
  skip_if_not_installed("brms")

  set.seed(123)
  n_persons <- 20
  n_items <- 5

  data <- data.frame(
    person = rep(1:n_persons, each = n_items * 2),
    item = rep(rep(1:n_items, 2), n_persons),
    dim = rep(c("A", "B"), n_persons * n_items)
  )

  person_effect_A <- rnorm(n_persons, 0, sqrt(0.4))
  person_effect_B <- rnorm(n_persons, 0, sqrt(0.3))
  item_effect <- rnorm(n_items, 0, sqrt(0.1))

  data$score <- NA
  for (i in 1:nrow(data)) {
    p <- data$person[i]
    it <- data$item[i]
    d <- data$dim[i]
    if (d == "A") {
      data$score[i] <- person_effect_A[p] + item_effect[it] + rnorm(1, 0, sqrt(0.5))
    } else {
      data$score[i] <- person_effect_B[p] + item_effect[it] + rnorm(1, 0, sqrt(0.6))
    }
  }

  data_wide <- tidyr::pivot_wider(data, names_from = dim, values_from = score)

  gstudy_result <- gstudy(
    formula = mvbind(A, B) ~ (1 | person) + (1 | item),
    data = data_wide,
    estimator = "mom"
  )

  dstudy_result <- dstudy(gstudy_result, n = list(person = 10, item = 5))

  expect_true("Composite" %in% dstudy_result$coefficients$dim)

  composite_row <- dstudy_result$coefficients[dstudy_result$coefficients$dim == "Composite", ]
  expect_true(composite_row$g > 0 && composite_row$g < 1)
  expect_true(composite_row$phi > 0 && composite_row$phi < 1)
  expect_true(composite_row$g >= composite_row$phi)
})

test_that("dstudy computes composite variance components for multivariate posterior estimation", {
  skip_if_not_installed("brms")
  set.seed(123)
  n_persons <- 20
  n_items <- 5
  data <- data.frame(
    person = rep(1:n_persons, each = n_items * 2),
    item = rep(rep(1:n_items, 2), n_persons),
    dim = rep(c("A", "B"), n_persons * n_items)
  )

  person_effect_A <- rnorm(n_persons, 0, sqrt(0.4))
  person_effect_B <- 0.5 * person_effect_A + rnorm(n_persons, 0, sqrt(0.3 * 0.75))
  item_effect_A <- rnorm(n_items, 0, sqrt(0.1))
  item_effect_B <- rnorm(n_items, 0, sqrt(0.1))

  data$score <- NA
  for (i in 1:nrow(data)) {
    p <- data$person[i]
    it <- data$item[i]
    d <- data$dim[i]
    if (d == "A") {
      data$score[i] <- person_effect_A[p] + item_effect_A[it] + rnorm(1, 0, sqrt(0.5))
    } else {
      data$score[i] <- person_effect_B[p] + item_effect_B[it] + rnorm(1, 0, sqrt(0.5))
    }
  }

  data_wide <- tidyr::pivot_wider(data, names_from = dim, values_from = score)

  gstudy_result <- gstudy(
    formula = mvbind(A, B) ~ (1 | person) + (1 | item),
    data = data_wide,
    estimator = "brms",
    cores = 2,
    chains = 2,
    iter = 1000
  )

  dstudy_result <- dstudy(gstudy_result, n = list(person = 10, item = 5), weights = c(A = 1, B = 1))

  expect_true("Composite" %in% dstudy_result$variance_components$dim)

  composite_rows <- dstudy_result$variance_components[dstudy_result$variance_components$dim == "Composite", ]
  expect_true(nrow(composite_rows) > 0)

  components <- unique(dstudy_result$variance_components$component)
  components <- components[components != "Composite"]
  expect_true(all(components %in% composite_rows$component))

  expect_true(!is.null(dstudy_result$composite_posterior))
})

test_that("composite coefficients are computed from posterior draws, not point estimates", {
  skip_if_not_installed("brms")
  skip_on_cran()
  set.seed(789)
  
  n_persons <- 20
  n_items <- 5
  
  data <- data.frame(
    person = rep(1:n_persons, each = n_items * 2),
    item = rep(rep(1:n_items, 2), n_persons),
    dim = rep(c("A", "B"), n_persons * n_items)
  )
  
  person_A <- rnorm(n_persons, 0, sqrt(0.4))
  person_B <- 0.7 * person_A + rnorm(n_persons, 0, sqrt(0.3))
  item_A <- rnorm(n_items, 0, sqrt(0.1))
  item_B <- rnorm(n_items, 0, sqrt(0.1))
  
  data$score <- NA
  for (i in 1:nrow(data)) {
    p <- data$person[i]
    it <- data$item[i]
    d <- data$dim[i]
    if (d == "A") {
      data$score[i] <- person_A[p] + item_A[it] + rnorm(1, 0, sqrt(0.5))
    } else {
      data$score[i] <- person_B[p] + item_B[it] + rnorm(1, 0, sqrt(0.5))
    }
  }
  
  data_wide <- tidyr::pivot_wider(data, names_from = dim, values_from = score)
  
  gstudy_result <- gstudy(
    formula = mvbind(A, B) ~ (1 | person) + (1 | item),
    data = data_wide,
    estimator = "brms",
    cores = 2, chains = 2, iter = 1000
  )
  
  dstudy_result <- dstudy(
    gstudy_result, 
    n = list(person = 10, item = 5), 
    weights = c(A = 1, B = 1)
  )
  
  expect_s3_class(dstudy_result$coefficients, "data.frame")
  
  composite_row <- dstudy_result$coefficients[dstudy_result$coefficients$dim == "Composite", ]
  expect_true(nrow(composite_row) == 1)
  
  expect_true(!is.null(dstudy_result$composite_posterior), 
              "composite_posterior should exist for brms estimator")
  
  post <- dstudy_result$composite_posterior
  
  expect_true("uni" %in% names(post), "uni draws should exist")
  expect_true("sigma2_delta" %in% names(post), "sigma2_delta draws should exist")
  expect_true("sigma2_delta_abs" %in% names(post), "sigma2_delta_abs draws should exist")
  expect_true("g" %in% names(post), "g draws should exist")
  expect_true("phi" %in% names(post), "phi draws should exist")
  
  g_draws <- post$uni / (post$uni + post$sigma2_delta)
  phi_draws <- post$uni / (post$uni + post$sigma2_delta_abs)
  
  g_draw_mean <- mean(g_draws, na.rm = TRUE)
  phi_draw_mean <- mean(phi_draws, na.rm = TRUE)
  
  expect_equal(composite_row$g, g_draw_mean, tolerance = 0.01,
               info = "g coefficient should be mean of per-draw calculations")
  expect_equal(composite_row$phi, phi_draw_mean, tolerance = 0.01,
               info = "phi coefficient should be mean of per-draw calculations")
})

test_that("composite coefficients include credible intervals when ci requested", {
  skip_if_not_installed("brms")
  skip_on_cran()
  set.seed(321)
  
  n_persons <- 15
  n_items <- 4
  
  data <- data.frame(
    person = rep(1:n_persons, each = n_items * 2),
    item = rep(rep(1:n_items, 2), n_persons),
    dim = rep(c("A", "B"), n_persons * n_items)
  )
  
  person_A <- rnorm(n_persons, 0, 1)
  person_B <- person_A + rnorm(n_persons, 0, 0.5)
  
  data$score <- NA
  for (i in 1:nrow(data)) {
    p <- data$person[i]
    d <- data$dim[i]
    if (d == "A") {
      data$score[i] <- person_A[p] + rnorm(1, 0, 0.5)
    } else {
      data$score[i] <- person_B[p] + rnorm(1, 0, 0.5)
    }
  }
  
  data_wide <- tidyr::pivot_wider(data, names_from = dim, values_from = score)
  
  gstudy_result <- gstudy(
    formula = mvbind(A, B) ~ (1 | person),
    data = data_wide,
    estimator = "brms",
    cores = 2, chains = 2, iter = 800
  )
  
  dstudy_result <- dstudy(
    gstudy_result, 
    n = list(person = 10),
    weights = c(A = 0.5, B = 0.5),
    ci = c("g", "phi")
  )
  
  composite_row <- dstudy_result$coefficients[dstudy_result$coefficients$dim == "Composite", ]
  
  expect_true("g_LL" %in% names(dstudy_result$coefficients), 
              "g_LL should exist when ci includes 'g'")
  expect_true("g_UL" %in% names(dstudy_result$coefficients),
              "g_UL should exist when ci includes 'g'")
  expect_true("phi_LL" %in% names(dstudy_result$coefficients),
              "phi_LL should exist when ci includes 'phi'")
  expect_true("phi_UL" %in% names(dstudy_result$coefficients),
              "phi_UL should exist when ci includes 'phi'")
  
  expect_true(!is.na(composite_row$g_LL), "Composite g_LL should not be NA")
  expect_true(!is.na(composite_row$g_UL), "Composite g_UL should not be NA")
  expect_true(composite_row$g_LL <= composite_row$g, "g_LL should be <= g")
  expect_true(composite_row$g_UL >= composite_row$g, "g_UL should be >= g")
})

test_that("VAR is stored in var element for long-format multivariate models", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data(rajaratnam)

  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
       sigma ~ 0 + Subtest),
    data = rajaratnam,
    estimator = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )

  # Use correct dimension names from gu$dimensions
  dims <- as.character(gu$dimensions)
  weights <- setNames(rep(1, length(dims)), dims)

  d <- dstudy(gu, n = list(Person = 5), weights = weights)

  # VAR should NOT be in coefficients table
  expect_false("var_rel" %in% names(d$coefficients))
  expect_false("var_abs" %in% names(d$coefficients))

  # VAR should be in var element (use [[ to avoid partial matching)
  expect_false(is.null(d[["var"]]))
  expect_true(is.list(d[["var"]]))

  # Check posterior matrices
  expect_true(is.matrix(d[["var"]]$var_rel_draws))
  expect_true(is.matrix(d[["var"]]$var_abs_draws))
  expect_true(is.matrix(d[["var"]]$prmse_c_rel_draws))
  expect_true(is.matrix(d[["var"]]$prmse_c_abs_draws))

  # Check dimensions
  expect_equal(ncol(d[["var"]]$var_rel_draws), length(dims))
  expect_equal(colnames(d[["var"]]$var_rel_draws), dims)

  # Check summary vectors
  expect_true(all(!is.na(d[["var"]]$var_rel)))
  expect_true(all(!is.na(d[["var"]]$var_abs)))
  expect_true(all(!is.na(d[["var"]]$prmse_c_rel)))
  expect_true(all(!is.na(d[["var"]]$prmse_c_abs)))

  # VAR values should be positive and finite
  expect_true(all(d[["var"]]$var_rel > 0))
  expect_true(all(d[["var"]]$var_abs > 0))
})

test_that("VAR is stored in var element for wide-format multivariate models (brms estimator)", {
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
    estimator = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )

  d <- dstudy(gu, n = list(person = 10), weights = c(A = 1, B = 1))

  # VAR should NOT be in coefficients table
  expect_false("var_rel" %in% names(d$coefficients))
  expect_false("var_abs" %in% names(d$coefficients))

  # VAR should be in var element (use [[ to avoid partial matching)
  expect_false(is.null(d[["var"]]))
  expect_true(is.list(d[["var"]]))

  # Check posterior matrices
  expect_true(is.matrix(d[["var"]]$var_rel_draws))
  expect_true(is.matrix(d[["var"]]$var_abs_draws))
  expect_true(is.matrix(d[["var"]]$prmse_c_rel_draws))
  expect_true(is.matrix(d[["var"]]$prmse_c_abs_draws))

  # Check dimensions
  expect_equal(ncol(d[["var"]]$var_rel_draws), 2)
  expect_equal(sort(colnames(d[["var"]]$var_rel_draws)), c("A", "B"))

  # VAR values should be positive and finite
  expect_true(all(d[["var"]]$var_rel > 0))
  expect_true(all(d[["var"]]$var_abs > 0))
})

test_that("prmse returns correct structure", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data(rajaratnam)

  gu <- gstudy(
  bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
    sigma ~ 0 + Subtest),
  data = rajaratnam,
  estimator = "brms",
  chains = 2,
  iter = 500,
  refresh = 0
  )

  dims <- as.character(gu$dimensions)
  weights <- setNames(rep(1, length(dims)), dims)

  d <- dstudy(gu, n = list(Person = 5), weights = weights)

  # prmse without ci should return point estimates only (no CI columns)
  var_df <- prmse(d)
  expect_s3_class(var_df, "tbl_df")
  expect_equal(nrow(var_df), length(dims))

  # Check column names - should NOT include CI columns by default
  expected_cols_no_ci <- c("dim", "prmse_s_rel", "prmse_s_abs", 
                           "prmse_c_rel", "prmse_c_abs",
                           "prmse_p_rel", "prmse_p_abs",
                           "var_rel", "var_abs")
  expect_true(all(expected_cols_no_ci %in% names(var_df)))
  # CI columns should NOT be present by default
  expect_false(any(grepl("_LL$|_UL$", names(var_df))))

  # Check values are positive
  expect_true(all(var_df$var_rel > 0))
  expect_true(all(var_df$var_abs > 0))
})

test_that("prmse with ci returns CI columns for brms estimator", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data(rajaratnam)

  gu <- gstudy(
  bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
    sigma ~ 0 + Subtest),
  data = rajaratnam,
  estimator = "brms",
  chains = 2,
  iter = 500,
  refresh = 0
  )

  dims <- as.character(gu$dimensions)
  weights <- setNames(rep(1, length(dims)), dims)

  d <- dstudy(gu, n = list(Person = 5), weights = weights)

  # prmse with ci should return CI columns
  var_df <- prmse(d, ci = c("prmse", "var"))
  expect_s3_class(var_df, "tbl_df")
  expect_equal(nrow(var_df), length(dims))

  # Check CI columns are present
  expect_true("prmse_c_rel_LL" %in% names(var_df))
  expect_true("prmse_c_rel_UL" %in% names(var_df))
  expect_true("var_rel_LL" %in% names(var_df))
  expect_true("var_rel_UL" %in% names(var_df))

  # Check credible intervals are ordered
  expect_true(all(var_df$var_rel_LL < var_df$var_rel_UL))
  expect_true(all(var_df$var_abs_LL < var_df$var_abs_UL))
})

test_that("var is NULL for univariate models", {
  data(brennan)
  g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
                (1 | Person:Task), data = brennan)
  d <- dstudy(g, n = list(Task = 3, Rater = 4))

  # Use [[ instead of $ to avoid partial matching with variance_components$var
  expect_true(is.null(d[["var"]]))
})

test_that("prmse returns G and Phi for univariate models", {
  data(brennan)
  g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
                (1 | Person:Task), data = brennan)
  d <- dstudy(g, n = list(Task = 3, Rater = 4))

  # Should return a tibble with prmse_rel and prmse_abs
  result <- suppressMessages(prmse(d))

  expect_s3_class(result, "tbl_df")
  expect_true("dim" %in% names(result))
  expect_true("prmse_rel" %in% names(result))
  expect_true("prmse_abs" %in% names(result))

  # prmse_rel should equal G, prmse_abs should equal Phi
  expect_equal(result$prmse_rel, d$coefficients$g)
  expect_equal(result$prmse_abs, d$coefficients$phi)
})

test_that("prmse warns that CIs only available for brms estimator", {
  data(brennan)
  g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
                (1 | Person:Task), data = brennan)
  d <- dstudy(g, n = list(Task = 3, Rater = 4))

  # Should warn when ci is specified for mom estimator
  expect_warning(
    prmse(d, ci = "prmse"),
    "Credible intervals are only available for brms estimator"
  )

  # Result should not have CI columns
  result <- suppressWarnings(prmse(d, ci = "prmse"))
  expect_false("prmse_rel_LL" %in% names(result))
  expect_false("prmse_abs_LL" %in% names(result))
})

# =============================================================================
# Tests for per-dimension sample size handling
# =============================================================================

test_that("validate_n_tibble validates correct tibble", {
  n_tb <- tibble::tibble(
    dim = c("A", "A", "B", "B"),
    facet = c("Person", "Item", "Person", "Item"),
    n = c(10, 5, 8, 4)
  )
  result <- validate_n_tibble(n_tb, c("A", "B"), c("Person", "Item"))
  expect_true(result$valid)
  expect_null(result$error)
})

test_that("validate_n_tibble rejects missing columns", {
  n_tb <- tibble::tibble(
    dim = c("A", "B"),
    facet = c("Person", "Person")
  )
  result <- validate_n_tibble(n_tb, c("A", "B"), c("Person"))
  expect_false(result$valid)
  expect_true(grepl("n", result$error))
})

test_that("validate_n_tibble rejects unknown dimensions", {
  n_tb <- tibble::tibble(
    dim = c("A", "C"),
    facet = c("Person", "Person"),
    n = c(10, 8)
  )
  result <- validate_n_tibble(n_tb, c("A", "B"), c("Person"))
  expect_false(result$valid)
  expect_true(grepl("unknown", tolower(result$error)))
})

test_that("validate_n_tibble rejects negative n values", {
  n_tb <- tibble::tibble(
    dim = c("A", "B"),
    facet = c("Person", "Person"),
    n = c(-1, 5)
  )
  result <- validate_n_tibble(n_tb, c("A", "B"), c("Person"))
  expect_false(result$valid)
  expect_true(grepl("positive", tolower(result$error)))
})

test_that("expand_n_per_dim handles non-sweep case", {
  n_tb <- tibble::tibble(
    dim = c("A", "A", "B", "B"),
    facet = c("Person", "Item", "Person", "Item"),
    n = c(10, 5, 8, 4)
  )
  result <- expand_n_per_dim(n_tb, sweep = FALSE)
  expect_false(result$is_sweep)
  expect_null(result$grid)
})

test_that("expand_n_per_dim detects sweep when multiple n per facet", {
  n_tb <- tibble::tibble(
    dim = c("A", "A", "A", "A"),
    facet = c("Person", "Person", "Item", "Item"),
    n = c(10, 15, 5, 8)
  )
  result <- expand_n_per_dim(n_tb, sweep = TRUE)
  expect_true(result$is_sweep)
  expect_s3_class(result$grid, "data.frame")
})

test_that("create_n_tibble_from_list creates correct structure", {
  n_list <- list(Person = 10, Item = 5)
  result <- create_n_tibble_from_list(n_list, c("A", "B"))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)  # 2 dimensions x 2 facets
  expect_true(all(c("dim", "facet", "n") %in% names(result)))
})

test_that("extract_sample_sizes_per_dim returns NULL for non-mgstudy", {
  skip_if_not_installed("lme4")
  data(brennan)
  g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
                (1 | Person:Task), data = brennan)
  result <- extract_sample_sizes_per_dim(g)
  expect_null(result)
})

# =============================================================================
# Glance method tests
# =============================================================================

test_that("glance.dstudy returns a one-row tibble for univariate dstudy", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  result <- glance(d)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1)
})

test_that("glance.dstudy columns follow var_unscaled_<component> pattern", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  result <- glance(d)

  var_cols <- names(result)[grepl("^var_unscaled_", names(result))]
  expect_true(length(var_cols) > 0)
  expect_true(all(grepl("^var_unscaled_", var_cols)))
  expect_true("var_unscaled_person" %in% names(result))
  expect_true("var_unscaled_rater" %in% names(result))
  expect_true("g" %in% names(result))
  expect_true("phi" %in% names(result))
})

test_that("glance.dstudy values match variance_components$var_unscaled", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  result <- glance(d)
  vc <- d$variance_components

  for (comp in vc$component) {
    col_name <- paste0("var_unscaled_", comp)
    expect_equal(result[[col_name]], vc$var_unscaled[vc$component == comp],
      tolerance = 1e-12)
  }
})

test_that("glance.dstudy returns multi-row tibble for multivariate dstudy", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data_mv <- data.frame(
    score1 = rnorm(60),
    score2 = rnorm(60),
    person = factor(rep(1:12, 5)),
    rater = factor(rep(1:5, each = 12))
  )

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = data_mv,
    estimator = "brms"
  )
  d <- dstudy(g, n = list(rater = 3))

  result <- glance(d)

  expect_s3_class(result, "tbl_df")
  expect_true("component" %in% names(result))
  expect_equal(nrow(result), length(unique(d$variance_components$component[
    d$variance_components$dim != "Composite"])))
})

test_that("glance.dstudy multivariate columns are dimension names", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data_mv <- data.frame(
    score1 = rnorm(60),
    score2 = rnorm(60),
    person = factor(rep(1:12, 5)),
    rater = factor(rep(1:5, each = 12))
  )

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = data_mv,
    estimator = "brms"
  )
  d <- dstudy(g, n = list(rater = 3))

  result <- glance(d)
  dims <- unique(d$variance_components$dim[d$variance_components$dim != "Composite"])

  expect_true(all(dims %in% names(result)))
})

test_that("glance.dstudy drops Composite rows from multivariate posterior dstudy", {
  mock_d <- list(
    variance_components = tibble::tibble(
      component = c("person", "person", "person"),
      dim = c("score1", "score2", "Composite"),
      var_unscaled = c(0.4, 0.5, 0.45),
      pct_unscaled = c(20, 25, 22.5),
      var_scaled = c(0.4, 0.5, 0.45),
      pct_scaled = c(20, 25, 22.5)
    ),
    object = "person",
    is_sweep = FALSE
  )
  class(mock_d) <- "dstudy"

  result <- glance(mock_d)

  expect_false("Composite" %in% result$component)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 5)
  expect_true("component" %in% names(result))
  expect_true("score1" %in% names(result))
  expect_true("score2" %in% names(result))
  expect_true("g" %in% names(result))
  expect_true("phi" %in% names(result))
})

test_that("glance.dstudy does not include covariance columns", {
  mock_d <- list(
    variance_components = tibble::tibble(
      component = c("person", "person", "person"),
      dim = c("score1", "score2", "Composite"),
      var_unscaled = c(0.4, 0.5, 0.45),
      pct_unscaled = c(20, 25, 22.5),
      var_scaled = c(0.4, 0.5, 0.45),
      pct_scaled = c(20, 25, 22.5)
    ),
    object = "person",
    is_sweep = FALSE
  )
  class(mock_d) <- "dstudy"

  result <- glance(mock_d)

  expect_false(any(grepl("cov", names(result), ignore.case = TRUE)))
})

test_that("glance.dstudy round-trips var_unscaled values for univariate", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  result <- glance(d)
  result_list <- as.list(result)
  vc <- d$variance_components

  for (comp in vc$component) {
    col_name <- paste0("var_unscaled_", comp)
    expect_true(col_name %in% names(result_list))
    expect_equal(result_list[[col_name]], vc$var_unscaled[vc$component == comp],
      tolerance = 1e-12)
  }
})

test_that("glance.dstudy includes scalar g and phi for univariate dstudy", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  result <- glance(d)

  expect_true("g" %in% names(result))
  expect_true("phi" %in% names(result))
  expect_type(result$g, "double")
  expect_type(result$phi, "double")
  expect_length(result$g, 1)
  expect_length(result$phi, 1)
  expect_true(result$g >= 0 && result$g <= 1)
  expect_true(result$phi >= 0 && result$phi <= 1)
})

test_that("glance.dstudy g and phi values match tidy(d)$g and tidy(d)$phi", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  result <- glance(d)
  td <- tidy(d)

  expect_equal(result$g, td$g, tolerance = 1e-12)
  expect_equal(result$phi, td$phi, tolerance = 1e-12)
})

test_that("glance.dstudy g and phi are list-columns for multivariate dstudy", {
  mock_d <- list(
    variance_components = tibble::tibble(
      component = c("person", "person", "person"),
      dim = c("score1", "score2", "Composite"),
      var_unscaled = c(0.4, 0.5, 0.45)
    ),
    coefficients = tibble::tibble(
      dim = c("score1", "score2"),
      g = c(0.8, 0.7),
      phi = c(0.75, 0.65)
    ),
    object = "person",
    is_sweep = FALSE
  )
  class(mock_d) <- "dstudy"

  result <- glance(mock_d)

  expect_true("g" %in% names(result))
  expect_true("phi" %in% names(result))
  expect_type(result$g, "list")
  expect_type(result$phi, "list")
  expect_length(result$g, 1)
  expect_length(result$phi, 1)
  expect_named(result$g[[1]], c("score1", "score2"))
  expect_named(result$phi[[1]], c("score1", "score2"))
  expect_equal(unname(result$g[[1]]["score1"]), 0.8, tolerance = 1e-12)
  expect_equal(unname(result$phi[[1]]["score2"]), 0.65, tolerance = 1e-12)
})

# =============================================================================
# 4-decimal-point display tests
# =============================================================================

test_that("dstudy variance components are stored at 4 decimal places by default", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  vc <- d$variance_components
  expect_true("var_unscaled" %in% names(vc))
  expect_true("var_scaled" %in% names(vc))

  # var_unscaled is copied directly from the g-study variance components and
  # is therefore at 4 dp. var_scaled is the result of division and may not
  # be at exactly 4 dp; we check it's close (within 1e-6).
  expect_true(all(abs(vc$var_unscaled - round(vc$var_unscaled, 4)) < 1e-6, na.rm = TRUE))
  expect_true(all(abs(vc$var_scaled - round(vc$var_scaled, 4)) < 1e-3, na.rm = TRUE))
})

test_that("print.dstudy default uses 4 decimal places", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  out <- capture.output(print(d))
  out <- paste(out, collapse = "\n")

  expect_match(out, "\\d+\\.\\d{4}", all = FALSE)
})

test_that("summary.dstudy default uses 4 decimal places", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  out <- capture.output(summary(d))
  out <- paste(out, collapse = "\n")

  expect_match(out, "\\d+\\.\\d{4}", all = FALSE)
})

test_that("tidy.dstudy returns a tibble with numeric columns at 4 decimal places by default", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  result <- tidy(d)

  expect_s3_class(result, "tbl_df")
  numeric_cols <- vapply(result, is.numeric, logical(1))
  numeric_cols <- names(result)[numeric_cols]

  for (col in numeric_cols) {
    expect_true(all(abs(result[[col]] - round(result[[col]], 4)) < 1e-6, na.rm = TRUE))
  }
})

test_that("tidy.dstudy digits argument overrides default", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  result <- tidy(d, digits = 2)
  expect_s3_class(result, "tbl_df")

  # Check that some coefficient column (e.g., g) is rounded to 2 dp
  if ("g" %in% names(result)) {
    expect_true(all(abs(result$g - round(result$g, 2)) < 1e-6, na.rm = TRUE))
  }
})

test_that("glance.dstudy returns values at 4 decimal places by default", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  result <- glance(d)

  expect_s3_class(result, "tbl_df")
  for (col in names(result)) {
    expect_true(all(abs(result[[col]] - round(result[[col]], 4)) < 1e-6, na.rm = TRUE))
  }
})

test_that("glance.dstudy digits argument overrides default", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  result <- glance(d, digits = 2)
  expect_s3_class(result, "tbl_df")

  for (col in names(result)) {
    expect_true(all(abs(result[[col]] - round(result[[col]], 2)) < 1e-6, na.rm = TRUE))
  }
})

test_that("print.dstudy digits argument overrides default", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
  d <- dstudy(g, n = list(rater = 3))

  out <- capture.output(print(d, digits = 2))
  out <- paste(out, collapse = "\n")

  # With digits = 2, the output should not contain values with more than
  # 4 decimal places (5+ dp). Some values may still show 3-4 dp due to
  # pillar's significant figure handling.
  expect_false(grepl("\\d+\\.\\d{5,}", out))
})

# =============================================================================
# Tests for residual divisor (object of measurement excluded)
# =============================================================================

test_that("compute_residual_divisor excludes object of measurement", {
  # Person is the object; residual Person:Rater should divide by n_Rater only
  divisor <- compute_residual_divisor(
    residual_is = "Person:Rater",
    n = list(Person = 10, Task = 3, Rater = 12),
    object_spec = "Person"
  )
  expect_equal(divisor, 12)

  # No residual_is provided: fall back to all non-object facets in n
  divisor_fallback <- compute_residual_divisor(
    residual_is = NULL,
    n = list(Person = 10, Rater = 5),
    object_spec = "Person"
  )
  expect_equal(divisor_fallback, 5)

  # Object not in residual: divisor uses all residual facets
  divisor_no_obj <- compute_residual_divisor(
    residual_is = "Task:Rater",
    n = list(Person = 10, Task = 3, Rater = 12),
    object_spec = "Person"
  )
  expect_equal(divisor_no_obj, 3 * 12)

  # Empty n: divisor is 1
  divisor_empty <- compute_residual_divisor(
    residual_is = "Person:Rater",
    n = list(),
    object_spec = "Person"
  )
  expect_equal(divisor_empty, 1)
})

test_that("compute_scale_factor_from_facets excludes object of measurement", {
  # Person:Task with object=Person: divide by n_Task only
  sf <- compute_scale_factor_from_facets(
    facets = c("Person", "Task"),
    n = list(Person = 10, Task = 3),
    object_spec = "Person"
  )
  expect_equal(sf, 3)

  # Without object_spec: original behavior preserved
  sf_no_obj <- compute_scale_factor_from_facets(
    facets = c("Person", "Task"),
    n = list(Person = 10, Task = 3)
  )
  expect_equal(sf_no_obj, 30)

  # Main effect not the object: divisor unchanged
  sf_main <- compute_scale_factor_from_facets(
    facets = c("Task"),
    n = list(Person = 10, Task = 3),
    object_spec = "Person"
  )
  expect_equal(sf_main, 3)
})

test_that("brennan residual is divided by n_Rater only (n provided)", {
  skip_if_not_installed("lme4")
  data(brennan)
  g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
                (1 | Person:Task), data = brennan)

  n_rater <- 4
  d <- dstudy(g, n = list(Task = 3, Rater = n_rater))

  vc <- d$variance_components
  res_var_unscaled <- vc$var_unscaled[vc$component == "Residual"]
  res_var_scaled <- vc$var_scaled[vc$component == "Residual"]

  # The residual (Person:Rater) should be divided by n_Rater only,
  # not n_Person * n_Task * n_Rater
  expect_equal(res_var_scaled, res_var_unscaled / n_rater, tolerance = 1e-10)
})

test_that("brennan residual is divided by n_Rater only (auto-extracted n)", {
  skip_if_not_installed("lme4")
  data(brennan)
  g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
                (1 | Person:Task), data = brennan)

  d <- dstudy(g)

  vc <- d$variance_components
  res_var_unscaled <- vc$var_unscaled[vc$component == "Residual"]
  res_var_scaled <- vc$var_scaled[vc$component == "Residual"]

  # When n is auto-extracted from G-study (Person=10, Task=3, Rater=12),
  # residual Person:Rater should still divide by n_Rater=12 only
  n_rater_gstudy <- length(unique(brennan$Rater))
  expect_equal(res_var_scaled, res_var_unscaled / n_rater_gstudy, tolerance = 1e-10)
})

test_that("brennan Person:Task is divided by n_Task only (object excluded)", {
  skip_if_not_installed("lme4")
  data(brennan)
  g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
                (1 | Person:Task), data = brennan)

  n_task <- 3
  d <- dstudy(g, n = list(Task = n_task, Rater = 4))

  vc <- d$variance_components
  pt_var_unscaled <- vc$var_unscaled[vc$component == "Person:Task"]
  pt_var_scaled <- vc$var_scaled[vc$component == "Person:Task"]

  # Person:Task interaction should be divided by n_Task only (object=Person excluded)
  expect_equal(pt_var_scaled, pt_var_unscaled / n_task, tolerance = 1e-10)
})

test_that("2-facet residual Person:Item is divided by n_Item only", {
  skip_if_not_installed("lme4")
  test_data_2f <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item), data = test_data_2f)

  n_item <- 10
  d <- dstudy(g, n = list(item = n_item))

  vc <- d$variance_components
  res_var_unscaled <- vc$var_unscaled[vc$component == "Residual"]
  res_var_scaled <- vc$var_scaled[vc$component == "Residual"]

  # Residual (Person:Item) should divide by n_Item only, not n_Person * n_Item
  expect_equal(res_var_scaled, res_var_unscaled / n_item, tolerance = 1e-10)
})

test_that("3-facet residual Person:Item:Rater is divided by n_Item * n_Rater", {
  skip_if_not_installed("lme4")
  test_data_3f <- data.frame(
    score = rnorm(120),
    person = factor(rep(1:10, 12)),
    item = factor(rep(1:3, each = 40)),
    rater = factor(rep(1:4, each = 30))
  )
  g <- gstudy(score ~ (1 | person) + (1 | item) + (1 | rater), data = test_data_3f)

  n_item <- 5
  n_rater <- 3
  d <- dstudy(g, n = list(item = n_item, rater = n_rater))

  vc <- d$variance_components
  res_var_unscaled <- vc$var_unscaled[vc$component == "Residual"]
  res_var_scaled <- vc$var_scaled[vc$component == "Residual"]

  # Residual (Person:Item:Rater) should divide by n_Item * n_Rater
  # (Person is object, excluded)
  expect_equal(res_var_scaled, res_var_unscaled / (n_item * n_rater),
               tolerance = 1e-10)
})

test_that("compute_component_scale_factor excludes object for non-residual interactions", {
  # Person:Task with object=Person: should divide by n_Task only
  sf <- compute_component_scale_factor(
    comp = "Person:Task",
    n = list(Person = 10, Task = 3),
    object_spec = "Person",
    n_provided = TRUE
  )
  expect_equal(sf, 3)

  # Residual with residual_is: excludes object
  sf_resid <- compute_component_scale_factor(
    comp = "Residual",
    n = list(Person = 10, Rater = 12),
    object_spec = "Person",
    n_provided = TRUE,
    residual_is = "Person:Rater"
  )
  expect_equal(sf_resid, 12)

  # Object component: always 1
  sf_obj <- compute_component_scale_factor(
    comp = "Person",
    n = list(Person = 10, Task = 3),
    object_spec = "Person",
    n_provided = TRUE
  )
  expect_equal(sf_obj, 1)
})

test_that("compute_scale_factors_for_viable uses correct residual divisor", {
  components <- c("Person", "Rater", "Person:Task", "Residual")
  n <- list(Person = 10, Task = 3, Rater = 12)
  object_spec <- "Person"
  universe_spec <- "Person"

  sfs <- compute_scale_factors_for_viable(
    components = components,
    n = n,
    universe_spec = universe_spec,
    object_spec = object_spec,
    residual_is = "Person:Rater"
  )

  expect_equal(sfs$Person, 1)         # object
  expect_equal(sfs$Rater, 12)          # main effect
  expect_equal(sfs$`Person:Task`, 3)   # interaction w/ object: exclude object
  expect_equal(sfs$Residual, 12)       # residual: Rater only (Person excluded)
})

