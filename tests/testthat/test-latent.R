# Tests for latent() function

# Test data for reuse
test_data_latent <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, 5)),
  rater = factor(rep(1:5, each = 20)),
  item = factor(rep(1:10, each = 10))
)

# =============================================================================
# Basic latent tests
# =============================================================================

test_that("latent returns a data frame", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data_latent)
  result <- latent(g)
  expect_s3_class(result, "data.frame")
})

test_that("latent validates input type", {
  expect_error(latent("not a gstudy"), "must be a gstudy or mgstudy object")
  expect_error(latent(NULL), "must be a gstudy or mgstudy object")
  expect_error(latent(list()), "must be a gstudy or mgstudy object")
})

test_that("latent works with default universe (object only)", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data_latent)
  result <- latent(g)

  expect_true("person" %in% names(result))
  expect_true("latent" %in% names(result))
  expect_equal(nrow(result), length(unique(test_data_latent$person)))
})

test_that("latent returns correct number of rows", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data_latent)
  result <- latent(g)

  n_persons <- length(unique(test_data_latent$person))
  expect_equal(nrow(result), n_persons)
})

# =============================================================================
# Universe specification tests
# =============================================================================

test_that("latent accepts character vector universe", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data_latent)
  result <- latent(g, universe = c("person"))

  expect_s3_class(result, "data.frame")
  expect_true("person" %in% names(result))
})

test_that("latent accepts formula universe", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data_latent)
  result <- latent(g, universe = ~ person)

  expect_s3_class(result, "data.frame")
  expect_true("person" %in% names(result))
})

test_that("latent adds object to universe if missing", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data_latent)

  expect_warning(
    result <- latent(g, universe = "rater"),
    "did not include the object of measurement"
  )

  expect_s3_class(result, "data.frame")
})

# =============================================================================
# Estimator tests
# =============================================================================

test_that("latent works with lme4 estimator", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater),
              data = test_data_latent,
              estimator = "lme4")
  result <- latent(g)

  expect_s3_class(result, "data.frame")
  expect_true("latent" %in% names(result))
})

test_that("latent works with mom estimator", {
  g <- gstudy(score ~ (1 | person) + (1 | rater),
    data = test_data_latent,
    estimator = "mom"
  )
  result <- latent(g)

  expect_s3_class(result, "data.frame")
  expect_true("latent" %in% names(result))
})

test_that("latent works with univariate brms estimator", {
  skip_if_not_installed("brms")
  skip_on_cran()

  g <- gstudy(score ~ (1 | person) + (1 | rater),
    data = test_data_latent,
    estimator = "brms",
    chains = 1,
    iter = 500
  )
  result <- latent(g)

  expect_s3_class(result, "data.frame")
  expect_true("person" %in% names(result))
  expect_true("latent" %in% names(result))
  expect_equal(nrow(result), length(unique(test_data_latent$person)))
})

# =============================================================================
# Multivariate tests
# =============================================================================
# Multivariate tests
# =============================================================================

test_that("latent works with mgstudy objects", {
  skip_if_not_installed("brms")
  skip_on_cran()

  test_data_mv <- data.frame(
    score1 = rnorm(100),
    score2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    estimator = "brms"
  )

  result <- latent(g)

  expect_s3_class(result, "data.frame")
  expect_true("person" %in% names(result))
  expect_true("latent_score1" %in% names(result))
  expect_true("latent_score2" %in% names(result))
})

test_that("latent returns correct structure for mom multivariate", {
  test_data_mv <- data.frame(
    score1 = rnorm(100),
    score2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )

  g <- gstudy(
    mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    estimator = "mom"
  )

  result <- latent(g)

  expect_s3_class(result, "data.frame")
  expect_true("person" %in% names(result))
  expect_true("latent_score1" %in% names(result))
  expect_true("latent_score2" %in% names(result))
})

test_that("latent returns wide format for mgstudy", {
  skip_if_not_installed("brms")
  skip_on_cran()

  test_data_mv <- data.frame(
    score1 = rnorm(100),
    score2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    estimator = "brms"
  )

  result <- latent(g)

  n_persons <- length(unique(test_data_mv$person))
  expect_equal(nrow(result), n_persons)
  expect_true(ncol(result) >= 3)
})

# =============================================================================
# Interaction tests
# =============================================================================

test_that("latent handles interaction in universe", {
  skip_if_not_installed("lme4")

  set.seed(42)
  test_data_interaction <- expand.grid(
    person = factor(1:10),
    rater = factor(1:5)
  )
  test_data_interaction <- rbind(test_data_interaction, test_data_interaction)
  test_data_interaction$score <- rnorm(nrow(test_data_interaction))

  g <- gstudy(
    score ~ (1 | person) + (1 | rater) + (1 | person:rater),
    data = test_data_interaction
  )

  result <- latent(g, universe = c("person", "person:rater"))

  expect_s3_class(result, "data.frame")
  expect_true("person" %in% names(result))
  expect_true("latent" %in% names(result))
})

test_that("latent handles interaction in universe with brms estimator", {
  skip_if_not_installed("brms")
  skip_on_cran()

  set.seed(42)
  test_data_interaction <- expand.grid(
    person = factor(1:10),
    rater = factor(1:5)
  )
  test_data_interaction <- rbind(test_data_interaction, test_data_interaction)
  test_data_interaction$score <- rnorm(nrow(test_data_interaction))

  g <- gstudy(
    score ~ (1 | person) + (1 | rater) + (1 | person:rater),
    data = test_data_interaction,
    estimator = "brms",
    chains = 1,
    iter = 500
  )

  result <- latent(g, universe = c("person", "person:rater"))

  expect_s3_class(result, "data.frame")
  expect_true("person" %in% names(result))
  expect_true("latent" %in% names(result))
})

# =============================================================================
# Edge case tests
# =============================================================================

test_that("latent handles empty results gracefully", {
  skip_if_not_installed("lme4")

  single_level_data <- data.frame(
    score = rnorm(10),
    person = factor(rep(1, 10)),
    rater = factor(rep(1:5, each = 2))
  )

  expect_error(
    gstudy(score ~ (1 | person) + (1 | rater), data = single_level_data),
    "sampled level"
  )
})

test_that("latent produces numeric latent scores", {
  skip_if_not_installed("lme4")
  g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data_latent)
  result <- latent(g)

  expect_true(is.numeric(result$latent))
})
