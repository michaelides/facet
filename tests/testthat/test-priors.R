# Tests for prior wrapper functions

test_data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, 5)),
  rater = factor(rep(1:5, each = 20)),
  item = factor(rep(1:10, each = 10))
)

# =============================================================================
# Prior wrapper function tests (require brms)
# =============================================================================

test_that("set_prior returns a brmsprior object", {
  skip_if_not_installed("brms")
  prior <- set_prior("normal(0, 1)", class = "sd")
  expect_s3_class(prior, "brmsprior")
  expect_true(nrow(prior) > 0)
})

test_that("set_prior with group argument works", {
  skip_if_not_installed("brms")
  prior <- set_prior("normal(0, 1)", class = "sd", group = "person")
  expect_s3_class(prior, "brmsprior")
  expect_equal(prior$group[1], "person")
})

test_that("default_prior works with formula and data", {
  skip_if_not_installed("brms")
  prior_info <- default_prior(score ~ (1 | person) + (1 | item), data = test_data)
  expect_s3_class(prior_info, "brmsprior")
  expect_true(nrow(prior_info) > 0)
})

test_that("prior returns a brmsprior object", {
  skip_if_not_installed("brms")
  prior <- prior(normal(0, 1), class = "sd")
  expect_s3_class(prior, "brmsprior")
})

test_that("prior_ returns a brmsprior object", {
  skip_if_not_installed("brms")
  prior <- prior_(~ normal(0, 1), class = "sd")
  expect_s3_class(prior, "brmsprior")
})

test_that("prior_string returns a brmsprior object", {
  skip_if_not_installed("brms")
  prior <- prior_string("normal(0, 1)", class = "sd")
  expect_s3_class(prior, "brmsprior")
})

test_that("empty_prior returns empty brmsprior object", {
  skip_if_not_installed("brms")
  prior <- empty_prior()
  expect_s3_class(prior, "brmsprior")
  expect_equal(nrow(prior), 0)
})

# =============================================================================
# gstudy prior integration tests (require brms)
# =============================================================================

test_that("gstudy accepts prior argument with brms backend", {
  skip_if_not_installed("brms")
  my_prior <- set_prior("normal(0, 1)", class = "sd", group = "person")
  
  result <- gstudy(
    score ~ (1 | person) + (1 | rater),
    data = test_data,
    prior = my_prior,
    backend = "brms",
    chains = 1,
    iter = 500,
    refresh = 0
  )
  
  expect_s3_class(result, "gstudy")
  expect_equal(result$backend, "brms")
})

test_that("gstudy errors when prior used with lme4 backend", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("brms")
  my_prior <- set_prior("normal(0, 1)", class = "sd")
  
  expect_error(
    gstudy(
      score ~ (1 | person) + (1 | rater),
      data = test_data,
      prior = my_prior,
      backend = "lme4"
    ),
    "prior is only supported with brms backend"
  )
})

test_that("gstudy errors when prior used with mom backend", {
  skip_if_not_installed("brms")
  my_prior <- set_prior("normal(0, 1)", class = "sd")
  
  expect_error(
    gstudy(
      score ~ (1 | person) + (1 | rater),
      data = test_data,
      prior = my_prior,
      backend = "mom"
    ),
    "prior is only supported with brms backend"
  )
})

test_that("gstudy works with prior for multiple random effects", {
  skip_if_not_installed("brms")
  my_prior <- c(
    set_prior("normal(0, 1)", class = "sd", group = "person"),
    set_prior("normal(0, 1)", class = "sd", group = "rater")
  )
  
  result <- gstudy(
    score ~ (1 | person) + (1 | rater),
    data = test_data,
    prior = my_prior,
    backend = "brms",
    chains = 1,
    iter = 500,
    refresh = 0
  )
  
  expect_s3_class(result, "gstudy")
})

test_that("gstudy works with empty_prior", {
  skip_if_not_installed("brms")
  my_prior <- empty_prior()
  
  result <- gstudy(
    score ~ (1 | person) + (1 | rater),
    data = test_data,
    prior = my_prior,
    backend = "brms",
    chains = 1,
    iter = 500,
    refresh = 0
  )
  
  expect_s3_class(result, "gstudy")
})
