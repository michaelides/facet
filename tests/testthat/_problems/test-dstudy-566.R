# Extracted from test-dstudy.R:566

# prequel ----------------------------------------------------------------------
test_data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, 5)),
  rater = factor(rep(1:5, each = 20))
)

# test -------------------------------------------------------------------------
skip_if_not_installed("lme4")
skip_if_not_installed("brms")
g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data, backend = "lme4")
expect_warning(
    result <- dstudy(g, n = list(rater = 3), estimation = "posterior"),
    "estimation = 'posterior' requires backend = 'brms'"
  )
