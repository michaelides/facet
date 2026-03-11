# Extracted from test-dstudy.R:377

# prequel ----------------------------------------------------------------------
test_data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, 5)),
  rater = factor(rep(1:5, each = 20))
)

# test -------------------------------------------------------------------------
skip_if_not_installed("lme4")
g <- gstudy(score ~ (1 | person) + (1 | rater), data = test_data)
expect_error(
    dstudy(g, n = list(rater = 3), object = "person", aggregation = "person"),
    "same facet cannot be both the object of measurement and the aggregation"
  )
