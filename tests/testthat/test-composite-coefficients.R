# Tests for composite coefficient helper functions

test_that("build_component_covariance_matrix creates correct diagonal matrix", {
  vc <- data.frame(
    component = c("Person", "Person", "Rater", "Rater", "Residual", "Residual"),
    dim = c("A", "B", "A", "B", "A", "B"),
    var = c(0.5, 0.3, 0.2, 0.15, 0.3, 0.25)
  )
  mat <- facet:::build_component_covariance_matrix(vc, "Person", c("A", "B"), NULL)
  expect_equal(mat["A", "A"], 0.5)
  expect_equal(mat["B", "B"], 0.3)
  expect_equal(mat["A", "B"], 0)
  expect_equal(mat["B", "A"], 0)
})

test_that("compute_weighted_variance computes correct weighted sum", {
  cov_mat <- matrix(c(0.5, 0.1, 0.1, 0.3), 2, 2)
  rownames(cov_mat) <- colnames(cov_mat) <- c("A", "B")
  weights <- c(1, 1)
  result <- facet:::compute_weighted_variance(cov_mat, weights)
  expect_equal(result, 1.0)
  weights2 <- c(2, 1)
  result2 <- facet:::compute_weighted_variance(cov_mat, weights2)
  expect_equal(result2, 2.7)
})
