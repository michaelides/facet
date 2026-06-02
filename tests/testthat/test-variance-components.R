test_that("extract_vc_brms_long_format parses long-format variance components", {
  skip_if_not_installed("brms")
  skip_if_not_installed("posterior")
  skip_on_cran()

  mock_draws <- matrix(
    c(0.5, 0.6, 0.55,
      0.3, 0.35, 0.32,
      0.0, 0.1, 0.05,
      0.2, 0.25, 0.22),
    nrow = 3
  )
  colnames(mock_draws) <- c(
    "sd_Person__SubtestA_Intercept",
    "sd_Person__SubtestB_Intercept",
    "b_sigma_SubtestA",
    "b_sigma_SubtestB"
  )

  result <- extract_vc_brms_long_format(
    draws = mock_draws,
    dimensions = c("SubtestA", "SubtestB"),
    facets = "Person"
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("component", "dim", "var", "pct") %in% names(result)))
  expect_true(all(unique(result$dim) %in% c("SubtestA", "SubtestB")))
})
