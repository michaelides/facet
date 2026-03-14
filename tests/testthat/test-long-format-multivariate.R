# Integration tests for long-format multivariate models
# These tests require brms and are skipped on CRAN

test_that("gstudy works with rajaratnam long-format model", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  data(rajaratnam)
  
  # Fit the model (with minimal iterations for speed)
  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||Item),
       sigma ~ 0 + Subtest),
    data = rajaratnam,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  
  # Check class
  expect_s3_class(gu, "mgstudy")
  expect_true(gu$long_format_multivariate)
  expect_equal(gu$dimension_var, "Subtest")
  
  # Check dimensions
  expect_true(length(gu$dimensions) == 3)
  
  # Check sample sizes per dimension
  expect_true("dim" %in% names(gu$sample_size_tibble))
  
  # Check variance components
  expect_true("dim" %in% names(gu$variance_components))
  expect_true(all(gu$dimensions %in% unique(gu$variance_components$dim)))
  
  # Check sigma is on original scale (not log scale)
  # The sigma values should be positive and reasonable
  residual_var <- gu$variance_components[
    gu$variance_components$component == "Residual",
  ]
  expect_true(all(residual_var$var > 0))
  expect_true(all(residual_var$var < 100))  # Reasonable upper bound
})

test_that("dstudy works with long-format multivariate", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  data(rajaratnam)
  
  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||Item),
       sigma ~ 0 + Subtest),
    data = rajaratnam,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  
  d <- dstudy(gu, n = list(Person = 5))
  
  expect_s3_class(d, "dstudy")
  expect_true("dim" %in% names(d$coefficients))
  expect_true(d$long_format_multivariate)
})

test_that("print and summary work for long-format multivariate", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  data(rajaratnam)
  
  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||Item),
       sigma ~ 0 + Subtest),
    data = rajaratnam,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )
  
  # Test print
  output <- capture.output(print(gu))
  expect_true(any(grepl("Long-Format", output)))
  expect_true(any(grepl("Dimension Variable", output)))
  
  # Test summary
  output <- capture.output(summary(gu))
  expect_true(any(grepl("Long-Format", output)))
})
