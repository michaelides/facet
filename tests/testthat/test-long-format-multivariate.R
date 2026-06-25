# Integration tests for long-format multivariate models
# These tests require brms and are skipped on CRAN

test_that("gstudy works with rajaratnam long-format model", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data(rajaratnam)

  # Fit the model (with minimal iterations for speed)
  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
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
  expect_true(all(residual_var$var < 100)) # Reasonable upper bound
})

test_that("dstudy works with long-format multivariate", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data(rajaratnam)

  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
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
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
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

test_that("dstudy posterior estimation works with long-format multivariate", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data(rajaratnam)

  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
       sigma ~ 0 + Subtest),
    data = rajaratnam,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )

  d <- dstudy(gu, n = list(Person = 5), estimation = "posterior")

  expect_s3_class(d, "dstudy")
  expect_true("dim" %in% names(d$coefficients))
  expect_true(d$long_format_multivariate)
  expect_true(!is.null(d$posterior))
})

test_that("dstudy sweep mode works with long-format multivariate", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data(rajaratnam)

  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
       sigma ~ 0 + Subtest),
    data = rajaratnam,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )

  d <- dstudy(gu, n = list(Person = c(2, 5, 10)))

  expect_s3_class(d, "dstudy")
  expect_true(d$is_sweep)
  expect_true("dim" %in% names(d$coefficients))
  expect_equal(length(unique(d$coefficients$dim)), 3)
})

test_that("dstudy composite coefficients work with long-format multivariate", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data(rajaratnam)

  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
       sigma ~ 0 + Subtest),
    data = rajaratnam,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )

  d <- dstudy(gu, n = list(Person = 5), weights = c(1, 1, 1))

  expect_s3_class(d, "dstudy")
  expect_true("dim" %in% names(d$coefficients))
  expect_true("Composite" %in% d$coefficients$dim)
})

test_that("extract_variance_draws works with long-format models", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_installed("posterior")

  data(rajaratnam)

  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
       sigma ~ 0 + Subtest),
    data = rajaratnam,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )

  draws <- brms::as_draws_matrix(gu$model)
  vc_draws <- facet:::extract_variance_draws(gu, draws)

  expect_true(is.list(vc_draws))
  expect_true(length(vc_draws) == 3)

  for (dim in names(vc_draws)) {
    expect_true(is.list(vc_draws[[dim]]))
    expect_true("Person" %in% names(vc_draws[[dim]]))
    expect_true("Residual" %in% names(vc_draws[[dim]]))

    for (comp in names(vc_draws[[dim]])) {
      expect_true(is.numeric(vc_draws[[dim]][[comp]]))
      expect_true(all(vc_draws[[dim]][[comp]] >= 0))
    }
  }
})

test_that("dstudy dimension_var is passed through from gstudy", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data(rajaratnam)

  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
       sigma ~ 0 + Subtest),
    data = rajaratnam,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )

  d <- dstudy(gu, n = list(Person = 5))

  expect_equal(d$dimension_var, "Subtest")
})

test_that("long-format multivariate coefficients are valid", {
  skip_if_not_installed("brms")
  skip_on_cran()

  data(rajaratnam)

  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
       sigma ~ 0 + Subtest),
    data = rajaratnam,
    backend = "brms",
    chains = 2,
    iter = 500,
    refresh = 0
  )

  d <- dstudy(gu, n = list(Person = 5))

  for (dim in unique(d$coefficients$dim)) {
    dim_coefs <- d$coefficients[d$coefficients$dim == dim, ]
    expect_true(dim_coefs$g >= 0 && dim_coefs$g <= 1)
    expect_true(dim_coefs$phi >= 0 && dim_coefs$phi <= 1)
    expect_true(dim_coefs$uni >= 0)
    expect_true(dim_coefs$sigma2_delta >= 0)
    expect_true(dim_coefs$sigma2_delta_abs >= 0)
  }
})

test_that("dstudy works with long-format multivariate", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  data(rajaratnam)
  
  gu <- gstudy(
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
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
    bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||ItemId),
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
