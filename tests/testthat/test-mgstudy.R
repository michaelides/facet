# Tests for mgstudy class (multivariate G-studies)

# Test data for multivariate studies
test_data_mv <- data.frame(
  score1 = rnorm(100),
  score2 = rnorm(100),
  person = factor(rep(1:20, 5)),
  rater = factor(rep(1:5, each = 20)),
  item = factor(rep(1:10, each = 10))
)

# =============================================================================
# Basic mgstudy class tests
# =============================================================================

test_that("is.mgstudy returns FALSE for non-mgstudy objects", {
  expect_false(is.mgstudy(list()))
  expect_false(is.mgstudy(NULL))
  expect_false(is.mgstudy("not an mgstudy"))
})

test_that("is.mgstudy returns TRUE for mgstudy objects", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )
  expect_true(is.mgstudy(result))
})

test_that("is.gstudy returns FALSE for mgstudy objects", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )
  expect_false(is.gstudy(result))
})

# =============================================================================
# mgstudy object structure tests
# =============================================================================

test_that("mgstudy object has correct structure", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_s3_class(result, "mgstudy")
  expect_true("model" %in% names(result))
  expect_true("variance_components" %in% names(result))
  expect_true("facets" %in% names(result))
  expect_true("object" %in% names(result))
  expect_true("dimensions" %in% names(result))
})

test_that("mgstudy stores dimensions correctly", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_equal(sort(result$dimensions), c("score1", "score2"))
})

# =============================================================================
# Variance components with dim column tests
# =============================================================================

test_that("mgstudy variance_components has dim column", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_true("dim" %in% names(result$variance_components))
})

test_that("mgstudy variance_components has rows for each dimension", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  vc <- result$variance_components
  unique_dims <- unique(vc$dim)

  expect_equal(sort(unique_dims), c("score1", "score2"))

  vc_score1 <- vc[vc$dim == "score1", ]
  vc_score2 <- vc[vc$dim == "score2", ]

  expect_equal(nrow(vc_score1), nrow(vc_score2))
})

test_that("mgstudy variance_components component names are consistent", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  vc <- result$variance_components
  components_score1 <- sort(unique(vc$component[vc$dim == "score1"]))
  components_score2 <- sort(unique(vc$component[vc$dim == "score2"]))

  expect_equal(components_score1, components_score2)
})

# =============================================================================
# Univariate gstudy dim column tests
# =============================================================================

test_that("univariate gstudy has dim column", {
  skip_if_not_installed("lme4")

  result <- gstudy(score1 ~ (1 | person) + (1 | rater), data = test_data_mv)

  expect_true("dim" %in% names(result$variance_components))
})

test_that("univariate gstudy dim column has response name", {
  skip_if_not_installed("lme4")

  result <- gstudy(score1 ~ (1 | person) + (1 | rater), data = test_data_mv)

  vc <- result$variance_components
  expect_true(all(vc$dim == "score1"))
})

# =============================================================================
# Print and summary methods
# =============================================================================

test_that("print.mgstudy works", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_output(print(result), "Multivariate Generalizability Study")
  expect_output(print(result), "Dimensions:")
})

test_that("summary.mgstudy works", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_output(summary(result), "Multivariate G-Study Summary")
  expect_output(summary(result), "Dimensions:")
})

test_that("print.mgstudy displays residual correlations in long format by default", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_output(print(result), "Correlated Residuals")
  expect_output(print(result), "dim1.*dim2.*estimate")
})

test_that("print.mgstudy can display correlations in matrix format", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_output(print(result, cor_format = "matrix"), "Residual Correlations:")
})

test_that("summary.mgstudy displays residual correlations in long format by default", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_output(summary(result), "Correlated Residuals")
})

test_that("summary.mgstudy can display correlations in matrix format", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_output(summary(result, cor_format = "matrix"), "Residual Correlations:")
})

test_that("mgstudy correlations include full statistics", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_true(!is.null(result$correlations$residual_cor))
  expect_true(all(c("dim1", "dim2", "estimate", "se", "lower", "upper") %in% names(result$correlations$residual_cor)))
})

test_that("mgstudy stores both tibble and matrix correlation formats", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_true(!is.null(result$correlations$residual_cor))
  expect_true(!is.null(result$correlations$residual_cor_matrix))
  expect_true(is.data.frame(result$correlations$residual_cor))
  expect_true(is.matrix(result$correlations$residual_cor_matrix))
})

# =============================================================================
# D-study with mgstudy
# =============================================================================

test_that("dstudy accepts mgstudy objects", {
  skip_if_not_installed("brms")
  skip_on_cran()

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_s3_class(g, "mgstudy")

  d <- dstudy(g, n = list(rater = 3))

  expect_s3_class(d, "dstudy")
})

test_that("dstudy coefficients have dim column for mgstudy", {
  skip_if_not_installed("brms")
  skip_on_cran()

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  d <- dstudy(g, n = list(rater = 3))

  expect_true("dim" %in% names(d$coefficients))
})

test_that("dstudy computes coefficients per dimension", {
  skip_if_not_installed("brms")
  skip_on_cran()

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  d <- dstudy(g, n = list(rater = 3))

  expect_equal(length(unique(d$coefficients$dim)), 2)
  expect_true("g" %in% names(d$coefficients))
  expect_true("phi" %in% names(d$coefficients))
})

test_that("dstudy is_multivariate field is set correctly", {
  skip_if_not_installed("brms")
  skip_on_cran()

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  d <- dstudy(g, n = list(rater = 3))

  expect_true(d$is_multivariate)
})

test_that("dstudy handles multivariate sweep mode", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )
  
  d <- dstudy(g, n = list(rater = c(2, 3, 4)))
  
  expect_true(d$is_sweep)
  expect_true(d$is_multivariate)
  expect_true("dim" %in% names(d$coefficients))
  expect_equal(length(unique(d$coefficients$dim)), 2)
  expect_equal(nrow(d$coefficients), 6) # 3 rater values * 2 dimensions
})

test_that("print.dstudy separates tables by dimension for sweep", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )
  
  d <- dstudy(g, n = list(rater = c(2, 3, 4)))
  
  output <- capture.output(print(d))
  
  expect_true(any(grepl("Dimension: score1", output)))
  expect_true(any(grepl("Dimension: score2", output)))
})

test_that("summary.dstudy separates tables by dimension for sweep", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )
  
  d <- dstudy(g, n = list(rater = c(2, 3, 4)))
  
  output <- capture.output(summary(d))
  
  expect_true(any(grepl("Dimension: score1", output)))
  expect_true(any(grepl("Dimension: score2", output)))
  expect_true(any(grepl("Highest G coefficient for", output)))
})

test_that("plot.dstudy facets by dimension for multivariate sweep", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_on_cran()
  
  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )
  
  d <- dstudy(g, n = list(rater = c(2, 3, 4)))
  
  p <- plot(d, type = "sweep", coefficient = "both")
  
  expect_s3_class(p, "ggplot")
  expect_true("FacetGrid" %in% class(p$facet))
})

test_that("plot.dstudy facets by dimension for single coefficient", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_on_cran()
  
  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )
  
  d <- dstudy(g, n = list(rater = c(2, 3, 4)))
  
  p_g <- plot(d, type = "sweep", coefficient = "g")
  p_phi <- plot(d, type = "sweep", coefficient = "phi")
  
  expect_s3_class(p_g, "ggplot")
  expect_s3_class(p_phi, "ggplot")
  expect_true("FacetWrap" %in% class(p_g$facet))
  expect_true("FacetWrap" %in% class(p_phi$facet))
})

test_that("posterior estimation produces separate results per dimension", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )
  
  d_post <- dstudy(g, n = list(rater = 3), estimation = "posterior")
  
  expect_true("dim" %in% names(d_post$coefficients))
  expect_equal(length(unique(d_post$coefficients$dim)), 2)
  expect_true(is.list(d_post$posterior))
  expect_true(all(c("score1", "score2") %in% names(d_post$posterior)))
})

test_that("posterior sweep produces separate results per dimension", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )
  
  d_post <- dstudy(g, n = list(rater = c(2, 3, 4)), estimation = "posterior")
  
  expect_true("dim" %in% names(d_post$coefficients))
  expect_equal(nrow(d_post$coefficients), 6) # 3 rater values * 2 dimensions
  expect_equal(length(unique(d_post$coefficients$dim)), 2)
})

test_that("posterior results are consistent with simple estimation", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms",
    chains = 2,
    iter = 1000
  )
  
  d_simple <- dstudy(g, n = list(rater = 3), estimation = "simple")
  d_post <- dstudy(g, n = list(rater = 3), estimation = "posterior")
  
  # Check that posterior means are close to simple estimates (within 0.15 tolerance)
  for (dim in c("score1", "score2")) {
    simple_g <- d_simple$coefficients$g[d_simple$coefficients$dim == dim]
    post_g <- d_post$coefficients$g[d_post$coefficients$dim == dim]
    expect_equal(simple_g, post_g, tolerance = 0.15)
    
    simple_phi <- d_simple$coefficients$phi[d_simple$coefficients$dim == dim]
    post_phi <- d_post$coefficients$phi[d_post$coefficients$dim == dim]
    expect_equal(simple_phi, post_phi, tolerance = 0.15)
  }
})

test_that("no deprecated posterior_samples warning", {
  skip_if_not_installed("brms")
  skip_on_cran()
  
  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )
  
  # Should not produce deprecated method warning
  expect_warning(
    d <- dstudy(g, n = list(rater = 3), estimation = "posterior"),
    regexp = "deprecated"
  ) |> expect_null()
})

test_that("extract_variance_draws returns correct structure for multivariate", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_installed("posterior")
  
  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )
  
  draws <- brms::as_draws_matrix(g$model)
  vc_draws <- facet:::extract_variance_draws(g, draws)
  
  # Should be nested list
  expect_true(is.list(vc_draws))
  expect_true(all(c("score1", "score2") %in% names(vc_draws)))
  
  # Each dimension should have variance components
  expect_true(is.list(vc_draws[["score1"]]))
  expect_true("person" %in% names(vc_draws[["score1"]]))
  expect_true("Residual" %in% names(vc_draws[["score1"]]))
  
  # Each component should be numeric vector
  expect_true(is.numeric(vc_draws[["score1"]][["person"]]))
  expect_true(is.numeric(vc_draws[["score1"]][["Residual"]]))
})

# =============================================================================
# Object of measurement
# =============================================================================

test_that("mgstudy object of measurement is single value", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_equal(result$object, "person")
})

test_that("mgstudy accepts explicit object of measurement", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms",
    object = "rater"
  )

  expect_equal(result$object, "rater")
})

# =============================================================================
# Plot method for mgstudy
# =============================================================================

test_that("plot.mgstudy returns ggplot for variance type", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  p <- plot(result, type = "variance")

  expect_s3_class(p, "ggplot")
  expect_true("FacetWrap" %in% class(p$facet))
})

test_that("plot.mgstudy returns ggplot for proportion type", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  p <- plot(result, type = "proportion")

  expect_s3_class(p, "ggplot")
  expect_true("FacetWrap" %in% class(p$facet))
})

test_that("plot.mgstudy returns ggplot for forest type with CIs", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  p <- plot(result, type = "forest")

  expect_s3_class(p, "ggplot")
  expect_true("FacetWrap" %in% class(p$facet))
})

test_that("plot.mgstudy falls back to variance for forest without CIs", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_warning(
    p <- plot(result, type = "forest"),
    "Confidence intervals not available"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot.mgstudy works with mom backend", {
  skip_if_not_installed("ggplot2")

  test_data_mv_mom <- data.frame(
    score1 = rnorm(100),
    score2 = rnorm(100),
    person = factor(rep(1:20, 5)),
    rater = factor(rep(1:5, each = 20))
  )

  result <- gstudy(
    mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv_mom,
    backend = "mom"
  )

  p <- plot(result, type = "variance")

  expect_s3_class(p, "ggplot")
  expect_true("FacetWrap" %in% class(p$facet))
})
