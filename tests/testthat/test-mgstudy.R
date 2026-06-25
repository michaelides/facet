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

  expect_output(summary(result), "Multivariate G Study Summary")
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

test_that("mgstudy object of measurement defaults to first facet", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_equal(result$object, "person")
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

# =============================================================================
# Long-format multivariate tests
# =============================================================================

test_that("print.mgstudy handles long_format_multivariate flag", {
  mock_mgstudy <- list(
    backend = "brms",
    long_format_multivariate = TRUE,
    dimension_var = "Subtest",
    dimensions = c("A", "B"),
    sample_size_tibble = tibble::tibble(
      dim = c("A", "A", "B", "B"),
      effect = c("Person", "Residual", "Person", "Residual"),
      type = c("main", "residual", "main", "residual"),
      n = c(10, 40, 8, 32)
    ),
    variance_components = tibble::tibble(
      component = c("Person", "Residual", "Person", "Residual"),
      dim = c("A", "A", "B", "B"),
      type = c("main", "residual", "main", "residual"),
      var = c(0.5, 0.3, 0.4, 0.35),
      sd = c(0.707, 0.548, 0.632, 0.592),
      pct = c(62.5, 37.5, 53.3, 46.7)
    ),
    object = "Person",
    facets = c("Person"),
    n_obs = 80,
    correlations = list(
      random_effect_cor = list(),
      residual_cor = NULL
    )
  )
  class(mock_mgstudy) <- "mgstudy"

  output <- capture.output(print(mock_mgstudy))

  expect_true(any(grepl("Long-Format", output)))
  expect_true(any(grepl("Dimension Variable", output)))
  expect_true(any(grepl("Subtest", output)))
})

test_that("print.mgstudy shows per-dimension sample sizes", {
  mock_mgstudy <- list(
    backend = "brms",
    long_format_multivariate = TRUE,
    dimension_var = "Subtest",
    dimensions = c("A", "B"),
    sample_size_tibble = tibble::tibble(
      dim = c("A", "B"),
      effect = c("Person", "Person"),
      type = c("main", "main"),
      n = c(10, 8)
    ),
    variance_components = tibble::tibble(
      component = character(),
      dim = character(),
      type = character(),
      var = numeric(),
      sd = numeric(),
      pct = numeric()
    ),
    object = "Person",
    facets = c("Person"),
    n_obs = 80,
    correlations = list()
  )
  class(mock_mgstudy) <- "mgstudy"

  output <- capture.output(print(mock_mgstudy))

  expect_true(any(grepl("dim", output)))
})

test_that("summary.mgstudy handles long_format_multivariate flag", {
  mock_mgstudy <- list(
    backend = "brms",
    long_format_multivariate = TRUE,
    dimension_var = "Subtest",
    dimensions = c("A", "B"),
    sample_size_tibble = tibble::tibble(
      dim = c("A", "B"),
      effect = c("Person", "Person"),
      type = c("main", "main"),
      n = c(10, 8)
    ),
    variance_components = tibble::tibble(
      component = c("Person", "Residual", "Person", "Residual"),
      dim = c("A", "A", "B", "B"),
      type = c("main", "residual", "main", "residual"),
      var = c(0.5, 0.3, 0.4, 0.35),
      pct = c(62.5, 37.5, 53.3, 46.7)
    ),
    object = "Person",
    facets = c("Person"),
    n_obs = 80,
    correlations = list()
  )
  class(mock_mgstudy) <- "mgstudy"

  output <- capture.output(summary(mock_mgstudy))

  expect_true(any(grepl("Long-Format", output)))
})

# =============================================================================
# Random effect correlations (correlated effects with "pr" syntax)
# =============================================================================

test_that("mgstudy stores both correlations and covariances for brms models", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | pr | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms",
    chains = 2,
    iter = 500
  )

  expect_true(!is.null(result$correlations$residual_cor))
  expect_true(!is.null(result$correlations$residual_cov))
  expect_true(!is.null(result$correlations$residual_cor_matrix))
  expect_true(!is.null(result$correlations$residual_cov_matrix))
})

test_that("mgstudy stores random effect correlations for correlated effects", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | pr | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms",
    chains = 2,
    iter = 500
  )

  expect_true(!is.null(result$correlations$random_effect_cor))
  expect_true(!is.null(result$correlations$random_effect_cor[["person"]]))
  expect_true(!is.null(result$correlations$random_effect_cov))
  expect_true(!is.null(result$correlations$random_effect_cov[["person"]]))
})

test_that("summary.mgstudy displays random effect correlations for correlated effects", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | pr | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms",
    chains = 2,
    iter = 500
  )

  output <- capture.output(summary(result))

  expect_true(any(grepl("Correlated Random Effects.*person", output, ignore.case = TRUE)))
})

test_that("summary.mgstudy displays random effect correlations with vc_format = dimension", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | pr | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms",
    chains = 2,
    iter = 500
  )

  output <- capture.output(summary(result, vc_format = "dimension"))

  expect_true(any(grepl("Correlated Random Effects", output)))
})

test_that("summary.mgstudy displays random effect correlations with vc_format = facet", {
  skip_if_not_installed("brms")
  skip_on_cran()

  result <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | pr | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms",
    chains = 2,
    iter = 500
  )

  output <- capture.output(summary(result, vc_format = "facet"))

  expect_true(any(grepl("Correlations.*person", output, ignore.case = TRUE)))
})

# =============================================================================
# Mom backend with unbalanced designs (Henderson's Method III)
# =============================================================================

test_data_mv_unbalanced <- data.frame(
  id   = factor(rep(1:20, 5)),
  item = factor(rep(1:5, each = 20)),
  A    = rnorm(100),
  B    = rnorm(100)
)
test_data_mv_unbalanced$B[sample(1:100, 20)] <- NA

test_that("gstudy with mom backend and unbalanced = TRUE sets is_unbalanced and n_per_dim", {
  result <- gstudy(
    brms::mvbind(A, B) ~ (1 | id) + (1 | item),
    data = test_data_mv_unbalanced,
    backend = "mom",
    unbalanced = TRUE
  )

  expect_true(isTRUE(result$is_unbalanced))
  expect_equal(result$backend, "mom")
  expect_equal(result$dimensions, c("A", "B"))
  expect_false(is.null(result$n_per_dim))
  expect_equal(result$n_per_dim$A, 100)
  expect_equal(result$n_per_dim$B, 80)
})

test_that("gstudy with mom balanced (default) does not set is_unbalanced", {
  balanced_data <- data.frame(
    id   = factor(rep(1:20, 5)),
    item = factor(rep(1:5, each = 20)),
    A    = rnorm(100),
    B    = rnorm(100)
  )

  result <- gstudy(
    brms::mvbind(A, B) ~ (1 | id) + (1 | item),
    data = balanced_data,
    backend = "mom"
  )

  expect_false(isTRUE(result$is_unbalanced))
  expect_null(result$n_per_dim)
})

test_that("print.mgstudy surfaces the unbalanced banner and per-dim totals", {
  result <- gstudy(
    brms::mvbind(A, B) ~ (1 | id) + (1 | item),
    data = test_data_mv_unbalanced,
    backend = "mom",
    unbalanced = TRUE
  )

  output <- capture.output(print(result))

  expect_true(any(grepl("Unbalanced", output)))
  expect_true(any(grepl("A.*100", output)))
  expect_true(any(grepl("B.*80", output)))
})

test_that("dstudy with per-dim n tibble uses per-dim sample sizes in non-posterior path", {
  result <- gstudy(
    brms::mvbind(A, B) ~ (1 | id) + (1 | item),
    data = test_data_mv_unbalanced,
    backend = "mom",
    unbalanced = TRUE
  )

  n_tb <- tibble::tibble(
    dim   = c("A", "A", "B", "B"),
    facet = c("id", "item", "id", "item"),
    n     = c(20, 5, 16, 4)
  )
  d <- dstudy(result, n = n_tb)

  expect_false(isTRUE(d$is_sweep))

  vc_A_resid <- result$variance_components$var[
    result$variance_components$dim == "A" &
      result$variance_components$component == "Residual"
  ]
  vc_B_resid <- result$variance_components$var[
    result$variance_components$dim == "B" &
      result$variance_components$component == "Residual"
  ]

  sigma_A <- d$coefficients$sigma2_delta[d$coefficients$dim == "A"]
  sigma_B <- d$coefficients$sigma2_delta[d$coefficients$dim == "B"]
  expect_equal(sigma_A, vc_A_resid / (20 * 5), tolerance = 1e-10)
  expect_equal(sigma_B, vc_B_resid / (16 * 4), tolerance = 1e-10)
})

test_that("dstudy with n as a per-dim tibble keeps the coefficients tibble free of facet columns", {
  result <- gstudy(
    brms::mvbind(A, B) ~ (1 | id) + (1 | item),
    data = test_data_mv_unbalanced,
    backend = "mom",
    unbalanced = TRUE
  )

  n_tb <- tibble::tibble(
    dim   = c("A", "A", "B", "B"),
    facet = c("id", "item", "id", "item"),
    n     = c(20, 5, 16, 4)
  )
  d <- dstudy(result, n = n_tb)

  forbidden <- c("id", "item", "person", "rater")
  leaked <- intersect(forbidden, names(d$coefficients))
  expect_equal(length(leaked), 0)
  expect_true("dim" %in% names(d$coefficients))
  expect_true("sigma2_delta" %in% names(d$coefficients))
})

test_that("expand_n_per_dim detects is_sweep per (dim, facet), not per facet", {
  n_tb <- tibble::tibble(
    dim   = c("A", "A", "B", "B"),
    facet = c("id", "item", "id", "item"),
    n     = c(20, 5, 16, 4)
  )
  out <- facet:::expand_n_per_dim(n_tb, sweep = TRUE)
  expect_false(out$is_sweep)

  n_tb_sweep <- tibble::tibble(
    dim   = c("A", "A", "A", "A", "B", "B"),
    facet = c("id", "id", "item", "item", "id", "item"),
    n     = c(20, 30, 5, 5, 16, 4)
  )
  out_sweep <- facet:::expand_n_per_dim(n_tb_sweep, sweep = TRUE)
  expect_true(out_sweep$is_sweep)
})

# =============================================================================
# 4-decimal-point display tests
# =============================================================================

test_that("mgstudy variance components are stored at 4 decimal places by default", {
  skip_if_not_installed("brms")
  skip_on_cran()

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  expect_s3_class(g, "mgstudy")
  vc <- g$variance_components
  expect_true("var" %in% names(vc))

  # Stored at 4 dp (allowing for floating point error of 1e-6)
  expect_true(all(abs(vc$var - round(vc$var, 4)) < 1e-6, na.rm = TRUE))
})

test_that("print.mgstudy default uses 4 decimal places", {
  skip_if_not_installed("brms")
  skip_on_cran()

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  out <- capture.output(print(g))
  out <- paste(out, collapse = "\n")

  expect_match(out, "\\d+\\.\\d{4}", all = FALSE)
})

test_that("summary.mgstudy default uses 4 decimal places", {
  skip_if_not_installed("brms")
  skip_on_cran()

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  out <- capture.output(summary(g))
  out <- paste(out, collapse = "\n")

  expect_match(out, "\\d+\\.\\d{4}", all = FALSE)
})

test_that("tidy.mgstudy returns a tibble with numeric columns at 4 decimal places by default", {
  skip_if_not_installed("brms")
  skip_on_cran()

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  result <- tidy(g)

  expect_s3_class(result, "tbl_df")
  numeric_cols <- vapply(result, is.numeric, logical(1))
  numeric_cols <- names(result)[numeric_cols]

  for (col in numeric_cols) {
    expect_true(all(abs(result[[col]] - round(result[[col]], 4)) < 1e-6, na.rm = TRUE))
  }
})

test_that("tidy.mgstudy digits argument overrides default", {
  skip_if_not_installed("brms")
  skip_on_cran()

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  result <- tidy(g, digits = 2)
  expect_s3_class(result, "tbl_df")

  if ("var" %in% names(result)) {
    expect_true(all(abs(result$var - round(result$var, 2)) < 1e-6, na.rm = TRUE))
  }
})

test_that("print.mgstudy digits argument overrides default", {
  skip_if_not_installed("brms")
  skip_on_cran()

  g <- gstudy(
    brms::mvbind(score1, score2) ~ (1 | person) + (1 | rater),
    data = test_data_mv,
    backend = "brms"
  )

  out <- capture.output(print(g, digits = 2))
  out <- paste(out, collapse = "\n")

  expect_false(grepl("\\d+\\.\\d{4}", out))
})


