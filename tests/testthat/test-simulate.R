# Tests for simulate_gtheory_data function

test_that("simulate_gtheory_data creates basic univariate data", {
  data <- simulate_gtheory_data(
    facets = list(person = 20, item = 10),
    vc = list(person = 1, item = 0.5),
    sd_residual = 1,
    seed = 123
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 200)
  expect_true("person" %in% names(data))
  expect_true("item" %in% names(data))
  expect_true("y" %in% names(data))
  expect_equal(nlevels(data$person), 20)
  expect_equal(nlevels(data$item), 10)
})

test_that("simulate_gtheory_data creates univariate data with interactions", {
  data <- simulate_gtheory_data(
    facets = list(person = 50, item = 10),
    vc = list(person = 1, item = 0.5, "person:item" = 0.3),
    sd_residual = 1,
    seed = 456
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 500)
  expect_true("y" %in% names(data))
})

test_that("simulate_gtheory_data creates wide-format multivariate data", {
  data <- simulate_gtheory_data(
    facets = list(person = 30, item = 5),
    vc = list(person = 1, item = 0.5),
    multivariate = TRUE,
    n_dims = 2,
    format = "wide",
    seed = 789
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 150)
  expect_true("person" %in% names(data))
  expect_true("item" %in% names(data))
  expect_true("y1" %in% names(data))
  expect_true("y2" %in% names(data))
  expect_equal(nlevels(data$person), 30)
  expect_equal(nlevels(data$item), 5)
})

test_that("simulate_gtheory_data creates long-format multivariate data", {
  data <- simulate_gtheory_data(
    facets = list(person = 25, item = 8),
    vc = list(person = 1, item = 0.5),
    multivariate = TRUE,
    n_dims = 2,
    format = "long",
    dim_var = "dimension",
    response_var = "score",
    seed = 101
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 400) # 25 * 8 * 2
  expect_true("person" %in% names(data))
  expect_true("item" %in% names(data))
  expect_true("dimension" %in% names(data))
  expect_true("score" %in% names(data))
  expect_equal(nlevels(data$person), 25)
  expect_equal(nlevels(data$item), 8)
  expect_equal(nlevels(data$dimension), 2)
})

test_that("simulate_gtheory_data respects custom dimension names", {
  data <- simulate_gtheory_data(
    facets = list(person = 10, item = 5),
    vc = list(person = 1, item = 0.5),
    multivariate = TRUE,
    n_dims = 3,
    dim_names = c("math", "reading", "writing"),
    format = "wide",
    seed = 202
  )

  expect_true("math" %in% names(data))
  expect_true("reading" %in% names(data))
  expect_true("writing" %in% names(data))
})

test_that("simulate_gtheory_data handles residual correlations", {
  cor_matrix <- matrix(c(1, 0.7, 0.7, 1), nrow = 2)

  data <- simulate_gtheory_data(
    facets = list(person = 50, item = 20),
    vc = list(person = 0, item = 0),
    multivariate = TRUE,
    n_dims = 2,
    format = "wide",
    residual_cor = cor_matrix,
    sd_residual = 1,
    seed = 303
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 1000)

  cor_est <- cor(data$y1, data$y2)
  expect_true(abs(cor_est - 0.7) < 0.1)
})

test_that("simulate_gtheory_data handles random effect correlations", {
  cor_matrix <- matrix(c(1, 0.6, 0.6, 1), nrow = 2)

  data <- simulate_gtheory_data(
    facets = list(person = 100, item = 10),
    vc = list(person = 1, item = 0.5),
    multivariate = TRUE,
    n_dims = 2,
    format = "wide",
    re_cor = list(person = cor_matrix),
    seed = 404
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 1000)
})

test_that("simulate_gtheory_data validates correlation matrices", {
  bad_cor <- matrix(c(1, 1.5, 1.5, 1), nrow = 2)

  expect_error(
    simulate_gtheory_data(
      facets = list(person = 10, item = 5),
      multivariate = TRUE,
      n_dims = 2,
      residual_cor = bad_cor
    ),
    "positive definite"
  )

  wrong_size_cor <- matrix(c(1, 0.5, 0.5, 1, 0.3, 0.3), nrow = 2)

  expect_error(
    simulate_gtheory_data(
      facets = list(person = 10, item = 5),
      multivariate = TRUE,
      n_dims = 2,
      residual_cor = wrong_size_cor
    ),
    "2x2"
  )
})

test_that("simulate_gtheory_data requires facets or formula", {
  expect_error(
    simulate_gtheory_data(vc = list(person = 1)),
    "Either 'facets' or 'formula'"
  )
})

test_that("simulate_gtheory_data extracts facets from formula", {
  data <- simulate_gtheory_data(
    formula = y ~ (1 | person) + (1 | item),
    sd_residual = 1,
    seed = 505
  )

  expect_s3_class(data, "data.frame")
  expect_true("person" %in% names(data))
  expect_true("item" %in% names(data))
  expect_true("y" %in% names(data))
})

test_that("simulate_gtheory_data produces reproducible results with seed", {
  data1 <- simulate_gtheory_data(
    facets = list(person = 20, item = 10),
    vc = list(person = 1, item = 0.5),
    sd_residual = 1,
    seed = 42
  )

  data2 <- simulate_gtheory_data(
    facets = list(person = 20, item = 10),
    vc = list(person = 1, item = 0.5),
    sd_residual = 1,
    seed = 42
  )

  expect_equal(data1$y, data2$y)
})

test_that("simulate_gtheory_legacy creates univariate data", {
  data <- simulate_gtheory_legacy(
    n_p = 100,
    n_i = 20,
    sd_res = 1,
    sd_p = 1,
    sd_i = 0.5,
    seed = 606
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 2000)
  expect_true("p" %in% names(data))
  expect_true("i" %in% names(data))
  expect_true("y" %in% names(data))
})

test_that("simulate_gtheory_legacy creates multivariate data", {
  data <- simulate_gtheory_legacy(
    n_p = 50,
    n_i = 10,
    n_d = 5,
    n_l = 10,
    multivariate = TRUE,
    n_dims = 2,
    seed = 707
  )

  expect_s3_class(data, "data.frame")
  expect_true("y1" %in% names(data))
  expect_true("y2" %in% names(data))
})

test_that("simulate_gtheory_legacy matches makedata structure", {
  data <- simulate_gtheory_legacy(
    n_p = 200,
    n_i = 25,
    n_d = 5,
    n_l = 20,
    sd_res = 1,
    sd_p = 1,
    sd_i = 0.5,
    sd_d = 0.3,
    sd_l = 0.5,
    sd_pd = 0.3,
    sd_li = 0.3,
    sd_ld = 0.3,
    seed = 808
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 5000)
  expect_equal(length(unique(data$p)), 200)
  expect_equal(length(unique(data$i)), 25)
  expect_equal(length(unique(data$d)), 5)
  expect_equal(length(unique(data$l)), 20)
})

test_that("simulate_gtheory_data works with nested designs", {
  nested_cor <- list(item = "person")

  data <- simulate_gtheory_data(
    facets = list(person = 20, item = 5),
    vc = list(person = 1, item = 0.5),
    nested = nested_cor,
    sd_residual = 1,
    seed = 909
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 100)
})

test_that("simulate_gtheory_data handles zero variance components", {
  data <- simulate_gtheory_data(
    facets = list(person = 20, item = 10),
    vc = list(person = 0, item = 0),
    sd_residual = 1,
    seed = 1010
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 200)
})

test_that("long format has correct structure for brms models", {
  data <- simulate_gtheory_data(
    facets = list(person = 30, rater = 5),
    vc = list(person = 1, rater = 0.5),
    multivariate = TRUE,
    n_dims = 3,
    format = "long",
    dim_var = "dimension",
    response_var = "score",
    seed = 1111
  )

  expect_true("dimension" %in% names(data))
  expect_true("score" %in% names(data))
  expect_true("person" %in% names(data))
  expect_true("rater" %in% names(data))

  expect_equal(length(unique(data$dimension)), 3)
  expect_true(all(c("dim1", "dim2", "dim3") %in% unique(data$dimension)))

  expect_equal(
    nrow(data),
    length(unique(data$person)) * length(unique(data$rater)) * length(unique(data$dimension))
  )
})

test_that("simulate_gtheory_data handles multiple interactions", {
  data <- simulate_gtheory_data(
    facets = list(person = 30, item = 10, rater = 3),
    vc = list(
      person = 1,
      item = 0.5,
      rater = 0.3,
      "person:item" = 0.2,
      "person:rater" = 0.15,
      "item:rater" = 0.1
    ),
    sd_residual = 1,
    seed = 1212
  )

  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 900)
  expect_true("y" %in% names(data))
})

test_that("simulate_gtheory_data variance components affect data variability", {
  data_high_vc <- simulate_gtheory_data(
    facets = list(person = 100, item = 10),
    vc = list(person = 2, item = 1),
    sd_residual = 0.5,
    seed = 1313
  )

  data_low_vc <- simulate_gtheory_data(
    facets = list(person = 100, item = 10),
    vc = list(person = 0.2, item = 0.1),
    sd_residual = 0.5,
    seed = 1414
  )

  var_person_high <- var(tapply(data_high_vc$y, data_high_vc$person, mean))
  var_person_low <- var(tapply(data_low_vc$y, data_low_vc$person, mean))

  expect_true(var_person_high > var_person_low)
})

test_that("prefix parameter affects level names", {
  data <- simulate_gtheory_data(
    facets = list(person = 10, item = 5),
    vc = list(person = 1, item = 0.5),
    prefix = "P",
    seed = 1515
  )

  expect_true(all(grepl("^P", levels(data$person))))
})
