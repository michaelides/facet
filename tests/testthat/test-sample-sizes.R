# Tests for sample size calculations in G-studies

test_that("calculate_sample_size_info returns correct structure for fully crossed design", {
  # Fully crossed: each person sees each item with each rater
  crossed_data <- expand.grid(
    person = factor(1:10),
    item = factor(1:5),
    rater = factor(1:2)
  )
  crossed_data$score <- rnorm(nrow(crossed_data))

  f <- score ~ (1 | person) + (1 | item) + (1 | rater)
  ssi <- calculate_sample_size_info(f, crossed_data)

  expect_type(ssi, "list")
  expect_true("main" %in% names(ssi))
  expect_true("interactions" %in% names(ssi))
  expect_true("residual" %in% names(ssi))
  expect_true("nested" %in% names(ssi))

  expect_equal(length(ssi$main), 3)
  expect_equal(ssi$main["person"], c(person = 10))
  expect_equal(ssi$main["item"], c(item = 5))
  expect_equal(ssi$main["rater"], c(rater = 2))

  expect_equal(length(ssi$interactions), 0)

  # For fully crossed, residual should be all facets
  expect_equal(ssi$residual$facets, "person:item:rater")
  expect_equal(ssi$residual$n, 100)  # 10 * 5 * 2 = 100

  expect_null(ssi$nested)
})

test_that("calculate_sample_size_info handles interaction terms", {
  # Interaction design with rater:task
  data_with_interaction <- expand.grid(
    person = factor(1:10),
    rater = factor(1:5),
    task = factor(1:2)
  )
  data_with_interaction$score <- rnorm(nrow(data_with_interaction))

  f <- score ~ (1 | person) + (1 | rater:task)
  ssi <- calculate_sample_size_info(f, data_with_interaction)

  expect_equal(length(ssi$main), 1)
  expect_equal(ssi$main["person"], c(person = 10))

  expect_equal(length(ssi$interactions), 1)
  expect_true("rater:task" %in% names(ssi$interactions))
  expect_equal(ssi$interactions["rater:task"], c(`rater:task` = 10))  # 5 * 2 = 10
})

test_that("calculate_sample_size_info detects nested effects", {
  # Nested design: task is nested in rater
  # Rater 1 sees tasks 1-5, Rater 2 sees tasks 6-10
  nested_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:10, 10)),
    task = factor(c(rep(1:5, 10), rep(6:10, 10))),
    rater = factor(rep(1:2, each = 50))
  )

  f <- score ~ (1 | person) + (1 | task) + (1 | rater)
  ssi <- calculate_sample_size_info(f, nested_data)

  expect_true("task_in_rater" %in% names(ssi$nested))

  nested_info <- ssi$nested$task_in_rater
  expect_equal(nested_info$nested_facet, "task")
  expect_equal(nested_info$nesting_facet, "rater")
  expect_equal(nested_info$n_groups, 2)
  expect_equal(nested_info$mean_per_group, 5)
  expect_equal(nested_info$harmonic_mean_per_group, 5)
})

test_that("calculate_sample_size_info handles unbalanced nested designs", {
  # Unbalanced nested design
  # Rater 1 has tasks 1-5 (50 obs), Rater 2 has tasks 6-8 (30 obs)
  n_persons <- 10
  n_rater1_tasks <- 5
  n_rater2_tasks <- 3

  data_rater1 <- data.frame(
    person = factor(rep(1:n_persons, n_rater1_tasks)),
    task = factor(rep(1:n_rater1_tasks, each = n_persons)),
    rater = factor(1),
    score = rnorm(n_rater1_tasks * n_persons)
  )

  data_rater2 <- data.frame(
    person = factor(rep(1:n_persons, n_rater2_tasks)),
    task = factor(rep((n_rater1_tasks + 1):(n_rater1_tasks + n_rater2_tasks), each = n_persons)),
    rater = factor(2),
    score = rnorm(n_rater2_tasks * n_persons)
  )

  unbalanced_data <- rbind(data_rater1, data_rater2)

  f <- score ~ (1 | person) + (1 | task) + (1 | rater)
  ssi <- calculate_sample_size_info(f, unbalanced_data)

  expect_true("task_in_rater" %in% names(ssi$nested))

  nested_info <- ssi$nested$task_in_rater
  expect_equal(nested_info$n_groups, 2)
  expect_true(nested_info$mean_per_group != nested_info$harmonic_mean_per_group)
  expect_equal(nested_info$mean_per_group, 4)  # (5 + 3) / 2
})

test_that("detect_nesting_patterns correctly identifies nested facets", {
  # Nested: task in rater
  nested_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:10, 10)),
    task = factor(c(rep(1:5, 10), rep(6:10, 10))),
    rater = factor(rep(1:2, each = 50))
  )

  patterns <- detect_nesting_patterns(nested_data, c("person", "task", "rater"))

  expect_true("task_in_rater" %in% names(patterns))
  expect_false("rater_in_task" %in% names(patterns))
  expect_false(any(grepl("person_in", names(patterns))))
})

test_that("detect_nesting_patterns returns empty list for crossed design", {
  # Fully crossed design
  crossed_data <- expand.grid(
    person = factor(1:10),
    item = factor(1:5),
    rater = factor(1:2)
  )
  crossed_data$score <- rnorm(nrow(crossed_data))

  patterns <- detect_nesting_patterns(crossed_data, c("person", "item", "rater"))

  expect_equal(length(patterns), 0)
})

test_that("calculate_nested_sample_sizes computes correct statistics", {
  # Balanced nested data
  data <- data.frame(
    score = rnorm(50),
    rater = factor(rep(1:5, 10)),
    task = factor(c(rep(1:5, 5), rep(6:10, 5)))
  )

  nesting_info <- list(
    task_in_rater = list(nested_facet = "task", nesting_facet = "rater")
  )

  result <- calculate_nested_sample_sizes(data, nesting_info)

  expect_true("task_in_rater" %in% names(result))
  expect_equal(result$task_in_rater$nested_facet, "task")
  expect_equal(result$task_in_rater$nesting_facet, "rater")
  expect_equal(result$task_in_rater$n_groups, 5)
  expect_true(result$task_in_rater$mean_per_group > 0)
  expect_true(result$task_in_rater$harmonic_mean_per_group > 0)
})

test_that("override_nesting_detection works correctly", {
  detected <- list(
    task_in_rater = list(nested_facet = "task", nesting_facet = "rater")
  )

  user_nested <- list(task = "person")

  result <- override_nesting_detection(detected, user_nested, c("person", "task", "rater"))

  expect_true("task_in_person" %in% names(result))
  expect_false("task_in_rater" %in% names(result))
})

test_that("calculate_residual_sample_size computes correctly", {
  # Fully crossed data
  data <- expand.grid(
    person = factor(1:10),
    item = factor(1:5),
    rater = factor(1:2)
  )
  data$score <- rnorm(nrow(data))

  n <- calculate_residual_sample_size(data, "person:item:rater")
  expect_equal(n, 100)  # 10 * 5 * 2

  n2 <- calculate_residual_sample_size(data, "person:item")
  expect_equal(n2, 50)  # 10 * 5
})

test_that("gstudy integrates sample_size_info correctly", {
  # Fully crossed design
  test_data <- expand.grid(
    person = factor(1:10),
    item = factor(1:5),
    rater = factor(1:2)
  )
  test_data$score <- rnorm(nrow(test_data))

  g <- gstudy(score ~ (1 | person) + (1 | item) + (1 | rater), data = test_data)

  expect_true("sample_size_info" %in% names(g))
  expect_equal(length(g$sample_size_info$main), 3)
  expect_equal(g$sample_size_info$residual$facets, "person:item:rater")
})

test_that("gstudy handles nested parameter", {
  # Nested design with user-specified nesting
  nested_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:10, 10)),
    task = factor(c(rep(1:5, 10), rep(6:10, 10))),
    rater = factor(rep(1:2, each = 50))
  )

  g <- gstudy(score ~ (1 | person) + (1 | task) + (1 | rater),
              data = nested_data)

  expect_true("sample_size_info" %in% names(g))
  expect_true("task_in_rater" %in% names(g$sample_size_info$nested))
})

test_that("format_nested_info produces correct output", {
  nested_info <- list(
    task_in_rater = list(
      nested_facet = "task",
      nesting_facet = "rater",
      n_groups = 2,
      mean_per_group = 5,
      harmonic_mean_per_group = 4.5,
      counts_per_group = c(rater1 = 5, rater2 = 4)
    )
  )

  lines <- format_nested_info(nested_info, indent = 1)

  expect_true(any(grepl("task nested in rater", lines)))
  expect_true(any(grepl("raters: 2", lines)))
  expect_true(any(grepl("mean=5", lines)))
  expect_true(any(grepl("harmonic_mean=4", lines)))
})

test_that("calculate_sample_size_info works with single facet", {
  single_data <- data.frame(
    score = rnorm(50),
    person = factor(rep(1:10, 5))
  )

  f <- score ~ (1 | person)
  ssi <- calculate_sample_size_info(f, single_data)

  expect_equal(length(ssi$main), 1)
  expect_equal(ssi$main["person"], c(person = 10))
  expect_equal(ssi$residual$facets, "person")
  expect_equal(ssi$residual$n, 10)
})

test_that("gstudy print method shows sample size info", {
  # Nested design
  nested_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:10, 10)),
    task = factor(c(rep(1:5, 10), rep(6:10, 10))),
    rater = factor(rep(1:2, each = 50))
  )

  g <- gstudy(score ~ (1 | person) + (1 | task) + (1 | rater),
              data = nested_data)

  output <- capture.output(print(g))

  # Should show sample sizes section
  expect_true(any(grepl("Sample Sizes", output)))

  # Should show tibble format
  expect_true(any(grepl("tibble", output)))

  # Should show nested details
  expect_true(any(grepl("Nested details", output)))

  # Should show task per rater
  expect_true(any(grepl("Task per Rater", output)))
})

test_that("calculate_single_sample_size_tibble works", {
  crossed_data <- expand.grid(
    person = factor(1:10),
    item = factor(1:5),
    rater = factor(1:2)
  )
  crossed_data$score <- rnorm(nrow(crossed_data))

  f <- score ~ (1 | person) + (1 | item) + (1 | rater)
  ssi <- calculate_sample_size_info(f, crossed_data)
  ss_tibble <- calculate_single_sample_size_tibble(ssi)

  expect_s3_class(ss_tibble, "tbl_df")
  expect_equal(names(ss_tibble), c("effect", "type", "n"))
  expect_equal(nrow(ss_tibble), 4)  # 3 main + 1 residual
  expect_true("person" %in% ss_tibble$effect)
  expect_true("person:item:rater (res)" %in% ss_tibble$effect)
  expect_true("main" %in% ss_tibble$type)
  expect_true("residual" %in% ss_tibble$type)
})

test_that("format_nested_footnote works", {
  nested_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:10, 10)),
    task = factor(c(rep(1:5, 10), rep(6:10, 10))),
    rater = factor(rep(1:2, each = 50))
  )

  f <- score ~ (1 | person) + (1 | task) + (1 | rater)
  ssi <- calculate_sample_size_info(f, nested_data)

  footnote <- format_nested_footnote(ssi$nested)
  expect_type(footnote, "character")
  expect_true(length(footnote) > 0)
  expect_true(any(grepl("task per rater", tolower(footnote))))
  expect_true(any(grepl("Task per Rater", footnote)))
})

test_that("gstudy print uses tibble format", {
  test_data <- expand.grid(
    person = factor(1:10),
    item = factor(1:5),
    rater = factor(1:2)
  )
  test_data$score <- rnorm(nrow(test_data))

  g <- gstudy(score ~ (1 | person) + (1 | item) + (1 | rater),
              data = test_data)

  output <- capture.output(print(g))

  # Should have effect, type, n columns
  expect_true(any(grepl("effect", output, ignore.case = TRUE)))
  expect_true(any(grepl("type", output, ignore.case = TRUE)))
  
  # Should show tibble format
  expect_true(any(grepl("A tibble", output)))
})

test_that("gstudy with nested effects shows footnote", {
  nested_data <- data.frame(
    score = rnorm(100),
    person = factor(rep(1:10, 10)),
    task = factor(c(rep(1:5, 10), rep(6:10, 10))),
    rater = factor(rep(1:2, each = 50))
  )

  g <- gstudy(score ~ (1 | person) + (1 | task) + (1 | rater),
              data = nested_data)

  output <- capture.output(print(g))
  
  # Should show nested details section
  expect_true(any(grepl("Nested details", output)))
  
  # Should show Task per Rater
  expect_true(any(grepl("Task per Rater", output)))
})
