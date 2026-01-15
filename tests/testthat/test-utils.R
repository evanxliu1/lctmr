# Tests for utility functions

test_that("assign_classes returns correct assignments", {
  prob_matrix <- matrix(c(
    0.9, 0.05, 0.05,
    0.1, 0.8, 0.1,
    0.1, 0.2, 0.7
  ), nrow = 3, byrow = TRUE)

  assignments <- assign_classes(prob_matrix)

  expect_equal(assignments, c(1L, 2L, 3L))
  expect_type(assignments, "integer")
})

test_that("assign_classes handles ties by choosing first", {
  prob_matrix <- matrix(c(
    0.5, 0.5,
    0.3, 0.7
  ), nrow = 2, byrow = TRUE)

  assignments <- assign_classes(prob_matrix)

  expect_equal(assignments[1], 1L)  # First column wins ties
  expect_equal(assignments[2], 2L)
})

test_that("assign_classes validates input", {
  expect_error(assign_classes("not numeric"), "must be numeric")
})

test_that("assign_classes handles data frames", {
  prob_df <- data.frame(
    prob1 = c(0.9, 0.1),
    prob2 = c(0.1, 0.9)
  )

  assignments <- assign_classes(prob_df)
  expect_equal(assignments, c(1L, 2L))
})

test_that("validate_lctm_data catches missing columns", {
  data <- data.frame(x = 1:5, y = 1:5)

  expect_error(
    validate_lctm_data(data, "weight", "time", "id"),
    "Missing required columns"
  )
})

test_that("validate_lctm_data catches non-numeric outcome", {
  data <- data.frame(
    id = 1:5,
    time = 1:5,
    weight = letters[1:5]
  )

  expect_error(
    validate_lctm_data(data, "weight", "time", "id"),
    "must be numeric"
  )
})

test_that("validate_lctm_data accepts valid data", {
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(0, 1, 0, 1),
    weight = c(5.0, 5.5, 6.0, 6.5)
  )

  result <- validate_lctm_data(data, "weight", "time", "id")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
})

test_that("validate_lctm_data warns about single observations", {
  data <- data.frame(
    id = c(1, 1, 2),  # Subject 2 has only 1 observation
    time = c(0, 1, 0),
    weight = c(5.0, 5.5, 6.0)
  )

  expect_warning(
    validate_lctm_data(data, "weight", "time", "id"),
    "fewer than 2 observations"
  )
})
