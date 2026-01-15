# Tests for adequacy metric calculations

test_that("calc_appa returns correct values for perfect separation", {
  # Perfect separation: each observation clearly belongs to one class
  prob_matrix <- matrix(c(
    1.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    0.0, 0.0, 1.0
  ), nrow = 6, byrow = TRUE)

  appa <- calc_appa(prob_matrix)

  expect_length(appa, 3)
  expect_equal(appa, c(Class1 = 1.0, Class2 = 1.0, Class3 = 1.0))
})

test_that("calc_appa returns correct values for imperfect separation", {
  prob_matrix <- matrix(c(
    0.9, 0.05, 0.05,
    0.8, 0.15, 0.05,
    0.1, 0.85, 0.05,
    0.05, 0.9, 0.05
  ), nrow = 4, byrow = TRUE)

  appa <- calc_appa(prob_matrix)

  expect_length(appa, 3)
  expect_equal(appa["Class1"], c(Class1 = 0.85))  # mean(0.9, 0.8)
  expect_equal(appa["Class2"], c(Class2 = 0.875)) # mean(0.85, 0.9)
  expect_true(is.na(appa["Class3"]))              # No assignments to class 3
})

test_that("calc_appa validates input", {
  expect_error(calc_appa("not a matrix"), "must be numeric")
  expect_error(calc_appa(matrix(c(-0.1, 1.1), nrow = 1)),
               "Probabilities must be between 0 and 1")
})

test_that("calc_occ returns values > 1 for good classification", {
  # Good separation
  prob_matrix <- matrix(c(
    0.95, 0.05,
    0.90, 0.10,
    0.10, 0.90,
    0.05, 0.95
  ), nrow = 4, byrow = TRUE)

  occ <- calc_occ(prob_matrix)

  expect_length(occ, 2)
  expect_true(all(occ > 5))  # Should be well above threshold
})

test_that("calc_occ handles custom class proportions", {
  prob_matrix <- matrix(c(
    0.9, 0.1,
    0.1, 0.9
  ), nrow = 2, byrow = TRUE)

  # Custom proportions
  occ <- calc_occ(prob_matrix, class_proportions = c(0.5, 0.5))

  expect_length(occ, 2)
  expect_true(all(is.finite(occ)))
})

test_that("calc_relative_entropy returns value in [0, 1]", {
  # High entropy (clear separation)
  prob_high <- matrix(c(
    0.99, 0.01,
    0.98, 0.02,
    0.02, 0.98,
    0.01, 0.99
  ), nrow = 4, byrow = TRUE)

  entropy_high <- calc_relative_entropy(prob_high)
  expect_true(entropy_high > 0.8)
  expect_true(entropy_high <= 1.0)

  # Lower entropy (less clear)
  prob_low <- matrix(c(
    0.6, 0.4,
    0.55, 0.45,
    0.45, 0.55,
    0.4, 0.6
  ), nrow = 4, byrow = TRUE)

  entropy_low <- calc_relative_entropy(prob_low)
  expect_true(entropy_low >= 0)
  expect_true(entropy_low < entropy_high)
})

test_that("calc_relative_entropy handles edge cases", {
  # Matrix with zeros
  prob_matrix <- matrix(c(
    1.0, 0.0,
    0.0, 1.0
  ), nrow = 2, byrow = TRUE)

  entropy <- calc_relative_entropy(prob_matrix)
  expect_true(is.finite(entropy))
  expect_equal(entropy, 1.0)
})

test_that("calc_relative_entropy requires at least 2 classes", {
  prob_matrix <- matrix(c(1, 1, 1), nrow = 3)
  expect_error(calc_relative_entropy(prob_matrix), "at least 2 classes")
})
