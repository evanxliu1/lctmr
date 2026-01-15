# Integration tests for lctm_setup

test_that("lctm_setup creates valid setup object", {
  skip_on_cran()

  # Load sample data
  data("sample_growth", package = "lctmr")

  # Run setup with minimal K range for speed
  setup <- lctm_setup(
    data = sample_growth,
    outcome = "weight_raw",
    time_var = "anthroage",
    id_var = "childid",
    k_range = 1:3,
    verbose = FALSE
  )

  # Check class
  expect_s3_class(setup, "lctm_setup")

  # Check components
  expect_equal(setup$outcome, "weight_raw")
  expect_equal(setup$time_var, "anthroage")
  expect_equal(setup$id_var, "childid")
  expect_equal(setup$degree, 2L)

  # Check BIC table
  expect_s3_class(setup$bic_table, "data.frame")
  expect_true("k" %in% names(setup$bic_table))
  expect_true("bic" %in% names(setup$bic_table))
  expect_equal(nrow(setup$bic_table), 3)

  # Check k_ranking
  expect_type(setup$k_ranking, "integer")
  expect_true(length(setup$k_ranking) <= 3)

  # Check base model
  expect_s3_class(setup$base_model, "hlme")
})

test_that("lctm_setup validates inputs", {
  data("sample_growth", package = "lctmr")

  # Missing column
  expect_error(
    lctm_setup(sample_growth, "nonexistent", "anthroage", "childid", verbose = FALSE),
    "Missing required columns"
  )

  # Invalid k_range
  expect_error(
    lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
               k_range = -1, verbose = FALSE),
    "positive integers"
  )

  # Invalid degree
  expect_error(
    lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
               degree = 3, verbose = FALSE),
    "degree must be 1"
  )
})

test_that("lctm_setup print method works", {
  skip_on_cran()
  data("sample_growth", package = "lctmr")

  setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
                      k_range = 1:2, verbose = FALSE)

  expect_output(print(setup), "LCTM Setup")
  expect_output(print(setup), "BIC Comparison")
})
