# End-to-end tests for lctm_auto

test_that("lctm_auto returns valid result object", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  # Run with minimal search for speed
  result <- lctm_auto(
    data = sample_growth,
    outcome = "weight_raw",
    time_var = "anthroage",
    id_var = "childid",
    k_range = 2:3,
    models = "F",
    try_splines = FALSE,
    try_linear = FALSE,
    verbose = FALSE
  )

  # Check class
 expect_s3_class(result, "lctm_result")

  # Check components exist
  expect_true("best_model" %in% names(result))
  expect_true("bic_table" %in% names(result))
  expect_true("search_history" %in% names(result))
  expect_true("all_models" %in% names(result))
})

test_that("lctm_auto finds adequate model when one exists", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  # Use lenient thresholds to ensure we find an adequate model
  result <- lctm_auto(
    data = sample_growth,
    outcome = "weight_raw",
    time_var = "anthroage",
    id_var = "childid",
    k_range = 2:3,
    models = "F",
    adequacy_thresholds = list(appa = 0.5, occ = 2.0, entropy = 0.3),
    try_splines = FALSE,
    try_linear = FALSE,
    verbose = FALSE
  )

  # Should find an adequate model with lenient thresholds
  expect_false(is.null(result$best_model))
  expect_true(result$adequacy$overall_pass)
})

test_that("lctm_auto records search history", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  result <- lctm_auto(
    data = sample_growth,
    outcome = "weight_raw",
    time_var = "anthroage",
    id_var = "childid",
    k_range = 2:3,
    models = c("F", "E"),
    try_splines = FALSE,
    try_linear = FALSE,
    verbose = FALSE
  )

  # Check search history
  expect_s3_class(result$search_history, "data.frame")
  expect_true(nrow(result$search_history) >= 1)
  expect_true("k" %in% names(result$search_history))
  expect_true("model_type" %in% names(result$search_history))
  expect_true("overall_pass" %in% names(result$search_history))
})

test_that("lctm_auto tries splines when enabled", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  result <- lctm_auto(
    data = sample_growth,
    outcome = "weight_raw",
    time_var = "anthroage",
    id_var = "childid",
    k_range = 2,
    models = "F",
    try_splines = TRUE,
    try_linear = FALSE,
    verbose = FALSE
  )

  # Should have tried splines
  expect_true(any(result$search_history$splines))
})

test_that("lctm_auto tries linear when enabled", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  result <- lctm_auto(
    data = sample_growth,
    outcome = "weight_raw",
    time_var = "anthroage",
    id_var = "childid",
    k_range = 2,
    models = "F",
    try_splines = FALSE,
    try_linear = TRUE,
    verbose = FALSE
  )

  # Should have tried linear
  expect_true(any(result$search_history$linear))
})

test_that("lctm_auto print method works", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  result <- lctm_auto(
    data = sample_growth,
    outcome = "weight_raw",
    time_var = "anthroage",
    id_var = "childid",
    k_range = 2:3,
    models = "F",
    try_splines = FALSE,
    try_linear = FALSE,
    verbose = FALSE
  )

  expect_output(print(result), "LCTM Auto Result")
})

test_that("lctm_auto returns class assignments when model found", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  result <- lctm_auto(
    data = sample_growth,
    outcome = "weight_raw",
    time_var = "anthroage",
    id_var = "childid",
    k_range = 2:3,
    models = "F",
    adequacy_thresholds = list(appa = 0.5, occ = 2.0, entropy = 0.3),
    try_splines = FALSE,
    try_linear = FALSE,
    verbose = FALSE
  )

  if (!is.null(result$best_model)) {
    expect_s3_class(result$class_assignments, "data.frame")
    expect_true("id" %in% names(result$class_assignments))
    expect_true("class" %in% names(result$class_assignments))
  }
})
