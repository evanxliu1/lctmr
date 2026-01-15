# Integration tests for lctm_fit

test_that("lctm_fit creates valid model object", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
                      k_range = 1:2, verbose = FALSE)

  # Fit Model F with K=2
  model <- lctm_fit(setup, k = 2, model = "F", verbose = FALSE)

  # Check class
  expect_s3_class(model, "lctm_model")

  # Check components
  expect_equal(model$k, 2L)
  expect_equal(model$model_type, "F")
  expect_type(model$bic, "double")
  expect_length(model$class_proportions, 2)

  # Check hlme model
  expect_s3_class(model$model, "hlme")
})

test_that("lctm_fit works with Model E", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
                      k_range = 1:2, verbose = FALSE)

  model_e <- lctm_fit(setup, k = 2, model = "E", verbose = FALSE)

  expect_equal(model_e$model_type, "E")
  expect_s3_class(model_e, "lctm_model")
})

test_that("lctm_fit works with linear option", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
                      k_range = 1:2, verbose = FALSE)

  model_lin <- lctm_fit(setup, k = 2, model = "F", linear = TRUE, verbose = FALSE)

  expect_true(model_lin$linear)
  expect_s3_class(model_lin, "lctm_model")
})

test_that("lctm_fit with K=1 works", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
                      k_range = 1:2, verbose = FALSE)

  model_k1 <- lctm_fit(setup, k = 1, model = "F", verbose = FALSE)

  expect_equal(model_k1$k, 1L)
  expect_equal(model_k1$class_proportions, c(Class1 = 1.0))
})

test_that("lctm_fit validates inputs", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
                      k_range = 1:2, verbose = FALSE)

  # Invalid model type
  expect_error(lctm_fit(setup, k = 2, model = "G"), "must be 'E' or 'F'")

  # Invalid k
  expect_error(lctm_fit(setup, k = 0), "positive integer")

  # Invalid setup
  expect_error(lctm_fit("not a setup", k = 2), "lctm_setup object")
})

test_that("lctm_fit print method works", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
                      k_range = 1:2, verbose = FALSE)
  model <- lctm_fit(setup, k = 2, model = "F", verbose = FALSE)

  expect_output(print(model), "LCTM Model")
  expect_output(print(model), "Class proportions")
})
