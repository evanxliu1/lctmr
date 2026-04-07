# Integration tests for lctm_adequacy

test_that("lctm_adequacy creates valid adequacy object", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
                      k_range = 1:2, verbose = FALSE)
  model <- lctm_fit(setup, k = 2, model = "B", verbose = FALSE)

  adequacy <- lctm_adequacy(model)

  # Check class
  expect_s3_class(adequacy, "lctm_adequacy")

  # Check components
  expect_type(adequacy$appa, "double")
  expect_type(adequacy$occ, "double")
  expect_type(adequacy$entropy, "double")
  expect_type(adequacy$overall_pass, "logical")

  # Check length matches K
  expect_length(adequacy$appa, 2)
  expect_length(adequacy$occ, 2)
})

test_that("lctm_adequacy works with K=1", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
                      k_range = 1:2, verbose = FALSE)
  model <- lctm_fit(setup, k = 1, model = "B", verbose = FALSE)

  adequacy <- lctm_adequacy(model)

  # K=1 should always pass
  expect_true(adequacy$overall_pass)
  expect_equal(adequacy$appa, c(Class1 = 1.0))
  expect_equal(adequacy$entropy, 1.0)
})

test_that("lctm_adequacy uses custom thresholds", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
                      k_range = 1:2, verbose = FALSE)
  model <- lctm_fit(setup, k = 2, model = "B", verbose = FALSE)

  # Very strict thresholds
  adequacy_strict <- lctm_adequacy(model,
                                   thresholds = list(appa = 0.99, occ = 100, entropy = 0.99))

  # Very lenient thresholds
  adequacy_lenient <- lctm_adequacy(model,
                                    thresholds = list(appa = 0.5, occ = 1, entropy = 0.1))

  # Lenient should be more likely to pass
  expect_equal(adequacy_strict$thresholds$appa, 0.99)
  expect_equal(adequacy_lenient$thresholds$appa, 0.5)
})

test_that("lctm_adequacy print method works", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid",
                      k_range = 1:2, verbose = FALSE)
  model <- lctm_fit(setup, k = 2, model = "B", verbose = FALSE)
  adequacy <- lctm_adequacy(model)

  expect_output(print(adequacy), "LCTM Model Adequacy")
  expect_output(print(adequacy), "APPA")
  expect_output(print(adequacy), "OCC")
  expect_output(print(adequacy), "Entropy")
})

test_that("lctm_adequacy validates input", {
  expect_error(lctm_adequacy("not a model"), "lctm_model object")
})
