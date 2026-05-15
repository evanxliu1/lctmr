# Integration tests for lctm_adequacy

# Helper: fit a small refine result for adequacy tests
.fit_test_result <- function(k_range = 2:2) {
  data("sample_growth", package = "lctmr")
  init <- lctm_initial(sample_growth, outcome = "weight_raw",
                       time_var = "anthroage", id_var = "childid",
                       k = 2, degree = 2, verbose = FALSE)
  lctm_refine(init, k_range = k_range, models = "B", verbose = FALSE)
}

test_that("lctm_adequacy creates valid adequacy object", {
  skip_on_cran()

  result <- .fit_test_result()
  skip_if(is.null(result$best_model), "no adequate model found in test data")

  adequacy <- lctm_adequacy(result$best_model)

  expect_s3_class(adequacy, "lctm_adequacy")
  expect_type(adequacy$appa, "double")
  expect_type(adequacy$occ, "double")
  expect_type(adequacy$entropy, "double")
  expect_type(adequacy$overall_pass, "logical")
  expect_length(adequacy$appa, 2)
  expect_length(adequacy$occ, 2)
})

test_that("lctm_adequacy uses custom thresholds", {
  skip_on_cran()

  result <- .fit_test_result()
  skip_if(is.null(result$best_model), "no adequate model found in test data")

  adequacy_strict <- lctm_adequacy(
    result$best_model,
    thresholds = list(appa = 0.99, occ = 100, entropy = 0.99)
  )
  adequacy_lenient <- lctm_adequacy(
    result$best_model,
    thresholds = list(appa = 0.5, occ = 1, entropy = 0.1)
  )

  expect_equal(adequacy_strict$thresholds$appa, 0.99)
  expect_equal(adequacy_lenient$thresholds$appa, 0.5)
})

test_that("lctm_adequacy print method works", {
  skip_on_cran()

  result <- .fit_test_result()
  skip_if(is.null(result$best_model), "no adequate model found in test data")
  adequacy <- lctm_adequacy(result$best_model)

  expect_output(print(adequacy), "LCTM Model Adequacy")
  expect_output(print(adequacy), "APPA")
  expect_output(print(adequacy), "OCC")
  expect_output(print(adequacy), "Entropy")
})

test_that("lctm_adequacy validates input", {
  expect_error(lctm_adequacy("not a model"), "lctm_model object")
})
