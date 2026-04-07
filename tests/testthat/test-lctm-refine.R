# Tests for lctm_refine

test_that("lctm_refine returns valid lctm_result object", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)
  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  result <- lctm_refine(init, random = ~ 1 + anthroage,
                         k_range = 2:3, verbose = FALSE)

  expect_s3_class(result, "lctm_result")
  expect_true("best_model" %in% names(result))
  expect_true("bic_table" %in% names(result))
  expect_true("search_history" %in% names(result))
  expect_true("all_models" %in% names(result))
})

test_that("lctm_refine records search history correctly", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)
  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  result <- lctm_refine(init, random = ~ 1 + anthroage,
                         k_range = 2:3, models = c("B", "A"),
                         verbose = FALSE)

  # 2 K values x 2 models = up to 4 entries (may stop early if adequate found)
  expect_s3_class(result$search_history, "data.frame")
  expect_true(nrow(result$search_history) >= 1)
  expect_true(all(c("k", "model_type", "bic", "overall_pass") %in%
                    names(result$search_history)))
})

test_that("lctm_refine finds adequate model with lenient thresholds", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)
  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  result <- lctm_refine(init, random = ~ 1 + anthroage,
                         k_range = 2:3,
                         adequacy_thresholds = list(appa = 0.4, occ = 1.0,
                                                    entropy = 0.2),
                         verbose = FALSE)

  expect_false(is.null(result$best_model))
  expect_true(result$adequacy$overall_pass)
  expect_s3_class(result$class_assignments, "data.frame")
})

test_that("lctm_refine inherits degree and knots from initial", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  init <- lctm_initial(cleaned, k = 2, degree = 1, verbose = FALSE)

  # Should inherit degree = 1 from init
  result <- lctm_refine(init, k_range = 2, verbose = FALSE)

  expect_s3_class(result, "lctm_result")
})

test_that("lctm_refine allows degree override", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  # Override to linear
  result <- lctm_refine(init, degree = 1,
                         random = ~ 1 + anthroage,
                         k_range = 2, verbose = FALSE)

  expect_s3_class(result, "lctm_result")
})

test_that("lctm_refine validates covariates", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)
  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  expect_error(
    lctm_refine(init, covariates = c("nonexistent_var"), k_range = 2,
                verbose = FALSE),
    "Covariates not found"
  )
})

test_that("lctm_refine works with covariates", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)
  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  result <- lctm_refine(init, random = ~ 1 + anthroage,
                         covariates = c("waz"),
                         k_range = 2, verbose = FALSE)

  expect_s3_class(result, "lctm_result")
})

test_that("lctm_refine start_simple option works", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)
  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  result <- lctm_refine(init, random = ~ 1 + anthroage,
                         k_range = 2, start_simple = TRUE,
                         verbose = FALSE)

  expect_s3_class(result, "lctm_result")
})

test_that("lctm_refine validates inputs", {
  expect_error(
    lctm_refine("not_an_initial"),
    "lctm_initial object"
  )
})

test_that("lctm_refine BIC table has correct structure", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)
  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  result <- lctm_refine(init, random = ~ 1 + anthroage,
                         k_range = 2:3, verbose = FALSE)

  expect_true("k" %in% names(result$bic_table))
  expect_true("bic" %in% names(result$bic_table))
  expect_true("converged" %in% names(result$bic_table))
})
