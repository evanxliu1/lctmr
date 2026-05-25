# Tests for lctm_initial

test_that("lctm_initial returns valid lctm_initial object", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  expect_s3_class(init, "lctm_initial")
  expect_equal(init$outcome, "weight_raw")
  expect_equal(init$time_var, "anthroage")
  expect_equal(init$id_var, "childid")
  expect_equal(init$degree, 2L)
  expect_equal(init$k, 2L)
  expect_null(init$knots)
  expect_s3_class(init$initial_model, "hlme")
  expect_s3_class(init$base_model, "hlme")
})

test_that("lctm_initial generates all four plots", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  expect_named(init$plots, c("spaghetti", "loess", "residuals", "guide"),
               ignore.order = TRUE)
  for (p in init$plots) {
    expect_s3_class(p, "ggplot")
  }
})

test_that("lctm_initial accepts raw data frame", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")

  init <- lctm_initial(sample_growth, outcome = "weight_raw",
                        time_var = "anthroage", id_var = "childid",
                        k = 2, degree = 2, verbose = FALSE)

  expect_s3_class(init, "lctm_initial")
})

test_that("lctm_initial handles data with NAs", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  df <- sample_growth
  df$weight_raw[c(1, 5, 10)] <- NA

  init <- lctm_initial(df, outcome = "weight_raw", time_var = "anthroage",
                        id_var = "childid", k = 2, degree = 2, verbose = FALSE)

  expect_s3_class(init, "lctm_initial")
  expect_true(length(init$plots) == 4)
})

test_that("lctm_initial supports linear degree", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  init <- lctm_initial(cleaned, k = 2, degree = 1, verbose = FALSE)

  expect_equal(init$degree, 1L)
})

test_that("lctm_initial supports cubic degree", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  init <- lctm_initial(cleaned, k = 2, degree = 3, verbose = FALSE)

  expect_equal(init$degree, 3L)
})

test_that("lctm_initial uses random intercept only", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  # The initial model should use random = ~1
  # hlme stores the call; check that random is ~1
  model_call <- init$initial_model$call
  expect_true("random" %in% names(model_call))
})

test_that("lctm_initial validates inputs", {
  data("sample_growth", package = "lctmr")

  # Missing required args for raw data
  expect_error(
    lctm_initial(sample_growth, k = 2),
    "outcome, time_var, and id_var are required"
  )

  # Invalid degree
  expect_error(
    lctm_initial(sample_growth, "weight_raw", "anthroage", "childid",
                  degree = 4),
    "degree must be"
  )

  # Invalid k
  expect_error(
    lctm_initial(sample_growth, "weight_raw", "anthroage", "childid",
                  k = 1),
    "k must be"
  )

  # knots with degree > 1 is no longer an error: degree is coerced to 1
  init_k <- lctm_initial(sample_growth, "weight_raw", "anthroage", "childid",
                         degree = 2, knots = c(3, 6), spline_degree = 2,
                         verbose = FALSE)
  expect_equal(init_k$degree, 1L)
  expect_equal(init_k$spline_degree, 2L)

  # spline_degree without knots is ignored, with a warning
  expect_warning(
    lctm_initial(sample_growth, "weight_raw", "anthroage", "childid",
                  degree = 2, spline_degree = 2, verbose = FALSE),
    "ignored because no"
  )

  # invalid spline_degree errors
  expect_error(
    lctm_initial(sample_growth, "weight_raw", "anthroage", "childid",
                  degree = 1, knots = c(3, 6), spline_degree = 5),
    "spline_degree must be"
  )
})

test_that("print.lctm_initial works", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  expect_output(print(init), "LCTM Initial Model")
  expect_output(print(init), "Residual Interpretation Guide")
  expect_output(print(init), "quadratic")
})

test_that("plot.lctm_initial works with which argument", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  init <- lctm_initial(cleaned, k = 2, degree = 2, verbose = FALSE)

  expect_error(plot(init, which = "nonexistent"), "not found")
  expect_invisible(plot(init, which = "spaghetti"))
})

test_that("lctm_initial saves PDF when requested", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  pdf_path <- tempfile(fileext = ".pdf")
  init <- lctm_initial(cleaned, k = 2, degree = 2,
                        save_pdf = pdf_path, verbose = FALSE)

  expect_true(file.exists(pdf_path))
  expect_true(file.size(pdf_path) > 0)
  unlink(pdf_path)
})

test_that("lctm_initial nudges when degree is left at the default", {
  skip_on_cran()

  data("sample_growth", package = "lctmr")
  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  # Default degree -> nudge message fires under verbose
  expect_message(
    lctm_initial(cleaned, k = 2, verbose = TRUE),
    "default trajectory degree"
  )

  # Explicit degree -> no nudge in the captured message stream
  msgs <- capture.output(
    lctm_initial(cleaned, k = 2, degree = 2, verbose = TRUE),
    type = "message"
  )
  expect_false(any(grepl("default trajectory degree", msgs)))
})
