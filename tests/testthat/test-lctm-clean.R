# Tests for lctm_clean

test_that("lctm_clean returns valid lctm_cleaned object", {
  data("sample_growth", package = "lctmr")

  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  expect_s3_class(cleaned, "lctm_cleaned")
  expect_true(is.data.frame(cleaned$data))
  expect_equal(cleaned$outcome, "weight_raw")
  expect_equal(cleaned$time_var, "anthroage")
  expect_equal(cleaned$id_var, "childid")
  expect_null(cleaned$sex_var)
})

test_that("lctm_clean removes NA rows", {
  data("sample_growth", package = "lctmr")

  df <- sample_growth
  df$weight_raw[1:5] <- NA

  cleaned <- lctm_clean(df, "weight_raw", "anthroage", "childid", verbose = FALSE)

  # At least 5 rows removed (NAs), possibly more if a subject dropped below 2 obs
  expect_true(cleaned$n_removed >= 5)
  expect_true(nrow(cleaned$data) < nrow(sample_growth))
})

test_that("lctm_clean removes subjects with < 2 observations", {
  df <- data.frame(
    id = c(1, 1, 1, 2, 3, 3),
    time = c(0, 3, 6, 0, 0, 3),
    outcome = c(5, 6, 7, 5, 4, 5)
  )

  expect_warning(
    cleaned <- lctm_clean(df, "outcome", "time", "id", verbose = FALSE),
    "1 subject"
  )
  expect_false(2 %in% cleaned$data$id)
  expect_equal(length(unique(cleaned$data$id)), 2)
})

test_that("lctm_clean applies WHO z-score cutoffs", {
  data("sample_growth", package = "lctmr")

  # Inject an extreme z-score
  df <- sample_growth
  df$waz[1] <- -7  # below WHO cutoff of -6

  cleaned <- lctm_clean(df, "weight_raw", "anthroage", "childid",
                         type = "weight", standards = "WHO",
                         zscore_col = "waz", verbose = FALSE)

  expect_true(nrow(cleaned$data) < nrow(df))
})

test_that("lctm_clean applies CDC z-score cutoffs", {
  data("sample_growth", package = "lctmr")

  df <- sample_growth
  df$waz[1] <- -6  # below CDC cutoff of -5 but within WHO

  cleaned_cdc <- lctm_clean(df, "weight_raw", "anthroage", "childid",
                             type = "weight", standards = "CDC",
                             zscore_col = "waz", verbose = FALSE)

  cleaned_who <- lctm_clean(df, "weight_raw", "anthroage", "childid",
                             type = "weight", standards = "WHO",
                             zscore_col = "waz", verbose = FALSE)

  # CDC should flag this (cutoff -5), WHO should not (cutoff -6)
  expect_true(nrow(cleaned_cdc$data) < nrow(cleaned_who$data))
})

test_that("lctm_clean applies birth weight criteria", {
  df <- data.frame(
    id = rep(1:5, each = 3),
    time = rep(c(0, 3, 6), 5),
    outcome = rnorm(15, 7, 1),
    bw = rep(c(400, 3000, 3500, 5100, 3200), each = 3)
  )

  cleaned <- lctm_clean(df, "outcome", "time", "id",
                         birth_weight_col = "bw", verbose = FALSE)

  # Subjects 1 (400g) and 4 (5100g) should have rows removed
  expect_true(nrow(cleaned$data) < nrow(df))
})

test_that("lctm_clean warns when check_decrease used without height/hc type", {
  data("sample_growth", package = "lctmr")

  expect_warning(
    lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
               check_decrease = TRUE, verbose = FALSE),
    "check_decrease"
  )
})

test_that("lctm_clean stores sex_var", {
  data("sample_growth", package = "lctmr")

  # sample_growth doesn't have a sex column, so create one
  df <- sample_growth
  df$sex <- sample(c("M", "F"), nrow(df), replace = TRUE)

  cleaned <- lctm_clean(df, "weight_raw", "anthroage", "childid",
                         sex_var = "sex", verbose = FALSE)

  expect_equal(cleaned$sex_var, "sex")
})

test_that("lctm_clean validates inputs", {
  data("sample_growth", package = "lctmr")

  expect_error(
    lctm_clean(sample_growth, "nonexistent", "anthroage", "childid"),
    "Missing required"
  )

  expect_error(
    lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
               sex_var = "nonexistent"),
    "not found"
  )

  expect_error(
    lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
               type = "weight"),
    "zscore_col is required"
  )
})

test_that("print.lctm_cleaned works", {
  data("sample_growth", package = "lctmr")

  cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
                         verbose = FALSE)

  expect_output(print(cleaned), "LCTM Cleaned Data")
  expect_output(print(cleaned), "Ready for lctm_initial")
})
