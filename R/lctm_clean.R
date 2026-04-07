#' Clean Data for LCTM Analysis
#'
#' Prepares and validates longitudinal data for the LCTM workflow by removing
#' missing values, applying population-level outlier criteria based on growth
#' standards, and removing subjects with insufficient observations.
#'
#' @param data A data frame containing longitudinal data.
#' @param outcome Character string naming the outcome variable.
#' @param time_var Character string naming the time variable.
#' @param id_var Character string naming the subject ID variable.
#' @param sex_var Character string naming the sex variable (optional).
#' @param type Character string specifying outcome type for outlier detection:
#'   `"weight"`, `"height"`, or `"hc"` (head circumference). If NULL (default),
#'   no growth-standard-based outlier detection is performed.
#' @param standards Character string specifying which growth chart cutoffs to
#'   apply: `"WHO"`, `"CDC"`, or `"both"` (default `"WHO"`). Only used when
#'   `type` is specified.
#' @param zscore_col Character string naming the column containing pre-computed
#'   z-scores for population-level outlier detection. Required when `type` is
#'   specified. Users must provide z-scores computed externally (e.g., using the
#'   `zscorer` or `anthro` packages).
#' @param zscore_col_cdc Character string naming a second z-score column for
#'   CDC standards. Required when `standards = "both"`. If `standards = "CDC"`,
#'   uses `zscore_col` for CDC z-scores.
#' @param birth_weight_col Character string naming the birth weight column (in
#'   grams) for applying birth weight criteria (<=500g or >=5000g). If NULL
#'   (default), birth weight outlier check is skipped.
#' @param check_decrease Logical; if TRUE and `type` is `"height"` or `"hc"`,
#'   flags observations where height or head circumference decreases by more
#'   than 3 cm between consecutive visits for the same subject. Default FALSE.
#' @param verbose Logical; if TRUE, prints cleaning summary.
#'
#' @return An `lctm_cleaned` object containing:
#' \describe{
#'   \item{data}{The cleaned data frame}
#'   \item{outcome}{Name of outcome variable}
#'   \item{time_var}{Name of time variable}
#'   \item{id_var}{Name of ID variable}
#'   \item{sex_var}{Name of sex variable (or NULL)}
#'   \item{n_removed}{Number of rows removed}
#'   \item{subjects_removed}{Number of subjects removed}
#' }
#'
#' @details
#' This is an optional first step of the LCTM workflow:
#'
#' \code{[lctm_clean()] -> lctm_initial() -> lctm_refine()}
#'
#' For special populations (e.g., extremely low birth weight or premature
#' infants), you may skip this function and pass raw data directly to
#' [lctm_initial()].
#'
#' **Cleaning steps:**
#' 1. Validates that required columns exist and are correct types
#' 2. Removes rows with NA in outcome, time, or ID variables
#' 3. If `type` is specified, applies population-level z-score cutoffs:
#'    - **WHO**: WAZ < -6 or > +5, LAZ < -6 or > +6, WHZ < -5 or > +5
#'    - **CDC**: WeightMZ < -5 or > +5, HeightMZ < -5 or > +3
#' 4. If `birth_weight_col` is specified, removes birth weights <=500g or >=5000g
#' 5. If `check_decrease = TRUE` (height/hc only), removes observations where
#'    the measure decreases by more than 3 cm
#' 6. Removes subjects with fewer than 2 remaining observations
#'
#' @examples
#' \dontrun{
#' data(sample_growth)
#'
#' # Basic cleaning (no growth standard cutoffs)
#' cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid")
#'
#' # With WHO z-score cutoffs for weight
#' cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid",
#'                        type = "weight", standards = "WHO", zscore_col = "waz")
#' }
#'
#' @export
lctm_clean <- function(data, outcome, time_var, id_var,
                       sex_var = NULL,
                       type = NULL,
                       standards = "WHO",
                       zscore_col = NULL,
                       zscore_col_cdc = NULL,
                       birth_weight_col = NULL,
                       check_decrease = FALSE,
                       verbose = TRUE) {

  # Validate basic inputs
  data <- validate_lctm_data(data, outcome, time_var, id_var)

  # Validate sex_var if provided
  if (!is.null(sex_var)) {
    if (!sex_var %in% names(data)) {
      stop("sex_var '", sex_var, "' not found in data", call. = FALSE)
    }
  }

  # Validate type
  if (!is.null(type)) {
    type <- match.arg(type, c("weight", "height", "hc"))
    if (!standards %in% c("WHO", "CDC", "both")) {
      stop("standards must be 'WHO', 'CDC', or 'both'", call. = FALSE)
    }
    if (is.null(zscore_col)) {
      stop("zscore_col is required when type is specified", call. = FALSE)
    }
    if (!zscore_col %in% names(data)) {
      stop("zscore_col '", zscore_col, "' not found in data", call. = FALSE)
    }
    if (standards == "both" && is.null(zscore_col_cdc)) {
      stop("zscore_col_cdc is required when standards = 'both'", call. = FALSE)
    }
    if (!is.null(zscore_col_cdc) && !zscore_col_cdc %in% names(data)) {
      stop("zscore_col_cdc '", zscore_col_cdc, "' not found in data", call. = FALSE)
    }
  }

  if (check_decrease && (is.null(type) || !type %in% c("height", "hc"))) {
    warning("check_decrease is only applied when type is 'height' or 'hc'. Ignoring.",
            call. = FALSE)
  }

  if (!is.null(birth_weight_col) && !birth_weight_col %in% names(data)) {
    stop("birth_weight_col '", birth_weight_col, "' not found in data",
         call. = FALSE)
  }

  n_orig <- nrow(data)
  n_subj_orig <- length(unique(data[[id_var]]))

  # Step 1: Remove rows with NA in key variables
  key_vars <- c(outcome, time_var, id_var)
  complete_rows <- complete.cases(data[, key_vars])
  data <- data[complete_rows, , drop = FALSE]
  n_na_removed <- n_orig - nrow(data)

  # Step 2: Population-level z-score outlier detection
  n_zscore_flagged <- 0
  if (!is.null(type)) {
    # Define cutoffs based on type and standards
    flag <- rep(FALSE, nrow(data))

    if (standards %in% c("WHO", "both")) {
      z <- data[[zscore_col]]
      who_flag <- switch(type,
        weight = !is.na(z) & (z < -6 | z > 5),
        height = !is.na(z) & (z < -6 | z > 6),
        hc     = !is.na(z) & (z < -6 | z > 6)
      )
      flag <- flag | who_flag
    }

    if (standards %in% c("CDC", "both")) {
      z_cdc <- if (standards == "both") data[[zscore_col_cdc]] else data[[zscore_col]]
      cdc_flag <- switch(type,
        weight = !is.na(z_cdc) & (z_cdc < -5 | z_cdc > 5),
        height = !is.na(z_cdc) & (z_cdc < -5 | z_cdc > 3),
        hc     = !is.na(z_cdc) & (z_cdc < -5 | z_cdc > 5)
      )
      flag <- flag | cdc_flag
    }

    n_zscore_flagged <- sum(flag)
    if (n_zscore_flagged > 0) {
      data <- data[!flag, , drop = FALSE]
      if (verbose) {
        message("Flagged ", n_zscore_flagged,
                " observations as population-level z-score outliers (",
                standards, " ", type, " criteria)")
      }
    }
  }

  # Step 3: Birth weight criteria
  n_bw_flagged <- 0
  if (!is.null(birth_weight_col)) {
    bw <- data[[birth_weight_col]]
    bw_flag <- !is.na(bw) & (bw <= 500 | bw >= 5000)
    n_bw_flagged <- sum(bw_flag)
    if (n_bw_flagged > 0) {
      data <- data[!bw_flag, , drop = FALSE]
      if (verbose) {
        message("Flagged ", n_bw_flagged,
                " observations with birth weight <=500g or >=5000g")
      }
    }
  }

  # Step 4: Height/HC decrease check (>3 cm decrease between consecutive visits)
  n_decrease_flagged <- 0
  if (check_decrease && !is.null(type) && type %in% c("height", "hc")) {
    # Sort by subject and time
    data <- data[order(data[[id_var]], data[[time_var]]), ]
    flag_decrease <- rep(FALSE, nrow(data))

    subjects <- unique(data[[id_var]])
    for (subj in subjects) {
      idx <- which(data[[id_var]] == subj)
      if (length(idx) < 2) next
      vals <- data[[outcome]][idx]
      diffs <- diff(vals)
      # Flag observations where decrease > 3 cm
      bad <- which(diffs < -3) + 1  # the second observation in each bad pair
      if (length(bad) > 0) {
        flag_decrease[idx[bad]] <- TRUE
      }
    }

    n_decrease_flagged <- sum(flag_decrease)
    if (n_decrease_flagged > 0) {
      data <- data[!flag_decrease, , drop = FALSE]
      if (verbose) {
        message("Flagged ", n_decrease_flagged,
                " observations with ", type, " decrease > 3 cm")
      }
    }
  }

  # Step 5: Remove subjects with < 2 observations
  obs_per_subject <- table(data[[id_var]])
  subjects_to_remove <- names(obs_per_subject[obs_per_subject < 2])

  if (length(subjects_to_remove) > 0) {
    warning(length(subjects_to_remove),
            " subject(s) removed for having fewer than 2 observations",
            call. = FALSE)
    data <- data[!data[[id_var]] %in% subjects_to_remove, , drop = FALSE]
  }

  n_removed <- n_orig - nrow(data)
  subjects_removed <- n_subj_orig - length(unique(data[[id_var]]))

  if (verbose) {
    message("=== Data Cleaning Summary ===")
    message("Original: ", n_orig, " rows, ", n_subj_orig, " subjects")
    if (n_na_removed > 0) {
      message("Removed ", n_na_removed, " rows with missing values")
    }
    if (n_zscore_flagged > 0) {
      message("Removed ", n_zscore_flagged, " z-score outliers")
    }
    if (n_bw_flagged > 0) {
      message("Removed ", n_bw_flagged, " birth weight outliers")
    }
    if (n_decrease_flagged > 0) {
      message("Removed ", n_decrease_flagged, " decrease outliers")
    }
    if (length(subjects_to_remove) > 0) {
      message("Removed ", length(subjects_to_remove),
              " subjects with < 2 observations")
    }
    message("Final: ", nrow(data), " rows, ",
            length(unique(data[[id_var]])), " subjects")
  }

  new_lctm_cleaned(
    data = data,
    outcome = outcome,
    time_var = time_var,
    id_var = id_var,
    sex_var = sex_var,
    n_removed = as.integer(n_removed),
    subjects_removed = as.integer(subjects_removed)
  )
}
