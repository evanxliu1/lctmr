#' Utility Functions for lctmr
#'
#' Helper functions for class assignment and model comparison.
#'
#' @name lctmr-utils
NULL

#' Assign Observations to Classes
#'
#' Assigns each observation to the class with the highest posterior probability.
#'
#' @param prob_matrix A matrix or data frame of posterior probabilities where
#'   rows are observations and columns are classes.
#'
#' @return An integer vector of class assignments (1, 2, ..., K).
#'
#' @examples
#' prob_matrix <- matrix(c(
#'   0.9, 0.05, 0.05,
#'   0.1, 0.8, 0.1,
#'   0.1, 0.1, 0.8
#' ), nrow = 3, byrow = TRUE)
#' assign_classes(prob_matrix)
#'
#' @export
assign_classes <- function(prob_matrix) {
  prob_matrix <- as.matrix(prob_matrix)

  if (!is.numeric(prob_matrix)) {
    stop("prob_matrix must be numeric", call. = FALSE)
  }

  # Find column with maximum probability for each row
  assignments <- apply(prob_matrix, 1, which.max)

  return(as.integer(assignments))
}

#' Compare Multiple LCTM Models
#'
#' Creates a comparison table for multiple fitted lctm_model objects.
#'
#' @param model_list A named list of lctm_model objects to compare.
#'
#' @return A data frame with columns for model name, K, model type, BIC,
#'   and other comparison metrics, sorted by BIC (ascending).
#'
#' @examples
#' \dontrun{
#' # After fitting multiple models
#' models <- list(
#'   "K3_E" = model_k3_e,
#'   "K3_F" = model_k3_f,
#'   "K4_E" = model_k4_e
#' )
#' compare_models(models)
#' }
#'
#' @export
compare_models <- function(model_list) {
  if (!is.list(model_list) || length(model_list) == 0) {
    stop("model_list must be a non-empty list", call. = FALSE)
  }

  # Extract info from each model
  comparisons <- lapply(names(model_list), function(name) {
    m <- model_list[[name]]

    # Handle both lctm_model objects and raw hlme objects
    if (inherits(m, "lctm_model")) {
      data.frame(
        model = name,
        k = m$k,
        model_type = m$model_type,
        bic = m$bic,
        splines = m$use_splines,
        linear = m$linear,
        stringsAsFactors = FALSE
      )
    } else if (inherits(m, "hlme")) {
      data.frame(
        model = name,
        k = m$ng,
        model_type = NA_character_,
        bic = m$BIC,
        splines = NA,
        linear = NA,
        stringsAsFactors = FALSE
      )
    } else {
      warning("Skipping '", name, "': not an lctm_model or hlme object")
      NULL
    }
  })

  # Combine and sort
  comparison_df <- do.call(rbind, comparisons)

  if (is.null(comparison_df) || nrow(comparison_df) == 0) {
    stop("No valid models to compare", call. = FALSE)
  }

  # Sort by BIC
  comparison_df <- comparison_df[order(comparison_df$bic), ]
  rownames(comparison_df) <- NULL

  # Add rank column
  comparison_df$rank <- seq_len(nrow(comparison_df))
  comparison_df <- comparison_df[, c("rank", "model", "k", "model_type",
                                     "bic", "splines", "linear")]

  return(comparison_df)
}

#' Extract Posterior Probabilities from Model
#'
#' Extracts the posterior probability matrix from an lctm_model or hlme object.
#'
#' @param model An lctm_model or hlme object.
#'
#' @return A matrix of posterior probabilities with rows as observations
#'   and columns as classes.
#'
#' @keywords internal
extract_pprob <- function(model) {
  # Handle lctm_model objects
  if (inherits(model, "lctm_model")) {
    hlme_model <- model$model
  } else if (inherits(model, "hlme")) {
    hlme_model <- model
  } else {
    stop("model must be an lctm_model or hlme object", call. = FALSE)
  }

  # Extract posterior probabilities from hlme object
  pprob <- hlme_model$pprob

  # Get probability columns (prob1, prob2, etc.)
  prob_cols <- grep("^prob", names(pprob), value = TRUE)

  if (length(prob_cols) == 0) {
    stop("No probability columns found in model", call. = FALSE)
  }

  prob_matrix <- as.matrix(pprob[, prob_cols])
  return(prob_matrix)
}

#' Get Class Proportions from Model
#'
#' Calculates the proportion of observations in each class.
#'
#' @param model An lctm_model or hlme object.
#'
#' @return A named numeric vector of class proportions.
#'
#' @keywords internal
get_class_proportions <- function(model) {
  prob_matrix <- extract_pprob(model)
  assignments <- assign_classes(prob_matrix)
  K <- ncol(prob_matrix)

  class_counts <- table(factor(assignments, levels = seq_len(K)))
  proportions <- as.numeric(prop.table(class_counts))

  names(proportions) <- paste0("Class", seq_len(K))
  return(proportions)
}

#' Validate Input Data for LCTM
#'
#' Checks that input data meets requirements for LCTM analysis.
#'
#' @param data Data frame to validate
#' @param outcome Name of outcome variable
#' @param time_var Name of time variable
#' @param id_var Name of ID variable
#'
#' @return The validated data frame (possibly with minor cleaning)
#'
#' @keywords internal
validate_lctm_data <- function(data, outcome, time_var, id_var) {
  # Check data is a data frame
  if (!is.data.frame(data)) {
    stop("data must be a data frame", call. = FALSE)
  }

  # Check required columns exist
  required_cols <- c(outcome, time_var, id_var)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  # Check outcome is numeric
  if (!is.numeric(data[[outcome]])) {
    stop("Outcome variable '", outcome, "' must be numeric", call. = FALSE)
  }

  # Check time is numeric
  if (!is.numeric(data[[time_var]])) {
    stop("Time variable '", time_var, "' must be numeric", call. = FALSE)
  }

  # Check for complete cases in key variables
  complete_rows <- complete.cases(data[, required_cols])
  n_missing <- sum(!complete_rows)
  if (n_missing > 0) {
    message("Note: ", n_missing, " rows with missing values in key variables")
  }

  # Check minimum observations per subject
  obs_per_subject <- table(data[[id_var]])
  if (any(obs_per_subject < 2)) {
    n_single <- sum(obs_per_subject < 2)
    warning(n_single, " subject(s) have fewer than 2 observations")
  }

  return(data)
}
