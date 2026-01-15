#' S3 Class Constructors and Validators for lctmr
#'
#' Internal functions for creating and validating custom S3 classes.
#'
#' @name lctmr-classes
#' @keywords internal
NULL

# -----------------------------------------------------------------------------
# lctm_setup class
# -----------------------------------------------------------------------------

#' Create a new lctm_setup object
#'
#' @param data Validated data frame
#' @param outcome Name of outcome variable
#' @param time_var Name of time variable
#' @param id_var Name of ID variable
#' @param degree Polynomial degree (1 = linear, 2 = quadratic)
#' @param bic_table Data frame with BIC comparisons
#' @param k_ranking Vector of K values ranked by BIC
#' @param base_model The base model (ng=1) used for initial values
#'
#' @return An lctm_setup object
#' @keywords internal
new_lctm_setup <- function(data,
                           outcome,
                           time_var,
                           id_var,
                           degree = 2L,
                           bic_table = NULL,
                           k_ranking = NULL,
                           base_model = NULL) {
  structure(
    list(
      data = data,
      outcome = outcome,
      time_var = time_var,
      id_var = id_var,
      degree = degree,
      bic_table = bic_table,
      k_ranking = k_ranking,
      base_model = base_model
    ),
    class = "lctm_setup"
  )
}

#' Validate an lctm_setup object
#'
#' @param x An lctm_setup object to validate
#' @return The validated object (invisibly) or throws an error
#' @keywords internal
validate_lctm_setup <- function(x) {
  if (!inherits(x, "lctm_setup")) {
    stop("Object must be of class 'lctm_setup'", call. = FALSE)
  }

  # Check required components
 if (!is.data.frame(x$data)) {
    stop("$data must be a data frame", call. = FALSE)
  }

  if (!is.character(x$outcome) || length(x$outcome) != 1) {
    stop("$outcome must be a single character string", call. = FALSE)
  }

  if (!is.character(x$time_var) || length(x$time_var) != 1) {
    stop("$time_var must be a single character string", call. = FALSE)
  }

  if (!is.character(x$id_var) || length(x$id_var) != 1) {
    stop("$id_var must be a single character string", call. = FALSE)
  }

  # Check that variables exist in data
  required_vars <- c(x$outcome, x$time_var, x$id_var)
  missing_vars <- setdiff(required_vars, names(x$data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "),
         call. = FALSE)
  }

  invisible(x)
}

# -----------------------------------------------------------------------------
# lctm_model class
# -----------------------------------------------------------------------------

#' Create a new lctm_model object
#'
#' @param model The lcmm hlme model object
#' @param k Number of classes
#' @param model_type Model type ("E" or "F")
#' @param bic BIC value
#' @param class_proportions Proportion in each class
#' @param setup The lctm_setup object used
#' @param use_splines Whether splines were used
#' @param linear Whether linear (vs quadratic) was used
#'
#' @return An lctm_model object
#' @keywords internal
new_lctm_model <- function(model,
                           k,
                           model_type,
                           bic,
                           class_proportions,
                           setup,
                           use_splines = FALSE,
                           linear = FALSE) {
  structure(
    list(
      model = model,
      k = k,
      model_type = model_type,
      bic = bic,
      class_proportions = class_proportions,
      setup = setup,
      use_splines = use_splines,
      linear = linear
    ),
    class = "lctm_model"
  )
}

#' Validate an lctm_model object
#'
#' @param x An lctm_model object to validate
#' @return The validated object (invisibly) or throws an error
#' @keywords internal
validate_lctm_model <- function(x) {
  if (!inherits(x, "lctm_model")) {
    stop("Object must be of class 'lctm_model'", call. = FALSE)
  }

  if (!inherits(x$model, "hlme")) {
    stop("$model must be an hlme object from lcmm package", call. = FALSE)
  }

  if (!is.numeric(x$k) || length(x$k) != 1 || x$k < 1) {
    stop("$k must be a positive integer", call. = FALSE)
  }

  if (!x$model_type %in% c("E", "F")) {
    stop("$model_type must be 'E' or 'F'", call. = FALSE)
  }

  invisible(x)
}

# -----------------------------------------------------------------------------
# lctm_adequacy class
# -----------------------------------------------------------------------------

#' Create a new lctm_adequacy object
#'
#' @param appa Named vector of APPA values per class
#' @param appa_pass Logical, whether APPA criterion is met
#' @param occ Named vector of OCC values per class
#' @param occ_pass Logical, whether OCC criterion is met
#' @param entropy Relative entropy value
#' @param entropy_pass Logical, whether entropy criterion is met
#' @param overall_pass Logical, whether all criteria are met
#' @param thresholds List of threshold values used
#'
#' @return An lctm_adequacy object
#' @keywords internal
new_lctm_adequacy <- function(appa,
                              appa_pass,
                              occ,
                              occ_pass,
                              entropy,
                              entropy_pass,
                              overall_pass,
                              thresholds) {
  structure(
    list(
      appa = appa,
      appa_pass = appa_pass,
      occ = occ,
      occ_pass = occ_pass,
      entropy = entropy,
      entropy_pass = entropy_pass,
      overall_pass = overall_pass,
      thresholds = thresholds
    ),
    class = "lctm_adequacy"
  )
}

#' Validate an lctm_adequacy object
#'
#' @param x An lctm_adequacy object to validate
#' @return The validated object (invisibly) or throws an error
#' @keywords internal
validate_lctm_adequacy <- function(x) {
  if (!inherits(x, "lctm_adequacy")) {
    stop("Object must be of class 'lctm_adequacy'", call. = FALSE)
  }

  if (!is.numeric(x$appa)) {
    stop("$appa must be numeric", call. = FALSE)
  }

  if (!is.logical(x$overall_pass) || length(x$overall_pass) != 1) {
    stop("$overall_pass must be a single logical value", call. = FALSE)
  }

  invisible(x)
}

# -----------------------------------------------------------------------------
# lctm_result class (for lctm_auto output)
# -----------------------------------------------------------------------------

#' Create a new lctm_result object
#'
#' @param best_model The best lctm_model object
#' @param best_k Optimal number of classes
#' @param best_model_type Which model type ("E" or "F")
#' @param adequacy lctm_adequacy object for best model
#' @param bic_table Data frame with all BIC comparisons
#' @param all_models List of all fitted lctm_model objects
#' @param class_assignments Data frame with individual class assignments
#' @param search_history Data frame documenting the search process
#'
#' @return An lctm_result object
#' @keywords internal
new_lctm_result <- function(best_model,
                            best_k,
                            best_model_type,
                            adequacy,
                            bic_table,
                            all_models = NULL,
                            class_assignments = NULL,
                            search_history = NULL) {
  structure(
    list(
      best_model = best_model,
      best_k = best_k,
      best_model_type = best_model_type,
      adequacy = adequacy,
      bic_table = bic_table,
      all_models = all_models,
      class_assignments = class_assignments,
      search_history = search_history
    ),
    class = "lctm_result"
  )
}

#' Validate an lctm_result object
#'
#' @param x An lctm_result object to validate
#' @return The validated object (invisibly) or throws an error
#' @keywords internal
validate_lctm_result <- function(x) {
 if (!inherits(x, "lctm_result")) {
    stop("Object must be of class 'lctm_result'", call. = FALSE)
  }

  if (!is.null(x$best_model) && !inherits(x$best_model, "lctm_model")) {
    stop("$best_model must be an lctm_model object or NULL", call. = FALSE)
  }

  invisible(x)
}
