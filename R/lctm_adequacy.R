#' Check Model Adequacy (Step 4)
#'
#' Evaluates whether a fitted LCTM model meets adequacy criteria for
#' class separation and classification accuracy.
#'
#' @param model An `lctm_model` object from [lctm_fit()].
#' @param thresholds A named list of threshold values:
#' \describe{
#'   \item{appa}{Minimum APPA value (default 0.70)}
#'   \item{occ}{Minimum OCC value (default 5.0)}
#'   \item{entropy}{Minimum relative entropy (default 0.5)}
#' }
#'
#' @return An `lctm_adequacy` object containing:
#' \describe{
#'   \item{appa}{Named vector of APPA values per class (NA for empty classes)}
#'   \item{appa_pass}{Logical; TRUE if all APPA values meet threshold}
#'   \item{occ}{Named vector of OCC values per class (NA/Inf for empty classes)}
#'   \item{occ_pass}{Logical; TRUE if all OCC values meet threshold}
#'   \item{entropy}{Relative entropy value}
#'   \item{entropy_pass}{Logical; TRUE if entropy meets threshold}
#'   \item{overall_pass}{Logical; TRUE if ALL criteria are met}
#'   \item{thresholds}{List of threshold values used}
#'   \item{is_degenerate}{Logical; TRUE if model has empty classes (degenerate)}
#' }
#'
#' @details
#' This function performs Step 4 of the LCTM workflow - checking whether
#' the fitted model has adequate class separation.
#'
#' **Adequacy Criteria:**
#' \describe{
#'   \item{APPA > 0.70}{Average Posterior Probability of Assignment should be
#'     at least 0.70 for each class, indicating good classification certainty.}
#'   \item{OCC > 5.0}{Odds of Correct Classification should be at least 5.0
#'     for each class, indicating the model classifies better than chance.}
#'   \item{Entropy > 0.5}{Relative entropy should be at least 0.5,
#'     indicating clear separation between classes.}
#' }
#'
#' If adequacy fails, consider:
#' - Trying a different number of classes (K)
#' - Using splines (`use_splines = TRUE` in [lctm_fit()])
#' - Using linear instead of quadratic model
#' - Trying the other model type (A vs B)
#'
#' @examples
#' \dontrun{
#' data(sample_growth)
#' setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid")
#' model <- lctm_fit(setup, k = 3, model = "B")
#'
#' adequacy <- lctm_adequacy(model)
#' print(adequacy)
#'
#' if (!adequacy$overall_pass) {
#'   # Try with splines
#'   model_splines <- lctm_fit(setup, k = 3, model = "B", use_splines = TRUE)
#'   adequacy_splines <- lctm_adequacy(model_splines)
#' }
#' }
#'
#' @export
lctm_adequacy <- function(model,
                          thresholds = list(appa = 0.70, occ = 5.0, entropy = 0.5)) {

  # Validate inputs
  if (!inherits(model, "lctm_model")) {
    stop("model must be an lctm_model object from lctm_fit()", call. = FALSE)
  }

  # Set default thresholds if not provided
  default_thresholds <- list(appa = 0.70, occ = 5.0, entropy = 0.5)
  thresholds <- modifyList(default_thresholds, thresholds)

  # Handle K=1 case (no class separation to evaluate)
  if (model$k == 1) {
    result <- new_lctm_adequacy(
      appa = c(Class1 = 1.0),
      appa_pass = TRUE,
      occ = c(Class1 = Inf),
      occ_pass = TRUE,
      entropy = 1.0,
      entropy_pass = TRUE,
      overall_pass = TRUE,
      thresholds = thresholds,
      is_degenerate = FALSE
    )
    return(result)
  }

  # Extract posterior probabilities
  prob_matrix <- extract_pprob(model)

  # Calculate adequacy metrics
  appa_values <- calc_appa(prob_matrix)
  occ_values <- calc_occ(prob_matrix, model$class_proportions)
  entropy_value <- calc_relative_entropy(prob_matrix)

  # Check for degenerate model (empty classes)
  # A model is degenerate if any class has no members (NA in APPA)
  # This matches LCTMtools behavior where NA indicates empty class
  n_empty_classes <- sum(is.na(appa_values))
  is_degenerate <- n_empty_classes > 0

  # Check against thresholds

  # For APPA: ALL classes must have valid values >= threshold
  # If any class is empty (NA), the model fails
  if (is_degenerate) {
    appa_pass <- FALSE
  } else {
    appa_pass <- all(appa_values >= thresholds$appa)
  }

  # For OCC: ALL classes must have valid (non-NA) values >= threshold
  # NA in OCC indicates empty class or class proportion issues — fail
  # Inf means APPA = 1 (perfect classification) — always passes
  if (any(is.na(occ_values))) {
    occ_pass <- FALSE
  } else {
    occ_pass <- all(occ_values >= thresholds$occ)
  }

  # Entropy check
  entropy_pass <- entropy_value >= thresholds$entropy

  # Overall: must pass ALL criteria AND not be degenerate
  overall_pass <- !is_degenerate && appa_pass && occ_pass && entropy_pass

  # Create and return lctm_adequacy object
  result <- new_lctm_adequacy(
    appa = appa_values,
    appa_pass = appa_pass,
    occ = occ_values,
    occ_pass = occ_pass,
    entropy = entropy_value,
    entropy_pass = entropy_pass,
    overall_pass = overall_pass,
    thresholds = thresholds,
    is_degenerate = is_degenerate
  )

  validate_lctm_adequacy(result)
  return(result)
}
