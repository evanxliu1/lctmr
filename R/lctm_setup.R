#' Set Up LCTM Analysis (Steps 1-2)
#'
#' Prepares data for latent class trajectory modeling by validating inputs,
#' fitting a base model, and comparing BIC across different numbers of classes.
#'
#' @param data A data frame containing longitudinal data.
#' @param outcome Character string naming the outcome variable (e.g., "weight").
#' @param time_var Character string naming the time variable (e.g., "age").
#' @param id_var Character string naming the subject ID variable.
#' @param k_range Integer vector of K values to compare (default 2:7).
#'   K=1 is always calculated as a reference but excluded from ranking
#'   since it represents no latent classes.
#' @param degree Polynomial degree for time: 1 = linear, 2 = quadratic (default).
#' @param verbose Logical; if TRUE, prints progress messages.
#'
#' @return An `lctm_setup` object containing:
#' \describe{
#'   \item{data}{The validated data frame}
#'   \item{outcome}{Name of outcome variable}
#'   \item{time_var}{Name of time variable}
#'   \item{id_var}{Name of ID variable}
#'   \item{degree}{Polynomial degree used}
#'   \item{bic_table}{Data frame with K and BIC values, sorted by BIC}
#'   \item{k_ranking}{Vector of K values ranked by BIC (K >= 2 only)}
#'   \item{base_model}{The base model (ng=1) for use as starting values}
#'   \item{k1_bic}{BIC for K=1 model (reference only)}
#' }
#'
#' @details
#' This function performs Steps 1-2 of the LCTM workflow:
#'
#' **Step 1:** Fits a base model with ng=1 to establish random effects structure
#' and provide starting values for multi-class models.
#'
#' **Step 2:** Fits models with K = 2, 3, ..., max(k_range) classes and compares
#' BIC values. Lower BIC indicates better fit.
#'
#' **Note on K=1:** The K=1 model (no latent classes) is always calculated as a
#' reference point, but it is excluded from the ranking because LCTM requires
#' at least 2 classes to identify distinct trajectory groups. K=1 often has the
#' lowest BIC due to fewer parameters, but it defeats the purpose of latent
#' class analysis.
#'
#' The output provides a ranked list of K values (K >= 2) to try in subsequent
#' model fitting with [lctm_fit()].
#'
#' @examples
#' \dontrun{
#' data(sample_growth)
#' setup <- lctm_setup(sample_growth,
#'                     outcome = "weight_raw",
#'                     time_var = "anthroage",
#'                     id_var = "childid")
#' print(setup)
#' }
#'
#' @export
lctm_setup <- function(data,
                       outcome,
                       time_var,
                       id_var,
                       k_range = 2:7,
                       degree = 2,
                       verbose = TRUE) {

  # Validate inputs
  data <- validate_lctm_data(data, outcome, time_var, id_var)

  if (!is.numeric(k_range) || any(k_range < 1)) {
    stop("k_range must be positive integers", call. = FALSE)
  }

  if (!degree %in% c(1, 2)) {
    stop("degree must be 1 (linear) or 2 (quadratic)", call. = FALSE)
  }

  # Build formula components
  if (degree == 2) {
    fixed_formula <- stats::as.formula(
      paste0(outcome, " ~ 1 + ", time_var, " + I(", time_var, "^2)")
    )
    random_formula <- stats::as.formula(
      paste0("~ 1 + ", time_var, " + I(", time_var, "^2)")
    )
    mixture_formula <- random_formula
  } else {
    fixed_formula <- stats::as.formula(
      paste0(outcome, " ~ 1 + ", time_var)
    )
    random_formula <- stats::as.formula(
      paste0("~ 1 + ", time_var)
    )
    mixture_formula <- random_formula
  }

  if (verbose) message("Fitting base model (K=1)...")

  # Fit base model (ng = 1) - this provides starting values AND K=1 reference
 base_model <- lcmm::hlme(
    fixed = fixed_formula,
    random = random_formula,
    ng = 1,
    idiag = FALSE,
    data = data,
    subject = id_var
  )

  # Store K=1 BIC as reference
  k1_bic <- base_model$BIC

  # Ensure k_range doesn't include 1 (we handle it separately)
  k_range <- k_range[k_range >= 2]

  if (length(k_range) == 0) {
    stop("k_range must include at least one value >= 2", call. = FALSE)
  }

  if (verbose) message("Comparing BIC across K = ", min(k_range), " to ", max(k_range), "...")

  # Compare BIC across K values (K >= 2 only)
  bic_results <- data.frame(
    k = integer(),
    bic = numeric(),
    converged = logical(),
    note = character(),
    stringsAsFactors = FALSE
  )

  for (k in k_range) {
    if (verbose) message("  Fitting K = ", k, "...")

    tryCatch({
      # Fit multi-class model using base model as starting values
      model_k <- lcmm::hlme(
        fixed = fixed_formula,
        mixture = mixture_formula,
        random = random_formula,
        ng = k,
        nwg = TRUE,  # Use Model F (proportional variance) for BIC comparison
        idiag = FALSE,
        data = data,
        subject = id_var,
        B = base_model
      )

      bic_results <- rbind(bic_results, data.frame(
        k = k,
        bic = model_k$BIC,
        converged = (model_k$conv == 1),
        note = "",
        stringsAsFactors = FALSE
      ))

    }, error = function(e) {
      if (verbose) message("    Warning: K = ", k, " failed to converge")
      bic_results <<- rbind(bic_results, data.frame(
        k = k,
        bic = NA_real_,
        converged = FALSE,
        note = "failed",
        stringsAsFactors = FALSE
      ))
    })
  }

  # Add K=1 as reference row
  k1_row <- data.frame(
    k = 1L,
    bic = k1_bic,
    converged = TRUE,
    note = "reference",
    stringsAsFactors = FALSE
  )
  bic_results <- rbind(k1_row, bic_results)

  # Sort by BIC (ascending - lower is better)
  bic_results <- bic_results[order(bic_results$bic, na.last = TRUE), ]
  rownames(bic_results) <- NULL

  # Add rank column
  bic_results$rank <- seq_len(nrow(bic_results))
  bic_results <- bic_results[, c("rank", "k", "bic", "converged", "note")]

  # Get K ranking (EXCLUDE K=1 since it's just reference)
  k_ranking <- bic_results$k[!is.na(bic_results$bic) & bic_results$k >= 2]

  if (verbose) {
    message("\nBIC Comparison (lower is better):")
    message("  K=1 (reference): BIC = ", round(k1_bic, 2), " (no latent classes)")
    message("  Best K for LCTM: K = ", k_ranking[1], " (BIC = ",
            round(bic_results$bic[bic_results$k == k_ranking[1]], 2), ")")
  }

  # Create and return lctm_setup object
  result <- new_lctm_setup(
    data = data,
    outcome = outcome,
    time_var = time_var,
    id_var = id_var,
    degree = as.integer(degree),
    bic_table = bic_results,
    k_ranking = k_ranking,
    base_model = base_model
  )

  # Add K=1 BIC as attribute
  attr(result, "k1_bic") <- k1_bic

  validate_lctm_setup(result)
  return(result)
}
