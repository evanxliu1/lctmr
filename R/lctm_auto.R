#' Automated LCTM Analysis
#'
#' Performs a complete latent class trajectory modeling analysis automatically,
#' selecting the best model that meets adequacy criteria.
#'
#' @note Consider using the new investigator-driven workflow instead:
#'   `lctm_clean() -> lctm_initial() -> lctm_refine()`.
#'   This function remains available for backward compatibility.
#'
#' @param data A data frame containing longitudinal data.
#' @param outcome Character string naming the outcome variable.
#' @param time_var Character string naming the time variable.
#' @param id_var Character string naming the subject ID variable.
#' @param k_range Integer vector of K values to search (default 2:7).
#'   K=1 is excluded because LCTM requires at least 2 classes.
#' @param models Character vector of model types to try (default c("B", "A")).
#' @param adequacy_thresholds List of adequacy thresholds:
#' \describe{
#'   \item{appa}{Minimum APPA (default 0.70)}
#'   \item{occ}{Minimum OCC (default 5.0)}
#'   \item{entropy}{Minimum relative entropy (default 0.5)}
#' }
#' @param try_splines Logical; if TRUE, tries adding splines when adequacy fails.
#' @param try_linear Logical; if TRUE, tries linear models when adequacy fails.
#' @param degree Polynomial degree: 1 = linear, 2 = quadratic (default).
#' @param verbose Logical; if TRUE, prints progress messages.
#'
#' @return An `lctm_result` object containing:
#' \describe{
#'   \item{best_model}{The best lctm_model object that passed adequacy, or NULL}
#'   \item{best_k}{Number of classes in best model}
#'   \item{best_model_type}{Model type of best model ("A" or "B")}
#'   \item{adequacy}{lctm_adequacy object for best model}
#'   \item{bic_table}{BIC comparison table from setup}
#'   \item{all_models}{List of all fitted models (optional)}
#'   \item{class_assignments}{Data frame with class assignments for best model}
#'   \item{search_history}{Data frame documenting the search process}
#' }
#'
#' @details
#' This function automates the 5-step LCTM workflow:
#'
#' 1. **Setup:** Validates data and compares BIC across K values
#' 2. **Search:** For each K (in BIC order):
#'    - Tries Model B, then Model A
#'    - If adequacy fails and `try_splines = TRUE`, tries with splines
#'    - If adequacy fails and `try_linear = TRUE`, tries linear model
#' 3. **Selection:** Returns first model that passes all adequacy criteria
#'
#' If no model passes adequacy, returns result with `best_model = NULL`.
#'
#' @examples
#' \dontrun{
#' data(sample_growth)
#'
#' # Fully automated analysis
#' result <- lctm_auto(sample_growth,
#'                     outcome = "weight_raw",
#'                     time_var = "anthroage",
#'                     id_var = "childid")
#'
#' # View results
#' print(result)
#' plot(result)
#'
#' # Access components
#' result$best_k
#' result$adequacy
#' }
#'
#' @export
lctm_auto <- function(data,
                      outcome,
                      time_var,
                      id_var,
                      k_range = 2:7,
                      models = c("B", "A"),
                      adequacy_thresholds = list(appa = 0.70, occ = 5.0, entropy = 0.5),
                      try_splines = TRUE,
                      try_linear = TRUE,
                      degree = 2,
                      verbose = TRUE) {

  # Set default thresholds
  default_thresholds <- list(appa = 0.70, occ = 5.0, entropy = 0.5)
  adequacy_thresholds <- modifyList(default_thresholds, adequacy_thresholds)

  if (verbose) message("=== LCTM Auto Analysis ===\n")

  # Step 1: Setup
  if (verbose) message("Step 1-2: Setting up and comparing BIC across K values...")

  setup <- lctm_setup(
    data = data,
    outcome = outcome,
    time_var = time_var,
    id_var = id_var,
    k_range = k_range,
    degree = degree,
    verbose = verbose
  )

  if (verbose) {
    message("\nBIC ranking: K = ", paste(setup$k_ranking, collapse = ", "))
    message("\n")
  }

  # Initialize search tracking
  search_history <- data.frame(
    k = integer(),
    model_type = character(),
    splines = logical(),
    linear = logical(),
    bic = numeric(),
    appa_pass = logical(),
    occ_pass = logical(),
    entropy_pass = logical(),
    overall_pass = logical(),
    stringsAsFactors = FALSE
  )

  all_models <- list()
  best_model <- NULL
  best_adequacy <- NULL

  # Step 2: Search for adequate model
  if (verbose) message("Step 3-4: Searching for model that meets adequacy criteria...\n")

  # Try K values in BIC order
  for (k in setup$k_ranking) {
    if (!is.null(best_model)) break  # Found adequate model

    if (verbose) message("--- Testing K = ", k, " ---")

    # Build list of model configurations to try
    configs <- list()

    # Base configurations: try each model type
    for (model_type in models) {
      configs[[length(configs) + 1]] <- list(
        model = model_type,
        splines = FALSE,
        linear = FALSE
      )
    }

    # Add spline configurations if enabled
    if (try_splines) {
      for (model_type in models) {
        configs[[length(configs) + 1]] <- list(
          model = model_type,
          splines = TRUE,
          linear = FALSE
        )
      }
    }

    # Add linear configurations if enabled
    if (try_linear) {
      for (model_type in models) {
        configs[[length(configs) + 1]] <- list(
          model = model_type,
          splines = FALSE,
          linear = TRUE
        )
      }
    }

    # Try each configuration
    for (config in configs) {
      if (!is.null(best_model)) break

      config_desc <- paste0(
        "Model ", config$model,
        if (config$splines) " + splines" else "",
        if (config$linear) " (linear)" else ""
      )

      if (verbose) message("  Trying ", config_desc, "...")

      # Fit model
      tryCatch({
        fitted_model <- lctm_fit(
          setup = setup,
          k = k,
          model = config$model,
          use_splines = config$splines,
          linear = config$linear,
          verbose = FALSE
        )

        # Check adequacy
        adequacy <- lctm_adequacy(fitted_model, thresholds = adequacy_thresholds)

        # Record in history
        search_history <- rbind(search_history, data.frame(
          k = k,
          model_type = config$model,
          splines = config$splines,
          linear = config$linear,
          bic = fitted_model$bic,
          appa_pass = adequacy$appa_pass,
          occ_pass = adequacy$occ_pass,
          entropy_pass = adequacy$entropy_pass,
          overall_pass = adequacy$overall_pass,
          stringsAsFactors = FALSE
        ))

        # Store model
        model_name <- paste0("K", k, "_", config$model,
                            if (config$splines) "_splines" else "",
                            if (config$linear) "_linear" else "")
        all_models[[model_name]] <- fitted_model

        # Check if adequate
        if (adequacy$overall_pass) {
          if (verbose) {
            message("    SUCCESS! Adequacy criteria met.")
            message("    APPA: ", paste(round(adequacy$appa, 2), collapse = ", "))
            message("    OCC: ", paste(round(adequacy$occ, 1), collapse = ", "))
            message("    Entropy: ", round(adequacy$entropy, 3))
          }
          best_model <- fitted_model
          best_adequacy <- adequacy
        } else {
          if (verbose) {
            fails <- c()
            if (!adequacy$appa_pass) fails <- c(fails, "APPA")
            if (!adequacy$occ_pass) fails <- c(fails, "OCC")
            if (!adequacy$entropy_pass) fails <- c(fails, "Entropy")
            message("    Failed: ", paste(fails, collapse = ", "))
          }
        }

      }, error = function(e) {
        if (verbose) message("    Error: ", e$message)
        search_history <<- rbind(search_history, data.frame(
          k = k,
          model_type = config$model,
          splines = config$splines,
          linear = config$linear,
          bic = NA_real_,
          appa_pass = NA,
          occ_pass = NA,
          entropy_pass = NA,
          overall_pass = FALSE,
          stringsAsFactors = FALSE
        ))
      })
    }
  }

  # Build result
  if (verbose) message("\n=== Results ===")

  if (!is.null(best_model)) {
    # Get class assignments from hlme's pprob table (one row per subject)
    pprob_table <- best_model$model$pprob
    class_assignments <- data.frame(
      id = pprob_table[[1]],  # First column is subject ID
      class = pprob_table$class
    )

    if (verbose) {
      message("Best model found: K = ", best_model$k,
              ", Model ", best_model$model_type)
      message("BIC: ", round(best_model$bic, 2))
    }

    result <- new_lctm_result(
      best_model = best_model,
      best_k = best_model$k,
      best_model_type = best_model$model_type,
      adequacy = best_adequacy,
      bic_table = setup$bic_table,
      all_models = all_models,
      class_assignments = class_assignments,
      search_history = search_history
    )
  } else {
    if (verbose) {
      message("No model met all adequacy criteria.")
      message("Consider:")
      message("  - Adjusting adequacy thresholds")
      message("  - Expanding k_range")
      message("  - Using lctm_setup() and lctm_fit() for manual exploration")
    }

    result <- new_lctm_result(
      best_model = NULL,
      best_k = NA_integer_,
      best_model_type = NA_character_,
      adequacy = NULL,
      bic_table = setup$bic_table,
      all_models = all_models,
      class_assignments = NULL,
      search_history = search_history
    )
  }

  validate_lctm_result(result)
  return(result)
}
