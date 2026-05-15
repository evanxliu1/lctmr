#' Refine LCTM Model with User-Specified Random Effects
#'
#' Searches across K values and model types (A/B) for the best model that
#' meets adequacy criteria, using the random effects structure chosen by the
#' investigator after reviewing [lctm_initial()] diagnostics.
#'
#' @param initial An `lctm_initial` object from [lctm_initial()].
#' @param random Formula for random effects (e.g., `~ 1 + age`). If NULL,
#'   defaults to the full polynomial matching `degree` from the initial model.
#' @param k_range Integer vector of K values to search (default 2:7).
#' @param degree Polynomial degree for the trajectory: 1 = linear, 2 = quadratic,
#'   3 = cubic. If NULL (default), inherits from the initial model.
#' @param knots Numeric vector of knot positions for natural splines. If NULL
#'   (default), inherits from the initial model. Overrides `degree` when specified.
#' @param covariates Character vector of covariate names to add to the fixed
#'   formula. Covariates are added to `fixed` only, not `mixture` (covariates
#'   adjust overall mean, not class-specific trajectory shape).
#' @param models Character vector of model types to compare (default `c("B", "A")`).
#'   Model B = proportional variance (nwg=TRUE), Model A = common variance (nwg=FALSE).
#' @param adequacy_thresholds List of adequacy thresholds:
#' \describe{
#'   \item{appa}{Minimum APPA (default 0.70)}
#'   \item{occ}{Minimum OCC (default 5.0)}
#'   \item{entropy}{Minimum relative entropy (default 0.5)}
#' }
#' @param start_simple Logical; if TRUE, uses a simple base model (1 class, no
#'   random effects) for starting values instead of the full base model. This
#'   can improve efficiency for large datasets. Default FALSE.
#' @param save_pdf Character string; file path to save result plots as PDF.
#' @param verbose Logical; if TRUE, prints progress messages.
#'
#' @return An `lctm_result` object containing:
#' \describe{
#'   \item{best_model}{The best lctm_model object that passed adequacy, or NULL}
#'   \item{best_k}{Number of classes in best model}
#'   \item{best_model_type}{Model type ("A" or "B")}
#'   \item{adequacy}{lctm_adequacy object for best model}
#'   \item{bic_table}{BIC comparison table}
#'   \item{all_models}{List of all fitted models}
#'   \item{class_assignments}{Data frame with class assignments for best model}
#'   \item{search_history}{Data frame documenting the search process}
#' }
#'
#' @details
#' This function performs the third step of the new LCTM workflow:
#'
#' \code{lctm_clean() -> lctm_initial() -> [user reviews] -> lctm_refine()}
#'
#' The user specifies the trajectory form (degree/knots), random effects
#' structure, and covariates based on their review of the initial model.
#' The same model specification is used across all K values and Model A/B
#' comparisons.
#'
#' The search proceeds as follows:
#' 1. Fit a base model (ng=1) with the refined formulas for starting values
#' 2. Compare BIC across K values in `k_range`
#' 3. For each K (in BIC order), try Model B then A
#' 4. Stop at the first model that passes all adequacy criteria
#'
#' @examples
#' \dontrun{
#' data(sample_growth)
#' cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid")
#' init <- lctm_initial(cleaned, k = 2, degree = 2)
#' plot(init)  # review residuals
#'
#' # Based on residual patterns, add linear random slope
#' result <- lctm_refine(init, random = ~ 1 + anthroage, k_range = 2:4)
#' print(result)
#' plot(result)
#'
#' # With covariates
#' result <- lctm_refine(init, random = ~ 1 + anthroage,
#'                        covariates = c("waz"), k_range = 2:4)
#'
#' # With cubic model and simpler starting values for large data
#' result <- lctm_refine(init, degree = 3,
#'                        random = ~ 1 + anthroage,
#'                        start_simple = TRUE, k_range = 2:4)
#'
#' # With spline knots
#' result <- lctm_refine(init, knots = c(6, 12, 24),
#'                        random = ~ 1 + anthroage, k_range = 2:4)
#' }
#'
#' @export
lctm_refine <- function(initial,
                        random = NULL,
                        k_range = 2:7,
                        degree = NULL,
                        knots = NULL,
                        covariates = NULL,
                        models = c("B", "A"),
                        adequacy_thresholds = list(appa = 0.70, occ = 5.0, entropy = 0.5),
                        start_simple = FALSE,
                        save_pdf = NULL,
                        verbose = TRUE) {

  # Validate inputs
  if (!inherits(initial, "lctm_initial")) {
    stop("initial must be an lctm_initial object from lctm_initial()", call. = FALSE)
  }

  # Extract from initial
  data <- initial$data
  outcome <- initial$outcome
  time_var <- initial$time_var
  id_var <- initial$id_var

  # Inherit degree and knots from initial if not specified
  if (is.null(degree)) degree <- initial$degree
  if (is.null(knots)) knots <- initial$knots

  if (!degree %in% c(1, 2, 3)) {
    stop("degree must be 1 (linear), 2 (quadratic), or 3 (cubic)", call. = FALSE)
  }

  if (!is.null(knots) && degree > 1) {
    stop("knots cannot be combined with degree > 1. Use degree = 1 with knots, or a polynomial degree without knots.",
         call. = FALSE)
  }

  # Set default thresholds
  default_thresholds <- list(appa = 0.70, occ = 5.0, entropy = 0.5)
  adequacy_thresholds <- modifyList(default_thresholds, adequacy_thresholds)

  # Build random formula (default: full polynomial matching degree)
  if (is.null(random)) {
    if (degree == 3) {
      random_formula <- stats::as.formula(
        paste0("~ 1 + ", time_var, " + I(", time_var, "^2) + I(", time_var, "^3)")
      )
    } else if (degree == 2) {
      random_formula <- stats::as.formula(
        paste0("~ 1 + ", time_var, " + I(", time_var, "^2)")
      )
    } else {
      random_formula <- stats::as.formula(
        paste0("~ 1 + ", time_var)
      )
    }
    if (verbose) message("Using default random effects: ", deparse(random_formula))
  } else {
    random_formula <- random
    if (verbose) message("Using specified random effects: ", deparse(random_formula))
  }

  # Build trajectory formulas using shared helper
  formula_parts <- .build_trajectory_formulas(outcome, time_var, degree, knots)
  trajectory_terms <- formula_parts$terms

  # Add covariates to fixed formula only
  if (!is.null(covariates)) {
    missing_covs <- setdiff(covariates, names(data))
    if (length(missing_covs) > 0) {
      stop("Covariates not found in data: ",
           paste(missing_covs, collapse = ", "), call. = FALSE)
    }
    fixed_terms <- paste0(trajectory_terms, " + ", paste(covariates, collapse = " + "))
    if (verbose) message("Covariates added to fixed formula: ", paste(covariates, collapse = ", "))
  } else {
    fixed_terms <- trajectory_terms
  }

  fixed_formula <- stats::as.formula(paste0(outcome, " ~ ", fixed_terms))
  mixture_formula <- stats::as.formula(paste0("~ ", trajectory_terms))

  if (verbose) message("\n=== LCTM Refine Analysis ===\n")

  # Step 1: Fit base model (ng=1)
  if (start_simple) {
    if (verbose) message("Fitting simple base model (no random effects) for starting values...")
    base_args <- list(
      fixed = fixed_formula,
      random = ~ 1,
      ng = 1,
      idiag = FALSE,
      data = data,
      subject = id_var
    )
  } else {
    if (verbose) message("Fitting base model with refined formulas...")
    base_args <- list(
      fixed = fixed_formula,
      random = random_formula,
      ng = 1,
      idiag = FALSE,
      data = data,
      subject = id_var
    )
  }
  base_model <- do.call(lcmm::hlme, base_args)

  # Step 2: BIC comparison across k_range
  k_range <- k_range[k_range >= 2]
  if (length(k_range) == 0) {
    stop("k_range must include at least one value >= 2", call. = FALSE)
  }

  if (verbose) message("Comparing BIC across K = ", min(k_range), " to ", max(k_range), "...")

  bic_results <- data.frame(
    k = integer(), bic = numeric(), converged = logical(),
    stringsAsFactors = FALSE
  )

  for (k in k_range) {
    if (verbose) message("  Fitting K = ", k, "...")
    tryCatch({
      model_args <- list(
        fixed = fixed_formula,
        mixture = mixture_formula,
        random = random_formula,
        ng = k,
        nwg = TRUE,
        idiag = FALSE,
        data = data,
        subject = id_var,
        B = base_model
      )
      model_k <- do.call(lcmm::hlme, model_args)

      bic_results <- rbind(bic_results, data.frame(
        k = k, bic = model_k$BIC, converged = (model_k$conv == 1),
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      if (verbose) message("    Warning: K = ", k, " failed: ", e$message)
      bic_results <<- rbind(bic_results, data.frame(
        k = k, bic = NA_real_, converged = FALSE,
        stringsAsFactors = FALSE
      ))
    })
  }

  # Sort by BIC and get ranking
  bic_results <- bic_results[order(bic_results$bic, na.last = TRUE), ]
  rownames(bic_results) <- NULL

  # Add rank and note columns to match lctm_setup's bic_table schema
  bic_results$rank <- seq_len(nrow(bic_results))
  bic_results$note <- NA_character_
  bic_results <- bic_results[, c("rank", "k", "bic", "converged", "note")]

  k_ranking <- bic_results$k[!is.na(bic_results$bic)]

  if (verbose) {
    message("\nBIC ranking: K = ", paste(k_ranking, collapse = ", "), "\n")
  }

  # Step 3: Search for adequate model
  if (verbose) message("Searching for model that meets adequacy criteria...\n")

  search_history <- data.frame(
    k = integer(), model_type = character(), bic = numeric(),
    appa_pass = logical(), occ_pass = logical(), entropy_pass = logical(),
    overall_pass = logical(), stringsAsFactors = FALSE
  )

  all_models <- list()
  best_model <- NULL
  best_adequacy <- NULL

  # Internal setup holder used as context for fitted models
  temp_setup <- new_lctm_setup(
    data = data,
    outcome = outcome,
    time_var = time_var,
    id_var = id_var,
    degree = as.integer(degree),
    bic_table = bic_results,
    k_ranking = k_ranking,
    base_model = base_model
  )

  for (k in k_ranking) {
    if (!is.null(best_model)) break

    if (verbose) message("--- Testing K = ", k, " ---")

    for (model_type in models) {
      if (!is.null(best_model)) break

      nwg <- (model_type == "B")

      if (verbose) message("  Trying Model ", model_type, "...")

      tryCatch({
        hlme_args <- list(
          fixed = fixed_formula,
          mixture = mixture_formula,
          random = random_formula,
          ng = k,
          nwg = nwg,
          idiag = FALSE,
          data = data,
          subject = id_var,
          B = base_model
        )
        hlme_model <- do.call(lcmm::hlme, hlme_args)

        # Check convergence before proceeding
        if (hlme_model$conv != 1) {
          stop("Did not converge (code: ", hlme_model$conv, ")")
        }

        class_counts <- table(factor(hlme_model$pprob$class, levels = seq_len(k)))
        class_proportions <- as.numeric(prop.table(class_counts))
        names(class_proportions) <- paste0("Class", seq_len(k))

        fitted_model <- new_lctm_model(
          model = hlme_model,
          k = k,
          model_type = model_type,
          bic = hlme_model$BIC,
          class_proportions = class_proportions,
          setup = temp_setup,
          use_splines = !is.null(knots),
          linear = (degree == 1 && is.null(knots))
        )

        adequacy <- lctm_adequacy(fitted_model, thresholds = adequacy_thresholds)

        search_history <- rbind(search_history, data.frame(
          k = k, model_type = model_type, bic = fitted_model$bic,
          appa_pass = adequacy$appa_pass, occ_pass = adequacy$occ_pass,
          entropy_pass = adequacy$entropy_pass, overall_pass = adequacy$overall_pass,
          stringsAsFactors = FALSE
        ))

        model_name <- paste0("K", k, "_", model_type)
        all_models[[model_name]] <- fitted_model

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
          k = k, model_type = model_type, bic = NA_real_,
          appa_pass = NA, occ_pass = NA, entropy_pass = NA, overall_pass = FALSE,
          stringsAsFactors = FALSE
        ))
      })
    }
  }

  # Build result
  if (verbose) message("\n=== Results ===")

  if (!is.null(best_model)) {
    pprob_table <- best_model$model$pprob
    class_assignments <- data.frame(
      id = pprob_table[[1]],
      class = pprob_table$class
    )

    if (verbose) {
      message("Best model: K = ", best_model$k, ", Model ", best_model$model_type)
      message("BIC: ", round(best_model$bic, 2))
    }

    result <- new_lctm_result(
      best_model = best_model,
      best_k = best_model$k,
      best_model_type = best_model$model_type,
      adequacy = best_adequacy,
      bic_table = bic_results,
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
      message("  - Trying different random effects structure")
    }

    result <- new_lctm_result(
      best_model = NULL,
      best_k = NA_integer_,
      best_model_type = NA_character_,
      adequacy = NULL,
      bic_table = bic_results,
      all_models = all_models,
      class_assignments = NULL,
      search_history = search_history
    )
  }

  # Save PDF if requested
  if (!is.null(save_pdf) && !is.null(best_model)) {
    if (verbose) message("Saving result plots to: ", save_pdf)
    grDevices::pdf(save_pdf, width = 10, height = 7)
    on.exit(grDevices::dev.off(), add = TRUE)
    tryCatch({
      # Use the public plot interface so legend labels include n and %
      print(lctm_plot_trajectories(best_model, type = "mean",  ci = TRUE))
      print(lctm_plot_trajectories(best_model, type = "spaghetti"))
    }, error = function(e) {
      if (verbose) message("Warning: Could not generate PDF plots: ", e$message)
    })
  }

  validate_lctm_result(result)
  return(result)
}
