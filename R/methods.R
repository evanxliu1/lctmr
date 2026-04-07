#' S3 Methods for lctmr Classes
#'
#' Print, summary, and plot methods for lctmr objects.
#'
#' @name lctmr-methods
NULL

# =============================================================================
# lctm_setup methods
# =============================================================================

#' @export
print.lctm_setup <- function(x, ...) {
  cat("LCTM Setup\n")
  cat("----------\n")
  cat("Outcome:", x$outcome, "\n")
  cat("Time variable:", x$time_var, "\n")
  cat("ID variable:", x$id_var, "\n")
  cat("Polynomial degree:", x$degree, ifelse(x$degree == 2, "(quadratic)", "(linear)"), "\n")
  cat("Observations:", nrow(x$data), "\n")
  cat("Subjects:", length(unique(x$data[[x$id_var]])), "\n")
  cat("\n")
  cat("BIC Comparison (ranked by BIC, lower is better):\n")
  print(x$bic_table, row.names = FALSE)
  cat("\nNote: K=1 is shown for reference only (no latent classes).\n")
  cat("      LCTM requires K >= 2 to identify trajectory groups.\n")
  cat("\nRecommended K values to try:", paste(head(x$k_ranking, 3), collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.lctm_setup <- function(object, ...) {
  cat("LCTM Setup Summary\n")
  cat("==================\n\n")

  cat("Data:\n")
  cat("  Outcome variable:", object$outcome, "\n")
  cat("  Time variable:", object$time_var, "\n")
  cat("  ID variable:", object$id_var, "\n")
  cat("  Total observations:", nrow(object$data), "\n")
  cat("  Unique subjects:", length(unique(object$data[[object$id_var]])), "\n")
  cat("  Time range:", paste(range(object$data[[object$time_var]], na.rm = TRUE), collapse = " - "), "\n")
  cat("  Outcome range:", paste(round(range(object$data[[object$outcome]], na.rm = TRUE), 2), collapse = " - "), "\n")
  cat("\n")

  cat("Model specification:\n")
  cat("  Polynomial degree:", object$degree, "\n")
  cat("\n")

  cat("BIC Comparison:\n")
  print(object$bic_table, row.names = FALSE)
  cat("\n")
  cat("Note: K=1 is reference only (no latent classes).\n\n")

  cat("Recommendation (K >= 2):\n")
  cat("  Best K for LCTM:", object$k_ranking[1], "\n")
  cat("  Next best:", paste(head(object$k_ranking[-1], 2), collapse = ", "), "\n")

  invisible(object)
}

# =============================================================================
# lctm_model methods
# =============================================================================

#' @export
print.lctm_model <- function(x, ...) {
  cat("LCTM Model\n")
  cat("----------\n")
  cat("Model type:", x$model_type,
      ifelse(x$model_type == "A", "(common variance)", "(proportional variance)"), "\n")
  cat("Number of classes (K):", x$k, "\n")
  cat("BIC:", round(x$bic, 2), "\n")

  if (x$use_splines) cat("Splines: Yes\n")
  if (x$linear) cat("Linear: Yes (not quadratic)\n")

  cat("\nClass proportions:\n")
  for (i in seq_along(x$class_proportions)) {
    cat("  Class", i, ":", round(x$class_proportions[i] * 100, 1), "%\n")
  }

  # Convergence status
  conv_status <- x$model$conv
  if (conv_status == 1) {
    cat("\nConvergence: Yes\n")
  } else {
    cat("\nConvergence: No (code:", conv_status, ")\n")
  }

  invisible(x)
}

#' @export
summary.lctm_model <- function(object, ...) {
  print(object)
  cat("\n")

  cat("Fixed Effects:\n")
  print(summary(object$model)$fixed)
  cat("\n")

  invisible(object)
}

#' @export
plot.lctm_model <- function(x, type = "mean", ...) {
  lctm_plot_trajectories(x, type = type, ...)
}

#' @export
coef.lctm_model <- function(object, ...) {
  coef(object$model)
}

# =============================================================================
# lctm_adequacy methods
# =============================================================================

#' @export
print.lctm_adequacy <- function(x, ...) {
  cat("LCTM Model Adequacy\n")
  cat("-------------------\n\n")

  # Check for degenerate model first
  if (isTRUE(x$is_degenerate)) {
    cat("WARNING: DEGENERATE MODEL DETECTED\n")
    cat("One or more classes have no members assigned.\n")
    cat("This model has effectively collapsed to fewer classes.\n\n")
  }

  # APPA
  cat("APPA (threshold >=", x$thresholds$appa, "):\n")
  for (i in seq_along(x$appa)) {
    status <- ifelse(is.na(x$appa[i]), "N/A (empty class)",
                     ifelse(x$appa[i] >= x$thresholds$appa, "PASS", "FAIL"))
    cat("  ", names(x$appa)[i], ": ", round(x$appa[i], 3), " [", status, "]\n", sep = "")
  }
  cat("  Overall:", ifelse(x$appa_pass, "PASS", "FAIL"), "\n\n")

  # OCC
  cat("OCC (threshold >=", x$thresholds$occ, "):\n")
  for (i in seq_along(x$occ)) {
    status <- ifelse(is.na(x$occ[i]) | is.infinite(x$occ[i]), "N/A (empty class)",
                     ifelse(x$occ[i] >= x$thresholds$occ, "PASS", "FAIL"))
    val <- ifelse(is.infinite(x$occ[i]), "Inf", round(x$occ[i], 2))
    cat("  ", names(x$occ)[i], ": ", val, " [", status, "]\n", sep = "")
  }
  cat("  Overall:", ifelse(x$occ_pass, "PASS", "FAIL"), "\n\n")

  # Entropy
  cat("Relative Entropy (threshold >=", x$thresholds$entropy, "):\n")
  cat("  Value:", round(x$entropy, 3), "[",
      ifelse(x$entropy_pass, "PASS", "FAIL"), "]\n\n")

  # Overall
  cat("========================================\n")
  if (x$overall_pass) {
    cat("OVERALL: PASS - Model meets all adequacy criteria\n")
  } else {
    cat("OVERALL: FAIL - Model does not meet all criteria\n")
    if (isTRUE(x$is_degenerate)) {
      cat("\nReason: Degenerate model (empty classes detected)\n")
    }
    cat("\nSuggestions:\n")
    cat("  - Try a different number of classes (K)\n")
    cat("  - Try adding splines (use_splines = TRUE)\n")
    cat("  - Try linear model (linear = TRUE)\n")
    cat("  - Try the other model type (A vs B)\n")
  }
  cat("========================================\n")

  invisible(x)
}

# =============================================================================
# lctm_result methods (for lctm_auto / lctm_refine output)
# =============================================================================

#' @export
print.lctm_result <- function(x, ...) {
  cat("LCTM Result\n")
  cat("===========\n\n")

  if (is.null(x$best_model)) {
    cat("No adequate model found.\n")
    cat("Consider:\n")
    cat("  - Expanding the search range\n")
    cat("  - Adjusting adequacy thresholds\n")
    cat("  - Manual exploration with lctm_setup() and lctm_fit()\n")
  } else {
    cat("Best Model Found:\n")
    cat("  Model type:", x$best_model_type, "\n")
    cat("  Number of classes (K):", x$best_k, "\n")
    cat("  BIC:", round(x$best_model$bic, 2), "\n")

    if (x$best_model$use_splines) cat("  Uses splines: Yes\n")
    if (x$best_model$linear) cat("  Linear model: Yes\n")

    cat("\nClass proportions:\n")
    props <- x$best_model$class_proportions
    for (i in seq_along(props)) {
      cat("  Class", i, ":", round(props[i] * 100, 1), "%\n")
    }

    cat("\nAdequacy: PASS\n")
    cat("  APPA:", paste(round(x$adequacy$appa, 2), collapse = ", "), "\n")
    cat("  OCC:", paste(round(x$adequacy$occ, 1), collapse = ", "), "\n")
    cat("  Entropy:", round(x$adequacy$entropy, 3), "\n")
  }

  invisible(x)
}

#' @export
summary.lctm_result <- function(object, ...) {
  print(object)

  if (!is.null(object$best_model)) {
    cat("\nBIC Comparison Table:\n")
    print(object$bic_table, row.names = FALSE)

    if (!is.null(object$search_history)) {
      cat("\nSearch History:\n")
      print(object$search_history, row.names = FALSE)
    }
  }

  invisible(object)
}

#' @export
plot.lctm_result <- function(x, type = "mean", ...) {
  if (is.null(x$best_model)) {
    stop("No model to plot - no adequate model was found", call. = FALSE)
  }
  lctm_plot_trajectories(x$best_model, type = type, ...)
}

# =============================================================================
# lctm_cleaned methods
# =============================================================================

#' @export
print.lctm_cleaned <- function(x, ...) {
  cat("LCTM Cleaned Data\n")
  cat("-----------------\n")
  cat("Outcome:", x$outcome, "\n")
  cat("Time variable:", x$time_var, "\n")
  cat("ID variable:", x$id_var, "\n")
  if (!is.null(x$sex_var)) cat("Sex variable:", x$sex_var, "\n")
  cat("Observations:", nrow(x$data), "\n")
  cat("Subjects:", length(unique(x$data[[x$id_var]])), "\n")
  if (x$n_removed > 0) {
    cat("Rows removed:", x$n_removed, "\n")
  }
  if (x$subjects_removed > 0) {
    cat("Subjects removed:", x$subjects_removed, "\n")
  }
  cat("\nReady for lctm_initial()\n")
  invisible(x)
}

# =============================================================================
# lctm_initial methods
# =============================================================================

#' @export
print.lctm_initial <- function(x, ...) {
  cat("LCTM Initial Model\n")
  cat("------------------\n")
  cat("Outcome:", x$outcome, "\n")
  cat("Time variable:", x$time_var, "\n")
  degree_label <- c("1" = "linear", "2" = "quadratic", "3" = "cubic")[as.character(x$degree)]
  cat("Degree:", x$degree, paste0("(", degree_label, ")"), "\n")
  if (!is.null(x$knots)) {
    cat("Knots:", paste(x$knots, collapse = ", "), "\n")
  }
  cat("Number of classes (K):", x$k, "\n")
  cat("Random effects: intercept only (~1)\n")

  # Convergence
  if (!is.null(x$initial_model)) {
    conv <- x$initial_model$conv
    cat("Convergence:", ifelse(conv == 1, "Yes", paste0("No (code: ", conv, ")")), "\n")
    cat("BIC:", round(x$initial_model$BIC, 2), "\n")

    # Class proportions
    k <- x$k
    class_counts <- table(factor(x$initial_model$pprob$class, levels = seq_len(k)))
    class_proportions <- as.numeric(prop.table(class_counts))
    cat("\nClass proportions:\n")
    for (i in seq_len(k)) {
      cat("  Class", i, ":", round(class_proportions[i] * 100, 1), "%\n")
    }
  }

  cat("\n=== Residual Interpretation Guide ===\n")
  cat("  Horizontal band       -> random intercept only (random = ~ 1)\n")
  cat("  Diagonal/linear trend  -> add linear random slope\n")
  cat("  Curved pattern        -> add quadratic random effect\n")
  cat("  S-shaped pattern      -> consider cubic terms or splines\n")
  cat("\nUse plot() to view diagnostic plots, then pass to lctm_refine().\n")

  invisible(x)
}

#' @export
plot.lctm_initial <- function(x, which = NULL, ...) {
  if (length(x$plots) == 0) {
    stop("No plots available in this lctm_initial object", call. = FALSE)
  }

  if (is.null(which)) {
    # Print all plots
    for (p in x$plots) {
      print(p)
    }
  } else {
    if (!which %in% names(x$plots)) {
      stop("Plot '", which, "' not found. Available: ",
           paste(names(x$plots), collapse = ", "), call. = FALSE)
    }
    print(x$plots[[which]])
  }
  invisible(x)
}
