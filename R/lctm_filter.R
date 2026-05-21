#' Remove subjects in below-floor classes and return a filtered dataset
#'
#' Investigator-driven companion to the `min_prop` floor reported by
#' [lctm_adequacy()]. When a refined model has one or more classes smaller than
#' the minimum proportion, this function identifies the subjects assigned to
#' those classes and returns a NEW data frame with them removed, along with a
#' full record of what was removed. It never modifies the original data and it
#' never refits anything itself: you decide whether to feed the filtered data
#' back into [lctm_initial()]/[lctm_refine()].
#'
#' @param result An `lctm_result` (from [lctm_refine()]) or an `lctm_model`
#'   (e.g. `result$best_model`).
#' @param min_prop Minimum class proportion. Classes below this share of
#'   subjects are flagged for removal. Default 0.05.
#' @param data Optional data frame to filter. Defaults to the data stored in
#'   the model's setup (the data the model was fit on).
#' @param verbose Logical; print a detailed report of what was removed.
#'
#' @return Invisibly, a list with:
#' \describe{
#'   \item{data}{The new filtered data frame (subjects in below-floor classes removed).}
#'   \item{removed_ids}{The subject IDs that were removed.}
#'   \item{removed_classes}{Data frame of the below-floor classes (class, n, proportion).}
#'   \item{n_removed}{Number of subjects removed.}
#'   \item{n_remaining}{Number of subjects remaining.}
#'   \item{id_var}{The id column name used.}
#'   \item{min_prop}{The floor that was applied.}
#'   \item{any_removed}{Logical; TRUE if any subjects were removed.}
#' }
#'
#' @details
#' Workflow for the investigator-driven loop:
#' 1. Fit with [lctm_refine()] and inspect `lctm_adequacy()` (it reports which
#'    classes fall below `min_prop`).
#' 2. If you decide to act, call `lctm_filter_small_classes(result)` to get the
#'    filtered data and the list of removed subjects.
#' 3. Re-run [lctm_initial()] + [lctm_refine()] on `$data` (typically at the
#'    same K first).
#' 4. Repeat as many rounds as you judge appropriate, or stop.
#'
#' Because every removed subject ID is returned, the exclusions can be reported
#' transparently (for example in a participant-flow diagram).
#'
#' @examples
#' \dontrun{
#' result <- lctm_refine(init, random = ~ 1 + anthroage, k_range = 2:5)
#' print(result$adequacy)                       # see which classes are below floor
#' filt <- lctm_filter_small_classes(result)    # remove them, get new data
#' init2 <- lctm_initial(filt$data, k = result$best_k, degree = 2)
#' result2 <- lctm_refine(init2, random = ~ 1 + anthroage,
#'                        k_range = result$best_k)   # refit at same K first
#' }
#'
#' @export
lctm_filter_small_classes <- function(result, min_prop = 0.05,
                                      data = NULL, verbose = TRUE) {
  # Accept either an lctm_result or an lctm_model
  if (inherits(result, "lctm_result")) {
    model <- result$best_model
  } else if (inherits(result, "lctm_model")) {
    model <- result
  } else {
    stop("`result` must be an lctm_result or lctm_model object", call. = FALSE)
  }
  if (is.null(model)) {
    stop("No fitted model found (best_model is NULL). Nothing to filter.",
         call. = FALSE)
  }

  setup   <- model$setup
  id_var  <- setup$id_var
  if (is.null(data)) data <- setup$data

  k <- model$k
  pprob <- model$model$pprob
  cls <- data.frame(id = pprob[[1]], class = pprob$class,
                    stringsAsFactors = FALSE)

  # Class proportions by subject
  counts <- table(factor(cls$class, levels = seq_len(k)))
  props  <- as.numeric(prop.table(counts))

  below <- which(props < min_prop)

  removed_classes <- data.frame(
    class      = below,
    n          = as.integer(counts[below]),
    proportion = props[below],
    stringsAsFactors = FALSE
  )

  removed_ids <- cls$id[cls$class %in% below]
  filtered    <- data[!(data[[id_var]] %in% removed_ids), , drop = FALSE]

  n_removed   <- length(removed_ids)
  n_remaining <- length(setdiff(unique(cls$id), removed_ids))
  any_removed <- n_removed > 0

  if (verbose) {
    cat("lctm_filter_small_classes()\n")
    cat("---------------------------\n")
    cat("Model: K =", k, ", type", model$model_type, "\n")
    cat("Floor (min_prop):", min_prop, "(", 100 * min_prop, "% )\n\n")
    cat("Class sizes (by subject):\n")
    for (j in seq_len(k)) {
      flag <- if (props[j] < min_prop) "  <-- BELOW FLOOR" else ""
      cat("  Class ", j, ": n = ", as.integer(counts[j]),
          " (", round(100 * props[j], 1), "%)", flag, "\n", sep = "")
    }
    cat("\n")
    if (any_removed) {
      cat("Below-floor classes: ", paste(below, collapse = ", "), "\n", sep = "")
      cat("Subjects to remove:  ", n_removed, "\n", sep = "")
      cat("Subjects remaining:  ", n_remaining, "\n", sep = "")
      cat("Removed ", id_var, " IDs:\n", sep = "")
      cat("  ", paste(removed_ids, collapse = ", "), "\n", sep = "")
      cat("\nReturned $data has these subjects removed. Re-run lctm_initial() +\n")
      cat("lctm_refine() on it to refit (same K first is the usual choice).\n")
    } else {
      cat("No classes are below the floor. Nothing removed.\n")
      cat("$data is unchanged from the input.\n")
    }
  }

  invisible(list(
    data            = filtered,
    removed_ids     = removed_ids,
    removed_classes = removed_classes,
    n_removed       = n_removed,
    n_remaining     = n_remaining,
    id_var          = id_var,
    min_prop        = min_prop,
    any_removed     = any_removed
  ))
}
