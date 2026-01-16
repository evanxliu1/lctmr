#' Plot LCTM Trajectories (Step 5)
#'
#' Creates visualizations of fitted latent class trajectories.
#'
#' @param model An `lctm_model` object from [lctm_fit()].
#' @param type Character; type of plot - "mean" for mean trajectories,
#'   "spaghetti" for individual trajectories, or "both" for both plots.
#' @param ci Logical; if TRUE, adds confidence intervals to mean trajectory plot.
#'   Currently uses approximate intervals based on class-specific variation.
#' @param time_range Numeric vector of length 2 specifying the time range to plot.
#'   If NULL (default), uses the range of observed time values.
#' @param n_points Integer; number of points for mean trajectory prediction (default 100).
#' @param colors Character vector of colors for each class. If NULL, uses default palette.
#' @param alpha Numeric; transparency for spaghetti plot lines (default 0.3).
#'
#' @return A ggplot object (if type is "mean" or "spaghetti") or a list of
#'   ggplot objects (if type is "both").
#'
#' @details
#' This function performs Step 5 of the LCTM workflow - visualizing the results.
#'
#' **Plot Types:**
#' \describe{
#'   \item{mean}{Shows the predicted mean trajectory for each class with optional
#'     confidence intervals.}
#'   \item{spaghetti}{Shows individual trajectories colored by assigned class,
#'     useful for seeing within-class variation.}
#'   \item{both}{Returns a list with both plots.}
#' }
#'
#' @examples
#' \dontrun{
#' data(sample_growth)
#' setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid")
#' model <- lctm_fit(setup, k = 3, model = "F")
#'
#' # Mean trajectory plot
#' lctm_plot_trajectories(model, type = "mean")
#'
#' # Spaghetti plot
#' lctm_plot_trajectories(model, type = "spaghetti")
#'
#' # Both plots
#' plots <- lctm_plot_trajectories(model, type = "both")
#' plots$mean
#' plots$spaghetti
#' }
#'
#' @export
lctm_plot_trajectories <- function(model,
                                   type = "mean",
                                   ci = TRUE,
                                   time_range = NULL,
                                   n_points = 100,
                                   colors = NULL,
                                   alpha = 0.3) {

  # Validate inputs
  if (!inherits(model, "lctm_model")) {
    stop("model must be an lctm_model object from lctm_fit()", call. = FALSE)
  }

  if (!type %in% c("mean", "spaghetti", "both")) {
    stop("type must be 'mean', 'spaghetti', or 'both'", call. = FALSE)
  }

  # Extract info from model
  setup <- model$setup
  data <- setup$data
  outcome <- setup$outcome
  time_var <- setup$time_var
  id_var <- setup$id_var
  k <- model$k

  # Set default colors
  if (is.null(colors)) {
    colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                "#FF7F00", "#FFFF33", "#A65628")[seq_len(k)]
  }

  # Set time range
  if (is.null(time_range)) {
    time_range <- range(data[[time_var]], na.rm = TRUE)
  }

  # Get class assignments and merge back to data
  if (k == 1) {
    data$class <- factor(1)
  } else {
    # Get class assignments from pprob table (one row per subject)
    pprob_table <- model$model$pprob
    class_df <- data.frame(
      id = pprob_table[[1]],
      class = factor(pprob_table$class)
    )
    names(class_df)[1] <- id_var

    # Merge back to full data
    data <- merge(data, class_df, by = id_var, all.x = TRUE)
  }

  # Generate plots based on type
  if (type == "mean" || type == "both") {
    mean_plot <- .create_mean_plot(model, time_range, n_points, ci, colors)
  }

  if (type == "spaghetti" || type == "both") {
    spaghetti_plot <- .create_spaghetti_plot(data, outcome, time_var, id_var,
                                             time_range, colors, alpha)
  }

  # Return based on type
  if (type == "mean") {
    return(mean_plot)
  } else if (type == "spaghetti") {
    return(spaghetti_plot)
  } else {
    return(list(mean = mean_plot, spaghetti = spaghetti_plot))
  }
}

#' Create Mean Trajectory Plot
#' @keywords internal
.create_mean_plot <- function(model, time_range, n_points, ci, colors) {
  setup <- model$setup
  time_var <- setup$time_var
  outcome <- setup$outcome
  k <- model$k
  hlme_model <- model$model

  # Create prediction data frame with time sequence
  time_seq <- seq(time_range[1], time_range[2], length.out = n_points)
  newdata <- data.frame(time_seq)
  names(newdata) <- time_var

  # Use lcmm::predictY to get model-based predictions
  # This uses the fitted model coefficients (not raw data smoothing)
  plotpred <- lcmm::predictY(hlme_model, newdata, var.time = time_var, draws = ci)

  # Extract time values
  x_vals <- plotpred$times[[time_var]]

  # predictY$pred is a matrix, use colnames() and matrix subsetting
  pred_colnames <- colnames(plotpred$pred)

  # Reshape predictions from wide to long format for ggplot
  pred_list <- list()
  for (i in seq_len(k)) {
    class_name <- paste0("Ypred_class", i)

    pred_list[[i]] <- data.frame(
      time = x_vals,
      pred = plotpred$pred[, class_name],
      class = factor(i)
    )

    # Add CI bounds if available (when draws = TRUE)
    lower_col <- paste0("lower.Ypred_class", i)
    upper_col <- paste0("upper.Ypred_class", i)
    if (ci && lower_col %in% pred_colnames) {
      pred_list[[i]]$lower <- plotpred$pred[, lower_col]
      pred_list[[i]]$upper <- plotpred$pred[, upper_col]
    }
  }

  pred_df <- do.call(rbind, pred_list)

  # Create the plot
  p <- ggplot2::ggplot(pred_df, ggplot2::aes(x = .data$time, y = .data$pred,
                                              color = .data$class)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = colors,
                                labels = paste("Class", seq_len(k))) +
    ggplot2::labs(
      title = "Mean Trajectories by Class",
      x = time_var,
      y = outcome,
      color = "Class"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  # Add confidence interval ribbons if available
  if (ci && "lower" %in% names(pred_df)) {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                                         fill = .data$class),
                           alpha = 0.2, color = NA) +
      ggplot2::scale_fill_manual(values = colors,
                                 labels = paste("Class", seq_len(k)),
                                 guide = "none")
  }

  return(p)
}

#' Create Spaghetti Plot
#' @keywords internal
.create_spaghetti_plot <- function(data, outcome, time_var, id_var,
                                   time_range, colors, alpha) {
  k <- length(unique(data$class))

  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[time_var]],
                                           y = .data[[outcome]],
                                           group = .data[[id_var]],
                                           color = .data$class)) +
    ggplot2::geom_line(alpha = alpha) +
    ggplot2::scale_color_manual(values = colors,
                                labels = paste("Class", seq_len(k))) +
    ggplot2::labs(
      title = "Individual Trajectories by Class",
      x = time_var,
      y = outcome,
      color = "Class"
    ) +
    ggplot2::coord_cartesian(xlim = time_range) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  return(p)
}

#' Plot Residuals for LCTM Model
#'
#' Creates diagnostic residual plots for a fitted LCTM model.
#'
#' @param model An `lctm_model` object from [lctm_fit()].
#' @param standardized Logical; if TRUE (default), plots standardized residuals.
#' @param by_class Logical; if TRUE (default), facets by class assignment.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' model <- lctm_fit(setup, k = 3, model = "F")
#' lctm_plot_residuals(model)
#' }
#'
#' @export
lctm_plot_residuals <- function(model,
                                standardized = TRUE,
                                by_class = TRUE) {

  if (!inherits(model, "lctm_model")) {
    stop("model must be an lctm_model object from lctm_fit()", call. = FALSE)
  }

  setup <- model$setup
  time_var <- setup$time_var
  hlme_model <- model$model

  # Extract residuals
  resid_type <- if (standardized) "residSS" else "resid_m"
  residuals_data <- data.frame(
    time = setup$data[[time_var]],
    residual = hlme_model[[resid_type]]
  )

  # Add class assignments
  if (model$k > 1 && by_class) {
    prob_matrix <- extract_pprob(model)
    residuals_data$class <- factor(assign_classes(prob_matrix))
  } else {
    residuals_data$class <- factor(1)
    by_class <- FALSE
  }

  # Remove NA residuals
  residuals_data <- residuals_data[!is.na(residuals_data$residual), ]

  # Create plot
  y_label <- if (standardized) "Standardized Residual" else "Residual"

  p <- ggplot2::ggplot(residuals_data, ggplot2::aes(x = .data$time,
                                                    y = .data$residual)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "blue",
                         formula = y ~ x) +
    ggplot2::labs(
      title = paste(y_label, "vs Time"),
      x = time_var,
      y = y_label
    ) +
    ggplot2::theme_minimal()

  if (by_class) {
    p <- p + ggplot2::facet_wrap(~ class, labeller = ggplot2::label_both)
  }

  return(p)
}
