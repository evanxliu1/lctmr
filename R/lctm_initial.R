#' Fit Initial LCTM Model with Diagnostic Plots
#'
#' Fits an initial latent class trajectory model using random intercept only
#' (no random slopes) and generates diagnostic plots to guide the choice of
#' random effects structure for refinement.
#'
#' @param data An `lctm_cleaned` object from [lctm_clean()], or a raw data frame.
#' @param outcome Character string naming the outcome variable. Required if
#'   `data` is a raw data frame; ignored if `data` is an `lctm_cleaned` object.
#' @param time_var Character string naming the time variable.
#' @param id_var Character string naming the subject ID variable.
#' @param k Integer; number of latent classes (default 2).
#' @param degree Polynomial degree: 1 = linear, 2 = quadratic (default), 3 = cubic.
#' @param knots Numeric vector of knot positions for splines. If provided,
#'   spline terms are added to the fixed and mixture formulas. The knots are
#'   placed at the specified ages/time values (e.g., `knots = c(6, 12, 24)` for
#'   knots at 6, 12, and 24 months). Cannot be combined with `degree > 1`.
#' @param spline_degree Polynomial degree of the spline pieces when `knots` is
#'   supplied: `3` (default) builds a natural cubic spline via `splines::ns()`
#'   (piecewise cubic, C2 continuity); `2` builds a quadratic B-spline via
#'   `splines::bs(degree = 2)` (piecewise quadratic, C1 continuity, fewer
#'   parameters per class). Ignored when `knots` is NULL.
#' @param sex_var Character string naming the sex variable (optional). Used for
#'   faceting spaghetti and loess plots.
#' @param save_pdf Character string; file path to save plots as multi-page PDF.
#'   If NULL (default), no PDF is saved. The PDF includes all diagnostic plots
#'   plus a reference page with example residual patterns.
#' @param verbose Logical; if TRUE, prints progress messages.
#'
#' @return An `lctm_initial` object containing:
#' \describe{
#'   \item{data}{The data frame used}
#'   \item{outcome}{Name of outcome variable}
#'   \item{time_var}{Name of time variable}
#'   \item{id_var}{Name of ID variable}
#'   \item{sex_var}{Name of sex variable (or NULL)}
#'   \item{degree}{Polynomial degree}
#'   \item{k}{Number of classes}
#'   \item{knots}{Knot positions (or NULL)}
#'   \item{spline_degree}{Spline piece degree used when knots are supplied}
#'   \item{initial_model}{The fitted hlme object (random = ~1)}
#'   \item{base_model}{The ng=1 base model for starting values}
#'   \item{plots}{Named list of ggplot objects (spaghetti, loess, residuals, guide)}
#' }
#'
#' @details
#' This function fits an initial model with random intercept only. The residual
#' plots help determine what random effects structure to use in [lctm_refine()]:
#'
#' \describe{
#'   \item{Horizontal band}{Random intercept only is sufficient}
#'   \item{Diagonal/linear trend}{Add linear random slope}
#'   \item{Curved pattern}{Add quadratic random effect}
#'   \item{S-shaped pattern}{Consider cubic terms or splines}
#' }
#'
#' This is the second step of the new LCTM workflow:
#'
#' \code{lctm_clean() -> lctm_initial() -> [user reviews] -> lctm_refine()}
#'
#' @examples
#' \dontrun{
#' data(sample_growth)
#' cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid")
#' init <- lctm_initial(cleaned, k = 2, degree = 2)
#' plot(init)
#'
#' # With cubic model
#' init_cubic <- lctm_initial(cleaned, k = 2, degree = 3)
#'
#' # With natural cubic spline knots at specific ages
#' init_spline <- lctm_initial(cleaned, k = 2, degree = 1, knots = c(6, 12, 24))
#'
#' # With a quadratic (piecewise) spline instead of cubic
#' init_q <- lctm_initial(cleaned, k = 2, degree = 1, knots = c(6, 12, 24),
#'                        spline_degree = 2)
#' }
#'
#' @export
lctm_initial <- function(data, outcome = NULL, time_var = NULL, id_var = NULL,
                         k = 2, degree = 2, knots = NULL, spline_degree = 3,
                         sex_var = NULL, save_pdf = NULL, verbose = TRUE) {

  # Handle lctm_cleaned input
  if (inherits(data, "lctm_cleaned")) {
    outcome <- data$outcome
    time_var <- data$time_var
    id_var <- data$id_var
    if (is.null(sex_var)) sex_var <- data$sex_var
    data <- data$data
  } else {
    if (is.null(outcome) || is.null(time_var) || is.null(id_var)) {
      stop("outcome, time_var, and id_var are required when data is a data frame",
           call. = FALSE)
    }
    data <- validate_lctm_data(data, outcome, time_var, id_var)
  }

  if (!degree %in% c(1, 2, 3)) {
    stop("degree must be 1 (linear), 2 (quadratic), or 3 (cubic)", call. = FALSE)
  }

  if (!is.null(knots) && degree > 1) {
    stop("knots cannot be combined with degree > 1. Use degree = 1 with knots, or a polynomial degree without knots.",
         call. = FALSE)
  }

  if (!spline_degree %in% c(2, 3)) {
    stop("spline_degree must be 2 (quadratic) or 3 (cubic)", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 1 || k < 2) {
    stop("k must be an integer >= 2", call. = FALSE)
  }
  k <- as.integer(k)

  # Build formulas
  formula_parts <- .build_trajectory_formulas(outcome, time_var, degree, knots,
                                              spline_degree)
  fixed_formula <- formula_parts$fixed
  mixture_formula <- formula_parts$mixture

  # Random intercept only -- no random slopes
  random_formula <- ~ 1

  # Step 1: Fit base model (ng=1) for starting values
  if (verbose) message("Fitting base model (K=1)...")

  base_args <- list(
    fixed = fixed_formula,
    random = random_formula,
    ng = 1,
    idiag = FALSE,
    data = data,
    subject = id_var
  )
  base_model <- do.call(lcmm::hlme, base_args)

  # Step 2: Fit initial K-class model with random intercept only
  if (verbose) {
    degree_label <- c("1" = "linear", "2" = "quadratic", "3" = "cubic")[as.character(degree)]
    if (!is.null(knots)) {
      spline_label <- if (spline_degree == 2) "quadratic spline" else "natural cubic spline"
      degree_label <- spline_label
      knot_label <- paste0(" with knots at ", paste(knots, collapse = ", "))
    } else {
      knot_label <- ""
    }
    message("Fitting initial ", degree_label, " model with K = ", k,
            " classes (random intercept only)", knot_label, "...")
  }

  model_args <- list(
    fixed = fixed_formula,
    mixture = mixture_formula,
    random = random_formula,
    ng = k,
    nwg = FALSE,
    idiag = FALSE,
    data = data,
    subject = id_var,
    B = base_model
  )
  initial_model <- do.call(lcmm::hlme, model_args)

  # Check convergence
  if (initial_model$conv != 1) {
    warning("Initial model did not converge (conv = ", initial_model$conv, ")",
            call. = FALSE)
  }

  # Class proportions
  class_counts <- table(factor(initial_model$pprob$class, levels = seq_len(k)))
  class_proportions <- as.numeric(prop.table(class_counts))

  if (verbose) {
    message("  BIC = ", round(initial_model$BIC, 2))
    message("  Class proportions: ",
            paste(paste0("Class", seq_len(k), "=",
                         round(class_proportions * 100, 1), "%"),
                  collapse = ", "))
  }

  # Step 3: Generate plots
  if (verbose) message("Generating diagnostic plots...")
  plots <- list()

  # Default colors
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
              "#FF7F00", "#FFFF33", "#A65628")[seq_len(k)]

  # 3a: Spaghetti plot (all individuals, colored by sex if available)
  if (!is.null(sex_var) && sex_var %in% names(data)) {
    plots$spaghetti <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = .data[[time_var]], y = .data[[outcome]],
                   group = .data[[id_var]], color = factor(.data[[sex_var]]))
    ) +
      ggplot2::geom_line(alpha = 0.3) +
      ggplot2::labs(title = "Individual Trajectories",
                    x = time_var, y = outcome, color = sex_var) +
      ggplot2::theme_minimal()
  } else {
    plots$spaghetti <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = .data[[time_var]], y = .data[[outcome]],
                   group = .data[[id_var]])
    ) +
      ggplot2::geom_line(alpha = 0.2, color = "grey50") +
      ggplot2::labs(title = "Individual Trajectories",
                    x = time_var, y = outcome) +
      ggplot2::theme_minimal()
  }

  # 3b: Loess plot (optionally faceted by sex)
  if (!is.null(sex_var) && sex_var %in% names(data)) {
    plots$loess <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = .data[[time_var]], y = .data[[outcome]])
    ) +
      ggplot2::geom_point(alpha = 0.1, size = 0.5) +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x,
                           color = "blue", se = TRUE) +
      ggplot2::facet_wrap(stats::as.formula(paste("~", sex_var))) +
      ggplot2::labs(title = "Loess Smoothed Trajectories",
                    x = time_var, y = outcome) +
      ggplot2::theme_minimal()
  } else {
    plots$loess <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = .data[[time_var]], y = .data[[outcome]])
    ) +
      ggplot2::geom_point(alpha = 0.1, size = 0.5) +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x,
                           color = "blue", se = TRUE) +
      ggplot2::labs(title = "Loess Smoothed Trajectory",
                    x = time_var, y = outcome) +
      ggplot2::theme_minimal()
  }

  # 3c: Residual plots per class
  # Residuals are in the $pred data frame (not standalone vectors).
  # pred_df may have fewer rows than data if hlme dropped NA rows,
  # so we must use pred_df for both residuals AND reconstruct time from it.
  pred_df <- initial_model$pred

  resid_vals <- if ("resid_ss" %in% names(pred_df)) {
    pred_df$resid_ss
  } else {
    # Fallback: marginal residuals
    pred_df$obs - pred_df$pred_m
  }

  # Reconstruct time: subset data to match pred_df rows (accounting for NA removal)
  if (!is.null(initial_model$na.action)) {
    used_data <- data[-initial_model$na.action, , drop = FALSE]
  } else {
    used_data <- data
  }

  resid_data <- data.frame(
    id = pred_df[[1]],  # subject ID is first column of pred
    time = used_data[[time_var]],
    residual = resid_vals
  )

  # Merge class assignments
  pprob_table <- initial_model$pprob
  class_df <- data.frame(
    id = pprob_table[[1]],
    class = factor(pprob_table$class)
  )
  names(class_df)[1] <- id_var
  names(resid_data)[1] <- id_var

  resid_data_full <- merge(
    resid_data,
    class_df,
    by = id_var,
    all.x = TRUE
  )
  resid_data_full <- resid_data_full[!is.na(resid_data_full$residual), ]

  plots$residuals <- ggplot2::ggplot(
    resid_data_full,
    ggplot2::aes(x = .data$time, y = .data$residual)
  ) +
    ggplot2::geom_point(alpha = 0.3, size = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "blue",
                         formula = y ~ x) +
    ggplot2::facet_wrap(~ class, labeller = ggplot2::label_both) +
    ggplot2::labs(title = "Standardized Residuals by Class",
                  subtitle = "Review pattern to choose random effects for lctm_refine()",
                  x = time_var, y = "Standardized Residual") +
    ggplot2::theme_minimal()

  # 3d: Reference guide — synthetic example residual patterns
  plots$guide <- .create_residual_guide()

  # Print residual interpretation guide
  if (verbose) {
    message("\n=== Residual Interpretation Guide ===")
    message("  Horizontal band       -> random intercept only (random = ~ 1)")
    message("  Diagonal/linear trend  -> add linear random slope (random = ~ 1 + ", time_var, ")")
    message("  Curved pattern        -> add quadratic random effect (random = ~ 1 + ", time_var, " + I(", time_var, "^2))")
    message("  S-shaped pattern      -> consider cubic terms or splines")
    message("\nUse plot() to view diagnostic plots, then pass to lctm_refine().")
  }

  # Step 4: Save PDF if requested
  if (!is.null(save_pdf)) {
    if (verbose) message("Saving plots to: ", save_pdf)
    grDevices::pdf(save_pdf, width = 10, height = 7)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (p in plots) {
      print(p)
    }
  }

  # Return lctm_initial object
  new_lctm_initial(
    data = data,
    outcome = outcome,
    time_var = time_var,
    id_var = id_var,
    sex_var = sex_var,
    degree = as.integer(degree),
    k = k,
    knots = knots,
    spline_degree = as.integer(spline_degree),
    initial_model = initial_model,
    base_model = base_model,
    plots = plots
  )
}

#' Build trajectory formulas from degree and knots
#' @keywords internal
.build_trajectory_formulas <- function(outcome, time_var, degree, knots = NULL,
                                        spline_degree = 3) {
  if (!is.null(knots)) {
    # Spline-based formula with knots. spline_degree controls the polynomial
    # order of the spline pieces:
    #   3 = natural cubic spline (splines::ns), piecewise cubic, C^2 continuity
    #   2 = quadratic B-spline (splines::bs, degree = 2), piecewise quadratic,
    #       C^1 continuity; fewer parameters per class than cubic.
    knot_str <- paste(knots, collapse = ", ")
    if (spline_degree == 2) {
      spline_term <- paste0("splines::bs(", time_var,
                            ", degree = 2, knots = c(", knot_str, "))")
      # The quadratic B-spline basis already spans the linear term, so we do
      # not add a separate time_var term (it would be collinear).
      fixed_terms <- paste0("1 + ", spline_term)
    } else {
      spline_term <- paste0("splines::ns(", time_var, ", knots = c(", knot_str, "))")
      fixed_terms <- paste0("1 + ", time_var, " + ", spline_term)
    }
    mixture_terms <- fixed_terms
  } else if (degree == 3) {
    fixed_terms <- paste0("1 + ", time_var, " + I(", time_var, "^2) + I(", time_var, "^3)")
    mixture_terms <- fixed_terms
  } else if (degree == 2) {
    fixed_terms <- paste0("1 + ", time_var, " + I(", time_var, "^2)")
    mixture_terms <- fixed_terms
  } else {
    fixed_terms <- paste0("1 + ", time_var)
    mixture_terms <- fixed_terms
  }

  list(
    fixed = stats::as.formula(paste0(outcome, " ~ ", fixed_terms)),
    mixture = stats::as.formula(paste0("~ ", mixture_terms)),
    terms = fixed_terms
  )
}

#' Create reference guide plot with example residual patterns
#' @keywords internal
.create_residual_guide <- function() {
  # Generate synthetic data for 4 example patterns
  set.seed(42)
  n <- 200
  x <- seq(0, 10, length.out = n)

  guide_data <- data.frame(
    x = rep(x, 4),
    y = c(
      # Horizontal band — random intercept only
      stats::rnorm(n, 0, 1),
      # Diagonal trend — needs linear random slope
      0.3 * x + stats::rnorm(n, 0, 0.5),
      # Curved — needs quadratic random effect
      0.1 * (x - 5)^2 - 2.5 + stats::rnorm(n, 0, 0.5),
      # S-shaped — needs cubic or splines
      2 * sin(x * 0.6) + stats::rnorm(n, 0, 0.5)
    ),
    pattern = factor(rep(
      c("A: Horizontal band\n-> Random intercept only",
        "B: Diagonal/linear trend\n-> Add linear random slope",
        "C: Curved pattern\n-> Add quadratic random effect",
        "D: S-shaped pattern\n-> Consider cubic or splines"),
      each = n
    ), levels = c(
      "A: Horizontal band\n-> Random intercept only",
      "B: Diagonal/linear trend\n-> Add linear random slope",
      "C: Curved pattern\n-> Add quadratic random effect",
      "D: S-shaped pattern\n-> Consider cubic or splines"
    ))
  )

  ggplot2::ggplot(guide_data, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.8, color = "grey40") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "blue",
                         formula = y ~ x) +
    ggplot2::facet_wrap(~ pattern, scales = "free_y", ncol = 2) +
    ggplot2::labs(
      title = "Residual Pattern Reference Guide",
      subtitle = "Compare your residual plots above to these example patterns",
      x = "Time", y = "Standardized Residual"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11)
    )
}
