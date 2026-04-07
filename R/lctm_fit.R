#' Fit LCTM Model (Step 3)
#'
#' Fits a latent class trajectory model with specified number of classes
#' and model type.
#'
#' @note Consider using the new investigator-driven workflow instead:
#'   `lctm_clean() -> lctm_initial() -> lctm_refine()`.
#'   This function remains available for backward compatibility.
#'
#' @param setup An `lctm_setup` object from [lctm_setup()].
#' @param k Integer; number of latent classes to fit.
#' @param model Character; model type - "A" (common variance) or "B" (proportional variance).
#'   Default is "B".
#' @param use_splines Logical; if TRUE, adds natural splines to the fixed and mixture
#'   effects. Can help capture non-polynomial trajectories.
#' @param spline_df Integer; degrees of freedom for natural splines (default 3).
#'   Only used if use_splines = TRUE.
#' @param linear Logical; if TRUE, uses linear model instead of quadratic.
#'   Overrides the degree setting in setup.
#' @param verbose Logical; if TRUE, prints progress messages.
#'
#' @return An `lctm_model` object containing:
#' \describe{
#'   \item{model}{The fitted lcmm hlme model object}
#'   \item{k}{Number of classes}
#'   \item{model_type}{Model type ("A" or "B")}
#'   \item{bic}{BIC value}
#'   \item{class_proportions}{Proportion of observations in each class}
#'   \item{setup}{The lctm_setup object used}
#'   \item{use_splines}{Whether splines were used}
#'   \item{linear}{Whether linear (vs quadratic) was used}
#' }
#'
#' @details
#' This function performs Step 3 of the LCTM workflow - fitting a model with
#' a specific number of classes.
#'
#' **Model Types:**
#' \describe{
#'   \item{Model A}{Random quadratic with common variance across classes (nwg = FALSE)}
#'   \item{Model B}{Random quadratic with proportional variance (nwg = TRUE)}
#' }
#'
#' Model B is generally preferred as it allows more flexibility in variance structure.
#'
#' **Splines:**
#' Adding natural splines (`use_splines = TRUE`) can help when trajectories
#' don't follow a simple polynomial pattern. The splines are added to both
#' the fixed and mixture effects, but not to random effects.
#'
#' @examples
#' \dontrun{
#' data(sample_growth)
#' setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid")
#'
#' # Fit Model B with 3 classes
#' model <- lctm_fit(setup, k = 3, model = "B")
#'
#' # Try with splines if adequacy fails
#' model_splines <- lctm_fit(setup, k = 3, model = "B", use_splines = TRUE)
#' }
#'
#' @export
lctm_fit <- function(setup,
                     k,
                     model = "B",
                     use_splines = FALSE,
                     spline_df = 3,
                     linear = FALSE,
                     verbose = TRUE) {

  # Validate inputs
  if (!inherits(setup, "lctm_setup")) {
    stop("setup must be an lctm_setup object from lctm_setup()", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("k must be a positive integer", call. = FALSE)
  }
  k <- as.integer(k)

  if (!model %in% c("A", "B")) {
    stop("model must be 'A' or 'B'", call. = FALSE)
  }

  # Determine nwg based on model type
  # Model A: nwg = FALSE (common variance)
  # Model B: nwg = TRUE (proportional variance)
  nwg <- (model == "B")

  # Extract info from setup
  data <- setup$data
  outcome <- setup$outcome
  time_var <- setup$time_var
  id_var <- setup$id_var
  base_model <- setup$base_model

  # Determine degree (linear overrides setup$degree)
  degree <- if (linear) 1 else setup$degree

  # Build formulas
  if (degree == 2) {
    base_terms <- paste0("1 + ", time_var, " + I(", time_var, "^2)")
    random_terms <- base_terms
  } else {
    base_terms <- paste0("1 + ", time_var)
    random_terms <- base_terms
  }

  # Add splines to fixed and mixture if requested
  if (use_splines) {
    spline_term <- paste0(" + splines::ns(", time_var, ", df = ", spline_df, ")")
    fixed_terms <- paste0(base_terms, spline_term)
    mixture_terms <- paste0(base_terms, spline_term)
  } else {
    fixed_terms <- base_terms
    mixture_terms <- base_terms
  }

  fixed_formula <- stats::as.formula(paste0(outcome, " ~ ", fixed_terms))
  mixture_formula <- stats::as.formula(paste0("~ ", mixture_terms))
  random_formula <- stats::as.formula(paste0("~ ", random_terms))

  # Determine if we need a new base model
  # (when linear or splines differ from setup's base model)
  needs_new_base <- linear || use_splines || (setup$degree != degree)

  if (needs_new_base && k > 1) {
    if (verbose) message("Fitting new base model for starting values...")
    # Use do.call to ensure formulas are stored in model$call (needed for predictY)
    base_args <- list(
      fixed = fixed_formula,
      random = random_formula,
      ng = 1,
      idiag = FALSE,
      data = data,
      subject = id_var
    )
    base_model <- do.call(lcmm::hlme, base_args)
  }

  if (verbose) {
    message("Fitting Model ", model, " with K = ", k, " classes",
            if (use_splines) " (with splines)" else "",
            if (linear) " (linear)" else "")
  }

  # Fit the model using do.call to ensure formulas are stored in model$call
  # This is required for lcmm::predictY to work correctly
  if (k == 1) {
    # For K=1, fit without mixture term
    hlme_args <- list(
      fixed = fixed_formula,
      random = random_formula,
      ng = 1,
      idiag = FALSE,
      data = data,
      subject = id_var
    )
  } else {
    # For K > 1, use mixture term and starting values from base model
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
  }

  hlme_model <- do.call(lcmm::hlme, hlme_args)

  # Check convergence
  if (hlme_model$conv != 1) {
    warning("Model did not converge (conv = ", hlme_model$conv, ")")
  }

 # Calculate class proportions
  if (k == 1) {
    class_proportions <- c(Class1 = 1.0)
  } else {
    class_counts <- table(factor(hlme_model$pprob$class, levels = seq_len(k)))
    class_proportions <- as.numeric(prop.table(class_counts))
    names(class_proportions) <- paste0("Class", seq_len(k))
  }

  if (verbose) {
    message("  BIC = ", round(hlme_model$BIC, 2))
    message("  Class proportions: ",
            paste(names(class_proportions), "=",
                  round(class_proportions * 100, 1), "%", collapse = ", "))
  }

  # Create and return lctm_model object
  result <- new_lctm_model(
    model = hlme_model,
    k = k,
    model_type = model,
    bic = hlme_model$BIC,
    class_proportions = class_proportions,
    setup = setup,
    use_splines = use_splines,
    linear = linear
  )

  validate_lctm_model(result)
  return(result)
}
