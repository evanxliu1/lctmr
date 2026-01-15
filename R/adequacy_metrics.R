#' Adequacy Metrics for Latent Class Models
#'
#' Low-level functions for calculating model adequacy metrics including
#' Average Posterior Probability of Assignment (APPA), Odds of Correct
#' Classification (OCC), and Relative Entropy.
#'
#' @name adequacy-metrics
NULL

#' Calculate Average Posterior Probability of Assignment (APPA)
#'
#' Computes the average posterior probability of assignment for each class.
#' APPA values > 0.70 indicate adequate class separation.
#'
#' @param prob_matrix A matrix or data frame of posterior probabilities where
#'   rows are observations and columns are classes. Each row should sum to 1.
#'
#' @return A named numeric vector with APPA values for each class.
#'
#' @details
#' For each class k, APPA is calculated as the mean posterior probability

#' among observations assigned to that class:
#'
#' \deqn{APPA_k = \frac{1}{n_k} \sum_{i: c_i = k} p_{ik}}
#'
#' where \eqn{n_k} is the number of observations assigned to class k,
#' \eqn{c_i} is the assigned class for observation i, and \eqn{p_{ik}}
#' is the posterior probability of observation i belonging to class k.
#'
#' @examples
#' # Example with 3 classes
#' prob_matrix <- matrix(c(
#'   0.9, 0.05, 0.05,
#'   0.85, 0.10, 0.05,
#'   0.1, 0.8, 0.1,
#'   0.05, 0.9, 0.05,
#'   0.1, 0.1, 0.8
#' ), nrow = 5, byrow = TRUE)
#' calc_appa(prob_matrix)
#'
#' @export
calc_appa <- function(prob_matrix) {
  # Convert to matrix if data frame
 prob_matrix <- as.matrix(prob_matrix)

  # Validate input
  if (!is.numeric(prob_matrix)) {
    stop("prob_matrix must be numeric", call. = FALSE)
  }

  if (any(prob_matrix < 0, na.rm = TRUE) || any(prob_matrix > 1, na.rm = TRUE)) {
    stop("Probabilities must be between 0 and 1", call. = FALSE)
  }

  K <- ncol(prob_matrix)

  # Assign each observation to the class with highest probability
  class_assignment <- apply(prob_matrix, 1, which.max)

  # Calculate APPA for each class
  appa_values <- numeric(K)
  for (k in seq_len(K)) {
    # Get probabilities for observations assigned to class k
    class_probs <- prob_matrix[class_assignment == k, k]
    if (length(class_probs) > 0) {
      appa_values[k] <- mean(class_probs, na.rm = TRUE)
    } else {
      appa_values[k] <- NA_real_
    }
  }

  names(appa_values) <- paste0("Class", seq_len(K))
  return(appa_values)
}

#' Calculate Odds of Correct Classification (OCC)
#'
#' Computes the odds of correct classification for each class.
#' OCC values > 5.0 indicate adequate classification accuracy.
#'
#' @param prob_matrix A matrix or data frame of posterior probabilities where
#'   rows are observations and columns are classes.
#' @param class_proportions A numeric vector of class proportions (pi values).
#'   If NULL, proportions are estimated from class assignments.
#'
#' @return A named numeric vector with OCC values for each class.
#'
#' @details
#' OCC compares the odds of correct classification based on posterior
#' probabilities to the odds based on random assignment:
#'
#' \deqn{OCC_k = \frac{APPA_k / (1 - APPA_k)}{\pi_k / (1 - \pi_k)}}
#'
#' where \eqn{APPA_k} is the average posterior probability for class k
#' and \eqn{\pi_k} is the proportion of observations in class k.
#'
#' @examples
#' prob_matrix <- matrix(c(
#'   0.9, 0.05, 0.05,
#'   0.85, 0.10, 0.05,
#'   0.1, 0.8, 0.1,
#'   0.05, 0.9, 0.05,
#'   0.1, 0.1, 0.8
#' ), nrow = 5, byrow = TRUE)
#' calc_occ(prob_matrix)
#'
#' @export
calc_occ <- function(prob_matrix, class_proportions = NULL) {
  prob_matrix <- as.matrix(prob_matrix)
  K <- ncol(prob_matrix)

  # Calculate APPA
  appa_values <- calc_appa(prob_matrix)

  # Calculate class proportions if not provided
  if (is.null(class_proportions)) {
    class_assignment <- apply(prob_matrix, 1, which.max)
    class_counts <- table(factor(class_assignment, levels = seq_len(K)))
    class_proportions <- as.numeric(prop.table(class_counts))
  }

  # Validate proportions
  if (length(class_proportions) != K) {
    stop("class_proportions must have length equal to number of classes",
         call. = FALSE)
  }

  # Calculate OCC: (APPA / (1 - APPA)) / (pi / (1 - pi))
  numerator <- appa_values / (1 - appa_values)
  denominator <- class_proportions / (1 - class_proportions)

  occ_values <- numerator / denominator

  # Handle edge cases
  occ_values[is.nan(occ_values)] <- NA_real_
  occ_values[is.infinite(occ_values)] <- NA_real_

  names(occ_values) <- paste0("Class", seq_len(K))
  return(occ_values)
}

#' Calculate Relative Entropy
#'
#' Computes the relative entropy (normalized entropy) of the posterior
#' probability distribution. Values > 0.5 indicate adequate class separation.
#'
#' @param prob_matrix A matrix or data frame of posterior probabilities where
#'   rows are observations and columns are classes.
#'
#' @return A single numeric value representing relative entropy (between 0 and 1).
#'
#' @details
#' Relative entropy measures the certainty of class assignments:
#'
#' \deqn{E_k = 1 + \frac{\sum_{i=1}^{n} \sum_{k=1}^{K} p_{ik} \log(p_{ik})}{n \log(K)}}
#'
#' Values close to 1 indicate high certainty (clear class separation),
#' while values close to 0 indicate high uncertainty.
#'
#' @examples
#' # High entropy (clear separation)
#' prob_matrix <- matrix(c(
#'   0.95, 0.05,
#'   0.90, 0.10,
#'   0.10, 0.90,
#'   0.05, 0.95
#' ), nrow = 4, byrow = TRUE)
#' calc_relative_entropy(prob_matrix)
#'
#' # Lower entropy (less clear)
#' prob_matrix <- matrix(c(
#'   0.6, 0.4,
#'   0.55, 0.45,
#'   0.45, 0.55,
#'   0.4, 0.6
#' ), nrow = 4, byrow = TRUE)
#' calc_relative_entropy(prob_matrix)
#'
#' @export
calc_relative_entropy <- function(prob_matrix) {
  prob_matrix <- as.matrix(prob_matrix)

  K <- ncol(prob_matrix)
  n <- nrow(prob_matrix)

  if (K < 2) {
    stop("Need at least 2 classes to calculate entropy", call. = FALSE)
  }

  # Handle zeros by replacing with small value to avoid log(0)
  prob_matrix[prob_matrix == 0] <- .Machine$double.eps

  # Calculate entropy: 1 + sum(p * log(p)) / (n * log(K))
  entropy_sum <- sum(prob_matrix * log(prob_matrix), na.rm = TRUE)
  relative_entropy <- 1 + entropy_sum / (n * log(K))

  # Ensure result is in [0, 1]
  relative_entropy <- max(0, min(1, relative_entropy))

  return(relative_entropy)
}
