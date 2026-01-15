#' Sample Pediatric Growth Data
#'
#' A dataset containing longitudinal weight measurements for pediatric growth
#' trajectory modeling. Contains data from 290 children measured at multiple
#' timepoints.
#'
#' @name sample_growth
#' @docType data
#' @format A data frame with 935 rows and 6 variables:
#' \describe{
#'   \item{childid}{Unique identifier for each child}
#'   \item{timepoint}{Measurement timepoint label (birth, 0, 3, 6)}
#'   \item{weight_raw}{Raw weight measurement in kilograms}
#'   \item{anthroage}{Age in months at measurement (anthropometric age)}
#'   \item{waz}{Weight-for-age z-score}
#'   \item{month}{Month of measurement}
#' }
#'
#' @source Mazira study data (anonymized)
#' @keywords datasets
#'
#' @examples
#' data(sample_growth)
#' head(sample_growth)
#'
#' # Basic summary
#' summary(sample_growth$weight_raw)
#'
#' # Number of children
#' length(unique(sample_growth$childid))
NULL
