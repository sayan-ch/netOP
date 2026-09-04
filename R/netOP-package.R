#' Network Operations (netOP): Efficient Tools and Overlapping Partitions based Methods for Large Networks
#'
#' netOP provides network generators, losses, graph utilities, spectral and
#' latent-space embeddings, SBM/DCBM estimators, SONNET scalable clustering,
#' and NETCROP, ECV, NCV, and DKEST model-selection workflows.
#'
#' The package is licensed GPL (>= 2). Its internal ECV implementation is
#' derived from CRAN `randnet` 1.0; see `LICENSE.note`, `inst/COPYRIGHTS`, and
#' `citation("netOP")` for licensing and scholarly attribution.
#'
#' @keywords internal
#' @useDynLib netOP, .registration = TRUE
#' @importFrom Matrix Matrix
#' @importFrom parallel mclapply
#' @importFrom Rcpp evalCpp
#' @importFrom stats kmeans
#' @name netOP-package
"_PACKAGE"

#' Matrix summary generics
#'
#' These generics are re-exported from Matrix so that common summary
#' operations dispatch correctly for sparse network matrices after loading
#' netOP.
#'
#' @param ... Objects passed to the selected Matrix-aware generic.
#' @param na.rm Whether missing values should be removed by [sum()].
#' @name matrix-generics
#' @aliases mean diag rowMeans rowSums colMeans colSums sum
#' @importFrom Matrix mean diag rowMeans rowSums colMeans colSums
#' @export mean
#' @export diag
#' @export rowMeans
#' @export rowSums
#' @export colMeans
#' @export colSums
#' @export sum
sum <- function(..., na.rm = FALSE) {
  base::sum(..., na.rm = na.rm)
}

utils::globalVariables(c(
  "algorithm", "average_loss", "d", "K", "lower", "lower_percent",
  "mean_accuracy_percent", "mean_loss", "method", "model", "plotted_loss",
  "repetition", "tau", "upper", "upper_percent"
))
