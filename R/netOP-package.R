#' Network Data Operations and Overlapping Partitions Based Methods for Large Networks
#'
#' Tools for generating, estimating, embedding, clustering, and
#' selecting statistical network models. Provides 'SONNET', a scalable
#' subsampling-based divide-and-conquer method for community detection
#' described by Chakrabarty, Sengupta and Chen (2025)
#' <doi:10.5705/ss.202022.0108>, and 'NETCROP', an
#' overlapping-partition framework for network cross-validation, model
#' selection, and regularization tuning described by Chakrabarty, Sengupta
#' and Chen (2026) <doi:10.48550/arXiv.2504.06903>. Also includes
#' network generators, estimators, spectral and latent-space methods,
#' embedders, loss functions, and supporting network-analysis utilities.
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
#' The functions `mean()`, `diag()`, `rowMeans()`, `rowSums()`, `colMeans()`,
#' and `colSums()` are re-exported from Matrix. The `sum()` wrapper delegates
#' to [base::sum()], which dispatches to Matrix methods for sparse inputs.
#' These functions support common summaries of sparse network matrices.
#'
#' @param ... Objects passed to the selected Matrix-aware generic.
#' @param na.rm Whether missing values should be removed by [sum()].
#' @returns `sum()` and `mean()` return a numeric scalar. `diag()` returns
#'   the diagonal vector of a matrix. `rowMeans()` and `rowSums()` return
#'   one value per row; `colMeans()` and `colSums()` return one per column.
#'   See the corresponding Matrix help for other supported inputs.
#' @examples
#' A <- Matrix::Matrix(matrix(c(0, 1, 1, 0), nrow = 2), sparse = TRUE)
#' sum(A)
#' mean(A)
#' diag(A)
#' rowSums(A)
#' colMeans(A)
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
