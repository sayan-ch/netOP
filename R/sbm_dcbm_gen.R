wrapped.dcsbm.gen <- function(n, K, B, psi = rep(1, n),
                              PI = rep(1/K, K), ncore = 1)
{
  # Ensure memory is freed after function execution by calling garbage collection
  on.exit(gc())

  # Generate community assignments for each node (g) using probabilities in PI
  g <- sample.int(K, n, TRUE, PI)

  # Initialize the scaling of psi
  # For each community, scale psi such that the maximum value for each community is 1
  psi.scale <- psi
  for(kk in 1:K) {
    # Scale psi values for each community separately
    psi.scale[g == kk] <- psi[g == kk] / max(psi[g == kk])
  }

  # Create an adjacency list using parallel processing to speed up calculations
  stor <- do.call('rbind',
                  parallel::mclapply(1:(n - 1), function(i) {
                    # Identify edges between nodes i and j based on a binomial probability
                    tmp <- which(rbinom(n - i, 1,
                                        B[g[i], g[(i + 1):n]] * psi.scale[i] * psi.scale[(i + 1):n]) == 1)

                    # If no edges are found, return NULL for this iteration
                    if (length(tmp) == 0)
                      return(NULL)
                    else
                      # Return a matrix of pairs (i, i+tmp) for each edge found
                      return(cbind(rep(i, length(tmp)), i + tmp))
                  }, mc.cores = ncore))  # Use multiple cores to speed up the loop

  # Create a sparse adjacency matrix from the list of edges (stor)
  A <- Matrix::sparseMatrix(stor[, 1], stor[, 2], x = 1, dims = c(n, n),
                            symmetric = TRUE)

  # Return a list containing the adjacency matrix (A), community memberships (g), and the scaled psi values
  return(list(A = A, g = g, psi = psi.scale))
}

#' @title SBM or DCBM Network Generator
#' @description
#' SBM or DCBM generator
#'   Generates network data from SBM or DCBM
#' @param n Number of nodes
#' @param K Number of communities
#' @param B K \eqn{\times} K community connectivity matrix
#' @param sparsity Sparsity of the network, used to compute B if B is not provided.
#' @param out_in_ratio Ratio of out-community to in-community connections, used to compute B if B is not provided. B is computed as sparsity * (diag(1 - out_in_ratio, K, K) + matrix(out_in_ratio, K, K)).
#' @param sbm Logical flag indicating whether the network is generated from SBM (TRUE) or DCBM (FALSE). If sbm = TRUE, psi is ignored and set to rep(1, n).
#' @param psi \eqn{\psi} Vector of length n containing the degree correction parameter for DCBM.
#' @param PI Vector of length K containing the community membership probabilities. Defaults to equal proportions.
#' @param ncore Number of processors to use for parallel processing. Defaults to half the number of detected cores.
#' @return A list containing the adjacency matrix (A), community memberships (g), and the scaled psi values.
#' @return A sparse adjacency matrix of size n x n. \eqn{A_{ij} \sim} Bernoulli\eqn{\left(B(g_i, g_j) \psi_i \psi_j\right)}.
#' @return g Vector of length n containing the community memberships for each node. Communities are sampled from 1:K with replacement using probabilities in PI.
#' @return psi Vector of length n containing the scaled degree correction parameter for each node. psi is scaled such that the maximum value for each community is 1.
#' @export
#' @examples
#' dcsbm_gen(n = 100, K = 5, sparsity = 0.5, out_in_ratio = 0.1, sbm = TRUE)
dcsbm_gen <- function(n, K, B = NULL, sparsity = NULL, out_in_ratio = NULL,
                      sbm = TRUE, psi = NULL,
                      PI = rep(1/K, K),
                      ncore = max(floor(parallel::detectCores()/2), 1)) {

  # Check if n is a positive integer
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("Error: n must be a positive integer.")
  }

  # Check if K is a positive integer
  if (!is.numeric(K) || length(K) != 1 || K <= 0 || K != floor(K)) {
    stop("Error: K must be a positive integer.")
  }

  # Warning if K == n
  if (K == n) {
    warning("Warning: K equals n. Community membership might not make sense in this case.")
  }

  # If B is not provided, compute it using sparsity and out_in_ratio
  if (is.null(B)) {
    if (is.null(sparsity) || is.null(out_in_ratio)) {
      stop("Error: If B is not provided, sparsity and out_in_ratio must be provided to compute B.")
    }
    B <- sparsity * (diag(1 - out_in_ratio, K, K) + matrix(out_in_ratio, K, K))
  } else {
    # Check if B is a matrix of size K x K
    if (!is.matrix(B) || nrow(B) != K || ncol(B) != K) {
      stop("Error: B must be a K x K matrix.")
    }
  }

  # Warning if sbm = TRUE and psi is provided
  if (sbm) {
    # Ignore psi and set it to rep(1, n) when sbm = TRUE
    psi <- rep(1, n)
  } else {
    # Check if psi is a numeric vector of length n (unless sbm = TRUE)
    if (!is.numeric(psi) || length(psi) != n) {
      stop("Error: psi must be a numeric vector of length n.")
    }
  }

  # Check if PI is a numeric vector of length K
  if (!is.numeric(PI) || length(PI) != K) {
    warning("Warning: PI must be a numeric vector of length K. Using default PI of rep(1/K, K).")
    PI <- rep(1/K, K)
  }

  # Check if ncore is a valid positive integer; use detectCores() - 2 if invalid
  if (!is.numeric(ncore) || length(ncore) != 1 || ncore <= 0 || ncore != floor(ncore)) {
    default_ncore <- max(floor(parallel::detectCores()/2), 1)  # Ensure ncore is at least 1
    warning(paste("Warning: ncore must be a positive integer. Using default ncore =", default_ncore, "."))
    ncore <- default_ncore
  }

  # Call the original function if all checks pass
  result <- wrapped.dcsbm.gen(n, K, B, psi, PI, ncore)

  return(result)
}
