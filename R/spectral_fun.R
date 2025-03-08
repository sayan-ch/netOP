#' Title
#'
#' @param A
#' @param K
#' @param symmetric
#' @param fast
#'
#' @returns
#' @export
#'
#' @examples
eigen_decomp <- function( A, K = min(dim(A)),
                          symmetric = Matrix::isSymmetric(A), fast = TRUE ) {

  # Check if K is a positive integer
  if ( !is.numeric(K) || length(K) != 1 || K <= 0 || K != floor(K) ) {
    warning("Warning: K must be a positive integer. Full decomposition is done.")
    K <- min(dim(A))
    fast <- FALSE
  }

  # Check if A is a square matrix
  if (dim(A)[1] != dim(A)[2] ) {
    warning("Warning: A must be a square matrix for eigen decomposition. SVD with top K singular values returned. Vectors return the left singular vectors")
    tmp.out <- sv_decomp(A, K = K, fast = fast)
    return( list(values = tmp.out$values, vectors = tmp.out$u,
                 u = tmp.out$u, v = tmp.out$v) )
  }

  # Compute the eigenvalues and eigenvectors of A
  if(!fast || K == min(dim(A))) {
    eig <- eigen(A, symmetric = symmetric)

    #sort the eigenvalues and eigenvectors in decreasing order
    idx <- order(eig$values, decreasing = TRUE)
    eig$values <- eig$values[idx]
    eig$vectors <- eig$vectors[, idx]

    # Return the top K eigenvalues and eigenvectors
    eig$values <- eig$values[]
    eig$vectors <- eig$vectors[]
  } else {
    if(symmetric) {
      eig <- RSpectra::eigs_sym(A, k = K, which = "LM")
    } else {
    eig <- RSpectra::eigs(A, k = K, which = "LM")
    }
  }

  return(list(values = eig$values, vectors = eig$vectors))
}

#' Title
#'
#' @param A
#' @param K
#' @param nu
#' @param nv
#' @param fast
#'
#' @returns
#' @export
#'
#' @examples
sv_decomp <- function( A, K = min(dim(A)), nu = K, nv = K,
                       fast = TRUE ){

  # Check if K is a positive integer
  if ( !is.numeric(K) || length(K) != 1 || K <= 0 || K != floor(K) ) {
    warning("Warning: K must be a positive integer. Full decomposition is done.")
    K <- min(dim(A))
    fast <- FALSE
  }

  # Compute the singular value decomposition of A
  if(!fast || K == min(dim(A))) {
    svd <- svd(A)
    svd$d <- svd$d
    svd$u <- svd$u
    svd$v <- svd$v

    idx <- order(svd$d, decreasing = TRUE)
    svd$d <- svd$d[idx]
    svd$u <- svd$u[, idx]
    svd$v <- svd$v[, idx]
  } else {
    svd <- RSpectra::svds(A, k = K)
  }

  return(list(values = svd$d, u = svd$u, v = svd$v))
}

################################################################################
################################################################################

#' Title
#'
#' @param A
#' @param K
#' @param laplacian
#' @param tau
#' @param row.norm
#' @param spherical
#' @param fast
#'
#' @returns
#' @export
#'
#' @examples
spectral.cluster <- function(A, K, laplacian = F, tau = 0, row.norm = F,
                             spherical = F, fast = T){
  n <- nrow(A)
  nc <- ncol(A)

  deg <- Matrix::rowSums(A)
  avg.deg <- mean(deg)

  A.tau <- A + tau*avg.deg/(n-1)

  L.tau <- A.tau
  if(laplacian){
    D.pow.neg12.tau <- Matrix::sparseMatrix(i = 1:n, j = 1:n,
                                  x = 1/sqrt(deg + tau*avg.deg*n/(n-1)),
                                  dims = c(n, n))
    L.tau <- Matrix::tcrossprod(Matrix::crossprod(D.pow.neg12.tau, A.tau),
                        D.pow.neg12.tau[1:nc, 1:nc])
  }

  con.nodes <- (deg != 0)

  L.tau.con <- L.tau[con.nodes, intersect(which(con.nodes), 1:nc)]
  deg.con <- deg[con.nodes]

  out.lab <- vector(length = n)
  out.lab[!con.nodes] <- sample(1:K, sum(!con.nodes), replace = T)

  K.safe <- min(K, ncol(L.tau.con))
  eig <- eigen_decomp(L.tau.con, K = K.safe, fast = fast)

  eivec.rn <- eig$vectors[, 1:K.safe]
  if(row.norm){
    eivec.rn <- eivec.rn / sqrt(Matrix::rowSums(eivec.rn^2))
  }

  # safety
  eivec.rn[is.na(eivec.rn)] <- 0

  if(!spherical){
    out.lab[con.nodes] <- as.integer(kmeans(eivec.rn, K, nstart = 100,
                                            iter.max = 10^7)$cluster)
  } else {
    out.lab[con.nodes] <- as.integer(cluster::pam(eivec.rn, K,
                                         metric = "manhattan",
                                         do.swap = F, cluster.only = T,
                                         pamonce = 6,
                                         nstart = 100))
  }

  return(out.lab)

}


################################################################################
################################################################################
























