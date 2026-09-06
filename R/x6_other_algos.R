# Other network model-selection algorithms
#
# The three low-level ECV routines below are preserved from the supplied
# randnet-derived code. The corresponding upstream functions are ECV.block(),
# holdout.evaluation.fast.all(), iter.SVD.core.fast.all(),
# ECV.undirected.Rank(), and missing.undirected.Rank.fast.all() in
# randnet/R/RCode.R. Their legacy dotted identifiers are external-boundary
# exceptions to CONVENTIONS.md. Only selectable loss names and package
# namespace integration were adapted during this conversion.

#' Edge cross-validation stability selection
#'
#' Selects a block-model community count or RDPG dimension by repeated edge
#' sampling. netOP provides self-contained wrappers around an ECV
#' implementation derived from the CRAN package `randnet`; `randnet` is not a
#' package dependency and its ECV-specific helpers remain internal.
#'
#' @references
#' Li, T., Levina, E., and Zhu, J. (2020). Network cross-validation by edge
#' sampling. *Biometrika*, 107(2), 257-276. \doi{10.1093/biomet/asaa006}
#' @name ecv_stability
NULL

#' Node cross-validation stability selection
#'
#' Repeated node cross-validation for selecting the number of block-model
#' communities. This is the netOP authors' implementation of Chen and Lei
#' (2018), with numerical-stability checks and explicit failure handling.
#'
#' @references
#' Chen, K. and Lei, J. (2018). Network Cross-Validation for Determining the
#' Number of Communities in Network Data. *Journal of the American Statistical
#' Association*, 113(521), 241-251. \doi{10.1080/01621459.2016.1246365}
#' @name ncv_stability
NULL

holdout.evaluation.fast.all <- function(holdout.index, A, max.K, tau = 0,
                                        dc.est = 2, p.sample = 1,
                                        loss = c('sse', 'bin_dev', 'auc_as_loss')) {
  n <- nrow(A)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  A.new <- matrix(0, n, n)
  A.new[upper.tri(A.new)] <- A[edge.index]
  A.new[edge.index[holdout.index]] <- NA
  A.new <- A.new + t(A.new)
  degrees <- colSums(A.new, na.rm = TRUE)
  no.edge <- 0
  no.edge <- sum(degrees == 0)
  
  Omega <- which(is.na(A.new))
  non.miss <- which(!is.na(A.new))
  
  SVD.result <- iter.SVD.core.fast.all(A.new, max.K, p.sample = p.sample)
  
  dc.block.sq.err <-  dc.loglike <- roc.auc <- bin.dev <-
    block.sq.err <- impute.sq.err <- loglike <- rep(0, max.K)
  sbm.auc <- dc.auc <- rep(0, max.K)
  
  for (k in 1:max.K) {
    tmp.est <- SVD.result[[k]]
    A.approx <- tmp.est$A.thr
    
    response <- A[edge.index[holdout.index]]#A[Omega]
    predictors <- A.approx[edge.index[holdout.index]]#A.approx[Omega]
    
    trunc.predictors <- predictors
    trunc.predictors <- pmin(trunc.predictors, 1 - 1e-6)
    trunc.predictors <- pmax(trunc.predictors, 1e-6)
    
    if (k == 1) {
      pb <- (sum(A.new, na.rm = TRUE) + 1) / (sum(!is.na(A.new)) - sum(!is.na(diag(A.new))) +
                                                1)
      if (pb < 1e-6)
        pb <- 1e-6
      if (pb > 1 - 1e-6)
        pb <- 1 - 1e-6
      A.Omega <- A[Omega]
      block.sq.err[k] <- sum((pb - A[Omega])^2)
      if('bin_dev' %in% loss)
        loglike[k] <- -sum(A.Omega*log(pb)) - sum((1-A.Omega)*log(1-pb))
    }
    
    if (k == 1) {
      U.approx <- matrix(tmp.est$SVD$v, ncol = k)
    } else{
      U.approx <- tmp.est$SVD$v[, 1:k, drop = F]
      if (tau > 0) {
        A.approx <- A.approx + tau * mean(colSums(A.approx)) / n
        d.approx <- colSums(A.approx)
        L.approx <- diag(1 / sqrt(d.approx)) %*% A.approx %*% diag(1 / sqrt(d.approx))
        A.approx.svd <- irlba::irlba(L.approx, nu = k, nv = k)
        U.approx <- A.approx.svd$v[, 1:k]
      }
    }
    
    km <- kmeans(
      U.approx,
      centers = k,
      nstart = 30,
      iter.max = 30
    )
    B <- matrix(0, k, k)
    Theta <- matrix(0, n, k)
    for (i in 1:k) {
      for (j in i:k) {
        N.i <- which(km$cluster == i)
        N.j <- which(km$cluster == j)
        if (i != j) {
          B[i, j] <- B[j, i] <- (sum(A.new[N.i, N.j], na.rm = TRUE) + 1) / (sum(!is.na(A.new[N.i, N.j])) +
                                                                              1)
        } else{
          B[i, j] <- B[j, i] <- (sum(A.new[N.i, N.j], na.rm = TRUE) + 1) /
            (sum(!is.na(A.new[N.i, N.j])) - sum(!is.na(diag(A.new[N.i, N.j]))) + 1)
        }
        
      }
      Theta[N.i, i] <- 1
    }
    P.hat <- Theta %*% B %*% t(Theta)
    diag(P.hat) <- 0
    block.sq.err[k] <- sum((P.hat[Omega] - A[Omega])^2)
    P.hat.Omega <- P.hat[Omega]
    A.Omega <- A[Omega]
    P.hat.Omega <- pmax(P.hat.Omega, 1e-6)
    P.hat.Omega <- pmin(P.hat.Omega, 1 - 1e-6)
    if('bin_dev' %in% loss)
      loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
    if('auc_as_loss' %in% loss)
      sbm.auc[k] <- auc_as_loss(A.Omega, P.hat.Omega) ##SC addition
    
    #### Degree correct model
    V <- U.approx
    
    ptm <- proc.time()
    
    if (k == 1) {
      V.norms <- as.numeric(abs(V))
    } else{
      V.norms <- apply(V, 1, function(x)
        sqrt(sum(x^2)))
    }
    
    iso.index <- which(V.norms == 0)
    Psi <- V.norms
    Psi <- Psi / max(V.norms)
    inv.V.norms <- 1 / V.norms
    inv.V.norms[iso.index] <- 1
    
    V.normalized <- diag(as.numeric(inv.V.norms)) %*% V
    
    if (k == 1) {
      if (dc.est > 1) {
        B <- sum(A.new, na.rm = TRUE) + 0.01
        
        partial.d <- colSums(A.new, na.rm = TRUE)
        partial.gd <- B
        phi <- rep(0, n)
        B.g <- partial.gd
        phi <- as.numeric(partial.d / B.g)
        B <- B / p.sample
        P.hat <- t(t(matrix(B, n, n) * phi) * phi)
        diag(P.hat) <- 0
      }
      dc.block.sq.err[k] <- sum((P.hat[Omega] - A[Omega])^2)
      P.hat.Omega <- P.hat[Omega]
      A.Omega <- A[Omega]
      P.hat.Omega <- pmax(P.hat.Omega, 1e-6)
      P.hat.Omega <- pmin(P.hat.Omega, 1 - 1e-6)
      if('bin_dev' %in% loss)
        dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
      if('auc_as_loss' %in% loss)
        dc.auc[k] <- auc_as_loss(A.Omega, P.hat.Omega)
    } else{
      km <- kmeans(
        V.normalized,
        centers = k,
        nstart = 30,
        iter.max = 30
      )
      if (dc.est > 1) {
        B <- matrix(0, k, k)
        Theta <- matrix(0, n, k)
        for (i in 1:k) {
          for (j in 1:k) {
            N.i <- which(km$cluster == i)
            N.j <- which(km$cluster == j)
            B[i, j] <- sum(A.new[N.i, N.j], na.rm = TRUE) + 0.01
          }
          Theta[N.i, i] <- 1
        }
        Theta <- Matrix(Theta, sparse = TRUE)
        partial.d <- colSums(A.new, na.rm = TRUE)
        partial.gd <- colSums(B)
        phi <- rep(0, n)
        B.g <- Theta %*% partial.gd
        phi <- as.numeric(partial.d / B.g)
        B <- B / p.sample
        tmp.int.mat <- Theta * phi
        P.hat <- as.matrix(tmp.int.mat %*% B %*% t(tmp.int.mat))
        diag(P.hat) <- 0
      }
      dc.block.sq.err[k] <- sum((P.hat[Omega] - A[Omega])^2)
      P.hat.Omega <- P.hat[Omega]
      A.Omega <- A[Omega]
      P.hat.Omega <- pmax(P.hat.Omega, 1e-6)
      P.hat.Omega <- pmin(P.hat.Omega, 1 - 1e-6)
      if('bin_dev' %in% loss)
        dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
      if('auc_as_loss' %in% loss)
        dc.auc[k] <- auc_as_loss(A.Omega, P.hat.Omega)
    }
  }
  return(
    list(
      impute.sq.err = impute.sq.err,
      block.sq.err = block.sq.err,
      loglike = loglike,
      roc.auc = roc.auc,
      no.edge = no.edge,
      dc.block.sq.err = dc.block.sq.err,
      dc.loglike = dc.loglike,
      bin.dev = bin.dev,
      sbm.auc = sbm.auc,
      dc.auc = dc.auc
    )
  )
}

iter.SVD.core.fast.all <- function(A, Kmax, tol = 1e-5, max.iter = 100,
                                   sparse = T, tau = 0, p.sample = 1) {
  if(sparse)
    A <- Matrix(A, sparse = TRUE)
  avg.p <- mean(as.numeric(A), na.rm = TRUE)
  cap <- 1#kappa*avg.p
  A[which(is.na(A))] <- 0
  A <- A / p.sample
  svd.new <- irlba::irlba(A, nu = Kmax, nv = Kmax)
  result <- list()
  for (K in 1:Kmax) {
    if (K == 1) {
      A.new <- svd.new$d[1] * matrix(svd.new$u[, 1], ncol = 1) %*% t(matrix(svd.new$v[, 1], ncol =
                                                                              1))
    } else{
      A.new <- A.new + svd.new$d[K] * matrix(svd.new$u[, K], ncol = 1) %*% t(matrix(svd.new$v[, K], ncol =
                                                                                      1))
    }
    A.new.thr <- A.new
    A.new.thr <- pmax(A.new.thr, 0 + tau)
    A.new.thr <- pmin(A.new.thr, cap)
    
    tmp.SVD <- list(u = svd.new$u[, 1:K],
                    v = svd.new$v[, 1:K],
                    d = svd.new$d[1:K])
    result[[K]] <- list(
      SVD = tmp.SVD,
      A = A.new,
      A.thr = A.new.thr
    )
  }
  
  return(result)
}

ECV.BM <- function (A, max.K, cv = 3, holdout.p = 0.1, tau = 0, dc.est = 2,
                    loss = c('sse', 'bin_dev', 'auc_as_loss'), ncore = 1,
                    seed = 100){
  set.seed(seed)
  n <- nrow(A)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  holdout.index.list <- list()
  holdout.n <- floor(holdout.p * edge.n)
  for (j in 1:cv) {
    holdout.index.list[[j]] <- sample(x = edge.n, size = holdout.n)
  }
  
  result <- mclapply(holdout.index.list,
                     holdout.evaluation.fast.all,
                     A = A,
                     max.K = max.K,
                     tau = tau,
                     dc.est = dc.est,
                     p.sample = 1 - holdout.p,
                     loss = loss,
                     mc.cores = min(cv, ncore)
  )
  
  dc.block.err.mat <- dc.loglike.mat <- bin.dev.mat <- roc.auc.mat <-
    impute.err.mat <- block.err.mat <- loglike.mat <-
    matrix(0, nrow = cv, ncol = max.K)
  sbm.auc.mat <- dc.auc.mat <- matrix(0, nrow = cv, ncol = max.K)
  
  no.edge.seq <- rep(0, cv)
  Omega.list <- A.list <- Imputed.A.list <- list()
  for (b in 1:cv) {
    impute.err.mat[b, ] <- result[[b]]$impute.sq.err
    block.err.mat[b, ] <- result[[b]]$block.sq.err
    loglike.mat[b, ] <- result[[b]]$loglike
    roc.auc.mat[b, ] <- result[[b]]$roc.auc
    bin.dev.mat[b, ] <- result[[b]]$bin.dev
    no.edge.seq[b] <- result[[b]]$no.edge
    dc.block.err.mat[b, ] <- result[[b]]$dc.block.sq.err
    dc.loglike.mat[b, ] <- result[[b]]$dc.loglike
    sbm.auc.mat[b, ] <- result[[b]]$sbm.auc
    dc.auc.mat[b, ] <- result[[b]]$dc.auc
  }
  
  output <- list(
    sbm.l2.mat = block.err.mat,
    sbm.bin.dev.mat = loglike.mat,
    sbm.auc.mat = sbm.auc.mat,
    dc.l2.mat = dc.block.err.mat,
    dcbm.bindev.mat = dc.loglike.mat,
    dc.auc.mat = dc.auc.mat,
    sbm.l2 = colMeans(block.err.mat),
    sbm.bin.dev = colSums(loglike.mat),
    sbm.auc = colMeans(sbm.auc.mat),
    dcbm.l2 = colMeans(dc.block.err.mat),
    dcbm.bin.dev = colSums(dc.loglike.mat),
    dcbm.auc = colMeans(dc.auc.mat)
  )
  
  if (min(output$sbm.bin.dev) > min(output$dcbm.bin.dev)) {
    bin.dev.model <- paste("DCBM", which.min(output$dcbm.bin.dev), sep = "-")
  }else {
    bin.dev.model <- paste("SBM", which.min(output$sbm.bin.dev), sep = "-")
  }
  
  if (min(output$sbm.l2) > min(output$dcbm.l2)) {
    l2.model <- paste("DCBM", which.min(output$dcbm.l2), sep = "-")
  }else {
    l2.model <- paste("SBM", which.min(output$sbm.l2), sep = "-")
  }
  
  if (min(output$sbm.auc) > min(output$dcbm.auc)) {
    auc.model <- paste("DCBM", which.min(output$dcbm.auc), sep = "-")
  }else {
    auc.model <- paste("SBM", which.min(output$sbm.auc), sep = "-")
  }
  
  output$l2.model <- l2.model
  output$bin.dev.model <- bin.dev.model
  output$auc.model <- auc.model
  return(output)
}

#' Select block-model size with edge cross-validation
#'
#' Runs repeated edge-sampling cross-validation over every community count from
#' one through `max_K` using the preserved ECV block-model implementation.
#'
#' @details
#' The wrapper validates binary undirected loop-free input, candidate and
#' sampling feasibility, dependencies, seeds, workers, and losses. Repetitions
#' run through [uni_mclapply()] while folds within each repetition remain
#' sequential. Raw ECV results are audited and converted to the shared tidy
#' loss schema; failed repetitions either stop or are explicitly omitted.
#'
#' @param A Binary, symmetric, loop-free adjacency matrix.
#' @param max_K Maximum candidate community count; candidates are `1:max_K`.
#' @param train_proportion Proportion of dyads retained for training.
#' @param cv Positive number of edge-sampling folds.
#' @param nrep Positive number of repeated ECV runs.
#' @param tau Nonnegative block-model regularization value.
#' @param losses Canonical loss names to evaluate.
#' @param ncores Positive outer worker count; folds run sequentially.
#' @param seed Optional nonnegative reproducibility seed.
#' @param verbose Print progress messages.
#' @param force_windows Use the Windows-compatible parallel backend.
#' @param ram_check Report conservative RAM demand.
#' @param failure_handling Stop or omit failed repetitions.
#' @param retain_intermediates Retain all or minimal legacy results.
#' @return A `netcrop_blockmodel` object with `algorithm = "ECV"`, tidy losses,
#'   selections, candidate metadata, realized seeds, failures, worker counts,
#'   timing, RAM diagnostics, and optional raw output.
#' @examples
#' \donttest{
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.5, beta = 0.1,
#'                   seed = 6, ncores = 1)
#' ecv_stability_blockmodel(A, max_K = 5, cv = 2, nrep = 1,
#'                          ncores = 1, seed = 7, verbose = FALSE)
#' }
#' @references
#' Li, T., Levina, E., and Zhu, J. (2020). Network cross-validation by edge
#' sampling. *Biometrika*, 107(2), 257-276. \doi{10.1093/biomet/asaa006}
#' @seealso [ecv_stability_rdpg()], [ncv_stability_blockmodel()],
#'   [netcrop_blockmodel()]
# Run repeated ECV block-model selection with validated inputs and stable output.
#' @rdname ecv_stability_blockmodel
#' @export
ecv_stability_blockmodel <- function(
    A,
    max_K,
    train_proportion = 0.9,
    cv = 3L,
    nrep = 20L,
    tau = 0,
    losses = "sse",
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    seed = NULL,
    verbose = TRUE,
    force_windows = FALSE,
    ram_check = FALSE,
    failure_handling = c("stop", "omit"),
    retain_intermediates = c("all", "minimal")) {
  call <- match.call()
  failure_handling <- match.arg(failure_handling)
  retain_intermediates <- match.arg(retain_intermediates)
  
  required_helpers <- c(
    "uni_mclapply", "modal", "is_symmetric_matrix",
    "warn_if_insufficient_ram", "format_bytes", "offset_generator_seed"
  )
  missing_helpers <- required_helpers[!vapply(
    required_helpers,
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE, envir = environment()
  )]
  if (length(missing_helpers) > 0L) {
    stop(
      "Required internal netOP helpers are unavailable; reinstall netOP. Missing: ",
      paste(missing_helpers, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  validate_count <- function(value, name, minimum = 1L) {
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) || value < minimum || value != floor(value)) {
      stop(name, " must be one integer at least ", minimum, ".",
           call. = FALSE)
    }
    as.integer(value)
  }
  validate_flag <- function(value, name) {
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
    value
  }
  
  if (is.null(dim(A)) || length(dim(A)) != 2L ||
      nrow(A) != ncol(A) || nrow(A) < 3L) {
    stop("A must be a square matrix-like object with at least three nodes.",
         call. = FALSE)
  }
  is_numeric_matrix <- is.numeric(A) || inherits(A, "Matrix")
  stored_values <- if (inherits(A, "sparseMatrix") &&
                       "x" %in% methods::slotNames(A)) {
    methods::slot(A, "x")
  } else if (inherits(A, "sparseMatrix")) {
    rep.int(1, length(methods::slot(A, "i")))
  } else {
    as.numeric(A)
  }
  if (!is_numeric_matrix || any(!is.finite(stored_values))) {
    stop("A must contain only finite numeric values.", call. = FALSE)
  }
  if (!is_symmetric_matrix(A)) {
    stop("A must be symmetric for ECV block-model selection.", call. = FALSE)
  }
  if (any(!stored_values %in% c(0, 1))) {
    stop("A must be binary.", call. = FALSE)
  }
  A_diagonal <- if (inherits(A, "Matrix")) {
    Matrix::diag(A)
  } else {
    base::diag(A)
  }
  if (any(A_diagonal != 0)) {
    stop("A must have a zero diagonal.", call. = FALSE)
  }
  
  n <- nrow(A)
  max_K <- validate_count(max_K, "max_K")
  cv <- validate_count(cv, "cv")
  nrep <- validate_count(nrep, "nrep")
  ncores <- validate_count(ncores, "ncores")
  if (max_K >= n) {
    stop("max_K must be smaller than nrow(A) for the truncated SVD.",
         call. = FALSE)
  }
  if (length(train_proportion) != 1L || !is.numeric(train_proportion) ||
      is.na(train_proportion) || !is.finite(train_proportion) ||
      train_proportion <= 0 || train_proportion >= 1) {
    stop("train_proportion must be strictly between zero and one.",
         call. = FALSE)
  }
  if (length(tau) != 1L || !is.numeric(tau) || is.na(tau) ||
      !is.finite(tau) || tau < 0) {
    stop("tau must be one finite non-negative number.", call. = FALSE)
  }
  if (tau > 1) {
    warning(
      "tau exceeds one; the preserved lower routine clips reconstructed ",
      "values after using tau as a lower bound.",
      call. = FALSE
    )
  }
  if (!is.character(losses) || length(losses) < 1L || anyNA(losses) ||
      any(!nzchar(losses))) {
    stop("losses must contain non-empty canonical loss names.",
         call. = FALSE)
  }
  losses <- unique(losses)
  supported_losses <- c("sse", "bin_dev", "auc_as_loss")
  unsupported_losses <- setdiff(losses, supported_losses)
  if (length(unsupported_losses) > 0L) {
    stop(
      "Unsupported ECV loss(es): ",
      paste(unsupported_losses, collapse = ", "), ". Supported losses: ",
      paste(supported_losses, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if ("auc_as_loss" %in% losses &&
      !exists("auc_as_loss", mode = "function", inherits = TRUE, envir = environment())) {
    stop("auc_as_loss() must be available when that loss is requested.",
         call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required by the preserved ECV routines.",
         call. = FALSE)
  }
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required by the preserved ECV routines.",
         call. = FALSE)
  }
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
        !is.finite(seed) || seed < 0 || seed != floor(seed) ||
        seed > .Machine$integer.max) {
      stop(
        "seed must be NULL or an integer from zero through ",
        ".Machine$integer.max.",
        call. = FALSE
      )
    }
    seed <- as.integer(seed)
  }
  verbose <- validate_flag(verbose, "verbose")
  force_windows <- validate_flag(force_windows, "force_windows")
  ram_check <- validate_flag(ram_check, "ram_check")
  
  edge_count <- choose(n, 2)
  holdout_count <- floor((1 - train_proportion) * edge_count)
  if (holdout_count < 1L || holdout_count >= edge_count) {
    stop(
      "train_proportion must leave at least one held-out and one observed ",
      "unordered node pair.",
      call. = FALSE
    )
  }
  upper_values <- A[upper.tri(A)]
  edge_ones <- sum(upper_values == 1)
  edge_zeros <- length(upper_values) - edge_ones
  if ("auc_as_loss" %in% losses && (edge_ones == 0L || edge_zeros == 0L)) {
    stop(
      "auc_as_loss requires A to contain at least one edge and one non-edge ",
      "above the diagonal.",
      call. = FALSE
    )
  }
  if ("auc_as_loss" %in% losses && holdout_count < 2L) {
    stop("auc_as_loss requires at least two held-out unordered pairs.",
         call. = FALSE)
  }
  
  outer_ncores <- min(ncores, nrep)
  matrix_bytes <- 8 * as.double(n) * as.double(n)
  retained_matrix_count <- 2 * max_K
  working_matrix_count <- if (tau > 0) 12 else 9
  per_worker_bytes <- matrix_bytes *
    (retained_matrix_count + working_matrix_count)
  estimated_bytes <- per_worker_bytes * outer_ncores + matrix_bytes
  ram_formula <- paste0(
    "(", format_bytes(matrix_bytes), " per dense n x n matrix x ",
    retained_matrix_count + working_matrix_count,
    " sequential retained/working matrices) x ", outer_ncores,
    " parallel repetitions + ", format_bytes(matrix_bytes),
    " parent adjacency estimate"
  )
  ram_status <- NULL
  if (ram_check) {
    ram_status <- warn_if_insufficient_ram(
      estimated_bytes = estimated_bytes,
      operation = paste0("ECV preflight: ", ram_formula)
    )
  }
  ram_report <- list(
    formula = ram_formula,
    estimated_bytes = estimated_bytes,
    per_worker_bytes = per_worker_bytes,
    parallel_repetitions = outer_ncores,
    status = ram_status
  )
  
  matrix_components <- list(
    sse = c(SBM = "sbm.l2.mat", DCBM = "dc.l2.mat"),
    bin_dev = c(SBM = "sbm.bin.dev.mat", DCBM = "dcbm.bindev.mat"),
    auc_as_loss = c(SBM = "sbm.auc.mat", DCBM = "dc.auc.mat")
  )
  validate_raw_output <- function(value) {
    if (!is.list(value)) {
      stop("ECV.BM returned a non-list result.", call. = FALSE)
    }
    for (loss_name in losses) {
      for (component_name in unname(matrix_components[[loss_name]])) {
        component <- value[[component_name]]
        if (!is.matrix(component) ||
            !identical(dim(component), c(cv, max_K))) {
          stop(
            "ECV.BM returned malformed component '", component_name,
            "'.", call. = FALSE
          )
        }
        if (any(!is.finite(component))) {
          stop(
            "ECV.BM returned non-finite requested loss values in '",
            component_name, "'. This commonly indicates a one-class AUC ",
            "holdout or a degenerate spectral fit.",
            call. = FALSE
          )
        }
      }
    }
    invisible(TRUE)
  }
  
  auc_function_template <- if ("auc_as_loss" %in% losses) {
    get("auc", mode = "function", inherits = TRUE, envir = environment())
  } else {
    NULL
  }
  auc_loss_function_template <- if ("auc_as_loss" %in% losses) {
    get("auc_as_loss", mode = "function", inherits = TRUE, envir = environment())
  } else {
    NULL
  }
  repetition_seeds <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, size = nrep, replace = TRUE)
  } else {
    vapply(
      seq_len(nrep),
      function(repetition) {
        offset_generator_seed(seed, 100 * repetition)
      },
      integer(1)
    )
  }
  
  run_repetition <- function(repetition) {
    repetition_seed <- repetition_seeds[[repetition]]
    tryCatch({
      runtime_environment <- new.env(parent = environment(ECV.BM))
      assign("mclapply", parallel::mclapply, envir = runtime_environment)
      assign("Matrix", Matrix::Matrix, envir = runtime_environment)
      assign("kmeans", stats::kmeans, envir = runtime_environment)
      assign(
        "t",
        function(x) {
          if (inherits(x, "Matrix")) {
            return(get("t", envir = asNamespace("Matrix"))(x))
          }
          base::t(x)
        },
        envir = runtime_environment
      )
      assign(
        "which",
        function(x, arr.ind = FALSE, useNames = TRUE) {
          if (inherits(x, "Matrix")) {
            return(base::which(
              as.matrix(x), arr.ind = arr.ind, useNames = useNames
            ))
          }
          base::which(x, arr.ind = arr.ind, useNames = useNames)
        },
        envir = runtime_environment
      )
      
      svd_function <- iter.SVD.core.fast.all
      evaluation_function <- holdout.evaluation.fast.all
      ecv_function <- ECV.BM
      environment(svd_function) <- runtime_environment
      environment(evaluation_function) <- runtime_environment
      environment(ecv_function) <- runtime_environment
      assign("iter.SVD.core.fast.all", svd_function,
             envir = runtime_environment)
      assign("holdout.evaluation.fast.all", evaluation_function,
             envir = runtime_environment)
      if ("auc_as_loss" %in% losses) {
        auc_function <- auc_function_template
        auc_loss_function <- auc_loss_function_template
        environment(auc_function) <- runtime_environment
        environment(auc_loss_function) <- runtime_environment
        assign("auc", auc_function, envir = runtime_environment)
        assign("auc_as_loss", auc_loss_function,
               envir = runtime_environment)
      }
      
      value <- ecv_function(
        A = A,
        max.K = max_K,
        cv = cv,
        holdout.p = 1 - train_proportion,
        tau = tau,
        dc.est = 2,
        loss = losses,
        ncore = 1L,
        seed = repetition_seed
      )
      validate_raw_output(value)
      list(
        success = TRUE,
        repetition = repetition,
        seed = repetition_seed,
        value = value,
        message = NA_character_
      )
    }, error = function(e) {
      error_call <- conditionCall(e)
      list(
        success = FALSE,
        repetition = repetition,
        seed = repetition_seed,
        value = NULL,
        message = paste0(
          conditionMessage(e),
          if (is.null(error_call)) "" else paste0(
            " [call: ", paste(deparse(error_call), collapse = " "), "]"
          )
        )
      )
    })
  }
  
  if (verbose) {
    message(
      "Running ", nrep, " ECV repetitions with ", outer_ncores,
      " outer worker(s); each lower ECV call is sequential."
    )
  }
  elapsed <- system.time({
    repetition_results <- uni_mclapply(
      seq_len(nrep),
      run_repetition,
      ncores = outer_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      future_packages = c("Matrix", "irlba")
    )
  })
  
  success_mask <- vapply(
    repetition_results, `[[`, logical(1), "success"
  )
  failure_diagnostics <- data.frame(
    repetition = vapply(
      repetition_results[!success_mask], `[[`, integer(1), "repetition"
    ),
    seed = vapply(
      repetition_results[!success_mask], `[[`, integer(1), "seed"
    ),
    message = vapply(
      repetition_results[!success_mask], `[[`, character(1), "message"
    ),
    stringsAsFactors = FALSE
  )
  if (any(!success_mask) && failure_handling == "stop") {
    stop(
      "ECV failed in repetition(s) ",
      paste(failure_diagnostics$repetition, collapse = ", "), ": ",
      paste(failure_diagnostics$message, collapse = " | "),
      call. = FALSE
    )
  }
  valid_results <- repetition_results[success_mask]
  if (length(valid_results) == 0L) {
    stop("No valid ECV repetitions remain.", call. = FALSE)
  }
  if (any(!success_mask) && failure_handling == "omit") {
    warning(
      "Omitting failed ECV repetition(s): ",
      paste(failure_diagnostics$repetition, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  cv_all_loss_parts <- list()
  part_id <- 0L
  for (result in valid_results) {
    for (loss_name in losses) {
      for (model_name in c("SBM", "DCBM")) {
        component_name <- matrix_components[[loss_name]][[model_name]]
        loss_matrix <- result$value[[component_name]]
        part_id <- part_id + 1L
        cv_all_loss_parts[[part_id]] <- data.frame(
          repetition = result$repetition,
          fold = rep(seq_len(cv), times = max_K),
          K = rep(seq_len(max_K), each = cv),
          model = model_name,
          loss = loss_name,
          loss_value = as.numeric(loss_matrix),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  cv_all_loss <- do.call(rbind, cv_all_loss_parts)
  rownames(cv_all_loss) <- NULL
  
  selection_diagnostics <- data.frame(
    loss = character(),
    model = character(),
    K = integer(),
    status = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  grouping <- interaction(
    cv_all_loss$repetition,
    cv_all_loss$K,
    cv_all_loss$model,
    cv_all_loss$loss,
    drop = TRUE,
    lex.order = TRUE
  )
  cv_loss <- do.call(rbind, lapply(split(cv_all_loss, grouping), function(x) {
    finite_values <- x$loss_value[is.finite(x$loss_value)]
    data.frame(
      repetition = x$repetition[1L],
      K = x$K[1L],
      model = x$model[1L],
      loss = x$loss[1L],
      average_loss = if (length(finite_values) == 0L) {
        NA_real_
      } else {
        sum(finite_values)
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(cv_loss) <- NULL
  cv_loss <- cv_loss[order(
    cv_loss$repetition,
    match(cv_loss$loss, losses),
    match(cv_loss$model, c("SBM", "DCBM")),
    cv_loss$K
  ), ]
  
  valid_repetitions <- vapply(
    valid_results, `[[`, integer(1), "repetition"
  )
  best_model_cv <- do.call(rbind, lapply(valid_repetitions, function(rep_id) {
    do.call(rbind, lapply(losses, function(loss_name) {
      candidates <- cv_loss[
        cv_loss$repetition == rep_id & cv_loss$loss == loss_name &
          is.finite(cv_loss$average_loss),
        ,
        drop = FALSE
      ]
      if (nrow(candidates) == 0L) {
        stop(
          "No finite candidates remain for repetition ", rep_id,
          " and loss '", loss_name, "'.",
          call. = FALSE
        )
      }
      best_id <- which.min(candidates$average_loss)
      data.frame(
        repetition = rep_id,
        loss = loss_name,
        model = candidates$model[best_id],
        K = candidates$K[best_id],
        average_loss = candidates$average_loss[best_id],
        best_model = paste0(
          candidates$model[best_id], "-", candidates$K[best_id]
        ),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(best_model_cv) <- NULL
  overall_best <- do.call(rbind, lapply(losses, function(loss_name) {
    selected <- best_model_cv[
      best_model_cv$loss == loss_name, , drop = FALSE
    ]
    data.frame(
      loss = loss_name,
      best_model = modal(selected$best_model),
      mean_average_loss = mean(selected$average_loss),
      stringsAsFactors = FALSE
    )
  }))
  rownames(overall_best) <- NULL
  
  candidate_models <- expand.grid(
    model = c("SBM", "DCBM"),
    K = seq_len(max_K),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  candidate_models <- candidate_models[order(
    match(candidate_models$model, c("SBM", "DCBM")),
    candidate_models$K
  ), ]
  
  raw_output <- lapply(valid_results, `[[`, "value")
  names(raw_output) <- paste0("repetition_", valid_repetitions)
  out <- list(
    algorithm = "ECV",
    cv_loss = cv_loss,
    cv_all_loss = cv_all_loss,
    best_model_cv = best_model_cv,
    overall_best = overall_best,
    candidate_models = candidate_models,
    nrep = nrep,
    valid_repetitions = valid_repetitions,
    completed_repetitions = length(valid_repetitions),
    cv = cv,
    train_proportion = train_proportion,
    holdout_proportion = 1 - train_proportion,
    holdout_count = holdout_count,
    losses = losses,
    model_candidates = c("SBM", "DCBM"),
    max_K = max_K,
    tau = tau,
    selection_diagnostics = selection_diagnostics,
    failure_diagnostics = failure_diagnostics,
    failure_handling = failure_handling,
    retain_intermediates = retain_intermediates,
    raw_output = if (retain_intermediates == "all") raw_output else NULL,
    ncores = list(
      requested = ncores,
      repetitions = outer_ncores,
      inner = 1L
    ),
    timing = c(total = unname(elapsed["elapsed"])),
    ram_preflight = ram_report,
    seed = seed,
    repetition_seeds = repetition_seeds,
    call = call
  )
  class(out) <- "netcrop_blockmodel"
  out
}

# Run one node-cross-validation repetition for SBM/DCBM selection.
ncv_bm <- function(
    A,
    max_K,
    cv = 3L,
    dc_est = c("spectral", "plugin"),
    tau = 0,
    use_laplacian = FALSE,
    losses = "sse",
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    seed = NULL,
    force_windows = FALSE,
    validate_inputs = TRUE) {
  dc_est <- match.arg(dc_est)
  required_helpers <- c(
    "uni_mclapply", "singular_decomp", "estimate_sbm", "estimate_dcbm",
    "sse", "bin_dev", "auc_as_loss", "offset_generator_seed",
    "is_symmetric_matrix"
  )
  missing_helpers <- required_helpers[!vapply(
    required_helpers, exists, logical(1), mode = "function", inherits = TRUE, envir = environment()
  )]
  if (length(missing_helpers) > 0L) {
    stop(
      "Required internal netOP helpers are unavailable; reinstall netOP. Missing: ",
      paste(missing_helpers, collapse = ", "), ".",
      call. = FALSE
    )
  }
  validate_count <- function(value, name, minimum = 1L) {
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) || value < minimum || value != floor(value)) {
      stop(name, " must be one integer at least ", minimum, ".",
           call. = FALSE)
    }
    as.integer(value)
  }
  validate_flag <- function(value, name) {
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
    value
  }
  max_K <- validate_count(max_K, "max_K", 2L)
  cv <- validate_count(cv, "cv", 2L)
  ncores <- validate_count(ncores, "ncores")
  use_laplacian <- validate_flag(use_laplacian, "use_laplacian")
  force_windows <- validate_flag(force_windows, "force_windows")
  validate_inputs <- validate_flag(validate_inputs, "validate_inputs")
  if (length(tau) != 1L || !is.numeric(tau) || is.na(tau) ||
      !is.finite(tau) || tau < 0) {
    stop("tau must be one finite non-negative number.", call. = FALSE)
  }
  supported_losses <- c("sse", "bin_dev", "auc_as_loss")
  if (!is.character(losses) || length(losses) < 1L || anyNA(losses) ||
      any(!nzchar(losses)) || any(!losses %in% supported_losses)) {
    stop(
      "losses must contain values from: ",
      paste(supported_losses, collapse = ", "), ".",
      call. = FALSE
    )
  }
  losses <- unique(losses)
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
        !is.finite(seed) || seed < 0 || seed != floor(seed) ||
        seed > .Machine$integer.max) {
      stop(
        "seed must be NULL or an integer from zero through ",
        ".Machine$integer.max.", call. = FALSE
      )
    }
    seed <- as.integer(seed)
  }
  if (validate_inputs) {
    if (is.null(dim(A)) || length(dim(A)) != 2L ||
        nrow(A) != ncol(A) || nrow(A) < 4L) {
      stop("A must be a square matrix-like object with at least four nodes.",
           call. = FALSE)
    }
    stored_values <- if (inherits(A, "sparseMatrix") &&
                         "x" %in% methods::slotNames(A)) {
      methods::slot(A, "x")
    } else if (inherits(A, "sparseMatrix")) {
      rep.int(1, length(methods::slot(A, "i")))
    } else {
      as.numeric(A)
    }
    if (!(is.numeric(A) || inherits(A, "Matrix")) ||
        any(!is.finite(stored_values)) ||
        any(!stored_values %in% c(0, 1))) {
      stop("A must contain only finite binary values.", call. = FALSE)
    }
    if (!is_symmetric_matrix(A)) {
      stop("A must be symmetric for NCV block-model selection.",
           call. = FALSE)
    }
    if (any(diag(A) != 0)) {
      stop("A must have a zero diagonal.", call. = FALSE)
    }
  }
  n <- nrow(A)
  fold_base_size <- floor(n / cv)
  if (fold_base_size < 2L) {
    stop("Every NCV fold must contain at least two nodes.", call. = FALSE)
  }
  fold_sizes <- rep.int(fold_base_size, cv)
  fold_sizes[cv] <- fold_sizes[cv] + n - sum(fold_sizes)
  minimum_train_size <- n - max(fold_sizes)
  if (max_K >= minimum_train_size) {
    stop(
      "max_K must be smaller than the number of training rows in every fold ",
      "(", minimum_train_size, ").",
      call. = FALSE
    )
  }
  
  if (!is.null(seed)) set.seed(seed)
  node_permutation <- sample(seq_len(n), size = n, replace = FALSE)
  fold_ends <- cumsum(fold_sizes)
  fold_starts <- c(1L, utils::head(fold_ends, -1L) + 1L)
  fold_seeds <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, size = cv, replace = TRUE)
  } else {
    vapply(
      seq_len(cv),
      function(fold_id) offset_generator_seed(seed, 100 * cv + 200 * fold_id),
      integer(1)
    )
  }
  
  evaluate_fold <- function(fold_id) {
    set.seed(fold_seeds[[fold_id]])
    fold_nodes <- node_permutation[
      seq.int(fold_starts[[fold_id]], fold_ends[[fold_id]])
    ]
    train_nodes <- setdiff(seq_len(n), fold_nodes)
    train_matrix <- A[train_nodes, , drop = FALSE]
    test_matrix <- A[fold_nodes, fold_nodes, drop = FALSE]
    if (tau > 0) {
      train_matrix <- train_matrix +
        tau * mean(rowSums(train_matrix)) / n
    }
    if (use_laplacian) {
      deg_row <- rowSums(train_matrix)
      deg_column <- colSums(train_matrix)
      if (any(!is.finite(deg_row)) || any(!is.finite(deg_column)) ||
          any(deg_row <= 0) || any(deg_column <= 0)) {
        stop(
          "Laplacian scaling encountered a non-positive row or column degree.",
          call. = FALSE
        )
      }
      train_matrix <- Matrix::Diagonal(x = 1 / sqrt(deg_row)) %*%
        train_matrix %*% Matrix::Diagonal(x = 1 / sqrt(deg_column))
    }
    decomposition <- singular_decomp(
      A = train_matrix,
      d = max_K,
      nu = max_K,
      nv = max_K,
      engine = "irlba",
      force_engine = TRUE,
      order_by = "magnitude",
      safe_d_multiplier = 1,
      validate_inputs = FALSE
    )
    train_eigenvectors <- decomposition$v[, seq_len(max_K), drop = FALSE]
    test_mask <- upper.tri(test_matrix, diag = FALSE)
    observed_test <- as.numeric(test_matrix[test_mask])
    if (length(observed_test) == 0L) {
      stop("The fold has no off-diagonal test dyads.", call. = FALSE)
    }
    
    records <- vector("list", 2L * max_K * length(losses))
    record_id <- 0L
    for (K in seq_len(max_K)) {
      if (K == 1L) {
        g_sbm <- rep.int(1L, n)
        g_dcbm <- rep.int(1L, n)
        row_norm <- NULL
      } else {
        representation <- train_eigenvectors[, seq_len(K), drop = FALSE]
        row_norm <- sqrt(rowSums(representation^2))
        safe_row_norm <- row_norm
        safe_row_norm[safe_row_norm == 0] <- 1e-6
        normalized_representation <- representation / safe_row_norm
        g_sbm <- stats::kmeans(
          representation, centers = K, nstart = 30, iter.max = 30
        )$cluster
        g_dcbm <- stats::kmeans(
          normalized_representation, centers = K,
          nstart = 30, iter.max = 30
        )$cluster
      }
      
      B_hat_sbm <- estimate_sbm(
        A = train_matrix,
        g = g_sbm,
        K = K,
        fold_nodes = fold_nodes,
        directed = FALSE,
        self_loops = FALSE,
        validate_inputs = FALSE
      )
      g_sbm_fold <- g_sbm[fold_nodes]
      P_hat_sbm <- B_hat_sbm[g_sbm_fold, g_sbm_fold, drop = FALSE]
      P_hat_sbm <- pmin(pmax(P_hat_sbm, 1e-6), 1 - 1e-6)
      diag(P_hat_sbm) <- 0
      
      dcbm_fit <- estimate_dcbm(
        A = train_matrix,
        g = g_dcbm,
        K = K,
        method = dc_est,
        fold_nodes = fold_nodes,
        row_norm = if (dc_est == "spectral" && K > 1L) row_norm else NULL,
        spectral_engine = "irlba",
        spectral_options = list(force_engine = TRUE),
        validate_inputs = FALSE
      )
      g_dcbm_fold <- g_dcbm[fold_nodes]
      psi_hat <- dcbm_fit$psi_hat[as.character(fold_nodes)]
      P_hat_dcbm <- dcbm_fit$B_hat[
        g_dcbm_fold, g_dcbm_fold, drop = FALSE
      ] * tcrossprod(psi_hat)
      P_hat_dcbm <- pmin(pmax(P_hat_dcbm, 1e-6), 1 - 1e-6)
      diag(P_hat_dcbm) <- 0
      
      predictions <- list(
        SBM = as.numeric(P_hat_sbm[test_mask]),
        DCBM = as.numeric(P_hat_dcbm[test_mask])
      )
      for (model_name in names(predictions)) {
        for (loss_name in losses) {
          loss_value <- switch(
            loss_name,
            sse = sse(observed_test, predictions[[model_name]]),
            bin_dev = bin_dev(observed_test, predictions[[model_name]]),
            auc_as_loss = auc_as_loss(
              observed_test, predictions[[model_name]]
            )
          )
          if (!is.finite(loss_value)) {
            stop(
              "Non-finite ", loss_name, " for ", model_name,
              ", K = ", K, ".",
              call. = FALSE
            )
          }
          record_id <- record_id + 1L
          records[[record_id]] <- data.frame(
            fold = fold_id,
            K = K,
            model = model_name,
            loss = loss_name,
            loss_value = loss_value,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    list(
      fold = fold_id,
      fold_nodes = fold_nodes,
      train_nodes = train_nodes,
      losses = do.call(rbind, records)
    )
  }
  
  fold_output <- uni_mclapply(
    seq_len(cv),
    evaluate_fold,
    ncores = min(ncores, cv),
    force_windows = force_windows,
    stop_on_error = TRUE,
    future_packages = "Matrix"
  )
  fold_loss <- do.call(rbind, lapply(fold_output, `[[`, "losses"))
  rownames(fold_loss) <- NULL
  grouping <- interaction(
    fold_loss$K, fold_loss$model, fold_loss$loss,
    drop = TRUE, lex.order = TRUE
  )
  cv_loss <- do.call(rbind, lapply(split(fold_loss, grouping), function(x) {
    data.frame(
      K = x$K[1L],
      model = x$model[1L],
      loss = x$loss[1L],
      average_loss = mean(x$loss_value),
      stringsAsFactors = FALSE
    )
  }))
  rownames(cv_loss) <- NULL
  cv_loss <- cv_loss[order(
    match(cv_loss$loss, losses),
    match(cv_loss$model, c("SBM", "DCBM")),
    cv_loss$K
  ), ]
  best <- do.call(rbind, lapply(losses, function(loss_name) {
    candidates <- cv_loss[cv_loss$loss == loss_name, , drop = FALSE]
    best_id <- which.min(candidates$average_loss)
    data.frame(
      loss = loss_name,
      model = candidates$model[best_id],
      K = candidates$K[best_id],
      average_loss = candidates$average_loss[best_id],
      best_model = paste0(candidates$model[best_id], "-", candidates$K[best_id]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(best) <- NULL
  list(
    fold_loss = fold_loss,
    cv_loss = cv_loss,
    best = best,
    folds = lapply(fold_output, `[[`, "fold_nodes"),
    fold_sizes = fold_sizes,
    node_permutation = node_permutation,
    fold_seeds = fold_seeds
  )
}

#' Select block-model size with node cross-validation
#'
#' Repeats node cross-validation over every community count from one through
#' `max_K` for both SBM and DCBM candidates.
#'
#' @details
#' Each repetition partitions nodes into folds, fits candidates on rectangular
#' training blocks, and evaluates held-out dyads. Repetitions are parallelized,
#' fold stages remain sequential, and probabilities are clipped consistently
#' for logarithmic losses. Failed repetitions can stop or be omitted.
#'
#' @param A Binary, symmetric, loop-free adjacency matrix.
#' @param max_K Maximum community count; candidates are `1:max_K`.
#' @param cv,nrep Positive numbers of folds and repetitions.
#' @param dc_est DCBM estimator, either `"spectral"` or `"plugin"`.
#' @param tau Nonnegative spectral regularization value.
#' @param use_laplacian Use Laplacian spectral representations.
#' @param losses Canonical loss names to evaluate.
#' @param ncores Positive outer worker count.
#' @param seed Optional nonnegative reproducibility seed.
#' @param verbose Print progress messages.
#' @param force_windows Use the Windows-compatible parallel backend.
#' @param ram_check Report conservative RAM demand.
#' @param failure_handling Stop or omit failed repetitions.
#' @param retain_intermediates Retain all or minimal repetition output.
#' @return A `netcrop_blockmodel` object with `algorithm = "NCV"`, tidy fold
#'   and repetition losses, selections, fold metadata, realized seeds,
#'   failures, worker counts, timing, and optional raw output.
#' @examples
#' \donttest{
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.55, beta = 0.08,
#'                   seed = 24, ncores = 1)
#' ncv_stability_blockmodel(A, max_K = 5, cv = 2, nrep = 1,
#'                          ncores = 1, seed = 26, verbose = FALSE)
#' }
#' @references
#' Chen, K. and Lei, J. (2018). Network Cross-Validation for Determining the
#' Number of Communities in Network Data. *Journal of the American Statistical
#' Association*, 113(521), 241-251. \doi{10.1080/01621459.2016.1246365}
#' @seealso [ecv_stability_blockmodel()], [netcrop_blockmodel()]
# Stabilize NCV block-model selection across repeated node partitions.
#' @rdname ncv_stability_blockmodel
#' @export
ncv_stability_blockmodel <- function(
    A,
    max_K,
    cv = 3L,
    nrep = 20L,
    dc_est = c("spectral", "plugin"),
    tau = 0,
    use_laplacian = FALSE,
    losses = "sse",
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    seed = NULL,
    verbose = TRUE,
    force_windows = FALSE,
    ram_check = FALSE,
    failure_handling = c("stop", "omit"),
    retain_intermediates = c("all", "minimal")) {
  call <- match.call()
  dc_est <- match.arg(dc_est)
  failure_handling <- match.arg(failure_handling)
  retain_intermediates <- match.arg(retain_intermediates)
  validate_count <- function(value, name, minimum = 1L) {
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) || value < minimum || value != floor(value)) {
      stop(name, " must be one integer at least ", minimum, ".",
           call. = FALSE)
    }
    as.integer(value)
  }
  validate_flag <- function(value, name) {
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
    value
  }
  nrep <- validate_count(nrep, "nrep")
  ncores <- validate_count(ncores, "ncores")
  verbose <- validate_flag(verbose, "verbose")
  force_windows <- validate_flag(force_windows, "force_windows")
  ram_check <- validate_flag(ram_check, "ram_check")
  
  required_helpers <- c(
    "uni_mclapply", "singular_decomp", "estimate_sbm", "estimate_dcbm",
    "sse", "bin_dev", "auc_as_loss", "modal", "offset_generator_seed",
    "is_symmetric_matrix", "format_bytes", "warn_if_insufficient_ram"
  )
  missing_helpers <- required_helpers[!vapply(
    required_helpers, exists, logical(1), mode = "function", inherits = TRUE, envir = environment()
  )]
  if (length(missing_helpers) > 0L) {
    stop(
      "Required internal netOP helpers are unavailable; reinstall netOP. Missing: ",
      paste(missing_helpers, collapse = ", "), ".",
      call. = FALSE
    )
  }
  use_laplacian <- validate_flag(use_laplacian, "use_laplacian")
  if (is.null(dim(A)) || length(dim(A)) != 2L ||
      nrow(A) != ncol(A) || nrow(A) < 4L) {
    stop("A must be a square matrix-like object with at least four nodes.",
         call. = FALSE)
  }
  stored_values <- if (inherits(A, "sparseMatrix") &&
                       "x" %in% methods::slotNames(A)) {
    methods::slot(A, "x")
  } else if (inherits(A, "sparseMatrix")) {
    rep.int(1, length(methods::slot(A, "i")))
  } else {
    as.numeric(A)
  }
  if (!(is.numeric(A) || inherits(A, "Matrix")) ||
      any(!is.finite(stored_values)) ||
      any(!stored_values %in% c(0, 1))) {
    stop("A must contain only finite binary values.", call. = FALSE)
  }
  if (!is_symmetric_matrix(A)) {
    stop("A must be symmetric for NCV block-model selection.",
         call. = FALSE)
  }
  A_diagonal <- if (inherits(A, "Matrix")) {
    Matrix::diag(A)
  } else {
    base::diag(A)
  }
  if (any(A_diagonal != 0)) {
    stop("A must have a zero diagonal.", call. = FALSE)
  }
  max_K <- validate_count(max_K, "max_K", 2L)
  cv <- validate_count(cv, "cv", 2L)
  if (floor(nrow(A) / cv) < 2L) {
    stop("Every NCV fold must contain at least two nodes.", call. = FALSE)
  }
  fold_sizes <- rep.int(floor(nrow(A) / cv), cv)
  fold_sizes[cv] <- fold_sizes[cv] + nrow(A) - sum(fold_sizes)
  if (max_K >= nrow(A) - max(fold_sizes)) {
    stop("max_K is too large for at least one NCV training matrix.",
         call. = FALSE)
  }
  if (length(tau) != 1L || !is.numeric(tau) || is.na(tau) ||
      !is.finite(tau) || tau < 0) {
    stop("tau must be one finite non-negative number.", call. = FALSE)
  }
  supported_losses <- c("sse", "bin_dev", "auc_as_loss")
  if (!is.character(losses) || length(losses) < 1L || anyNA(losses) ||
      any(!nzchar(losses)) || any(!losses %in% supported_losses)) {
    stop(
      "losses must contain values from: ",
      paste(supported_losses, collapse = ", "), ".",
      call. = FALSE
    )
  }
  losses <- unique(losses)
  if ("auc_as_loss" %in% losses) {
    edge_count <- choose(nrow(A), 2)
    observed_edges <- as.numeric(sum(A)) / 2
    if (observed_edges <= 0 || observed_edges >= edge_count) {
      stop(
        "auc_as_loss requires at least one edge and one non-edge above ",
        "the diagonal.",
        call. = FALSE
      )
    }
  }
  if (inherits(A, "sparseMatrix") && tau > 0) {
    warning(
      "Positive tau densifies each rectangular NCV training matrix.",
      call. = FALSE
    )
  }
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
       !is.finite(seed) || seed < 0 || seed != floor(seed) ||
       seed > .Machine$integer.max)) {
    stop("seed must be NULL or a valid non-negative integer.",
         call. = FALSE)
  }
  if (!is.null(seed)) seed <- as.integer(seed)
  
  outer_ncores <- min(ncores, nrep)
  n <- nrow(A)
  largest_fold <- max(fold_sizes)
  training_rows <- n - min(fold_sizes)
  training_matrix_bytes <- 8 * as.double(training_rows) * n
  test_matrix_bytes <- 8 * as.double(largest_fold)^2
  per_worker_bytes <- 5 * training_matrix_bytes + 6 * test_matrix_bytes +
    8 * as.double(n) * max_K
  estimated_bytes <- per_worker_bytes * outer_ncores
  ram_formula <- paste0(
    "(", format_bytes(training_matrix_bytes),
    " per rectangular training matrix x 5 sequential operations + ",
    format_bytes(test_matrix_bytes),
    " per fold probability/test matrix x 6 sequential operations) x ",
    outer_ncores, " parallel repetitions"
  )
  ram_status <- NULL
  if (ram_check) {
    ram_status <- warn_if_insufficient_ram(
      estimated_bytes,
      operation = paste0("NCV preflight: ", ram_formula)
    )
  }
  ram_report <- list(
    formula = ram_formula,
    estimated_bytes = estimated_bytes,
    per_worker_bytes = per_worker_bytes,
    parallel_repetitions = outer_ncores,
    status = ram_status
  )
  repetition_seeds <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, size = nrep, replace = TRUE)
  } else {
    vapply(
      seq_len(nrep),
      function(repetition) offset_generator_seed(seed, 1000 * repetition),
      integer(1)
    )
  }
  run_repetition <- function(repetition) {
    tryCatch(
      list(
        success = TRUE,
        repetition = repetition,
        seed = repetition_seeds[[repetition]],
        value = ncv_bm(
          A = A,
          max_K = max_K,
          cv = cv,
          dc_est = dc_est,
          tau = tau,
          use_laplacian = use_laplacian,
          losses = losses,
          ncores = 1L,
          seed = repetition_seeds[[repetition]],
          force_windows = FALSE,
          validate_inputs = FALSE
        ),
        message = NA_character_
      ),
      error = function(e) {
        list(
          success = FALSE,
          repetition = repetition,
          seed = repetition_seeds[[repetition]],
          value = NULL,
          message = conditionMessage(e)
        )
      }
    )
  }
  if (verbose) {
    message(
      "Running ", nrep, " NCV repetitions with ", outer_ncores,
      " outer worker(s); each fold stage is sequential."
    )
  }
  elapsed <- system.time({
    repetition_output <- uni_mclapply(
      seq_len(nrep),
      run_repetition,
      ncores = outer_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      future_packages = "Matrix"
    )
  })
  success_mask <- vapply(repetition_output, `[[`, logical(1), "success")
  failures <- repetition_output[!success_mask]
  failure_diagnostics <- data.frame(
    repetition = vapply(failures, `[[`, integer(1), "repetition"),
    seed = vapply(failures, `[[`, integer(1), "seed"),
    message = vapply(failures, `[[`, character(1), "message"),
    stringsAsFactors = FALSE
  )
  if (length(failures) > 0L && failure_handling == "stop") {
    stop(
      "NCV failed in repetition(s) ",
      paste(failure_diagnostics$repetition, collapse = ", "), ": ",
      paste(failure_diagnostics$message, collapse = " | "),
      call. = FALSE
    )
  }
  valid_output <- repetition_output[success_mask]
  if (length(valid_output) == 0L) {
    stop("No valid NCV repetitions remain.", call. = FALSE)
  }
  if (length(failures) > 0L) {
    warning(
      "Omitting failed NCV repetition(s): ",
      paste(failure_diagnostics$repetition, collapse = ", "), ".",
      call. = FALSE
    )
  }
  cv_all_loss <- do.call(rbind, lapply(valid_output, function(result) {
    transform(result$value$fold_loss, repetition = result$repetition)
  }))
  cv_all_loss <- cv_all_loss[
    , c("repetition", "fold", "K", "model", "loss", "loss_value")
  ]
  cv_loss <- do.call(rbind, lapply(valid_output, function(result) {
    transform(result$value$cv_loss, repetition = result$repetition)
  }))
  cv_loss <- cv_loss[
    , c("repetition", "K", "model", "loss", "average_loss")
  ]
  best_model_cv <- do.call(rbind, lapply(valid_output, function(result) {
    transform(result$value$best, repetition = result$repetition)
  }))
  best_model_cv <- best_model_cv[
    , c("repetition", "loss", "model", "K", "average_loss", "best_model")
  ]
  overall_best <- do.call(rbind, lapply(losses, function(loss_name) {
    selected <- best_model_cv[
      best_model_cv$loss == loss_name, , drop = FALSE
    ]
    data.frame(
      loss = loss_name,
      best_model = modal(selected$best_model),
      mean_average_loss = mean(selected$average_loss),
      stringsAsFactors = FALSE
    )
  }))
  candidate_models <- expand.grid(
    model = c("SBM", "DCBM"),
    K = seq_len(max_K),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  valid_repetitions <- vapply(valid_output, `[[`, integer(1), "repetition")
  raw_output <- lapply(valid_output, `[[`, "value")
  names(raw_output) <- paste0("repetition_", valid_repetitions)
  out <- list(
    algorithm = "NCV",
    cv_loss = cv_loss,
    cv_all_loss = cv_all_loss,
    best_model_cv = best_model_cv,
    overall_best = overall_best,
    candidate_models = candidate_models,
    nrep = nrep,
    valid_repetitions = valid_repetitions,
    completed_repetitions = length(valid_repetitions),
    cv = cv,
    fold_sizes = fold_sizes,
    losses = losses,
    model_candidates = c("SBM", "DCBM"),
    max_K = max_K,
    dc_est = dc_est,
    tau = tau,
    use_laplacian = use_laplacian,
    failure_diagnostics = failure_diagnostics,
    failure_handling = failure_handling,
    retain_intermediates = retain_intermediates,
    raw_output = if (retain_intermediates == "all") raw_output else NULL,
    ncores = list(requested = ncores, repetitions = outer_ncores, inner = 1L),
    timing = c(total = unname(elapsed["elapsed"])),
    ram_preflight = ram_report,
    seed = seed,
    repetition_seeds = repetition_seeds,
    call = call
  )
  class(out) <- "netcrop_blockmodel"
  out
}

# Preserve the pasted RDPG ECV helpers in an isolated namespace. This avoids a
# collision with the block-model helper bearing the same legacy name.
ecv_rdpg_legacy_environment <- local({
  legacy_environment <- new.env(parent = parent.frame())
  evalq({
    iter.SVD.core.fast.all <- function(A,Kmax,tol=1e-5,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE,tau=0,fast=FALSE,p.sample=1){
      if(sparse) A <- Matrix(A,sparse=TRUE)
      avg.p <- mean(as.numeric(A),na.rm=TRUE)
      cap <- 1#kappa*avg.p
      A[which(is.na(A))] <- 0
      A <- A/p.sample
      #svd.new <- svd(A,nu=K,nv=K)
      ##print("begin SVD")
      svd.new <- irlba::irlba(A,nu=Kmax,nv=Kmax)
      ##print("end SVD")
      result <- list()
      for(K in 1:Kmax){
        #print(K)
        if(K==1){
          A.new <- svd.new$d[1]*matrix(svd.new$u[,1],ncol=1)%*%t(matrix(svd.new$v[,1],ncol=1))
        }else{
          A.new <- A.new + svd.new$d[K]*matrix(svd.new$u[,K],ncol=1)%*%t(matrix(svd.new$v[,K],ncol=1))
        }
        A.new.thr <- A.new
        A.new.thr <- pmax(A.new.thr, 0 + tau)
        A.new.thr <- pmin(A.new.thr, cap)
        
        tmp.SVD <- list(u=svd.new$u[,1:K],v=svd.new$v[,1:K],d=svd.new$d[1:K])
        result[[K]] <- list(iter=NA,SVD=tmp.SVD,A=A.new,err.seq=NA,A.thr=A.new.thr)
      }
      return(result)
    }
    
    missing.undirected.Rank.weighted.fast.all <-
      function(holdout.index,A,max.K,soft=FALSE,p.sample=1, loss = loss){
        n <- nrow(A)
        #A.new <- A
        #A.new[holdout.index] <- NA
        edge.index <- which(upper.tri(A))
        edge.n <- length(edge.index)
        A.new <- matrix(0,n,n)
        A.new[upper.tri(A.new)] <- A[edge.index]
        A.new[edge.index[holdout.index]] <- NA
        A.new <- A.new + t(A.new)
        diag(A.new) <- diag(A)
        degrees <- colSums(A.new,na.rm=TRUE)
        no.edge <- 0
        no.edge <- sum(degrees==0)
        
        Omega <- which(is.na(A.new))
        imputed.A <- list()
        sse <- roc.auc <- dev <- rep(0,max.K)
        SVD.result <- iter.SVD.core.fast.all(A.new,max.K,fast=TRUE,p.sample=p.sample)
        for(k in 1:max.K){
          # print(k)
          tmp.est <- SVD.result[[k]]
          #if(k==1){
          #A.approx <- matrix(tmp.est$SVD$u,ncol=1)%*%t(matrix(tmp.est$SVD$v,ncol=1))*tmp.est$SVD$d[1]
          #}else{
          #   A.approx <- tmp.est$SVD$u%*%t(tmp.est$SVD$v*tmp.est$SVD$d)
          #}
          A.approx <- tmp.est$A
          response <- A[Omega]
          predictors <- A.approx[Omega]
          #aa <- AUC::roc(predictions=predictors,labels=factor(response))
          if('AUC' %in% loss)
            roc.auc[k] <- AUC(response, predictors)
          if('l2' %in% loss)
            sse[k] <- mean((response-predictors)^2)
          
          predictors <- pmax(predictors, 1e-6)
          predictors <- pmin(predictors, 1 - 1e-6)
          if('bin.dev' %in% loss)
            dev[k] <- bin.dev(matrix(response, ncol = k), matrix(predictors, ncol = k))
          
          imputed.A[[k]] <- A.approx
        }
        return(list(imputed.A=imputed.A,Omega=Omega, roc.auc = roc.auc, sse=sse, dev = dev))
      }
    
    ECV.undirected.Rank <- function(A,max.K,B=3,holdout.p=0.1,soft=FALSE,
                                    loss = c('l2', 'bin.dev', 'AUC'),
                                    ncore = 1, seed = NULL){
      if(!is.null(seed))
        set.seed(seed)
      
      n <- nrow(A)
      #edge.index <- 1:n^2
      #edge.n <- length(edge.index)
      edge.index <- which(upper.tri(A))
      edge.n <- length(edge.index)
      
      holdout.index.list <- list()
      
      holdout.n <- floor(holdout.p*edge.n)
      
      for(j in 1:B){
        holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
      }
      result <- mclapply(holdout.index.list,
                         missing.undirected.Rank.weighted.fast.all,
                         A=A,max.K=max.K,soft=soft,p.sample=1-holdout.p,
                         loss = loss, mc.cores = ncore)
      
      sse.mat <- roc.auc.mat <- dev.mat <- matrix(0,nrow=B,ncol=max.K)
      
      for(b in 1:B){
        roc.auc.mat[b,] <- result[[b]]$roc.auc
        sse.mat[b,] <- result[[b]]$sse
        dev.mat[b,] <- result[[b]]$dev
      }
      
      auc.seq <- colMeans(roc.auc.mat)
      #auc.sd <- apply(roc.auc.mat,2,sd)/sqrt(B)
      sse.seq <- colMeans(sse.mat)
      #sse.sd <- apply(sse.mat,2,sd)/sqrt(B)
      dev.seq <- colMeans(dev.mat)
      # return(list(sse=sse.seq,sse.sd=sse.sd))
      return(list(rank.sse=which.min(sse.seq), sse=sse.seq,
                  rank.dev = which.min(dev.seq), dev = dev.seq,
                  rank.auc=which.min(auc.seq),auc=auc.seq
      ))
    }
  }, envir = legacy_environment)
  legacy_environment
})

#' Select RDPG dimension with edge cross-validation
#'
#' Runs repeated edge-sampling cross-validation over every symmetric-RDPG
#' dimension from one through `max_d`.
#'
#' @details
#' The preserved RDPG ECV reconstruction routines run inside an isolated
#' internal environment. The wrapper validates graph and holdout feasibility,
#' supplies namespace-safe dependencies, audits all results, converts them to
#' tidy losses, and retains independently realized repetition seeds.
#'
#' @param A Binary, symmetric, loop-free adjacency matrix.
#' @param max_d Maximum candidate dimension; candidates are `1:max_d`.
#' @param cv Positive number of edge-sampling folds.
#' @param nrep Positive number of repeated ECV runs.
#' @param train_proportion Proportion of dyads retained for training.
#' @param losses Canonical loss names to evaluate: `"mse"`, `"sse"`,
#'   `"bin_dev"`, and `"auc_as_loss"`.
#' @param ncores Positive outer worker count; folds run sequentially.
#' @param seed Optional nonnegative reproducibility seed.
#' @param verbose Print progress messages.
#' @param force_windows Use the Windows-compatible parallel backend.
#' @param ram_check Report conservative RAM demand.
#' @param failure_handling Stop or omit failed repetitions.
#' @param retain_intermediates Retain all or minimal legacy results.
#' @return A `netcrop_rdpg` object with `algorithm = "ECV"`, tidy loss curves,
#'   selections, holdout metadata, realized seeds, failures, worker counts,
#'   timings, RAM diagnostics, and optional legacy output.
#' @examples
#' \donttest{
#' A <- as.matrix(generate_rdpg(n = 200, d = 3, average_degree = 8,
#'                              seed = 33, ncores = 1))
#' ecv_stability_rdpg(A, max_d = 5, cv = 2, nrep = 1,
#'                    ncores = 1, seed = 34, verbose = FALSE)
#' }
#' @references
#' Li, T., Levina, E., and Zhu, J. (2020). Network cross-validation by edge
#' sampling. *Biometrika*, 107(2), 257-276. \doi{10.1093/biomet/asaa006}
#' @seealso [ecv_stability_blockmodel()], [netcrop_rdpg()]
# Stabilize ECV dimension selection for a symmetric RDPG.
#' @rdname ecv_stability_rdpg
#' @export
ecv_stability_rdpg <- function(
    A,
    max_d,
    cv = 3L,
    nrep = 1L,
    train_proportion = 0.9,
    losses = c("mse", "bin_dev", "auc_as_loss"),
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    seed = NULL,
    verbose = TRUE,
    force_windows = FALSE,
    ram_check = FALSE,
    failure_handling = c("stop", "omit"),
    retain_intermediates = c("all", "minimal")) {
  call <- match.call()
  failure_handling <- match.arg(failure_handling)
  retain_intermediates <- match.arg(retain_intermediates)
  validate_count <- function(value, name, minimum = 1L) {
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) || value < minimum || value != floor(value)) {
      stop(name, " must be one integer at least ", minimum, ".",
           call. = FALSE)
    }
    as.integer(value)
  }
  validate_flag <- function(value, name) {
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
    value
  }
  required_helpers <- c(
    "uni_mclapply", "is_symmetric_matrix", "modal",
    "offset_generator_seed", "warn_if_insufficient_ram", "format_bytes",
    "bin_dev", "auc_as_loss"
  )
  missing_helpers <- required_helpers[!vapply(
    required_helpers, exists, logical(1), mode = "function", inherits = TRUE, envir = environment()
  )]
  if (length(missing_helpers) > 0L) {
    stop(
      "Required internal netOP helpers are unavailable; reinstall netOP. Missing: ",
      paste(missing_helpers, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!requireNamespace("Matrix", quietly = TRUE) ||
      !requireNamespace("irlba", quietly = TRUE)) {
    stop("ECV RDPG requires the Matrix and irlba packages.", call. = FALSE)
  }
  if (is.null(dim(A)) || length(dim(A)) != 2L ||
      nrow(A) != ncol(A) || nrow(A) < 3L) {
    stop("A must be a square matrix-like object with at least three rows.",
         call. = FALSE)
  }
  if (!(is.numeric(A) || inherits(A, "Matrix"))) {
    stop("A must contain only finite numeric values.", call. = FALSE)
  }
  A_values <- if (inherits(A, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(A)) {
      methods::slot(A, "x")
    } else {
      rep.int(1, length(methods::slot(A, "i")))
    }
  } else {
    as.numeric(A)
  }
  if (any(!is.finite(A_values))) {
    stop("A must contain only finite numeric values.", call. = FALSE)
  }
  if (!is_symmetric_matrix(A)) {
    stop("A must be symmetric for symmetric RDPG ECV.", call. = FALSE)
  }
  A_diagonal <- if (inherits(A, "Matrix")) Matrix::diag(A) else diag(A)
  if (any(A_diagonal != 0)) {
    stop("A must have a zero diagonal (no self-loops).", call. = FALSE)
  }
  if (any(A_values < 0)) {
    stop("A must be non-negative.", call. = FALSE)
  }
  if (sum(A_values) == 0) {
    stop("A is an empty graph; RDPG ECV is undefined.", call. = FALSE)
  }
  if (inherits(A, "symmetricMatrix")) {
    A <- methods::as(A, "generalMatrix")
  }
  max_d <- validate_count(max_d, "max_d")
  cv <- validate_count(cv, "cv")
  nrep <- validate_count(nrep, "nrep")
  ncores <- validate_count(ncores, "ncores")
  verbose <- validate_flag(verbose, "verbose")
  force_windows <- validate_flag(force_windows, "force_windows")
  ram_check <- validate_flag(ram_check, "ram_check")
  n <- nrow(A)
  if (max_d >= n) {
    stop("max_d must be smaller than nrow(A) for irlba.", call. = FALSE)
  }
  if (length(train_proportion) != 1L || !is.numeric(train_proportion) ||
      is.na(train_proportion) || !is.finite(train_proportion) ||
      train_proportion <= 0 || train_proportion >= 1) {
    stop("train_proportion must be strictly between zero and one.",
         call. = FALSE)
  }
  supported_losses <- c("mse", "sse", "bin_dev", "auc_as_loss")
  if (!is.character(losses) || length(losses) < 1L || anyNA(losses) ||
      any(!nzchar(losses)) || any(!losses %in% supported_losses)) {
    stop(
      "losses must contain only mse, sse, bin_dev, and auc_as_loss.",
      call. = FALSE
    )
  }
  losses <- unique(losses)
  binary_losses <- intersect(losses, c("bin_dev", "auc_as_loss"))
  if (length(binary_losses) > 0L && any(!A_values %in% c(0, 1))) {
    stop(
      paste(binary_losses, collapse = " and "),
      " require a binary adjacency matrix.",
      call. = FALSE
    )
  }
  edge_count <- choose(n, 2)
  holdout_proportion <- 1 - train_proportion
  holdout_count <- floor(holdout_proportion * edge_count)
  if (holdout_count < 1L || holdout_count >= edge_count) {
    stop(
      "train_proportion must leave at least one held-out and one observed ",
      "unordered pair per fold.",
      call. = FALSE
    )
  }
  if ("auc_as_loss" %in% losses) {
    observed_edges <- as.numeric(sum(A)) / 2
    if (observed_edges <= 0 || observed_edges >= edge_count) {
      stop(
        "auc_as_loss requires at least one edge and one non-edge above ",
        "the diagonal.",
        call. = FALSE
      )
    }
  }
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
       !is.finite(seed) || seed < 0 || seed != floor(seed) ||
       seed > .Machine$integer.max)) {
    stop("seed must be NULL or a valid non-negative integer.",
         call. = FALSE)
  }
  if (!is.null(seed)) seed <- as.integer(seed)
  if (inherits(A, "sparseMatrix")) {
    warning(
      "The immutable RDPG ECV helpers densify A and all rank estimates.",
      call. = FALSE
    )
  }
  
  outer_ncores <- min(ncores, nrep)
  dense_matrix_bytes <- 8 * as.double(n) * n
  per_worker_bytes <- dense_matrix_bytes * (2 * cv * max_d + 6)
  estimated_bytes <- per_worker_bytes * outer_ncores
  ram_formula <- paste0(
    format_bytes(dense_matrix_bytes), " per dense n-by-n matrix x (2 x ",
    cv, " folds x ", max_d, " ranks + 6 working matrices) x ",
    outer_ncores, " parallel repetitions"
  )
  ram_status <- NULL
  if (ram_check) {
    ram_status <- warn_if_insufficient_ram(
      estimated_bytes,
      operation = paste0("RDPG ECV preflight: ", ram_formula)
    )
  }
  ram_report <- list(
    formula = ram_formula,
    estimated_bytes = estimated_bytes,
    per_worker_bytes = per_worker_bytes,
    parallel_repetitions = outer_ncores,
    status = ram_status
  )
  legacy_losses <- unname(c(
    mse = "l2", sse = "l2", bin_dev = "bin.dev", auc_as_loss = "AUC"
  )[losses])
  legacy_losses <- unique(legacy_losses)
  repetition_seeds <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, size = nrep, replace = TRUE)
  } else {
    vapply(
      seq_len(nrep),
      function(repetition) offset_generator_seed(seed, 1000 * repetition),
      integer(1)
    )
  }
  legacy_ecv <- get(
    "ECV.undirected.Rank",
    envir = ecv_rdpg_legacy_environment,
    inherits = FALSE
  )
  run_repetition <- function(repetition) {
    tryCatch({
      runtime_environment <- new.env(parent = environment(legacy_ecv))
      assign("mclapply", parallel::mclapply, envir = runtime_environment)
      assign("Matrix", Matrix::Matrix, envir = runtime_environment)
      assign(
        "which",
        function(x, arr.ind = FALSE, useNames = TRUE) {
          base::which(as.matrix(x), arr.ind = arr.ind, useNames = useNames)
        },
        envir = runtime_environment
      )
      assign(
        "t",
        function(x) {
          if (inherits(x, "Matrix")) Matrix::t(x) else base::t(x)
        },
        envir = runtime_environment
      )
      assign(
        "diag",
        function(x) {
          if (inherits(x, "Matrix")) Matrix::diag(x) else base::diag(x)
        },
        envir = runtime_environment
      )
      assign(
        "AUC",
        function(response, predictors) auc_as_loss(response, predictors),
        envir = runtime_environment
      )
      assign(
        "bin.dev",
        function(response, predictors) {
          original_length <- if (exists(
            "response", envir = parent.frame(), inherits = FALSE
          )) {
            length(get("response", envir = parent.frame(), inherits = FALSE))
          } else {
            length(response)
          }
          bin_dev(
            as.numeric(response)[seq_len(original_length)],
            as.numeric(predictors)[seq_len(original_length)]
          )
        },
        envir = runtime_environment
      )
      legacy_function_names <- c(
        "iter.SVD.core.fast.all",
        "missing.undirected.Rank.weighted.fast.all",
        "ECV.undirected.Rank"
      )
      for (function_name in legacy_function_names) {
        legacy_function <- get(
          function_name,
          envir = ecv_rdpg_legacy_environment,
          inherits = FALSE
        )
        environment(legacy_function) <- runtime_environment
        assign(function_name, legacy_function, envir = runtime_environment)
      }
      repetition_function <- get(
        "ECV.undirected.Rank", envir = runtime_environment, inherits = FALSE
      )
      value <- repetition_function(
        A = A,
        max.K = max_d,
        B = cv,
        holdout.p = holdout_proportion,
        loss = legacy_losses,
        ncore = 1L,
        seed = repetition_seeds[[repetition]]
      )
      if ("sse" %in% losses) value$sse_sum <- value$sse * holdout_count
      component_names <- c(
        mse = "sse", sse = "sse_sum", bin_dev = "dev", auc_as_loss = "auc"
      )
      for (loss_name in losses) {
        component <- value[[component_names[[loss_name]]]]
        if (!is.numeric(component) || length(component) != max_d ||
            any(!is.finite(component))) {
          stop(
            "Legacy RDPG ECV returned malformed or non-finite ", loss_name,
            " values. This commonly indicates a one-class AUC holdout or ",
            "a degenerate spectral fit.",
            call. = FALSE
          )
        }
      }
      list(
        success = TRUE,
        repetition = repetition,
        seed = repetition_seeds[[repetition]],
        value = value,
        message = NA_character_
      )
    }, error = function(e) {
      list(
        success = FALSE,
        repetition = repetition,
        seed = repetition_seeds[[repetition]],
        value = NULL,
        message = conditionMessage(e)
      )
    })
  }
  if (verbose) {
    message(
      "Running ", nrep, " RDPG ECV repetitions with ", outer_ncores,
      " outer worker(s); each fold stage is sequential."
    )
  }
  elapsed <- system.time({
    repetition_output <- uni_mclapply(
      seq_len(nrep),
      run_repetition,
      ncores = outer_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      future_packages = c("Matrix", "irlba")
    )
  })
  success_mask <- vapply(repetition_output, `[[`, logical(1), "success")
  failures <- repetition_output[!success_mask]
  failure_diagnostics <- data.frame(
    repetition = vapply(failures, `[[`, integer(1), "repetition"),
    seed = vapply(failures, `[[`, integer(1), "seed"),
    message = vapply(failures, `[[`, character(1), "message"),
    stringsAsFactors = FALSE
  )
  if (length(failures) > 0L && failure_handling == "stop") {
    stop(
      "RDPG ECV failed in repetition(s) ",
      paste(failure_diagnostics$repetition, collapse = ", "), ": ",
      paste(failure_diagnostics$message, collapse = " | "),
      call. = FALSE
    )
  }
  valid_output <- repetition_output[success_mask]
  if (length(valid_output) == 0L) {
    stop("No valid RDPG ECV repetitions remain.", call. = FALSE)
  }
  if (length(failures) > 0L) {
    warning(
      "Omitting failed RDPG ECV repetition(s): ",
      paste(failure_diagnostics$repetition, collapse = ", "), ".",
      call. = FALSE
    )
  }
  component_names <- c(
    mse = "sse", sse = "sse_sum", bin_dev = "dev", auc_as_loss = "auc"
  )
  cv_loss <- do.call(rbind, lapply(valid_output, function(result) {
    do.call(rbind, lapply(losses, function(loss_name) {
      data.frame(
        repetition = result$repetition,
        d = seq_len(max_d),
        loss = loss_name,
        average_loss = as.numeric(
          result$value[[component_names[[loss_name]]]]
        ),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(cv_loss) <- NULL
  best_dimension_cv <- do.call(rbind, lapply(valid_output, function(result) {
    do.call(rbind, lapply(losses, function(loss_name) {
      candidates <- cv_loss[
        cv_loss$repetition == result$repetition & cv_loss$loss == loss_name,
        , drop = FALSE
      ]
      best_id <- which.min(candidates$average_loss)
      data.frame(
        repetition = result$repetition,
        loss = loss_name,
        d_hat = candidates$d[best_id],
        average_loss = candidates$average_loss[best_id],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(best_dimension_cv) <- NULL
  overall_best <- do.call(rbind, lapply(losses, function(loss_name) {
    selected <- best_dimension_cv[
      best_dimension_cv$loss == loss_name, , drop = FALSE
    ]
    data.frame(
      loss = loss_name,
      d_hat = modal(selected$d_hat),
      mean_average_loss = mean(selected$average_loss),
      stringsAsFactors = FALSE
    )
  }))
  rownames(overall_best) <- NULL
  cv_all_loss <- transform(
    cv_loss,
    fold = NA_integer_,
    loss_value = average_loss
  )[, c("repetition", "fold", "d", "loss", "loss_value")]
  valid_repetitions <- vapply(valid_output, `[[`, integer(1), "repetition")
  raw_output <- lapply(valid_output, `[[`, "value")
  names(raw_output) <- paste0("repetition_", valid_repetitions)
  out <- list(
    algorithm = "ECV",
    cv_loss = cv_loss,
    cv_all_loss = cv_all_loss,
    best_dimension_cv = best_dimension_cv,
    overall_best = overall_best,
    d_candidates = seq_len(max_d),
    nrep = nrep,
    completed_repetitions = length(valid_repetitions),
    valid_repetitions = valid_repetitions,
    cv = cv,
    train_proportion = train_proportion,
    holdout_proportion = holdout_proportion,
    holdout_count = holdout_count,
    losses = losses,
    failure_diagnostics = failure_diagnostics,
    failure_handling = failure_handling,
    retain_intermediates = retain_intermediates,
    raw_output = if (retain_intermediates == "all") raw_output else NULL,
    ncores = list(requested = ncores, repetitions = outer_ncores, inner = 1L),
    timing = c(total = unname(elapsed["elapsed"])),
    ram_preflight = ram_report,
    seed = seed,
    repetition_seeds = repetition_seeds,
    call = call
  )
  class(out) <- "netcrop_rdpg"
  out
}

# Build a comparison plot for NETCROP and ECV RDPG-selection results.
plot_rdpg_comparison <- function(..., loss_scale = c("relative", "raw")) {
  outputs <- list(...)
  loss_scale <- match.arg(loss_scale)
  if (length(outputs) != 2L ||
      any(!vapply(outputs, inherits, logical(1), "netcrop_rdpg"))) {
    stop("Supply two objects inheriting from netcrop_rdpg.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting.", call. = FALSE)
  }
  algorithms <- vapply(outputs, function(output) {
    if (is.null(output$algorithm)) "NETCROP" else toupper(output$algorithm)
  }, character(1))
  if (!setequal(algorithms, c("NETCROP", "ECV"))) {
    stop("Comparison requires one NETCROP and one ECV result.", call. = FALSE)
  }
  loss_family <- function(loss_name) {
    switch(
      loss_name,
      sse = "squared_error",
      mse = "squared_error",
      bin_dev = "binary_deviance",
      bin_dev_mean = "binary_deviance",
      loss_name
    )
  }
  comparison_names <- lapply(outputs, function(output) {
    if (loss_scale == "raw") output$losses else {
      vapply(output$losses, loss_family, character(1))
    }
  })
  common_losses <- Reduce(intersect, comparison_names)
  common_d <- Reduce(
    intersect, lapply(outputs, function(output) unique(output$cv_loss$d))
  )
  if (length(common_losses) == 0L) {
    stop(
      if (loss_scale == "raw") {
        "Raw comparison requires at least one identical loss name."
      } else {
        "The results have no comparable loss family."
      },
      call. = FALSE
    )
  }
  if (length(common_d) == 0L) {
    stop("The results have no common d candidates.", call. = FALSE)
  }
  prepare_plot_data <- function(output, algorithm) {
    data <- output$cv_loss
    data$comparison_loss <- if (loss_scale == "raw") {
      data$loss
    } else {
      vapply(data$loss, loss_family, character(1))
    }
    data <- data[
      data$comparison_loss %in% common_losses & data$d %in% common_d &
        is.finite(data$average_loss),
      , drop = FALSE
    ]
    grouping <- interaction(
      data$d, data$comparison_loss, drop = TRUE, lex.order = TRUE
    )
    result <- do.call(rbind, lapply(split(data, grouping), function(x) {
      data.frame(
        d = x$d[1L],
        loss = x$comparison_loss[1L],
        mean_loss = mean(x$average_loss),
        stringsAsFactors = FALSE
      )
    }))
    result$algorithm <- algorithm
    result
  }
  plot_data <- do.call(rbind, Map(prepare_plot_data, outputs, algorithms))
  plot_data$plotted_loss <- plot_data$mean_loss
  y_label <- "Mean CV loss"
  plot_caption <- NULL
  if (loss_scale == "relative") {
    normalization_group <- interaction(
      plot_data$algorithm, plot_data$loss, drop = TRUE, lex.order = TRUE
    )
    plot_data$plotted_loss <- unsplit(
      lapply(split(plot_data$mean_loss, normalization_group), function(x) {
        value_range <- range(x)
        if (diff(value_range) == 0) return(rep(0, length(x)))
        (x - value_range[1L]) / diff(value_range)
      }),
      normalization_group
    )
    y_label <- "Within-algorithm relative CV loss"
    plot_caption <- paste0(
      "Relative mode compares compatible families; SSE and MSE are ",
      "normalized within each algorithm."
    )
  }
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = d,
      y = plotted_loss,
      color = algorithm,
      linetype = algorithm,
      group = algorithm
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = sort(common_d)) +
    ggplot2::facet_wrap(
      ~loss, scales = if (loss_scale == "raw") "free_y" else "fixed"
    ) +
    ggplot2::labs(
      title = paste(algorithms, collapse = " vs "),
      subtitle = "Symmetric-RDPG cross-validation comparison",
      caption = plot_caption,
      x = "Latent dimension (d)",
      y = y_label,
      color = "Algorithm",
      linetype = "Algorithm"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' Tune a degree-regularized spectral model
#'
#' Selects a spectral regularizer using the DKEST eigenspectral-ratio
#' criterion.
#'
#' @details
#' Each candidate regularizes the observed matrix, optionally constructs its
#' normalized Laplacian, clusters the fixed-`K` embedding, estimates SBM or
#' DCBM probabilities, and computes the ratio between selected residual and
#' fitted eigenvalues. Candidate failures are audited explicitly.
#'
#' @param A Finite symmetric, loop-free adjacency matrix.
#' @param K Fixed positive number of communities.
#' @param tau_candidates Unique finite nonnegative regularization candidates.
#' @param use_laplacian Use a normalized graph Laplacian.
#' @param use_dcbm Fit a degree-corrected rather than ordinary block model.
#' @param dcbm_est_method DCBM estimator, either `"plugin"` or `"spectral"`.
#' @param ncores Positive worker count.
#' @param seed Optional nonnegative reproducibility seed.
#' @param verbose Print progress messages.
#' @param force_windows Use the Windows-compatible parallel backend.
#' @param ram_check Report conservative RAM demand.
#' @param failure_handling Stop or omit failed candidates.
#' @param retain_intermediates Retain all or minimal candidate fits.
#' @return A `netcrop_regularizer` object containing `tau_hat`, the selected DK
#'   statistic, per-candidate numerator, denominator, and ratio, diagnostics,
#'   seeds, worker counts, timing, RAM information, and optional raw fits.
#' @examples
#' \donttest{
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.5, beta = 0.1,
#'                   seed = 35, ncores = 1)
#' dkest_tune_regularizer(A, K = 3, tau_candidates = c(0, 0.1),
#'                        ncores = 1, seed = 36, verbose = FALSE)
#' }
#' @seealso [netcrop_tune_regularizer()], [mult_reg_spectral_cluster()]
# Select a spectral regularization parameter using the DK statistic.
#' @rdname dkest_tune_regularizer
#' @export
dkest_tune_regularizer <- function(
    A,
    K,
    tau_candidates,
    use_laplacian = TRUE,
    use_dcbm = TRUE,
    dcbm_est_method = c("plugin", "spectral"),
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    seed = NULL,
    verbose = TRUE,
    force_windows = FALSE,
    ram_check = FALSE,
    failure_handling = c("stop", "omit"),
    retain_intermediates = c("all", "minimal")) {
  call <- match.call()
  dcbm_est_method <- match.arg(dcbm_est_method)
  failure_handling <- match.arg(failure_handling)
  retain_intermediates <- match.arg(retain_intermediates)
  validate_count <- function(value, name, minimum = 1L) {
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) || value < minimum || value != floor(value)) {
      stop(name, " must be one integer at least ", minimum, ".",
           call. = FALSE)
    }
    as.integer(value)
  }
  validate_flag <- function(value, name) {
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
    value
  }
  required_helpers <- c(
    "uni_mclapply", "is_symmetric_matrix", "estimate_sbm",
    "estimate_dcbm", "offset_generator_seed", "warn_if_insufficient_ram",
    "format_bytes"
  )
  missing_helpers <- required_helpers[!vapply(
    required_helpers, exists, logical(1), mode = "function", inherits = TRUE, envir = environment()
  )]
  if (length(missing_helpers) > 0L) {
    stop(
      "Required internal netOP helpers are unavailable; reinstall netOP. Missing: ",
      paste(missing_helpers, collapse = ", "), ".",
      call. = FALSE
    )
  }
  base_packages <- c("Matrix", "RSpectra")
  missing_packages <- base_packages[!vapply(
    base_packages, requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing_packages) > 0L) {
    stop(
      "Required package(s) are not installed: ",
      paste(missing_packages, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (is.null(dim(A)) || length(dim(A)) != 2L ||
      nrow(A) != ncol(A) || nrow(A) < 2L ||
      !(is.numeric(A) || inherits(A, "Matrix"))) {
    stop("A must be a numeric square matrix with at least two nodes.",
         call. = FALSE)
  }
  A_values <- if (inherits(A, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(A)) {
      methods::slot(A, "x")
    } else {
      rep.int(1, length(methods::slot(A, "i")))
    }
  } else {
    as.numeric(A)
  }
  if (anyNA(A_values) || any(!is.finite(A_values)) || any(A_values < 0)) {
    stop("A must contain only finite non-negative values.", call. = FALSE)
  }
  if (!is_symmetric_matrix(A)) {
    stop("A must be symmetric for DK regularizer selection.", call. = FALSE)
  }
  A_diagonal <- if (inherits(A, "Matrix")) Matrix::diag(A) else diag(A)
  if (any(A_diagonal != 0)) {
    stop("A must have a zero diagonal (no self-loops).", call. = FALSE)
  }
  if (sum(A_values) == 0) {
    stop("A is an empty graph; the DK statistic is undefined.",
         call. = FALSE)
  }
  if (inherits(A, "symmetricMatrix")) {
    A <- methods::as(A, "generalMatrix")
  }
  K <- validate_count(K, "K")
  ncores <- validate_count(ncores, "ncores")
  use_laplacian <- validate_flag(use_laplacian, "use_laplacian")
  use_dcbm <- validate_flag(use_dcbm, "use_dcbm")
  verbose <- validate_flag(verbose, "verbose")
  force_windows <- validate_flag(force_windows, "force_windows")
  ram_check <- validate_flag(ram_check, "ram_check")
  n <- nrow(A)
  if (K >= n) {
    stop("K must be smaller than nrow(A) for RSpectra.", call. = FALSE)
  }
  if (!is.numeric(tau_candidates) || length(tau_candidates) < 1L ||
      anyNA(tau_candidates) || any(!is.finite(tau_candidates)) ||
      any(tau_candidates < 0) || anyDuplicated(tau_candidates)) {
    stop("tau_candidates must contain unique finite non-negative numbers.",
         call. = FALSE)
  }
  tau_candidates <- as.numeric(tau_candidates)
  worker_packages <- base_packages
  if (K > 1L && (use_dcbm || any(tau_candidates > 0))) {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("The cluster package is required for this DK configuration.",
           call. = FALSE)
    }
    worker_packages <- c(worker_packages, "cluster")
  }
  deg <- if (inherits(A, "Matrix")) Matrix::rowSums(A) else rowSums(A)
  avg_deg <- mean(deg)
  good_node_count <- sum(deg > 0)
  if (K > 1L && any(tau_candidates == 0) && good_node_count <= K) {
    stop(
      "For tau = 0 and K > 1, more than K positive-degree nodes are ",
      "required by RSpectra.",
      call. = FALSE
    )
  }
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
       !is.finite(seed) || seed < 0 || seed != floor(seed) ||
       seed > .Machine$integer.max)) {
    stop("seed must be NULL or a valid non-negative integer.",
         call. = FALSE)
  }
  if (!is.null(seed)) seed <- as.integer(seed)
  if (inherits(A, "sparseMatrix") && any(tau_candidates > 0)) {
    warning(
      "Positive tau candidates densify A during DK regularization.",
      call. = FALSE
    )
  }
  
  candidate_ncores <- min(ncores, length(tau_candidates))
  dense_matrix_bytes <- 8 * as.double(n) * n
  per_worker_bytes <- 12 * dense_matrix_bytes +
    8 * as.double(n) * max(K, 1L)
  estimated_bytes <- per_worker_bytes * candidate_ncores
  ram_formula <- paste0(
    format_bytes(dense_matrix_bytes),
    " per dense n-by-n matrix x 12 working matrices x ",
    candidate_ncores, " parallel tau candidates"
  )
  ram_status <- NULL
  if (ram_check) {
    ram_status <- warn_if_insufficient_ram(
      estimated_bytes,
      operation = paste0("DKEST preflight: ", ram_formula)
    )
  }
  ram_report <- list(
    formula = ram_formula,
    estimated_bytes = estimated_bytes,
    per_worker_bytes = per_worker_bytes,
    parallel_candidates = candidate_ncores,
    status = ram_status
  )
  candidate_seeds <- if (is.null(seed)) {
    sample.int(
      .Machine$integer.max,
      size = length(tau_candidates),
      replace = TRUE
    )
  } else {
    vapply(
      seq_along(tau_candidates),
      function(candidate_id) {
        offset_generator_seed(seed, 1000 * candidate_id)
      },
      integer(1)
    )
  }
  evaluate_candidate <- function(candidate_id) {
    tryCatch({
      set.seed(candidate_seeds[[candidate_id]])
      tau <- tau_candidates[[candidate_id]]
      A_tau <- A + tau * avg_deg / n
      d_tau <- Matrix::sparseMatrix(
        i = seq_len(n),
        j = seq_len(n),
        x = 1 / sqrt(deg + tau * avg_deg)
      )
      L_tau <- A_tau
      if (use_laplacian) {
        L_tau <- tcrossprod(crossprod(d_tau, A_tau), d_tau)
        L_tau[is.na(L_tau)] <- 0
      }
      
      if (K == 1L) {
        out_comm <- rep.int(1L, n)
      } else if (tau > 0) {
        eig_max <- RSpectra::eigs_sym(
          L_tau,
          K,
          which = "LM",
          opts = list(v0 = rep(1, n))
        )$vectors
        rn_eig <- eig_max
        if (use_dcbm) {
          row_norm <- sqrt(rowSums(eig_max^2))
          row_norm[row_norm == 0] <- 1
          rn_eig <- eig_max / row_norm
          out_comm <- as.integer(cluster::pam(
            rn_eig,
            K,
            metric = "euclidean",
            do.swap = FALSE,
            cluster.only = TRUE,
            pamonce = 6
          ))
        } else {
          out_comm <- as.integer(cluster::clara(
            eig_max,
            K,
            metric = "euclidean",
            cluster.only = TRUE
          ))
        }
      } else {
        out_comm <- integer(n)
        bad_node <- which(deg == 0)
        if (length(bad_node) > 0L) {
          out_comm[bad_node] <- sample(
            seq_len(K), length(bad_node), replace = TRUE
          )
        }
        good_node <- which(deg > 0)
        L_good <- L_tau[good_node, good_node, drop = FALSE]
        eig_max_good <- RSpectra::eigs_sym(
          L_good,
          K,
          which = "LM",
          opts = list(v0 = rep(1, length(good_node)))
        )$vectors
        rn_eig_good <- eig_max_good
        if (use_dcbm) {
          row_norm <- sqrt(rowSums(eig_max_good^2))
          row_norm[row_norm == 0] <- 1
          rn_eig_good <- eig_max_good / row_norm
          out_comm[good_node] <- as.integer(cluster::pam(
            rn_eig_good,
            K,
            metric = "euclidean",
            do.swap = FALSE,
            cluster.only = TRUE,
            pamonce = 6
          ))
        } else {
          out_comm[good_node] <- as.integer(stats::kmeans(
            eig_max_good,
            K,
            nstart = 100,
            iter.max = 10^7
          )$cluster)
        }
      }
      
      if (!is.numeric(out_comm) || length(out_comm) != n ||
          anyNA(out_comm) || any(out_comm < 1L) || any(out_comm > K)) {
        stop("Community estimation returned invalid labels.", call. = FALSE)
      }
      if (!use_dcbm) {
        B_hat <- estimate_sbm(
          A,
          g = out_comm,
          K = K,
          validate_inputs = FALSE
        )
        P_hat <- B_hat[out_comm, out_comm]
      } else {
        parameter_hat <- estimate_dcbm(
          A,
          g = out_comm,
          K = K,
          method = dcbm_est_method,
          validate_inputs = FALSE
        )
        P_hat <- parameter_hat$B_hat[out_comm, out_comm] *
          tcrossprod(parameter_hat$psi_hat)
      }
      P_hat <- pmin(P_hat, 1)
      P_hat <- pmax(P_hat, 0)
      deg_hat <- rowSums(P_hat)
      P_hat_tau <- P_hat + tau * mean(deg) / n
      d_tau_hat <- Matrix::sparseMatrix(
        i = seq_len(n),
        j = seq_len(n),
        x = 1 / sqrt(deg_hat + tau * mean(deg) / n)
      )
      L_hat_tau <- P_hat_tau + 1 - 1
      if (use_laplacian) {
        L_hat_tau <- tcrossprod(
          crossprod(d_tau_hat, P_hat_tau), d_tau_hat
        )
        L_hat_tau[is.na(L_hat_tau)] <- 0
      }
      numerator <- RSpectra::eigs_sym(
        L_tau - L_hat_tau,
        K,
        which = "LM",
        opts = list(v0 = rep(1, n))
      )$values[1L]
      denominator <- RSpectra::eigs_sym(
        L_hat_tau,
        K,
        which = "LM",
        opts = list(v0 = rep(1, n))
      )$values[K]
      if (!is.finite(numerator) || !is.finite(denominator) ||
          abs(denominator) <= sqrt(.Machine$double.eps)) {
        stop(
          "The DK ratio has a non-finite numerator or a zero denominator.",
          call. = FALSE
        )
      }
      dk_stat <- abs(numerator / denominator)
      if (!is.finite(dk_stat)) {
        stop("The DK statistic is non-finite.", call. = FALSE)
      }
      list(
        success = TRUE,
        candidate_id = candidate_id,
        seed = candidate_seeds[[candidate_id]],
        tau = tau,
        dk_stat = dk_stat,
        communities = out_comm,
        numerator = as.numeric(numerator),
        denominator = as.numeric(denominator),
        message = NA_character_
      )
    }, error = function(e) {
      list(
        success = FALSE,
        candidate_id = candidate_id,
        seed = candidate_seeds[[candidate_id]],
        tau = tau_candidates[[candidate_id]],
        dk_stat = NA_real_,
        communities = NULL,
        numerator = NA_real_,
        denominator = NA_real_,
        message = conditionMessage(e)
      )
    })
  }
  if (verbose) {
    message(
      "Evaluating ", length(tau_candidates), " DK regularizer candidate(s) ",
      "with ", candidate_ncores, " worker(s)."
    )
  }
  elapsed <- system.time({
    candidate_output <- uni_mclapply(
      seq_along(tau_candidates),
      evaluate_candidate,
      ncores = candidate_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      future_packages = worker_packages
    )
  })
  success_mask <- vapply(candidate_output, `[[`, logical(1), "success")
  failures <- candidate_output[!success_mask]
  diagnostics <- data.frame(
    tau = vapply(failures, `[[`, numeric(1), "tau"),
    seed = vapply(failures, `[[`, integer(1), "seed"),
    message = vapply(failures, `[[`, character(1), "message"),
    stringsAsFactors = FALSE
  )
  if (length(failures) > 0L && failure_handling == "stop") {
    stop(
      "DKEST failed for tau candidate(s) ",
      paste(diagnostics$tau, collapse = ", "), ": ",
      paste(diagnostics$message, collapse = " | "),
      call. = FALSE
    )
  }
  if (!any(success_mask)) {
    stop("No valid DKEST candidates remain.", call. = FALSE)
  }
  if (length(failures) > 0L) {
    warning(
      "Omitting failed DKEST tau candidate(s): ",
      paste(diagnostics$tau, collapse = ", "), ".",
      call. = FALSE
    )
  }
  dk_values <- vapply(candidate_output, `[[`, numeric(1), "dk_stat")
  numerator_values <- vapply(candidate_output, `[[`, numeric(1), "numerator")
  denominator_values <- vapply(
    candidate_output, `[[`, numeric(1), "denominator"
  )
  dk_statistic <- data.frame(
    tau = tau_candidates,
    dk_stat = dk_values,
    numerator = numerator_values,
    denominator = denominator_values,
    valid = success_mask,
    stringsAsFactors = FALSE
  )
  valid_ids <- which(success_mask)
  minimum_statistic <- min(dk_values[valid_ids])
  tied_ids <- valid_ids[dk_values[valid_ids] == minimum_statistic]
  best_id <- tied_ids[which.min(tau_candidates[tied_ids])]
  tau_hat <- tau_candidates[[best_id]]
  selected_dk_stat <- dk_values[[best_id]]
  overall_best <- data.frame(
    criterion = "dk_stat",
    tau_hat = tau_hat,
    statistic = selected_dk_stat,
    stringsAsFactors = FALSE
  )
  selection <- list(
    tau_hat = tau_hat,
    dk_stat = selected_dk_stat
  )
  out <- list(
    algorithm = "DKEST",
    dk_statistic = dk_statistic,
    tau_hat = tau_hat,
    selected_dk_stat = selected_dk_stat,
    selection = selection,
    overall_best = overall_best,
    tau_candidates = tau_candidates,
    K = K,
    model = if (use_dcbm) "DCBM" else "SBM",
    diagnostics = diagnostics,
    candidate_seeds = candidate_seeds,
    raw_output = if (retain_intermediates == "all") {
      candidate_output
    } else {
      NULL
    },
    retain_intermediates = retain_intermediates,
    options = list(
      use_laplacian = use_laplacian,
      dcbm_est_method = dcbm_est_method
    ),
    ncores = list(requested = ncores, tuning = candidate_ncores),
    timing = c(total = unname(elapsed["elapsed"])),
    ram_preflight = ram_report,
    seed = seed,
    call = call
  )
  class(out) <- "netcrop_regularizer"
  out
}

#' Fit spectral clusterings over multiple regularizers
#'
#' Fits spectral clustering across a requested regularization grid for one
#' network or a matched list of networks.
#'
#' @details
#' Network-by-candidate tasks are parallelized. The same realized seed is used
#' across candidates within a network so stochastic clustering comparisons are
#' fair. Results can retain complete fitted objects or labels only.
#'
#' @param A One adjacency matrix or a nonempty list of equal-sized matrices.
#' @param K Fixed positive number of communities.
#' @param tau_candidates Unique finite nonnegative regularization candidates.
#' @param laplacian Convert adjacency matrices to graph Laplacians.
#' @param normalize_laplacian Use symmetric normalized Laplacians.
#' @param handle_zero_degree_nodes Zero-degree node policy.
#' @param row_normalize Normalize embedding rows before clustering.
#' @param spectral_method Use an eigen- or singular-vector representation.
#' @param spectral_engine Decomposition backend.
#' @param spectral_options Named decomposition options.
#' @param cluster_engine Clustering backend.
#' @param cluster_options Named clustering options.
#' @param ncores Positive task-worker count.
#' @param seed Optional nonnegative reproducibility seed.
#' @param verbose Print progress messages.
#' @param force_windows Use the Windows-compatible parallel backend.
#' @param ram_check Report conservative RAM demand.
#' @param failure_handling Stop or omit failed candidate fits.
#' @param retain_fits Retain fitted `spectral_cluster` objects.
#' @return A `mult_reg_clustering` object containing labels and optional fits by
#'   network and candidate, resolved parameters, task metadata, failures,
#'   seeds, workers, timing, and RAM diagnostics.
#' @examples
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.5, beta = 0.1,
#'                   seed = 37, ncores = 1)
#' fits <- mult_reg_spectral_cluster(
#'   A, K = 3, tau_candidates = c(0, 0.1), ncores = 1,
#'   seed = 38, verbose = FALSE, spectral_engine = "base",
#'   cluster_engine = "kmeans"
#' )
#' @seealso [dkest_tune_regularizer()], [netcrop_tune_regularizer()]
# Fit spectral clustering over multiple regularization values and networks.
#' @rdname mult_reg_spectral_cluster
#' @export
mult_reg_spectral_cluster <- function(
    A,
    K,
    tau_candidates,
    laplacian = FALSE,
    normalize_laplacian = TRUE,
    handle_zero_degree_nodes = c("none", "random_label", "remove"),
    row_normalize = FALSE,
    spectral_method = c("eigen", "svd"),
    spectral_engine = c("RSpectra", "irlba", "base"),
    spectral_options = list(),
    cluster_engine = c("clara", "kmeans", "pam"),
    cluster_options = list(),
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    seed = NULL,
    verbose = TRUE,
    force_windows = FALSE,
    ram_check = FALSE,
    failure_handling = c("stop", "omit"),
    retain_fits = TRUE) {
  call <- match.call()
  handle_zero_degree_nodes <- match.arg(handle_zero_degree_nodes)
  spectral_method <- match.arg(spectral_method)
  spectral_engine <- match.arg(spectral_engine)
  cluster_engine <- match.arg(cluster_engine)
  failure_handling <- match.arg(failure_handling)
  validate_count <- function(value, name) {
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) || value < 1 || value != floor(value)) {
      stop(name, " must be one positive integer.", call. = FALSE)
    }
    as.integer(value)
  }
  validate_flag <- function(value, name) {
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
    value
  }
  required_helpers <- c(
    "uni_mclapply", "spectral_cluster", "is_symmetric_matrix",
    "offset_generator_seed", "warn_if_insufficient_ram", "format_bytes"
  )
  missing_helpers <- required_helpers[!vapply(
    required_helpers, exists, logical(1), mode = "function", inherits = TRUE, envir = environment()
  )]
  if (length(missing_helpers) > 0L) {
    stop(
      "Required internal netOP helpers are unavailable; reinstall netOP. Missing: ",
      paste(missing_helpers, collapse = ", "), ".",
      call. = FALSE
    )
  }
  A_is_single <- !is.list(A) || inherits(A, "Matrix")
  A_list <- if (A_is_single) list(A) else A
  if (length(A_list) < 1L) {
    stop("A must be a matrix-like object or a non-empty list of them.",
         call. = FALSE)
  }
  validate_network <- function(network, network_id) {
    if (is.null(dim(network)) || length(dim(network)) != 2L ||
        nrow(network) != ncol(network) || nrow(network) < 2L ||
        !(is.numeric(network) || inherits(network, "Matrix"))) {
      stop("A[[", network_id, "]] must be a numeric square matrix.",
           call. = FALSE)
    }
    values <- if (inherits(network, "sparseMatrix") &&
                  "x" %in% methods::slotNames(network)) {
      methods::slot(network, "x")
    } else {
      as.numeric(network)
    }
    if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
      stop("A[[", network_id,
           "]] must contain finite non-negative values.", call. = FALSE)
    }
    if (!is_symmetric_matrix(network)) {
      stop("A[[", network_id, "]] must be symmetric.", call. = FALSE)
    }
    diagonal <- if (inherits(network, "Matrix")) {
      Matrix::diag(network)
    } else {
      diag(network)
    }
    if (any(diagonal != 0)) {
      stop("A[[", network_id, "]] must have a zero diagonal.",
           call. = FALSE)
    }
    invisible(TRUE)
  }
  invisible(Map(validate_network, A_list, seq_along(A_list)))
  dimensions <- vapply(A_list, nrow, integer(1))
  if (length(unique(dimensions)) != 1L) {
    stop("All matrices in A must have the same dimensions.", call. = FALSE)
  }
  K <- validate_count(K, "K")
  ncores <- validate_count(ncores, "ncores")
  if (K >= dimensions[[1L]]) {
    stop("K must be smaller than the shared network dimension.",
         call. = FALSE)
  }
  if (!is.numeric(tau_candidates) || length(tau_candidates) < 1L ||
      anyNA(tau_candidates) || any(!is.finite(tau_candidates)) ||
      any(tau_candidates < 0) || anyDuplicated(tau_candidates)) {
    stop("tau_candidates must contain unique finite non-negative numbers.",
         call. = FALSE)
  }
  tau_candidates <- as.numeric(tau_candidates)
  option_lists <- list(
    spectral_options = spectral_options,
    cluster_options = cluster_options
  )
  invalid_options <- vapply(option_lists, function(options) {
    !is.list(options) || (length(options) > 0L &&
                            (is.null(names(options)) || any(!nzchar(names(options)))))
  }, logical(1))
  if (any(invalid_options)) {
    stop(
      paste(names(option_lists)[invalid_options], collapse = ", "),
      " must be named lists.", call. = FALSE
    )
  }
  logical_inputs <- list(
    laplacian = laplacian,
    normalize_laplacian = normalize_laplacian,
    row_normalize = row_normalize,
    verbose = verbose,
    force_windows = force_windows,
    ram_check = ram_check,
    retain_fits = retain_fits
  )
  for (name in names(logical_inputs)) {
    logical_inputs[[name]] <- validate_flag(logical_inputs[[name]], name)
  }
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
       !is.finite(seed) || seed < 0 || seed != floor(seed) ||
       seed > .Machine$integer.max)) {
    stop("seed must be NULL or a valid non-negative integer.",
         call. = FALSE)
  }
  if (!is.null(seed)) seed <- as.integer(seed)
  network_seeds <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, length(A_list), replace = TRUE)
  } else {
    vapply(
      seq_along(A_list),
      function(network_id) offset_generator_seed(seed, 10000 * network_id),
      integer(1)
    )
  }
  task_grid <- expand.grid(
    network_id = seq_along(A_list),
    tau_id = seq_along(tau_candidates),
    KEEP.OUT.ATTRS = FALSE
  )
  task_ncores <- min(ncores, nrow(task_grid))
  dense_matrix_bytes <- 8 * as.double(dimensions[[1L]])^2
  per_worker_bytes <- 5 * dense_matrix_bytes +
    8 * as.double(dimensions[[1L]]) * K
  estimated_bytes <- per_worker_bytes * task_ncores
  ram_status <- NULL
  if (ram_check) {
    ram_status <- warn_if_insufficient_ram(
      estimated_bytes,
      operation = paste0(
        "multi-tau spectral clustering preflight: ",
        format_bytes(per_worker_bytes), " per worker x ", task_ncores
      )
    )
  }
  run_task <- function(task_id) {
    network_id <- task_grid$network_id[[task_id]]
    tau_id <- task_grid$tau_id[[task_id]]
    task_seed <- network_seeds[[network_id]]
    tryCatch({
      set.seed(task_seed)
      fit <- spectral_cluster(
        A = A_list[[network_id]],
        K = K,
        laplacian = laplacian,
        normalize_laplacian = normalize_laplacian,
        regularize_tau = tau_candidates[[tau_id]],
        handle_zero_degree_nodes = handle_zero_degree_nodes,
        row_normalize = row_normalize,
        spectral_method = spectral_method,
        spectral_engine = spectral_engine,
        spectral_options = spectral_options,
        cluster_engine = cluster_engine,
        cluster_options = cluster_options,
        ram_check = FALSE,
        validate_inputs = TRUE
      )
      labels <- fit$g_hat
      if (!is.numeric(labels) || length(labels) != dimensions[[network_id]] ||
          anyNA(labels) || any(!is.finite(labels)) || any(labels < 1L) ||
          any(labels > K) || any(labels != floor(labels))) {
        stop("spectral_cluster() returned invalid labels.", call. = FALSE)
      }
      list(
        success = TRUE,
        network_id = network_id,
        tau_id = tau_id,
        tau = tau_candidates[[tau_id]],
        seed = task_seed,
        labels = as.integer(labels),
        fit = if (retain_fits) fit else NULL,
        message = NA_character_
      )
    }, error = function(e) {
      list(
        success = FALSE,
        network_id = network_id,
        tau_id = tau_id,
        tau = tau_candidates[[tau_id]],
        seed = task_seed,
        labels = NULL,
        fit = NULL,
        message = conditionMessage(e)
      )
    })
  }
  if (verbose) {
    message(
      "Running ", nrow(task_grid), " network/tau spectral fit(s) with ",
      task_ncores, " worker(s)."
    )
  }
  elapsed <- system.time({
    task_output <- uni_mclapply(
      seq_len(nrow(task_grid)),
      run_task,
      ncores = task_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      future_packages = if (any(vapply(
        A_list, inherits, logical(1), "Matrix"
      ))) "Matrix" else NULL
    )
  })
  success_mask <- vapply(task_output, `[[`, logical(1), "success")
  failures <- task_output[!success_mask]
  diagnostics <- data.frame(
    network_id = vapply(failures, `[[`, integer(1), "network_id"),
    tau = vapply(failures, `[[`, numeric(1), "tau"),
    seed = vapply(failures, `[[`, integer(1), "seed"),
    message = vapply(failures, `[[`, character(1), "message"),
    stringsAsFactors = FALSE
  )
  if (length(failures) > 0L && failure_handling == "stop") {
    stop(
      "Spectral clustering failed for ", length(failures), " task(s): ",
      paste(diagnostics$message, collapse = " | "), call. = FALSE
    )
  }
  if (length(failures) > 0L) {
    warning(length(failures), " spectral task(s) were omitted.",
            call. = FALSE)
  }
  labels <- lapply(seq_along(A_list), function(network_id) {
    lapply(seq_along(tau_candidates), function(tau_id) {
      task_id <- which(
        task_grid$network_id == network_id & task_grid$tau_id == tau_id
      )
      task_output[[task_id]]$labels
    })
  })
  names(labels) <- paste0("network_", seq_along(A_list))
  for (network_id in seq_along(labels)) {
    names(labels[[network_id]]) <- as.character(tau_candidates)
  }
  fits <- if (retain_fits) {
    lapply(seq_along(A_list), function(network_id) {
      result <- lapply(seq_along(tau_candidates), function(tau_id) {
        task_id <- which(
          task_grid$network_id == network_id & task_grid$tau_id == tau_id
        )
        task_output[[task_id]]$fit
      })
      names(result) <- as.character(tau_candidates)
      result
    })
  } else {
    NULL
  }
  result <- list(
    algorithm = "spectral_cluster",
    labels = labels,
    fits = fits,
    tau_candidates = tau_candidates,
    K = K,
    network_count = length(A_list),
    shared_dimension = dimensions[[1L]],
    parameters = list(
      laplacian = laplacian,
      normalize_laplacian = normalize_laplacian,
      handle_zero_degree_nodes = handle_zero_degree_nodes,
      row_normalize = row_normalize,
      spectral_method = spectral_method,
      spectral_engine = spectral_engine,
      spectral_options = spectral_options,
      cluster_engine = cluster_engine,
      cluster_options = cluster_options
    ),
    task_grid = task_grid,
    diagnostics = diagnostics,
    failure_handling = failure_handling,
    retain_fits = retain_fits,
    ncores = list(requested = ncores, tasks = task_ncores),
    seed = seed,
    network_seeds = network_seeds,
    timing = c(total = unname(elapsed["elapsed"])),
    ram_preflight = list(
      estimated_bytes = estimated_bytes,
      per_worker_bytes = per_worker_bytes,
      workers = task_ncores,
      status = ram_status
    ),
    call = call
  )
  class(result) <- "mult_reg_clustering"
  result
}

#' Fit SONNET over multiple regularizers
#'
#' Fits SONNET across a requested regularization grid for one network or a
#' matched list.
#'
#' @details
#' Candidate fits run sequentially so `ncores` remains available to SONNET's
#' internal subnetwork computations. The same realized seed is reused across
#' regularization candidates within each network.
#'
#' @param A One adjacency matrix or a nonempty list of equal-sized matrices.
#' @param K Fixed positive number of communities.
#' @param tau_candidates Unique finite nonnegative regularization candidates.
#' @param num_subnetworks Number of overlapping subnetworks; `NULL` selects it
#'   automatically.
#' @param overlap_size Number of overlap nodes; `NULL` selects it automatically.
#' @param extra_nrep Number of additional SONNET repetitions.
#' @param ncores Positive SONNET worker count.
#' @param seed Optional nonnegative reproducibility seed.
#' @param matching_method Label-alignment method.
#' @param confirm_large Confirm potentially factorial brute-force matching.
#' @param verbose Print progress messages.
#' @param force_windows Use the Windows-compatible parallel backend.
#' @param ram_check Report conservative RAM demand.
#' @param share_overlap Reuse one overlap across repetitions.
#' @param parameter_select_options Named automatic-partition options.
#' @param failure_handling Stop or omit failed candidate fits.
#' @param retain_fits Retain fitted `sonnet` objects.
#' @param ... Named spectral-clustering arguments forwarded to [sonnet()].
#' @return A `mult_reg_clustering` object containing SONNET labels and optional
#'   fits by network and regularizer, parameters, tasks, failures, seeds, and
#'   timing.
#' @examples
#' \donttest{
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.5, beta = 0.1,
#'                   seed = 39, ncores = 1)
#' mult_reg_sonnet(
#'   A, K = 3, tau_candidates = c(0, 0.1),
#'   num_subnetworks = 2, overlap_size = 20, ncores = 1,
#'   seed = 40, verbose = FALSE, spectral_engine = "base"
#' )
#' }
#' @seealso [sonnet()], [mult_reg_spectral_cluster()]
# Fit SONNET over multiple regularization values and networks.
#' @rdname mult_reg_sonnet
#' @export
mult_reg_sonnet <- function(
    A,
    K,
    tau_candidates,
    num_subnetworks = NULL,
    overlap_size = NULL,
    extra_nrep = 0L,
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    seed = NULL,
    matching_method = c("greedy", "hungarian", "brute_force"),
    confirm_large = TRUE,
    verbose = TRUE,
    force_windows = FALSE,
    ram_check = FALSE,
    share_overlap = FALSE,
    parameter_select_options = list(),
    failure_handling = c("stop", "omit"),
    retain_fits = TRUE,
    ...) {
  call <- match.call()
  matching_method <- match.arg(matching_method)
  failure_handling <- match.arg(failure_handling)
  required_helpers <- c("sonnet", "offset_generator_seed")
  missing_helpers <- required_helpers[!vapply(
    required_helpers, exists, logical(1), mode = "function", inherits = TRUE, envir = environment()
  )]
  if (length(missing_helpers) > 0L) {
    stop(
      "Required internal netOP helpers are unavailable; reinstall netOP. Missing: ",
      paste(missing_helpers, collapse = ", "), ".", call. = FALSE
    )
  }
  A_is_single <- !is.list(A) || inherits(A, "Matrix")
  A_list <- if (A_is_single) list(A) else A
  if (length(A_list) < 1L) {
    stop("A must be a matrix-like object or a non-empty list of them.",
         call. = FALSE)
  }
  dimensions <- vapply(A_list, function(network) {
    if (is.null(dim(network)) || length(dim(network)) != 2L ||
        nrow(network) != ncol(network) ||
        !(is.numeric(network) || inherits(network, "Matrix")) ||
        any(!is.finite(network))) {
      stop("Every A entry must be a finite numeric square matrix.",
           call. = FALSE)
    }
    nrow(network)
  }, integer(1))
  if (length(unique(dimensions)) != 1L) {
    stop("All matrices in A must have the same dimensions.", call. = FALSE)
  }
  validate_count <- function(value, name, minimum = 0L) {
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) || value < minimum || value != floor(value)) {
      stop(name, " must be one integer at least ", minimum, ".",
           call. = FALSE)
    }
    as.integer(value)
  }
  K <- validate_count(K, "K", 1L)
  extra_nrep <- validate_count(extra_nrep, "extra_nrep", 0L)
  ncores <- validate_count(ncores, "ncores", 1L)
  if (K > dimensions[[1L]]) {
    stop("K cannot exceed the shared network dimension.", call. = FALSE)
  }
  if (!is.numeric(tau_candidates) || length(tau_candidates) < 1L ||
      anyNA(tau_candidates) || any(!is.finite(tau_candidates)) ||
      any(tau_candidates < 0) || anyDuplicated(tau_candidates)) {
    stop("tau_candidates must contain unique finite non-negative numbers.",
         call. = FALSE)
  }
  tau_candidates <- as.numeric(tau_candidates)
  dots <- list(...)
  if (length(dots) > 0L &&
      (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    stop("All spectral-clustering arguments in ... must be named.",
         call. = FALSE)
  }
  protected_dots <- intersect(
    names(dots), c("A", "U", "K", "regularize_tau", "ram_check")
  )
  if (length(protected_dots) > 0L) {
    stop(
      "... cannot override ", paste(protected_dots, collapse = ", "), ".",
      call. = FALSE
    )
  }
  logical_inputs <- list(
    confirm_large = confirm_large,
    verbose = verbose,
    force_windows = force_windows,
    ram_check = ram_check,
    share_overlap = share_overlap,
    retain_fits = retain_fits
  )
  invalid_logical <- vapply(logical_inputs, function(value) {
    length(value) != 1L || !is.logical(value) || is.na(value)
  }, logical(1))
  if (any(invalid_logical)) {
    stop(
      paste(names(logical_inputs)[invalid_logical], collapse = ", "),
      " must be TRUE or FALSE.", call. = FALSE
    )
  }
  if (!is.list(parameter_select_options) ||
      (length(parameter_select_options) > 0L &&
       (is.null(names(parameter_select_options)) ||
        any(!nzchar(names(parameter_select_options)))))) {
    stop("parameter_select_options must be a named list.", call. = FALSE)
  }
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
       !is.finite(seed) || seed < 0 || seed != floor(seed) ||
       seed > .Machine$integer.max)) {
    stop("seed must be NULL or a valid non-negative integer.",
         call. = FALSE)
  }
  if (!is.null(seed)) seed <- as.integer(seed)
  network_seeds <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, length(A_list), replace = TRUE)
  } else {
    vapply(
      seq_along(A_list),
      function(network_id) offset_generator_seed(seed, 10000 * network_id),
      integer(1)
    )
  }
  task_grid <- expand.grid(
    network_id = seq_along(A_list),
    tau_id = seq_along(tau_candidates),
    KEEP.OUT.ATTRS = FALSE
  )
  task_output <- vector("list", nrow(task_grid))
  if (verbose) {
    message(
      "Running ", nrow(task_grid),
      " SONNET fit(s) sequentially with parallelism inside SONNET."
    )
  }
  elapsed <- system.time({
    for (task_id in seq_len(nrow(task_grid))) {
      network_id <- task_grid$network_id[[task_id]]
      tau_id <- task_grid$tau_id[[task_id]]
      task_seed <- network_seeds[[network_id]]
      task_output[[task_id]] <- tryCatch({
        arguments <- c(
          list(
            A = A_list[[network_id]],
            K = K,
            num_subnetworks = num_subnetworks,
            overlap_size = overlap_size,
            extra_nrep = extra_nrep,
            ncores = ncores,
            seed = task_seed,
            matching_method = matching_method,
            confirm_large = confirm_large,
            verbose = verbose,
            force_windows = force_windows,
            ram_check = ram_check,
            share_overlap = share_overlap,
            parameter_select_options = parameter_select_options,
            regularize_tau = tau_candidates[[tau_id]]
          ),
          dots
        )
        fit <- do.call(sonnet, arguments)
        labels <- fit$labels
        if (!is.numeric(labels) || length(labels) != dimensions[[network_id]] ||
            anyNA(labels) || any(!is.finite(labels)) || any(labels < 1L) ||
            any(labels > K) || any(labels != floor(labels))) {
          stop("sonnet() returned invalid labels.", call. = FALSE)
        }
        list(
          success = TRUE,
          network_id = network_id,
          tau_id = tau_id,
          tau = tau_candidates[[tau_id]],
          seed = task_seed,
          labels = as.integer(labels),
          fit = if (retain_fits) fit else NULL,
          message = NA_character_
        )
      }, error = function(e) {
        list(
          success = FALSE,
          network_id = network_id,
          tau_id = tau_id,
          tau = tau_candidates[[tau_id]],
          seed = task_seed,
          labels = NULL,
          fit = NULL,
          message = conditionMessage(e)
        )
      })
    }
  })
  success_mask <- vapply(task_output, `[[`, logical(1), "success")
  failures <- task_output[!success_mask]
  diagnostics <- data.frame(
    network_id = vapply(failures, `[[`, integer(1), "network_id"),
    tau = vapply(failures, `[[`, numeric(1), "tau"),
    seed = vapply(failures, `[[`, integer(1), "seed"),
    message = vapply(failures, `[[`, character(1), "message"),
    stringsAsFactors = FALSE
  )
  if (length(failures) > 0L && failure_handling == "stop") {
    stop(
      "SONNET failed for ", length(failures), " task(s): ",
      paste(diagnostics$message, collapse = " | "), call. = FALSE
    )
  }
  if (length(failures) > 0L) {
    warning(length(failures), " SONNET task(s) were omitted.",
            call. = FALSE)
  }
  labels <- lapply(seq_along(A_list), function(network_id) {
    result <- lapply(seq_along(tau_candidates), function(tau_id) {
      task_id <- which(
        task_grid$network_id == network_id & task_grid$tau_id == tau_id
      )
      task_output[[task_id]]$labels
    })
    names(result) <- as.character(tau_candidates)
    result
  })
  names(labels) <- paste0("network_", seq_along(A_list))
  fits <- if (retain_fits) {
    lapply(seq_along(A_list), function(network_id) {
      result <- lapply(seq_along(tau_candidates), function(tau_id) {
        task_id <- which(
          task_grid$network_id == network_id & task_grid$tau_id == tau_id
        )
        task_output[[task_id]]$fit
      })
      names(result) <- as.character(tau_candidates)
      result
    })
  } else {
    NULL
  }
  result <- list(
    algorithm = "sonnet",
    labels = labels,
    fits = fits,
    tau_candidates = tau_candidates,
    K = K,
    network_count = length(A_list),
    shared_dimension = dimensions[[1L]],
    parameters = list(
      num_subnetworks = num_subnetworks,
      overlap_size = overlap_size,
      extra_nrep = extra_nrep,
      ncores = ncores,
      matching_method = matching_method,
      confirm_large = confirm_large,
      force_windows = force_windows,
      ram_check = ram_check,
      share_overlap = share_overlap,
      parameter_select_options = parameter_select_options,
      spectral_arguments = dots
    ),
    task_grid = task_grid,
    diagnostics = diagnostics,
    failure_handling = failure_handling,
    retain_fits = retain_fits,
    seed = seed,
    network_seeds = network_seeds,
    timing = c(total = unname(elapsed["elapsed"])),
    call = call
  )
  class(result) <- "mult_reg_clustering"
  result
}

#' Compare clustering methods with known labels
#'
#' Evaluates SONNET, spectral clustering, or both against known community
#' labels over candidate regularizers.
#'
#' @details
#' Optional NETCROP and DKEST results add their selected regularizers to the
#' comparison, including off-grid values. Accuracy is one minus the validated
#' label-matching mismatch rate; multiple networks are summarized with means
#' and standard deviations.
#'
#' @param A One adjacency matrix or a list matching `g_true`.
#' @param g_true Ground-truth labels or a list matching `A`.
#' @param tau_candidates Unique finite nonnegative regularization candidates.
#' @param K Optional fixed number of communities; inferred from `g_true` when
#'   omitted.
#' @param netcrop_outcomes,dkest_outcomes Optional matching tuner results.
#' @param include_netcrop_mean,include_netcrop_mode Include NETCROP repetition
#'   summaries in addition to its first repetition.
#' @param losses Optional NETCROP loss names to include.
#' @param engines Fit `"sonnet"`, `"spectral_cluster"`, or both.
#' @param matching_method Label-alignment method used for accuracy.
#' @param confirm_large Confirm potentially factorial brute-force matching.
#' @param sonnet_options,spectral_cluster_options Named fitter option lists.
#' @param ncores Positive worker count.
#' @param seed Optional nonnegative reproducibility seed.
#' @param verbose Print progress messages.
#' @param force_windows Use the Windows-compatible parallel backend.
#' @param ram_check Report conservative RAM demand.
#' @return Invisibly returns the rendered `ggplot`. Raw accuracy data,
#'   aggregated plotting data, effective fit metadata, and diagnostics are
#'   attached as attributes.
#' @examples
#' \donttest{
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.5, beta = 0.1,
#'                   seed = 41, ncores = 1)
#' truth <- get_generator_parameters(A)$g_true
#' oracle_plotter(
#'   A, g_true = truth, tau_candidates = c(0, 0.1), K = 3,
#'   engines = "spectral_cluster", ncores = 1, seed = 42,
#'   verbose = FALSE,
#'   spectral_cluster_options = list(spectral_engine = "base")
#' )
#' }
#' @seealso [spectral_cluster()], [sonnet()]
# Plot oracle clustering accuracy and optional regularizer selections.
#' @rdname oracle_plotter
#' @export
oracle_plotter <- function(
    A,
    g_true,
    tau_candidates,
    K = NULL,
    netcrop_outcomes = NULL,
    dkest_outcomes = NULL,
    include_netcrop_mean = TRUE,
    include_netcrop_mode = TRUE,
    losses = NULL,
    engines = c("sonnet", "spectral_cluster"),
    matching_method = c("greedy", "hungarian", "brute_force"),
    confirm_large = NULL,
    sonnet_options = list(),
    spectral_cluster_options = list(),
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    seed = NULL,
    verbose = TRUE,
    force_windows = FALSE,
    ram_check = FALSE) {
  matching_method <- match.arg(matching_method)
  engines <- match.arg(
    engines,
    choices = c("sonnet", "spectral_cluster"),
    several.ok = TRUE
  )
  engines <- unique(engines)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for oracle plotting.",
         call. = FALSE)
  }
  required_helpers <- c(
    "mult_reg_sonnet", "mult_reg_spectral_cluster",
    "label_match_greedy", "label_match_brute_force", "modal",
    "offset_generator_seed"
  )
  missing_helpers <- required_helpers[!vapply(
    required_helpers, exists, logical(1), mode = "function", inherits = TRUE, envir = environment()
  )]
  if (length(missing_helpers) > 0L) {
    stop(
      "Missing oracle helper(s): ", paste(missing_helpers, collapse = ", "),
      ".", call. = FALSE
    )
  }
  A_is_single <- !is.list(A) || inherits(A, "Matrix")
  A_list <- if (A_is_single) list(A) else A
  truth_is_single <- !is.list(g_true)
  truth_list <- if (truth_is_single) list(g_true) else g_true
  if (A_is_single && truth_is_single) {
    if (!is.null(netcrop_outcomes) &&
        inherits(netcrop_outcomes, "netcrop_regularizer")) {
      netcrop_outcomes <- list(netcrop_outcomes)
    }
    if (!is.null(dkest_outcomes) &&
        inherits(dkest_outcomes, "netcrop_regularizer")) {
      dkest_outcomes <- list(dkest_outcomes)
    }
  }
  if (length(A_list) < 1L || length(A_list) != length(truth_list)) {
    stop("A and g_true must contain the same positive number of entries.",
         call. = FALSE)
  }
  dimensions <- vapply(A_list, function(network) {
    if (is.null(dim(network)) || length(dim(network)) != 2L ||
        nrow(network) != ncol(network) ||
        !(is.numeric(network) || inherits(network, "Matrix")) ||
        any(!is.finite(network))) {
      stop("Every A entry must be a finite numeric square matrix.",
           call. = FALSE)
    }
    nrow(network)
  }, integer(1))
  if (length(unique(dimensions)) != 1L) {
    stop("All A entries must share the same dimensions.", call. = FALSE)
  }
  truth_valid <- vapply(seq_along(truth_list), function(network_id) {
    truth <- truth_list[[network_id]]
    is.numeric(truth) && length(truth) == dimensions[[network_id]] &&
      !anyNA(truth) && all(is.finite(truth)) && all(truth >= 1) &&
      all(truth == floor(truth))
  }, logical(1))
  if (!all(truth_valid)) {
    stop(
      "Every g_true entry must be a positive-integer vector matching A.",
      call. = FALSE
    )
  }
  inferred_K <- max(unlist(truth_list, use.names = FALSE))
  if (is.null(K)) K <- inferred_K
  if (length(K) != 1L || !is.numeric(K) || is.na(K) ||
      !is.finite(K) || K < inferred_K || K != floor(K)) {
    stop("K must be one positive integer covering every g_true label.",
         call. = FALSE)
  }
  K <- as.integer(K)
  if (K >= dimensions[[1L]]) {
    stop("K must be smaller than the shared network dimension.",
         call. = FALSE)
  }
  if (!is.numeric(tau_candidates) || length(tau_candidates) < 1L ||
      anyNA(tau_candidates) || any(!is.finite(tau_candidates)) ||
      any(tau_candidates < 0) || anyDuplicated(tau_candidates)) {
    stop("tau_candidates must contain unique finite non-negative numbers.",
         call. = FALSE)
  }
  tau_candidates <- as.numeric(tau_candidates)
  logical_inputs <- list(
    include_netcrop_mean = include_netcrop_mean,
    include_netcrop_mode = include_netcrop_mode,
    verbose = verbose,
    force_windows = force_windows,
    ram_check = ram_check
  )
  invalid_logical <- vapply(logical_inputs, function(value) {
    length(value) != 1L || !is.logical(value) || is.na(value)
  }, logical(1))
  if (any(invalid_logical)) {
    stop(
      paste(names(logical_inputs)[invalid_logical], collapse = ", "),
      " must be TRUE or FALSE.", call. = FALSE
    )
  }
  if (length(ncores) != 1L || !is.numeric(ncores) || is.na(ncores) ||
      !is.finite(ncores) || ncores < 1 || ncores != floor(ncores)) {
    stop("ncores must be one positive integer.", call. = FALSE)
  }
  ncores <- as.integer(ncores)
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
       !is.finite(seed) || seed < 0 || seed != floor(seed) ||
       seed > .Machine$integer.max)) {
    stop("seed must be NULL or a valid non-negative integer.",
         call. = FALSE)
  }
  if (!is.null(seed)) seed <- as.integer(seed)
  option_lists <- list(
    sonnet_options = sonnet_options,
    spectral_cluster_options = spectral_cluster_options
  )
  invalid_options <- vapply(option_lists, function(options) {
    !is.list(options) || (length(options) > 0L &&
                            (is.null(names(options)) || any(!nzchar(names(options)))))
  }, logical(1))
  if (any(invalid_options)) {
    stop(
      paste(names(option_lists)[invalid_options], collapse = ", "),
      " must be named lists.", call. = FALSE
    )
  }
  protected_options <- c("A", "K", "tau_candidates", "g_true")
  option_conflicts <- c(
    intersect(names(sonnet_options), protected_options),
    intersect(names(spectral_cluster_options), protected_options)
  )
  if (length(option_conflicts) > 0L) {
    stop(
      "Option lists cannot override ",
      paste(unique(option_conflicts), collapse = ", "), ".",
      call. = FALSE
    )
  }
  validate_outcomes <- function(outcomes, algorithm) {
    if (is.null(outcomes)) return(invisible(NULL))
    if (!is.list(outcomes) || length(outcomes) != length(A_list)) {
      stop(
        algorithm, " outcomes must be a list with one entry per network.",
        call. = FALSE
      )
    }
    for (network_id in seq_along(outcomes)) {
      outcome <- outcomes[[network_id]]
      if (!inherits(outcome, "netcrop_regularizer")) {
        stop(algorithm, " outcome ", network_id,
             " must inherit from netcrop_regularizer.", call. = FALSE)
      }
      observed_algorithm <- if (is.null(outcome$algorithm)) {
        "NETCROP"
      } else {
        toupper(outcome$algorithm)
      }
      if (observed_algorithm != algorithm) {
        stop("Expected ", algorithm, " at outcome ", network_id, ".",
             call. = FALSE)
      }
      if (is.null(outcome$K) || outcome$K != K) {
        stop(algorithm, " outcome ", network_id,
             " does not use the requested K.", call. = FALSE)
      }
    }
    invisible(NULL)
  }
  validate_outcomes(netcrop_outcomes, "NETCROP")
  validate_outcomes(dkest_outcomes, "DKEST")
  if (!is.null(losses) &&
      (!is.character(losses) || length(losses) < 1L || anyNA(losses) ||
       any(!nzchar(losses)))) {
    stop("losses must be NULL or non-empty loss names.", call. = FALSE)
  }
  if (matching_method == "brute_force" && K > 8L &&
      !identical(confirm_large, TRUE)) {
    stop(
      "Brute-force matching above K = 8 requires confirm_large = TRUE.",
      call. = FALSE
    )
  }
  if (!is.null(confirm_large) &&
      (length(confirm_large) != 1L || !is.logical(confirm_large) ||
       is.na(confirm_large))) {
    stop("confirm_large must be NULL, TRUE, or FALSE.", call. = FALSE)
  }
  match_accuracy <- function(labels, truth) {
    if (matching_method == "brute_force") {
      match <- label_match_brute_force(
        labels,
        truth,
        K = K,
        confirm_large = confirm_large
      )
    } else {
      match <- label_match_greedy(
        labels,
        truth,
        K = K,
        algorithm = matching_method
      )
    }
    accuracy <- 1 - match$mismatch_rate
    if (length(accuracy) != 1L || !is.finite(accuracy) ||
        accuracy < 0 || accuracy > 1) {
      stop("Label matching returned an invalid accuracy.", call. = FALSE)
    }
    accuracy
  }
  tau_equal <- function(x, y) {
    abs(x - y) <= sqrt(.Machine$double.eps) * max(1, abs(x), abs(y))
  }
  plot_seed <- function(network_id, engine_offset) {
    if (is.null(seed)) return(NULL)
    offset_generator_seed(seed, engine_offset + 10000 * network_id)
  }
  accuracy_records <- list()
  record_id <- 0L
  fit_metadata <- list()
  fit_metadata_id <- 0L
  diagnostic_records <- list()
  diagnostic_id <- 0L
  register_fit <- function(fit, network_id, engine, tau, parameter_source) {
    fit_metadata_id <<- fit_metadata_id + 1L
    fit_metadata[[fit_metadata_id]] <<- list(
      network_id = network_id,
      engine = engine,
      tau = tau,
      parameter_source = parameter_source,
      seed = fit$network_seeds[[1L]],
      parameters = fit$parameters,
      timing = fit$timing
    )
    if (!is.null(fit$diagnostics) && nrow(fit$diagnostics) > 0L) {
      diagnostics <- fit$diagnostics
      diagnostics$network_id <- network_id
      diagnostics$engine <- engine
      diagnostics$parameter_source <- parameter_source
      diagnostic_id <<- diagnostic_id + 1L
      diagnostic_records[[diagnostic_id]] <<- diagnostics
    }
    invisible(NULL)
  }
  add_record <- function(network_id, engine, method, tau, accuracy) {
    record_id <<- record_id + 1L
    accuracy_records[[record_id]] <<- data.frame(
      network_id = network_id,
      engine = engine,
      method = method,
      tau = tau,
      accuracy = accuracy,
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }
  extract_label <- function(fit, tau) {
    tau_id <- which(vapply(fit$tau_candidates, tau_equal, logical(1), tau))
    if (length(tau_id) != 1L || is.null(fit$labels[[1L]][[tau_id]])) {
      stop("No valid fitted labels are available for tau = ", tau, ".",
           call. = FALSE)
    }
    fit$labels[[1L]][[tau_id]]
  }
  prepare_netcrop_selections <- function(outcome) {
    selections <- outcome$best_tau_each_rep
    if (!is.data.frame(selections) ||
        !all(c("repetition", "loss", "tau_hat") %in% names(selections))) {
      stop("NETCROP outcome lacks valid repetition-level tau selections.",
           call. = FALSE)
    }
    selected_losses <- unique(selections$loss)
    if (!is.null(losses)) selected_losses <- intersect(selected_losses, losses)
    if (length(selected_losses) == 0L) {
      stop("No requested NETCROP losses are available.", call. = FALSE)
    }
    do.call(rbind, lapply(selected_losses, function(loss_name) {
      tau_values <- selections$tau_hat[selections$loss == loss_name]
      if (length(tau_values) < 1L || any(!is.finite(tau_values))) {
        stop("NETCROP has invalid tau selections for ", loss_name, ".",
             call. = FALSE)
      }
      rows <- list(data.frame(
        method = paste0("NETCROP (", loss_name, ") - one repetition"),
        tau = tau_values[[1L]],
        stringsAsFactors = FALSE
      ))
      if (length(tau_values) > 1L && include_netcrop_mean) {
        rows <- c(rows, list(data.frame(
          method = paste0("NETCROP (", loss_name, ") - mean"),
          tau = mean(tau_values),
          stringsAsFactors = FALSE
        )))
      }
      if (length(tau_values) > 1L && include_netcrop_mode) {
        rows <- c(rows, list(data.frame(
          method = paste0("NETCROP (", loss_name, ") - mode"),
          tau = modal(tau_values),
          stringsAsFactors = FALSE
        )))
      }
      do.call(rbind, rows)
    }))
  }
  run_sonnet_values <- function(network_id, tau_values, outcome = NULL) {
    arguments <- sonnet_options
    arguments[c(
      "A", "K", "tau_candidates", "ncores", "seed", "verbose",
      "force_windows", "ram_check", "failure_handling", "retain_fits"
    )] <- NULL
    outcome_algorithm <- NULL
    if (!is.null(outcome)) {
      outcome_algorithm <- if (is.null(outcome$algorithm)) {
        "NETCROP"
      } else {
        toupper(outcome$algorithm)
      }
      derived <- list(
        laplacian = outcome$options$use_laplacian,
        normalize_laplacian = TRUE,
        row_normalize = identical(outcome$model, "DCBM"),
        spectral_method = "eigen",
        spectral_engine = "RSpectra"
      )
      if (outcome_algorithm == "NETCROP") {
        derived <- c(derived, list(
          num_subnetworks = outcome$num_subnetworks,
          overlap_size = outcome$effective_overlap_size,
          extra_nrep = 0L,
          matching_method = outcome$options$matching_method,
          spectral_options = outcome$options$spectral_options,
          cluster_engine = if (identical(outcome$model, "DCBM")) {
            "clara"
          } else {
            "kmeans"
          },
          cluster_options = utils::modifyList(
            if (identical(outcome$model, "DCBM")) {
              list(metric = "manhattan", cluster.only = TRUE, samples = 5L)
            } else {
              list(nstart = 100L, iter.max = 10^7)
            },
            outcome$options$cluster_options,
            keep.null = TRUE
          )
        ))
      }
      derived <- derived[!vapply(derived, is.null, logical(1))]
      arguments <- utils::modifyList(derived, arguments, keep.null = TRUE)
    }
    fit_by_tau <- vector("list", length(tau_values))
    for (tau_id in seq_along(tau_values)) {
      tau <- tau_values[[tau_id]]
      task_arguments <- arguments
      if (!is.null(outcome)) {
        task_defaults <- list(
          handle_zero_degree_nodes = if (tau == 0) {
            "random_label"
          } else {
            "none"
          }
        )
        if (outcome_algorithm == "DKEST") {
          use_dcbm <- identical(outcome$model, "DCBM")
          task_defaults$cluster_engine <- if (use_dcbm) {
            "pam"
          } else if (tau == 0) {
            "kmeans"
          } else {
            "clara"
          }
          task_defaults$cluster_options <- if (use_dcbm) {
            list(
              metric = "euclidean", do.swap = FALSE,
              cluster.only = TRUE, pamonce = 6
            )
          } else if (tau == 0) {
            list(nstart = 100, iter.max = 10^7)
          } else {
            list(metric = "euclidean", cluster.only = TRUE)
          }
        }
        task_arguments <- utils::modifyList(
          task_defaults, task_arguments, keep.null = TRUE
        )
      }
      fit_by_tau[[tau_id]] <- do.call(
        mult_reg_sonnet,
        c(list(
          A = A_list[[network_id]],
          K = K,
          tau_candidates = tau,
          ncores = ncores,
          seed = if (!is.null(seed)) {
            plot_seed(network_id, 100000L)
          } else if (!is.null(outcome$seed)) {
            outcome$seed
          } else {
            NULL
          },
          verbose = verbose,
          force_windows = force_windows,
          ram_check = ram_check,
          failure_handling = "stop",
          retain_fits = FALSE
        ), task_arguments)
      )
      register_fit(
        fit_by_tau[[tau_id]],
        network_id = network_id,
        engine = "sonnet",
        tau = tau,
        parameter_source = if (is.null(outcome_algorithm)) {
          "current_defaults"
        } else {
          outcome_algorithm
        }
      )
    }
    names(fit_by_tau) <- as.character(tau_values)
    fit_by_tau
  }
  run_spectral_values <- function(network_id, tau_values, outcome = NULL) {
    arguments <- spectral_cluster_options
    arguments[c(
      "A", "K", "tau_candidates", "ncores", "seed", "verbose",
      "force_windows", "ram_check", "failure_handling", "retain_fits"
    )] <- NULL
    fit_by_tau <- vector("list", length(tau_values))
    for (tau_id in seq_along(tau_values)) {
      tau <- tau_values[[tau_id]]
      task_arguments <- arguments
      if (!is.null(outcome)) {
        outcome_algorithm <- if (is.null(outcome$algorithm)) {
          "NETCROP"
        } else {
          toupper(outcome$algorithm)
        }
        use_dcbm <- identical(outcome$model, "DCBM")
        cluster_engine <- if (outcome_algorithm == "NETCROP") {
          if (use_dcbm) "clara" else "kmeans"
        } else if (use_dcbm) {
          "pam"
        } else if (tau == 0) {
          "kmeans"
        } else {
          "clara"
        }
        cluster_options <- if (outcome_algorithm == "NETCROP") {
          utils::modifyList(
            if (use_dcbm) {
              list(metric = "manhattan", cluster.only = TRUE, samples = 5L)
            } else {
              list(nstart = 100L, iter.max = 10^7)
            },
            outcome$options$cluster_options,
            keep.null = TRUE
          )
        } else if (use_dcbm) {
          list(
            metric = "euclidean", do.swap = FALSE,
            cluster.only = TRUE, pamonce = 6
          )
        } else if (tau == 0) {
          list(nstart = 100, iter.max = 10^7)
        } else {
          list(metric = "euclidean", cluster.only = TRUE)
        }
        task_arguments <- utils::modifyList(
          list(
            laplacian = outcome$options$use_laplacian,
            normalize_laplacian = TRUE,
            handle_zero_degree_nodes = if (tau == 0) {
              "random_label"
            } else {
              "none"
            },
            row_normalize = use_dcbm,
            spectral_method = "eigen",
            spectral_engine = "RSpectra",
            spectral_options = if (outcome_algorithm == "NETCROP") {
              if (is.null(outcome$options$spectral_options)) {
                list()
              } else {
                outcome$options$spectral_options
              }
            } else {
              if (is.null(task_arguments$spectral_options)) {
                list()
              } else {
                task_arguments$spectral_options
              }
            },
            cluster_engine = cluster_engine,
            cluster_options = cluster_options
          ),
          task_arguments,
          keep.null = TRUE
        )
      }
      fit_by_tau[[tau_id]] <- do.call(
        mult_reg_spectral_cluster,
        c(list(
          A = A_list[[network_id]],
          K = K,
          tau_candidates = tau,
          ncores = ncores,
          seed = if (!is.null(seed)) {
            plot_seed(network_id, 200000L)
          } else if (!is.null(outcome$seed)) {
            outcome$seed
          } else {
            NULL
          },
          verbose = verbose,
          force_windows = force_windows,
          ram_check = ram_check,
          failure_handling = "stop",
          retain_fits = FALSE
        ), task_arguments)
      )
      register_fit(
        fit_by_tau[[tau_id]],
        network_id = network_id,
        engine = "spectral_cluster",
        tau = tau,
        parameter_source = if (is.null(outcome)) {
          "current_defaults"
        } else {
          outcome_algorithm
        }
      )
    }
    names(fit_by_tau) <- as.character(tau_values)
    fit_by_tau
  }
  use_selection_plot <- !is.null(netcrop_outcomes) ||
    !is.null(dkest_outcomes)
  for (network_id in seq_along(A_list)) {
    truth <- as.integer(truth_list[[network_id]])
    netcrop_outcome <- if (is.null(netcrop_outcomes)) {
      NULL
    } else {
      netcrop_outcomes[[network_id]]
    }
    dkest_outcome <- if (is.null(dkest_outcomes)) {
      NULL
    } else {
      dkest_outcomes[[network_id]]
    }
    selections <- NULL
    if (!is.null(netcrop_outcome)) {
      selections <- prepare_netcrop_selections(netcrop_outcome)
    }
    if (!is.null(dkest_outcome)) {
      if (is.null(dkest_outcome$tau_hat) ||
          length(dkest_outcome$tau_hat) != 1L ||
          !is.numeric(dkest_outcome$tau_hat) ||
          !is.finite(dkest_outcome$tau_hat) ||
          dkest_outcome$tau_hat < 0) {
        stop("DKEST outcome ", network_id, " lacks a valid tau_hat.",
             call. = FALSE)
      }
      selections <- rbind(
        selections,
        data.frame(
          method = "DKEST",
          tau = as.numeric(dkest_outcome$tau_hat),
          stringsAsFactors = FALSE
        )
      )
    }
    required_tau <- if (use_selection_plot) {
      unique(c(tau_candidates, selections$tau))
    } else {
      tau_candidates
    }
    
    for (engine in engines) {
      engine_label <- if (engine == "sonnet") {
        "SONNET"
      } else {
        "Spectral clustering"
      }
      parameter_outcome <- if (engine == "sonnet") {
        if (!is.null(netcrop_outcome)) netcrop_outcome else dkest_outcome
      } else {
        if (!is.null(dkest_outcome)) dkest_outcome else netcrop_outcome
      }
      fits <- if (engine == "sonnet") {
        run_sonnet_values(network_id, required_tau, parameter_outcome)
      } else {
        run_spectral_values(network_id, required_tau, parameter_outcome)
      }
      candidate_accuracy <- numeric(length(tau_candidates))
      for (tau_id in seq_along(tau_candidates)) {
        fit_id <- which(vapply(
          required_tau, tau_equal, logical(1), tau_candidates[[tau_id]]
        ))
        candidate_accuracy[[tau_id]] <- match_accuracy(
          extract_label(fits[[fit_id]], tau_candidates[[tau_id]]),
          truth
        )
      }
      if (!use_selection_plot) {
        for (tau_id in seq_along(tau_candidates)) {
          add_record(
            network_id,
            engine_label,
            paste0("tau = ", tau_candidates[[tau_id]]),
            tau_candidates[[tau_id]],
            candidate_accuracy[[tau_id]]
          )
        }
        next
      }
      
      oracle_ids <- which(candidate_accuracy == max(candidate_accuracy))
      oracle_id <- oracle_ids[which.min(tau_candidates[oracle_ids])]
      add_record(
        network_id,
        engine_label,
        "Oracle",
        tau_candidates[[oracle_id]],
        candidate_accuracy[[oracle_id]]
      )
      zero_id <- which(vapply(tau_candidates, tau_equal, logical(1), 0))
      if (length(zero_id) == 1L) {
        add_record(
          network_id,
          engine_label,
          "tau = 0",
          0,
          candidate_accuracy[[zero_id]]
        )
      }
      for (selection_id in seq_len(nrow(selections))) {
        fit_id <- which(vapply(
          required_tau,
          tau_equal,
          logical(1),
          selections$tau[[selection_id]]
        ))
        add_record(
          network_id,
          engine_label,
          selections$method[[selection_id]],
          selections$tau[[selection_id]],
          match_accuracy(
            extract_label(fits[[fit_id]], selections$tau[[selection_id]]),
            truth
          )
        )
      }
    }
  }
  accuracy_data <- do.call(rbind, accuracy_records)
  grouping <- interaction(
    accuracy_data$engine,
    accuracy_data$method,
    drop = TRUE,
    lex.order = TRUE
  )
  plot_data <- do.call(rbind, lapply(split(accuracy_data, grouping), function(x) {
    accuracy_sd <- stats::sd(x$accuracy)
    if (is.na(accuracy_sd)) accuracy_sd <- 0
    data.frame(
      engine = x$engine[[1L]],
      method = x$method[[1L]],
      mean_accuracy = mean(x$accuracy),
      sd_accuracy = accuracy_sd,
      lower = max(0, mean(x$accuracy) - accuracy_sd),
      upper = min(1, mean(x$accuracy) + accuracy_sd),
      simulations = nrow(x),
      stringsAsFactors = FALSE
    )
  }))
  plot_data$mean_accuracy_percent <- 100 * plot_data$mean_accuracy
  plot_data$lower_percent <- 100 * plot_data$lower
  plot_data$upper_percent <- 100 * plot_data$upper
  method_levels <- unique(accuracy_data$method)
  plot_data$method <- factor(plot_data$method, levels = rev(method_levels))
  errorbar_data <- plot_data[
    plot_data$simulations > 1L,
    ,
    drop = FALSE
  ]
  oracle_lines <- if (use_selection_plot) {
    plot_data[
      as.character(plot_data$method) == "Oracle",
      c("engine", "mean_accuracy_percent"),
      drop = FALSE
    ]
  } else {
    NULL
  }
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = mean_accuracy_percent, y = method)
  ) +
    ggplot2::geom_point()
  if (nrow(errorbar_data) > 0L) {
    plot <- plot + ggplot2::geom_errorbar(
      data = errorbar_data,
      ggplot2::aes(xmin = lower_percent, xmax = upper_percent),
      orientation = "y"
    )
  }
  if (!is.null(oracle_lines) && nrow(oracle_lines) > 0L) {
    plot <- plot + ggplot2::geom_vline(
      data = oracle_lines,
      ggplot2::aes(xintercept = mean_accuracy_percent),
      linetype = "dashed",
      inherit.aes = FALSE
    )
  }
  plot <- plot +
    ggplot2::facet_wrap(~engine, scales = "free_y") +
    ggplot2::labs(
      title = if (use_selection_plot) {
        "Oracle and selected regularizer clustering accuracy"
      } else {
        "Clustering accuracy by regularization candidate"
      },
      x = if (nrow(errorbar_data) > 0L) {
        "Clustering accuracy (%) [mean plus or minus SD]"
      } else {
        "Clustering accuracy (%)"
      },
      y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  attr(plot, "accuracy_data") <- accuracy_data
  attr(plot, "plot_data") <- plot_data
  attr(plot, "metadata") <- list(
    network_count = length(A_list),
    shared_dimension = dimensions[[1L]],
    K = K,
    tau_candidates = tau_candidates,
    engines = engines,
    matching_method = matching_method,
    include_netcrop_mean = include_netcrop_mean,
    include_netcrop_mode = include_netcrop_mode,
    requested_losses = losses,
    netcrop_outcomes_provided = !is.null(netcrop_outcomes),
    dkest_outcomes_provided = !is.null(dkest_outcomes),
    ncores = ncores,
    seed = seed,
    fit_parameters = fit_metadata
  )
  attr(plot, "diagnostics") <- if (length(diagnostic_records) > 0L) {
    do.call(rbind, diagnostic_records)
  } else {
    data.frame(
      network_id = integer(),
      tau = numeric(),
      seed = integer(),
      message = character(),
      engine = character(),
      parameter_source = character(),
      stringsAsFactors = FALSE
    )
  }
  print(plot)
  invisible(plot)
}


# Build a comparison plot for two or three block-model selection results.
plot_blockmodel_comparison <- function(
    ...,
    loss_scale = c("relative", "raw")) {
  outputs <- list(...)
  loss_scale <- match.arg(loss_scale)
  if (length(outputs) < 2L || length(outputs) > 3L ||
      any(!vapply(outputs, inherits, logical(1), "netcrop_blockmodel"))) {
    stop(
      "Supply two or three objects inheriting from netcrop_blockmodel.",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting.", call. = FALSE)
  }
  algorithms <- vapply(outputs, function(output) {
    if (is.null(output$algorithm)) "NETCROP" else toupper(output$algorithm)
  }, character(1))
  if (any(!algorithms %in% c("NETCROP", "ECV", "NCV")) ||
      anyDuplicated(algorithms)) {
    stop(
      "Comparison inputs must represent distinct NETCROP, ECV, or NCV results.",
      call. = FALSE
    )
  }
  loss_family <- function(loss_name) {
    switch(
      loss_name,
      sse = "squared_error",
      mse = "squared_error",
      sae = "absolute_error",
      mae = "absolute_error",
      bin_dev = "binary_deviance",
      bin_dev_mean = "binary_deviance",
      loss_name
    )
  }
  comparison_names <- lapply(outputs, function(output) {
    if (loss_scale == "raw") {
      output$losses
    } else {
      vapply(output$losses, loss_family, character(1))
    }
  })
  common_losses <- Reduce(intersect, comparison_names)
  common_K <- Reduce(
    intersect, lapply(outputs, function(output) unique(output$cv_loss$K))
  )
  common_models <- Reduce(
    intersect, lapply(outputs, `[[`, "model_candidates")
  )
  if (length(common_losses) == 0L) {
    stop(
      if (loss_scale == "raw") {
        "Raw comparison requires at least one identical loss name."
      } else {
        "The results have no comparable loss family."
      },
      call. = FALSE
    )
  }
  if (length(common_K) == 0L) {
    stop("The results have no common K candidates.", call. = FALSE)
  }
  if (length(common_models) == 0L) {
    stop("The results have no common model candidates.", call. = FALSE)
  }
  prepare_plot_data <- function(output, algorithm) {
    data <- output$cv_loss
    data$comparison_loss <- if (loss_scale == "raw") {
      data$loss
    } else {
      vapply(data$loss, loss_family, character(1))
    }
    data <- data[
      data$comparison_loss %in% common_losses & data$K %in% common_K &
        data$model %in% common_models & is.finite(data$average_loss),
      , drop = FALSE
    ]
    grouping <- interaction(
      data$K, data$model, data$comparison_loss,
      drop = TRUE, lex.order = TRUE
    )
    result <- do.call(rbind, lapply(split(data, grouping), function(x) {
      data.frame(
        K = x$K[1L],
        model = x$model[1L],
        loss = x$comparison_loss[1L],
        mean_loss = mean(x$average_loss),
        stringsAsFactors = FALSE
      )
    }))
    result$algorithm <- algorithm
    result
  }
  plot_data <- do.call(rbind, Map(prepare_plot_data, outputs, algorithms))
  plot_data$plotted_loss <- plot_data$mean_loss
  y_label <- "Mean CV loss"
  plot_caption <- NULL
  if (loss_scale == "relative") {
    normalization_group <- interaction(
      plot_data$algorithm, plot_data$loss,
      drop = TRUE, lex.order = TRUE
    )
    plot_data$plotted_loss <- unsplit(
      lapply(split(plot_data$mean_loss, normalization_group), function(x) {
        value_range <- range(x)
        if (diff(value_range) == 0) return(rep(0, length(x)))
        (x - value_range[1L]) / diff(value_range)
      }),
      normalization_group
    )
    y_label <- "Within-algorithm relative CV loss"
    plot_caption <- paste0(
      "Relative mode compares compatible families; SSE/MSE and SAE/MAE ",
      "are normalized within each algorithm."
    )
  }
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = K,
      y = plotted_loss,
      color = model,
      linetype = algorithm,
      group = interaction(model, algorithm)
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = sort(common_K)) +
    ggplot2::facet_wrap(
      ~loss, scales = if (loss_scale == "raw") "free_y" else "fixed"
    ) +
    ggplot2::labs(
      title = paste(algorithms, collapse = ", "),
      subtitle = "Block-model cross-validation comparison",
      caption = plot_caption,
      x = "Number of communities (K)",
      y = y_label,
      color = "Model",
      linetype = "Algorithm"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

################################################################################
# library(Matrix)
# system.time(
#   net <- generate_dcbm(n = 3000, K = 5, ncores = 5, average_degree = 80)
# )
# 
# # system.time(
# #   net_all <- DCBM.gen(
# #     n = 10000, K = 10, B = diag(0.08, 10, 10) + matrix(0.02, 10, 10),
# #     avg.deg = 80, ncore = 8
# #   )
# # )
# 
# # net <- net_all$A
# # g_true <- net_all$member
# g_true <- attributes(net)$generator_parameters$g_true
# mean(rowSums(net))
# 
# system.time(
#   pram <- peakRAM::peakRAM({sc.out <- spectral_cluster(
#     A = net,
#     K = 5,
#     laplacian = FALSE,
#     regularize_tau = 0,
#     row_normalize = TRUE,
#     spectral_method = "eigen",
#     # spectral_engine = "RSpectra",
#     cluster_engine = "clara"
#     # cluster_options = list("metric" = "manhattan")
#   )})
# )
# 
# pram$Peak_RAM_Used_MiB / 1024^2
# label_match_greedy(sc.out$g_hat, g_true)$mismatch_rate
# 
# system.time(
#   nc.out <- netcrop_blockmodel(
#     A = net,
#     K_candidates = 1:10,
#     # num_subnetworks = 3,
#     # overlap_size = 8002,
#     nrep = 1,
#     losses = "sse",
#     model_candidates = c("SBM", "DCBM"),
#     ncores = 8,
#     verbose = TRUE,
#     force_windows = FALSE,
#     ram_check = FALSE,
#     sbm_est_options = list(
#       spectral_cluster = list(
#         cluster_engine = "clara",
#         # cluster_options = list(nstart = 30)
#         cluster_options = list(samples = 5L)
#         )
#     ),
#     dcbm_est_options = list(
#       spectral_cluster = list(
#         cluster_engine = "clara",
#         # cluster_options = list(nstart = 30)
#         cluster_options = list(samples = 5L)
#         )
#     )
#   )
# )
# nc.out
# summary(nc.out)
# plot(nc.out)
# 
# system.time(
#   ec.out <- ecv_stability_blockmodel(
#     A = net,
#     max_K = 10,
#     nrep = 1L,
#     losses = "sse",
#     ncores = 8L,
#     seed = NULL,
#     force_windows = FALSE
#   )
# )
# ec.out
# summary(ec.out)
# plot(ec.out)
# 
# plot(nc.out, ec.out)
