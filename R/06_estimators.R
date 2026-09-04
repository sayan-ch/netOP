# Network-model estimators
#
# Package loading resolves all internal helper relationships.
# This file currently uses base R and supports dense or Matrix-package inputs.

# Estimate the block-mean matrix of an SBM from a full or NCV adjacency matrix.
#
# With fold_nodes = NULL, A must be the full length(g)-by-length(g) matrix.
# With fold_nodes supplied, A may be any of the common NCV layouts:
#
#   * full n-by-n A (the non-fold rows are selected internally),
#   * non-fold rows by all columns,
#   * all rows by fold columns, or
#   * non-fold rows by fold columns.
#
# The rectangular estimator counts the dyads actually represented by A. It
# therefore fixes the historical denominator based only on matrix dimensions.
# For undirected networks, sufficient statistics from both available block
# orientations are pooled before division. This also handles folds whose
# community composition is unbalanced. Genuine self-pairs are removed using
# inferred original node indices rather than by assuming a square matrix.
#' @rdname embedding_estimating
#' @export
estimate_sbm <- function(
    A,
    g,
    K = max(g),
    fold_nodes = NULL,
    directed = FALSE,
    self_loops = FALSE,
    validate_inputs = TRUE) {
  if (length(validate_inputs) != 1L || !is.logical(validate_inputs) ||
      is.na(validate_inputs)) {
    stop("validate_inputs must be TRUE or FALSE.", call. = FALSE)
  }
  if (validate_inputs &&
      (is.null(dim(A)) || length(dim(A)) != 2L || any(dim(A) < 1L))) {
    stop("A must be a non-empty matrix-like object.", call. = FALSE)
  }
  if (validate_inputs &&
      (!(is.numeric(A) || inherits(A, "Matrix")) || any(!is.finite(A)))) {
    stop("A must contain only finite numeric values.", call. = FALSE)
  }
  if (validate_inputs &&
      (!is.numeric(g) || is.null(g) || anyNA(g) || any(!is.finite(g)) ||
       any(g < 1) || any(g != floor(g)))) {
    stop("g must contain positive integer community labels.", call. = FALSE)
  }
  n <- length(g)
  if (length(K) != 1L ||
      !is.numeric(K) ||
      is.na(K) ||
      !is.finite(K) ||
      K < 1 ||
      K != floor(K) ||
      any(g > K)) {
    stop("K must be a positive integer at least max(g).", call. = FALSE)
  }
  logical_inputs <- list(directed = directed, self_loops = self_loops)
  invalid_logical <- vapply(
    logical_inputs,
    function(value) {
      length(value) != 1L || !is.logical(value) || is.na(value)
    },
    logical(1)
  )
  if (any(invalid_logical)) {
    stop(
      paste(names(logical_inputs)[invalid_logical], collapse = ", "),
      " must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  K <- as.integer(K)

  if (is.null(fold_nodes)) {
    if (!identical(dim(A), c(n, n))) {
      stop(
        "Without fold_nodes, A must have length(g) rows and columns.",
        call. = FALSE
      )
    }
    row_nodes <- seq_len(n)
    col_nodes <- seq_len(n)
  } else {
    if (!is.numeric(fold_nodes) ||
        length(fold_nodes) < 1L ||
        length(fold_nodes) >= n ||
        anyNA(fold_nodes) ||
        any(!is.finite(fold_nodes)) ||
        any(fold_nodes != floor(fold_nodes)) ||
        any(fold_nodes < 1 | fold_nodes > n) ||
        anyDuplicated(fold_nodes)) {
      stop(
        paste0(
          "fold_nodes must contain between 1 and length(g) - 1 unique ",
          "node indices."
        ),
        call. = FALSE
      )
    }
    fold_nodes <- sort(as.integer(fold_nodes))
    nonfold_nodes <- setdiff(seq_len(n), fold_nodes)
    nfold <- length(fold_nodes)
    nnonfold <- length(nonfold_nodes)

    if (identical(dim(A), c(n, n))) {
      A <- A[nonfold_nodes, , drop = FALSE]
      row_nodes <- nonfold_nodes
      col_nodes <- seq_len(n)
    } else if (identical(dim(A), c(nnonfold, n))) {
      row_nodes <- nonfold_nodes
      col_nodes <- seq_len(n)
    } else if (identical(dim(A), c(n, nfold))) {
      row_nodes <- seq_len(n)
      col_nodes <- fold_nodes
    } else if (identical(dim(A), c(nnonfold, nfold))) {
      row_nodes <- nonfold_nodes
      col_nodes <- fold_nodes
    } else {
      stop(
        "A dimensions do not match a supported full or NCV layout.",
        call. = FALSE
      )
    }
  }

  g_row <- g[row_nodes]
  g_col <- g[col_nodes]
  row_groups <- lapply(seq_len(K), function(k) which(g_row == k))
  col_groups <- lapply(seq_len(K), function(k) which(g_col == k))
  numerator <- matrix(0, nrow = K, ncol = K)
  denominator <- matrix(0, nrow = K, ncol = K)
  full_undirected <- !directed && is.null(fold_nodes) &&
    identical(row_nodes, col_nodes)
  column_sequence <- if (full_undirected) {
    function(k) k:K
  } else {
    function(k) seq_len(K)
  }
  for (k in seq_len(K)) {
    row_indices <- row_groups[[k]]
    if (length(row_indices) == 0L) next
    for (l in column_sequence(k)) {
      col_indices <- col_groups[[l]]
      if (length(col_indices) == 0L) next
      block_sum <- sum(A[row_indices, col_indices, drop = FALSE])
      block_count <- length(row_indices) * length(col_indices)
      if (!self_loops && k == l) {
        col_lookup <- match(row_nodes[row_indices], col_nodes[col_indices])
        matched_rows <- which(!is.na(col_lookup))
        if (length(matched_rows) > 0L) {
          block_sum <- block_sum - sum(A[cbind(
            row_indices[matched_rows],
            col_indices[col_lookup[matched_rows]]
          )])
          block_count <- block_count - length(matched_rows)
        }
      }
      numerator[k, l] <- block_sum
      denominator[k, l] <- block_count
      if (full_undirected && l != k) {
        numerator[l, k] <- block_sum
        denominator[l, k] <- block_count
      }
    }
  }
  if (!directed && !full_undirected) {
    numerator <- numerator + t(numerator)
    denominator <- denominator + t(denominator)
  }

  B_hat <- numerator / denominator
  B_hat[denominator == 0] <- NA_real_
  dimnames(B_hat) <- list(
    paste0("community_", seq_len(K)),
    paste0("community_", seq_len(K))
  )
  B_hat
}

# Estimate degree-corrected SBM block and node parameters.
#
# method = "plugin" implements the degree/block-sum plug-in estimator and
# retains the historical small block-sum stabilizer. method = "spectral" uses
# row norms of a spectral representation. A full square network uses
# eig_decomp(); rectangular NCV layouts use singular_decomp() so both row and
# column representations are available.
#
# row_norm may be either a full length(g) vector, a vector corresponding to the
# columns of A when all row nodes also occur among those columns, or a named
# list with numeric components row and col matching nrow(A) and ncol(A).
#' @rdname embedding_estimating
#' @export
estimate_dcbm <- function(
    A,
    g,
    K = max(g),
    method = c("plugin", "spectral"),
    fold_nodes = NULL,
    row_norm = NULL,
    psi_omit = 0L,
    stabilizer = 0.01,
    spectral_engine = c("rspectra", "irlba", "base"),
    spectral_options = list(),
    validate_inputs = TRUE) {
  method <- match.arg(method)
  spectral_engine <- match.arg(spectral_engine)
  if (length(validate_inputs) != 1L || !is.logical(validate_inputs) ||
      is.na(validate_inputs)) {
    stop("validate_inputs must be TRUE or FALSE.", call. = FALSE)
  }
  if (validate_inputs &&
      (is.null(dim(A)) || length(dim(A)) != 2L || any(dim(A) < 1L))) {
    stop("A must be a non-empty matrix-like object.", call. = FALSE)
  }
  if (validate_inputs &&
      (!(is.numeric(A) || inherits(A, "Matrix")) || any(!is.finite(A)))) {
    stop("A must contain only finite numeric values.", call. = FALSE)
  }
  if (validate_inputs &&
      (!is.numeric(g) || is.null(g) || anyNA(g) || any(!is.finite(g)) ||
       any(g < 1) || any(g != floor(g)))) {
    stop("g must contain positive integer community labels.", call. = FALSE)
  }
  n <- length(g)
  node_names <- names(g)
  if (is.null(node_names) && nrow(A) == n) {
    node_names <- rownames(A)
  }
  if (is.null(node_names)) {
    node_names <- as.character(seq_len(n))
  }
  if (length(K) != 1L ||
      !is.numeric(K) ||
      is.na(K) ||
      !is.finite(K) ||
      K < 1 ||
      K != floor(K) ||
      any(g > K)) {
    stop("K must be a positive integer at least max(g).", call. = FALSE)
  }
  if (length(psi_omit) != 1L ||
      !is.numeric(psi_omit) ||
      is.na(psi_omit) ||
      !is.finite(psi_omit) ||
      psi_omit < 0 ||
      psi_omit != floor(psi_omit) ||
      psi_omit >= n) {
    stop("psi_omit must be an integer between 0 and length(g) - 1.",
         call. = FALSE)
  }
  if (length(stabilizer) != 1L ||
      !is.numeric(stabilizer) ||
      is.na(stabilizer) ||
      !is.finite(stabilizer) ||
      stabilizer < 0) {
    stop("stabilizer must be one finite non-negative number.",
         call. = FALSE)
  }
  invalid_spectral_options <- !is.list(spectral_options) ||
    (length(spectral_options) > 0L &&
     (is.null(names(spectral_options)) ||
      any(!nzchar(names(spectral_options)))))
  if (invalid_spectral_options) {
    stop("spectral_options must be a named list.", call. = FALSE)
  }
  K <- as.integer(K)
  psi_omit <- as.integer(psi_omit)

  if (is.null(fold_nodes)) {
    if (!identical(dim(A), c(n, n))) {
      stop(
        "Without fold_nodes, A must have length(g) rows and columns.",
        call. = FALSE
      )
    }
    row_nodes <- seq_len(n)
    col_nodes <- seq_len(n)
  } else {
    if (psi_omit > 0L) {
      stop("psi_omit and fold_nodes cannot be used together.", call. = FALSE)
    }
    if (!is.numeric(fold_nodes) ||
        length(fold_nodes) < 1L ||
        length(fold_nodes) >= n ||
        anyNA(fold_nodes) ||
        any(!is.finite(fold_nodes)) ||
        any(fold_nodes != floor(fold_nodes)) ||
        any(fold_nodes < 1 | fold_nodes > n) ||
        anyDuplicated(fold_nodes)) {
      stop(
        paste0(
          "fold_nodes must contain between 1 and length(g) - 1 unique ",
          "node indices."
        ),
        call. = FALSE
      )
    }
    fold_nodes <- sort(as.integer(fold_nodes))
    nonfold_nodes <- setdiff(seq_len(n), fold_nodes)
    nfold <- length(fold_nodes)
    nnonfold <- length(nonfold_nodes)

    if (identical(dim(A), c(n, n))) {
      A <- A[nonfold_nodes, , drop = FALSE]
      row_nodes <- nonfold_nodes
      col_nodes <- seq_len(n)
    } else if (identical(dim(A), c(nnonfold, n))) {
      row_nodes <- nonfold_nodes
      col_nodes <- seq_len(n)
    } else if (identical(dim(A), c(n, nfold))) {
      row_nodes <- seq_len(n)
      col_nodes <- fold_nodes
    } else if (identical(dim(A), c(nnonfold, nfold))) {
      row_nodes <- nonfold_nodes
      col_nodes <- fold_nodes
    } else {
      stop(
        "A dimensions do not match a supported full or NCV layout.",
        call. = FALSE
      )
    }
  }

  g_row <- g[row_nodes]
  g_col <- g[col_nodes]
  row_groups <- lapply(seq_len(K), function(k) which(g_row == k))
  col_groups <- lapply(seq_len(K), function(k) which(g_col == k))

  if (method == "plugin") {
    B_hat <- matrix(0, nrow = K, ncol = K)
    for (k in seq_len(K)) {
      for (l in k:K) {
        B_hat[k, l] <- B_hat[l, k] <-
          sum(A[row_groups[[k]], col_groups[[l]], drop = FALSE]) + stabilizer
      }
    }

    if (is.null(fold_nodes)) {
      degree <- if (inherits(A, "Matrix")) {
        Matrix::rowSums(A)
      } else {
        rowSums(A)
      }
      psi_denominator <- rowSums(B_hat)[g]
      psi_hat <- as.numeric(degree / psi_denominator)
      psi_hat[psi_denominator <= 0] <- NA_real_
      if (psi_omit > 0L) {
        psi_hat <- psi_hat[-seq_len(psi_omit)]
      }
    } else {
      degree <- if (inherits(A, "Matrix")) {
        Matrix::colSums(A)
      } else {
        colSums(A)
      }
      psi_denominator <- colSums(B_hat)[g_col]
      psi_col <- as.numeric(degree / psi_denominator)
      psi_col[psi_denominator <= 0] <- NA_real_
      psi_hat <- psi_col[match(fold_nodes, col_nodes)]
    }
  } else {
    required_helpers <- c("eig_decomp", "singular_decomp")
    missing_helpers <- required_helpers[!vapply(
      required_helpers,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE, envir = environment()
    )]
    if (length(missing_helpers) > 0L) {
      stop(
        "Required decomposition helpers are unavailable; reinstall netOP.",
        call. = FALSE
      )
    }

    if (is.null(row_norm)) {
      protected_options <- c(
        "A", "d", "only_values", "nu", "nv", "engine", "use_laplacian"
      )
      conflicts <- intersect(names(spectral_options), protected_options)
      if (length(conflicts) > 0L) {
        stop(
          "spectral_options cannot override ",
          paste(conflicts, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      if (identical(row_nodes, col_nodes) && nrow(A) == ncol(A)) {
        decomposition_arguments <- utils::modifyList(
          list(
            A = A,
            d = K,
            only_values = FALSE,
            use_laplacian = FALSE,
            engine = spectral_engine,
            order_by = "magnitude"
          ),
          spectral_options,
          keep.null = TRUE
        )
        decomposition <- do.call(eig_decomp, decomposition_arguments)
        psi_full <- sqrt(rowSums(decomposition$vectors^2))
        psi_row <- psi_full
        psi_col <- psi_full
      } else {
        decomposition_arguments <- utils::modifyList(
          list(
            A = A,
            d = K,
            only_values = FALSE,
            nu = K,
            nv = K,
            use_laplacian = FALSE,
            engine = spectral_engine,
            order_by = "value"
          ),
          spectral_options,
          keep.null = TRUE
        )
        decomposition <- do.call(singular_decomp, decomposition_arguments)
        psi_row <- sqrt(rowSums(decomposition$u^2))
        psi_col <- sqrt(rowSums(decomposition$v^2))
      }
    } else if (is.list(row_norm)) {
      if (!identical(sort(names(row_norm)), c("col", "row")) ||
          !is.numeric(row_norm$row) ||
          !is.numeric(row_norm$col) ||
          length(row_norm$row) != nrow(A) ||
          length(row_norm$col) != ncol(A) ||
          any(!is.finite(row_norm$row)) ||
          any(!is.finite(row_norm$col)) ||
          any(row_norm$row < 0) ||
          any(row_norm$col < 0)) {
        stop(
          paste0(
            "A list row_norm must contain finite non-negative row and col ",
            "vectors matching A."
          ),
          call. = FALSE
        )
      }
      psi_row <- as.numeric(row_norm$row)
      psi_col <- as.numeric(row_norm$col)
    } else {
      if (!is.numeric(row_norm) ||
          anyNA(row_norm) ||
          any(!is.finite(row_norm)) ||
          any(row_norm < 0)) {
        stop("row_norm must contain finite non-negative values.",
             call. = FALSE)
      }
      if (length(row_norm) == n) {
        psi_row <- as.numeric(row_norm[row_nodes])
        psi_col <- as.numeric(row_norm[col_nodes])
      } else if (length(row_norm) == ncol(A) &&
                 all(row_nodes %in% col_nodes)) {
        psi_col <- as.numeric(row_norm)
        psi_row <- psi_col[match(row_nodes, col_nodes)]
      } else {
        stop(
          paste0(
            "A vector row_norm must have length(g) entries, or one per ",
            "column when all row nodes occur among the columns."
          ),
          call. = FALSE
        )
      }
    }

    B_hat <- matrix(NA_real_, nrow = K, ncol = K)
    for (k in seq_len(K)) {
      for (l in k:K) {
        denominator <- sum(psi_row[row_groups[[k]]]) *
          sum(psi_col[col_groups[[l]]])
        if (denominator > 0) {
          B_hat[k, l] <- B_hat[l, k] <-
            sum(A[row_groups[[k]], col_groups[[l]], drop = FALSE]) /
            denominator
        }
      }
    }

    if (is.null(fold_nodes)) {
      psi_hat <- as.numeric(psi_col)
      if (psi_omit > 0L) {
        psi_hat <- psi_hat[-seq_len(psi_omit)]
      }
    } else {
      psi_hat <- as.numeric(psi_col[match(fold_nodes, col_nodes)])
    }
  }

  dimnames(B_hat) <- list(
    paste0("community_", seq_len(K)),
    paste0("community_", seq_len(K))
  )
  psi_nodes <- if (is.null(fold_nodes)) {
    if (psi_omit > 0L) seq.int(psi_omit + 1L, n) else seq_len(n)
  } else {
    fold_nodes
  }
  names(psi_hat) <- node_names[psi_nodes]
  list(B_hat = B_hat, psi_hat = psi_hat)
}

# Estimate an SBM probability matrix directly.
#
# In NCV mode, return the fold-by-fold probability matrix implied by the block
# estimate. Otherwise return probabilities for all nodes. The function name is
# retained exactly as requested; the returned object follows the P_hat naming
# convention.
#' @rdname embedding_estimating
#' @export
estimate_sbm_P_hat <- function(
    A,
    g,
    K = max(g),
    fold_nodes = NULL,
    directed = FALSE,
    self_loops = FALSE,
    lower_clip = 0,
    upper_clip = 1) {
  if (length(lower_clip) != 1L ||
      !is.numeric(lower_clip) ||
      is.na(lower_clip) ||
      !is.finite(lower_clip) ||
      length(upper_clip) != 1L ||
      !is.numeric(upper_clip) ||
      is.na(upper_clip) ||
      !is.finite(upper_clip) ||
      lower_clip > upper_clip) {
    stop("Clipping bounds must be finite scalars with lower_clip <= upper_clip.",
         call. = FALSE)
  }

  B_hat <- estimate_sbm(
    A = A,
    g = g,
    K = K,
    fold_nodes = fold_nodes,
    directed = directed,
    self_loops = self_loops
  )
  target_nodes <- if (is.null(fold_nodes)) {
    seq_along(g)
  } else {
    sort(as.integer(fold_nodes))
  }
  g_target <- g[target_nodes]
  P_hat <- B_hat[g_target, g_target, drop = FALSE]
  P_hat <- pmin(pmax(P_hat, lower_clip), upper_clip)
  if (!self_loops) {
    diag(P_hat) <- 0
  }
  node_names <- names(g)
  if (is.null(node_names)) {
    node_names <- as.character(seq_along(g))
  }
  dimnames(P_hat) <- list(node_names[target_nodes], node_names[target_nodes])
  P_hat
}

# Estimate a DCBM probability matrix directly.
#
# P_hat[i, j] = psi_hat[i] * psi_hat[j] * B_hat[g[i], g[j]]. In
# NCV mode psi_hat is estimated for fold nodes, so the result is fold-by-fold.
# With psi_omit, the result covers the retained trailing nodes.
#' @rdname embedding_estimating
#' @export
estimate_dcbm_P_hat <- function(
    A,
    g,
    K = max(g),
    method = c("plugin", "spectral"),
    fold_nodes = NULL,
    row_norm = NULL,
    psi_omit = 0L,
    stabilizer = 0.01,
    spectral_engine = c("rspectra", "irlba", "base"),
    spectral_options = list(),
    self_loops = FALSE,
    lower_clip = 0,
    upper_clip = 1) {
  method <- match.arg(method)
  spectral_engine <- match.arg(spectral_engine)
  if (length(self_loops) != 1L ||
      !is.logical(self_loops) ||
      is.na(self_loops)) {
    stop("self_loops must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(lower_clip) != 1L ||
      !is.numeric(lower_clip) ||
      is.na(lower_clip) ||
      !is.finite(lower_clip) ||
      length(upper_clip) != 1L ||
      !is.numeric(upper_clip) ||
      is.na(upper_clip) ||
      !is.finite(upper_clip) ||
      lower_clip > upper_clip) {
    stop("Clipping bounds must be finite scalars with lower_clip <= upper_clip.",
         call. = FALSE)
  }

  fit <- estimate_dcbm(
    A = A,
    g = g,
    K = K,
    method = method,
    fold_nodes = fold_nodes,
    row_norm = row_norm,
    psi_omit = psi_omit,
    stabilizer = stabilizer,
    spectral_engine = spectral_engine,
    spectral_options = spectral_options
  )
  target_nodes <- if (!is.null(fold_nodes)) {
    sort(as.integer(fold_nodes))
  } else if (psi_omit > 0L) {
    seq.int(as.integer(psi_omit) + 1L, length(g))
  } else {
    seq_along(g)
  }
  g_target <- g[target_nodes]
  block_probabilities <- fit$B_hat[g_target, g_target, drop = FALSE]
  P_hat <- tcrossprod(fit$psi_hat) * block_probabilities
  P_hat <- pmin(pmax(P_hat, lower_clip), upper_clip)
  if (!self_loops) {
    diag(P_hat) <- 0
  }
  node_names <- names(g)
  if (is.null(node_names)) {
    node_names <- as.character(seq_along(g))
  }
  dimnames(P_hat) <- list(node_names[target_nodes], node_names[target_nodes])
  P_hat
}

# Relabel communities to maximize agreement with a fixed labeling.
#
# The default uses inexpensive exact shortcuts for identical and two-community
# inputs, followed by the historical greedy maximum-overlap assignment for
# larger K. algorithm = "hungarian" uses a dependency-free O(K^3) assignment
# implementation and therefore maximizes total agreement globally.
#' @rdname embedding_estimating
#' @export
label_match_greedy <- function(
    match_this,
    standard,
    K = max(c(match_this, standard)),
    algorithm = c("greedy", "hungarian"),
    return_mapping = FALSE) {
  algorithm <- match.arg(algorithm)
  if (!is.numeric(match_this) ||
      !is.numeric(standard) ||
      length(match_this) != length(standard) ||
      length(match_this) < 1L ||
      anyNA(match_this) ||
      anyNA(standard) ||
      any(!is.finite(match_this)) ||
      any(!is.finite(standard)) ||
      any(match_this < 1) ||
      any(standard < 1) ||
      any(match_this != floor(match_this)) ||
      any(standard != floor(standard))) {
    stop(
      paste0(
        "match_this and standard must be equal-length, non-empty vectors ",
        "of positive integer labels."
      ),
      call. = FALSE
    )
  }
  if (length(K) != 1L ||
      !is.numeric(K) ||
      is.na(K) ||
      !is.finite(K) ||
      K < 1 ||
      K != floor(K) ||
      any(match_this > K) ||
      any(standard > K)) {
    stop(
      "K must be a positive integer covering both label vectors.",
      call. = FALSE
    )
  }
  if (length(return_mapping) != 1L ||
      !is.logical(return_mapping) ||
      is.na(return_mapping)) {
    stop("return_mapping must be TRUE or FALSE.", call. = FALSE)
  }
  K <- as.integer(K)
  match_this <- as.integer(match_this)
  standard <- as.integer(standard)
  overlap <- table(
    factor(match_this, levels = seq_len(K)),
    factor(standard, levels = seq_len(K))
  )
  overlap <- unclass(overlap)

  if (identical(match_this, standard) || K == 1L) {
    mapping <- seq_len(K)
  } else if (K == 2L) {
    identity_agreement <- overlap[1L, 1L] + overlap[2L, 2L]
    swapped_agreement <- overlap[1L, 2L] + overlap[2L, 1L]
    mapping <- if (identity_agreement >= swapped_agreement) {
      1:2
    } else {
      2:1
    }
  } else if (algorithm == "greedy") {
    remaining_overlap <- overlap
    mapping <- integer(K)
    for (step in seq_len(K)) {
      maximum <- max(remaining_overlap)
      selected <- which(
        remaining_overlap == maximum,
        arr.ind = TRUE
      )[1L, ]
      source_label <- selected[1L]
      fixed_label <- selected[2L]
      mapping[source_label] <- fixed_label
      remaining_overlap[source_label, ] <- -Inf
      remaining_overlap[, fixed_label] <- -Inf
    }
  } else {
    # Hungarian algorithm for a square minimization problem. Converting
    # overlap to max(overlap) - overlap makes its minimum assignment maximize
    # agreement. Slot 1 in the work vectors is the dummy zero-indexed column.
    cost <- max(overlap) - overlap
    potential_row <- numeric(K + 1L)
    potential_col <- numeric(K + 1L)
    matched_row <- integer(K + 1L)
    predecessor <- integer(K + 1L)

    for (source_label in seq_len(K)) {
      matched_row[1L] <- source_label
      current_col <- 1L
      minimum_reduced_cost <- rep(Inf, K + 1L)
      used_col <- rep(FALSE, K + 1L)

      repeat {
        used_col[current_col] <- TRUE
        current_row <- matched_row[current_col]
        delta <- Inf
        next_col <- 1L
        for (candidate_col in 2:(K + 1L)) {
          if (!used_col[candidate_col]) {
            reduced_cost <- cost[current_row, candidate_col - 1L] -
              potential_row[current_row + 1L] -
              potential_col[candidate_col]
            if (reduced_cost < minimum_reduced_cost[candidate_col]) {
              minimum_reduced_cost[candidate_col] <- reduced_cost
              predecessor[candidate_col] <- current_col
            }
            if (minimum_reduced_cost[candidate_col] < delta) {
              delta <- minimum_reduced_cost[candidate_col]
              next_col <- candidate_col
            }
          }
        }
        for (candidate_col in seq_len(K + 1L)) {
          if (used_col[candidate_col]) {
            matched_index <- matched_row[candidate_col] + 1L
            potential_row[matched_index] <-
              potential_row[matched_index] + delta
            potential_col[candidate_col] <-
              potential_col[candidate_col] - delta
          } else {
            minimum_reduced_cost[candidate_col] <-
              minimum_reduced_cost[candidate_col] - delta
          }
        }
        current_col <- next_col
        if (matched_row[current_col] == 0L) {
          break
        }
      }

      repeat {
        previous_col <- predecessor[current_col]
        matched_row[current_col] <- matched_row[previous_col]
        current_col <- previous_col
        if (current_col == 1L) {
          break
        }
      }
    }

    mapping <- integer(K)
    for (fixed_label in seq_len(K)) {
      source_label <- matched_row[fixed_label + 1L]
      mapping[source_label] <- fixed_label
    }
  }

  matched_labels <- unname(mapping[match_this])
  mismatch_rate <- mean(matched_labels != standard)
  if (!return_mapping) {
    return(list(
      matched_labels = matched_labels,
      mismatch_rate = mismatch_rate
    ))
  }
  list(
    matched_labels = matched_labels,
    mismatch_rate = mismatch_rate,
    mapping = stats::setNames(mapping, seq_len(K)),
    agreement = sum(matched_labels == standard),
    overlap = overlap
  )
}

# Relabel communities by checking every one of the K! permutations.
#
# Permutations are generated recursively and scored one at a time, avoiding a
# K!-by-K permutation matrix. Runtime remains factorial. For K > 8, an
# interactive call requires explicit confirmation; non-interactive callers
# must opt in with confirm_large = TRUE.
#' @rdname embedding_estimating
#' @export
label_match_brute_force <- function(
    match_this,
    standard,
    K = max(c(match_this, standard)),
    return_mapping = FALSE,
    confirm_large = NULL) {
  if (!is.numeric(match_this) ||
      !is.numeric(standard) ||
      length(match_this) != length(standard) ||
      length(match_this) < 1L ||
      anyNA(match_this) ||
      anyNA(standard) ||
      any(!is.finite(match_this)) ||
      any(!is.finite(standard)) ||
      any(match_this < 1) ||
      any(standard < 1) ||
      any(match_this != floor(match_this)) ||
      any(standard != floor(standard))) {
    stop(
      paste0(
        "match_this and standard must be equal-length, non-empty vectors ",
        "of positive integer labels."
      ),
      call. = FALSE
    )
  }
  if (length(K) != 1L ||
      !is.numeric(K) ||
      is.na(K) ||
      !is.finite(K) ||
      K < 1 ||
      K != floor(K) ||
      any(match_this > K) ||
      any(standard > K)) {
    stop(
      "K must be a positive integer covering both label vectors.",
      call. = FALSE
    )
  }
  if (length(return_mapping) != 1L ||
      !is.logical(return_mapping) ||
      is.na(return_mapping)) {
    stop("return_mapping must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(confirm_large) &&
      (length(confirm_large) != 1L ||
       !is.logical(confirm_large) ||
       is.na(confirm_large))) {
    stop("confirm_large must be NULL, TRUE, or FALSE.", call. = FALSE)
  }
  K <- as.integer(K)
  match_this <- as.integer(match_this)
  standard <- as.integer(standard)

  if (K > 8L) {
    warning_text <- paste0(
      "Brute-force label matching checks all ",
      format(factorial(K), scientific = FALSE, trim = TRUE),
      " permutations and can be very slow and memory consuming. ",
      "label_match_greedy() is preferred."
    )
    if (is.null(confirm_large)) {
      if (!interactive()) {
        stop(
          warning_text,
          " Set confirm_large = TRUE to continue non-interactively.",
          call. = FALSE
        )
      }
      choice <- utils::menu(
        choices = c("No, cancel", "Yes, continue"),
        title = paste0(warning_text, " Do you want to continue?")
      )
      if (choice != 2L) {
        stop("Brute-force label matching was cancelled.", call. = FALSE)
      }
    } else if (!confirm_large) {
      stop("Brute-force label matching was not confirmed.", call. = FALSE)
    }
  }

  overlap <- table(
    factor(match_this, levels = seq_len(K)),
    factor(standard, levels = seq_len(K))
  )
  overlap <- unclass(overlap)
  current_mapping <- integer(K)
  best_mapping <- integer(K)
  best_agreement <- -Inf
  permutations_evaluated <- 0

  enumerate_permutations <- function(source_label, available_fixed, score) {
    if (source_label > K) {
      permutations_evaluated <<- permutations_evaluated + 1
      if (score > best_agreement) {
        best_agreement <<- score
        best_mapping <<- current_mapping
      }
      return(invisible(NULL))
    }
    for (position in seq_along(available_fixed)) {
      fixed_label <- available_fixed[position]
      current_mapping[source_label] <<- fixed_label
      enumerate_permutations(
        source_label = source_label + 1L,
        available_fixed = available_fixed[-position],
        score = score + overlap[source_label, fixed_label]
      )
    }
    invisible(NULL)
  }
  enumerate_permutations(
    source_label = 1L,
    available_fixed = seq_len(K),
    score = 0
  )

  matched_labels <- unname(best_mapping[match_this])
  mismatch_rate <- mean(matched_labels != standard)
  if (!return_mapping) {
    return(list(
      matched_labels = matched_labels,
      mismatch_rate = mismatch_rate
    ))
  }
  list(
    matched_labels = matched_labels,
    mismatch_rate = mismatch_rate,
    mapping = stats::setNames(best_mapping, seq_len(K)),
    agreement = as.numeric(best_agreement),
    overlap = overlap,
    permutations_evaluated = permutations_evaluated
  )
}
