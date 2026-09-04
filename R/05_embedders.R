# Network embedding and latent-space fitting helpers
#
# Package loading resolves helper order and registers the optional C++
# accelerator for lsm_pgd().

# Compute an adjacency spectral embedding (ASE).
#
# For an undirected A, use the d eigenvalues of largest magnitude and form
# Z_hat = U |Lambda|^(1/2). The signed eigenvalues are retained so an optional
# rank-d reconstruction uses U Lambda U^T rather than Z_hat Z_hat^T.
#
# For a directed A, use the rank-d SVD and return left/right embeddings
# U D^(1/2) and V D^(1/2). If align_with is supplied, the left embedding is
# aligned to it and the same orthogonal rotation is applied to the right
# embedding, preserving their reconstructed matrix.
#' @rdname embedding_estimating
#' @export
ase <- function(
    A,
    d,
    directed = FALSE,
    reconstruct = FALSE,
    align_with = NULL,
    ram_check = FALSE) {
  if (!exists("eig_decomp", mode = "function", inherits = TRUE, envir = environment()) ||
      !exists("singular_decomp", mode = "function", inherits = TRUE, envir = environment()) ||
      !exists("procrustes", mode = "function", inherits = TRUE, envir = environment())) {
    stop("Required mathematical helpers are unavailable; reinstall netOP before calling ase().", call. = FALSE)
  }
  if (is.null(dim(A)) || length(dim(A)) != 2L || nrow(A) != ncol(A)) {
    stop("A must be a square matrix-like object.", call. = FALSE)
  }
  if (!(is.numeric(A) || inherits(A, "Matrix")) || any(!is.finite(A))) {
    stop("A must contain only finite numeric values.", call. = FALSE)
  }
  if (length(d) != 1L ||
      !is.numeric(d) ||
      is.na(d) ||
      !is.finite(d) ||
      d < 1 ||
      d != floor(d) ||
      d > nrow(A)) {
    stop("d must be one positive integer no larger than nrow(A).",
         call. = FALSE)
  }
  logical_inputs <- list(
    directed = directed,
    reconstruct = reconstruct,
    ram_check = ram_check
  )
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
  d <- as.integer(d)
  n <- nrow(A)
  if (!is.null(align_with)) {
    if (is.null(dim(align_with)) ||
        length(dim(align_with)) != 2L ||
        !is.numeric(align_with) ||
        !identical(dim(align_with), c(n, d)) ||
        any(!is.finite(align_with))) {
      stop("align_with must be a finite nrow(A)-by-d numeric matrix.",
           call. = FALSE)
    }
    align_with <- as.matrix(align_with)
  }

  if (ram_check) {
    ram_helpers <- c(
      "estimate_spectral_decomp_ram",
      "estimate_matrix_product_ram",
      "report_ram_preflight",
      "report_ram_formula"
    )
    missing_ram_helpers <- ram_helpers[!vapply(
      ram_helpers,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE, envir = environment()
    )]
    if (length(missing_ram_helpers) > 0L) {
      stop(
        "Required RAM helpers are unavailable; reinstall netOP. Missing: ",
        paste(missing_ram_helpers, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    decomposition_ram <- estimate_spectral_decomp_ram(
      n = n,
      p = n,
      K = d,
      method = if (directed) "svd" else "eigen",
      engine = "rspectra",
      nu = d,
      nv = d,
      dense_input = FALSE
    )
    report_ram_preflight(
      estimated_bytes = decomposition_ram$estimated_bytes,
      operation = paste0(
        "ASE ",
        if (directed) "singular" else "eigen",
        " decomposition"
      ),
      detail = paste0(n, " x ", n)
    )
    multiplication_estimates <- numeric(0)
    multiplication_details <- character(0)
    multiplication_count <- 0L
    if (!is.null(align_with)) {
      multiplication_estimates <- c(
        multiplication_estimates,
        estimate_matrix_product_ram(d, n, d),
        estimate_matrix_product_ram(d, d, d),
        estimate_matrix_product_ram(n, d, d)
      )
      multiplication_details <- c(
        multiplication_details,
        paste0(d, " x ", n, " by ", n, " x ", d),
        paste0(d, " x ", d, " by ", d, " x ", d),
        paste0(n, " x ", d, " by ", d, " x ", d)
      )
      multiplication_count <- multiplication_count + 3L
      if (directed) {
        multiplication_estimates <- c(
          multiplication_estimates,
          estimate_matrix_product_ram(n, d, d)
        )
        multiplication_details <- c(
          multiplication_details,
          paste0(n, " x ", d, " by ", d, " x ", d)
        )
        multiplication_count <- multiplication_count + 1L
      }
    }
    if (reconstruct) {
      multiplication_estimates <- c(
        multiplication_estimates,
        estimate_matrix_product_ram(n, d, n)
      )
      multiplication_details <- c(
        multiplication_details,
        paste0(n, " x ", d, " by ", d, " x ", n)
      )
      multiplication_count <- multiplication_count + 1L
    }
    if (multiplication_count > 0L) {
      largest_id <- which.max(multiplication_estimates)
      report_ram_preflight(
        estimated_bytes = multiplication_estimates[largest_id],
        operation = "ASE largest matrix multiplication",
        operation_count = multiplication_count,
        detail = multiplication_details[largest_id]
      )
    }
    ase_formula_terms <- list(list(
      estimated_bytes = decomposition_ram$estimated_bytes,
      sequential_count = 1L,
      parallel_count = 1L,
      label = "decomposition"
    ))
    if (multiplication_count > 0L) {
      ase_formula_terms <- c(ase_formula_terms, list(list(
        estimated_bytes = multiplication_estimates[largest_id],
        sequential_count = multiplication_count,
        parallel_count = 1L,
        label = "largest matrix multiplication"
      )))
    }
    report_ram_formula(
      terms = ase_formula_terms,
      operation = "ASE conservative combined preflight",
      detail = "decomposition plus counted matrix multiplications"
    )
  }

  if (!directed) {
    decomposition <- eig_decomp(
      A = A,
      d = d,
      only_values = FALSE,
      scale_by = "none",
      use_laplacian = FALSE,
      engine = "rspectra",
      force_engine = FALSE,
      order_by = "magnitude"
    )
    Z_hat <- sweep(
      decomposition$vectors,
      2L,
      sqrt(abs(decomposition$values)),
      FUN = "*"
    )
    rownames(Z_hat) <- rownames(A)
    alignment <- NULL
    if (!is.null(align_with)) {
      alignment <- procrustes(
        X = Z_hat,
        X_star = align_with,
        translate = FALSE,
        dilate = FALSE,
        sumsq = TRUE
      )
      Z_hat <- alignment$X_new
    }
    result <- list(
      Z_hat = Z_hat,
      eigenvalues = decomposition$values
    )
    if (!is.null(alignment)) {
      result$alignment <- list(
        rotation = alignment$rotation,
        sse = alignment$sse
      )
    }
    if (reconstruct) {
      scaled_vectors <- sweep(
        decomposition$vectors,
        2L,
        decomposition$values,
        FUN = "*"
      )
      A_hat <- tcrossprod(scaled_vectors, decomposition$vectors)
      dimnames(A_hat) <- dimnames(A)
      result$A_hat <- A_hat
    }
    return(result)
  }

  decomposition <- singular_decomp(
    A = A,
    d = d,
    only_values = FALSE,
    nu = d,
    nv = d,
    scale_by = "none",
    use_laplacian = FALSE,
    engine = "rspectra",
    force_engine = FALSE,
    order_by = "value"
  )
  sqrt_values <- sqrt(decomposition$values)
  Z_hat_left <- sweep(decomposition$u, 2L, sqrt_values, FUN = "*")
  Z_hat_right <- sweep(decomposition$v, 2L, sqrt_values, FUN = "*")
  rownames(Z_hat_left) <- rownames(A)
  rownames(Z_hat_right) <- colnames(A)
  alignment <- NULL
  if (!is.null(align_with)) {
    alignment <- procrustes(
      X = Z_hat_left,
      X_star = align_with,
      translate = FALSE,
      dilate = FALSE,
      sumsq = TRUE
    )
    rotation <- alignment$rotation
    Z_hat_left <- alignment$X_new
    Z_hat_right <- tcrossprod(Z_hat_right, t(rotation))
    rownames(Z_hat_right) <- colnames(A)
  }
  result <- list(
    Z_hat_left = Z_hat_left,
    Z_hat_right = Z_hat_right,
    singular_values = decomposition$values
  )
  if (!is.null(alignment)) {
    result$alignment <- list(
      rotation = alignment$rotation,
      sse = alignment$sse
    )
  }
  if (reconstruct) {
    A_hat <- tcrossprod(Z_hat_left, Z_hat_right)
    dimnames(A_hat) <- dimnames(A)
    result$A_hat <- A_hat
  }
  result
}

# Cluster an undirected network through a spectral representation.
#
# If U is supplied, it is clustered directly and A (including its degree
# information) is ignored. Otherwise asymmetric A is converted exactly as
# A + t(A), following the package convention for undirected networks.
# regularize_tau uses A_tau = A + tau * mean(deg) / n before either direct
# decomposition or Laplacian construction.
#' @rdname embedding_estimating
#' @export
spectral_cluster <- function(
    A = NULL,
    U = NULL,
    K,
    laplacian = FALSE,
    normalize_laplacian = TRUE,
    regularize_tau = 0,
    handle_zero_degree_nodes = c("none", "random_label", "remove"),
    row_normalize = FALSE,
    spectral_method = c("eigen", "svd"),
    spectral_engine = c("RSpectra", "irlba", "base"),
    spectral_options = list(),
    cluster_engine = c("clara", "kmeans", "pam"),
    cluster_options = list(),
    ram_check = FALSE,
    validate_inputs = TRUE) {
  ram_check_missing <- missing(ram_check)
  handle_zero_degree_nodes <- match.arg(handle_zero_degree_nodes)
  spectral_method <- match.arg(spectral_method)
  spectral_engine <- match.arg(spectral_engine)
  cluster_engine <- match.arg(cluster_engine)
  if (length(validate_inputs) != 1L || !is.logical(validate_inputs) ||
      is.na(validate_inputs)) {
    stop("validate_inputs must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.list(spectral_options) ||
      (length(spectral_options) > 0L &&
       (is.null(names(spectral_options)) ||
        any(!nzchar(names(spectral_options)))))) {
    stop("spectral_options must be a named list.", call. = FALSE)
  }
  if ("ram_check" %in% names(spectral_options)) {
    if (!ram_check_missing) {
      stop(
        "Specify ram_check directly, not both directly and in spectral_options.",
        call. = FALSE
      )
    }
    ram_check <- spectral_options$ram_check
    spectral_options$ram_check <- NULL
  }

  logical_inputs <- list(
    laplacian = laplacian,
    normalize_laplacian = normalize_laplacian,
    row_normalize = row_normalize,
    ram_check = ram_check
  )
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
  if (length(K) != 1L ||
      !is.numeric(K) ||
      is.na(K) ||
      !is.finite(K) ||
      K < 1 ||
      K != floor(K)) {
    stop("K must be one positive integer.", call. = FALSE)
  }
  K <- as.integer(K)
  if (length(regularize_tau) != 1L ||
      !is.numeric(regularize_tau) ||
      is.na(regularize_tau) ||
      !is.finite(regularize_tau) ||
      regularize_tau < 0) {
    stop("regularize_tau must be one finite non-negative number.",
         call. = FALSE)
  }
  option_lists <- list(
    spectral_options = spectral_options,
    cluster_options = cluster_options
  )
  invalid_options <- vapply(
    option_lists,
    function(options) {
      !is.list(options) ||
        (length(options) > 0L &&
         (is.null(names(options)) || any(!nzchar(names(options)))))
    },
    logical(1)
  )
  if (any(invalid_options)) {
    stop(
      paste(names(option_lists)[invalid_options], collapse = ", "),
      " must be a named list.",
      call. = FALSE
    )
  }

  merge_options <- function(defaults, options, protected, option_name) {
    conflicts <- intersect(names(options), protected)
    if (length(conflicts) > 0L) {
      stop(
        option_name,
        " cannot override ",
        paste(conflicts, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    utils::modifyList(defaults, options, keep.null = TRUE)
  }

  supplied_U <- !is.null(U)
  ram_formula_terms <- list()
  zero_degree_nodes <- integer(0)
  retained_nodes <- NULL
  spectral_values <- NULL
  decomposition <- NULL

  if (supplied_U) {
    if (validate_inputs &&
        (is.null(dim(U)) || length(dim(U)) != 2L || nrow(U) < 1L ||
        ncol(U) < 1L || !(is.numeric(U) || inherits(U, "Matrix")) ||
        any(!is.finite(U)))) {
      stop("U must be a finite numeric matrix-like object.", call. = FALSE)
    }
    if (!is.null(A) || handle_zero_degree_nodes != "none") {
      message(
        "U was supplied, so A and zero-degree handling are ignored."
      )
    }
    U_hat <- if (is.matrix(U)) U else as.matrix(U)
    retained_nodes <- seq_len(nrow(U_hat))
  } else {
    if (is.null(A)) {
      stop("At least one of A or U must be provided.", call. = FALSE)
    }
    required_helpers <- c(
      "is_symmetric_matrix",
      "graph_laplacian",
      "eig_decomp",
      "singular_decomp"
    )
    missing_helpers <- required_helpers[!vapply(
      required_helpers,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE, envir = environment()
    )]
    if (length(missing_helpers) > 0L) {
      stop("Required mathematical helpers are unavailable; reinstall netOP before calling spectral_cluster().",
           call. = FALSE)
    }
    if (is.null(dim(A)) || length(dim(A)) != 2L ||
        nrow(A) != ncol(A) || nrow(A) < 1L) {
      stop("A must be a non-empty square matrix-like object.", call. = FALSE)
    }
    if (!(is.numeric(A) || inherits(A, "Matrix")) ||
        any(!is.finite(A)) || any(A < 0)) {
      stop("A must contain only finite non-negative numeric values.",
           call. = FALSE)
    }
    if (!is_symmetric_matrix(A)) {
      message("Warning: A is not symmetric. Using A + t(A) instead.")
      A <- if (inherits(A, "Matrix")) A + Matrix::t(A) else A + t(A)
    }

    n <- nrow(A)
    deg <- if (inherits(A, "Matrix")) Matrix::rowSums(A) else rowSums(A)
    zero_degree_nodes <- which(deg == 0)
    retained_nodes <- if (handle_zero_degree_nodes == "none") {
      seq_len(n)
    } else {
      which(deg != 0)
    }
    if (length(retained_nodes) < K) {
      stop(
        "Fewer than K nonzero-degree nodes remain for clustering.",
        call. = FALSE
      )
    }
    A_spectral <- A[retained_nodes, retained_nodes, drop = FALSE]

    if (ram_check) {
      ram_helpers <- c(
        "estimate_spectral_decomp_ram",
        "estimate_matrix_product_ram",
        "report_ram_preflight",
        "report_ram_formula"
      )
      missing_ram_helpers <- ram_helpers[!vapply(
        ram_helpers,
        exists,
        logical(1),
        mode = "function",
        inherits = TRUE, envir = environment()
      )]
      if (length(missing_ram_helpers) > 0L) {
        stop(
          "Required RAM helpers are unavailable; reinstall netOP. Missing: ",
          paste(missing_ram_helpers, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      spectral_n <- nrow(A_spectral)
      safe_d_multiplier <- if (is.null(spectral_options$safe_d_multiplier)) {
        1
      } else {
        spectral_options$safe_d_multiplier
      }
      force_engine <- isTRUE(spectral_options$force_engine)
      dense_input <- laplacian || regularize_tau > 0
      decomposition_ram <- estimate_spectral_decomp_ram(
        n = spectral_n,
        p = spectral_n,
        K = K,
        method = spectral_method,
        engine = spectral_engine,
        force_engine = force_engine,
        safe_d_multiplier = safe_d_multiplier,
        nu = K,
        nv = if (spectral_method == "svd") 0L else K,
        dense_input = dense_input
      )
      report_ram_preflight(
        estimated_bytes = decomposition_ram$estimated_bytes,
        operation = "Spectral clustering decomposition",
        detail = paste0(
          spectral_n, " x ", spectral_n,
          if (laplacian) ", dense Laplacian input" else ""
        )
      )
      ram_formula_terms <- c(ram_formula_terms, list(list(
        estimated_bytes = decomposition_ram$estimated_bytes,
        sequential_count = 1L,
        parallel_count = 1L,
        label = "decomposition"
      )))
      if (laplacian) {
        laplacian_operation_count <- if (normalize_laplacian) 2L else 1L
        laplacian_operation_ram <- estimate_matrix_product_ram(
          spectral_n,
          spectral_n,
          spectral_n
        )
        report_ram_preflight(
          estimated_bytes = laplacian_operation_ram,
          operation = "Spectral clustering largest matrix operation",
          operation_count = laplacian_operation_count,
          detail = paste0(
            spectral_n, " x ", spectral_n,
            " dense Laplacian construction"
          )
        )
        ram_formula_terms <- c(ram_formula_terms, list(list(
          estimated_bytes = laplacian_operation_ram,
          sequential_count = laplacian_operation_count,
          parallel_count = 1L,
          label = "largest Laplacian matrix operation"
        )))
      } else if (regularize_tau > 0) {
        regularization_ram <- 1.25 * 8 * as.double(spectral_n)^2
        report_ram_preflight(
          estimated_bytes = regularization_ram,
          operation = "Spectral clustering dense regularization",
          detail = paste0(spectral_n, " x ", spectral_n)
        )
        ram_formula_terms <- c(ram_formula_terms, list(list(
          estimated_bytes = regularization_ram,
          sequential_count = 1L,
          parallel_count = 1L,
          label = "dense regularization"
        )))
      }
    }

    if (laplacian) {
      spectral_matrix <- graph_laplacian(
        A_spectral,
        normalized = normalize_laplacian,
        tau = regularize_tau
      )
      # The smallest Laplacian eigenvectors are the leading eigenvectors of
      # -L. For SVD, shifting by a spectral upper bound gives the same vectors
      # while keeping the transformed eigenvalues non-negative.
      if (spectral_method == "eigen") {
        decomposition_matrix <- -spectral_matrix
        spectral_order <- "value"
        spectral_shift <- 0
      } else {
        absolute_row_sums <- if (inherits(spectral_matrix, "Matrix")) {
          Matrix::rowSums(abs(spectral_matrix))
        } else {
          rowSums(abs(spectral_matrix))
        }
        spectral_shift <- max(absolute_row_sums)
        decomposition_matrix <- -spectral_matrix
        if (inherits(spectral_matrix, "Matrix")) {
          Matrix::diag(decomposition_matrix) <-
            Matrix::diag(decomposition_matrix) + spectral_shift
        } else {
          diag(decomposition_matrix) <-
            diag(decomposition_matrix) + spectral_shift
        }
        spectral_order <- "value"
      }
    } else {
      spectral_matrix <- A_spectral
      if (regularize_tau > 0) {
        deg_spectral <- if (inherits(A_spectral, "Matrix")) {
          Matrix::rowSums(A_spectral)
        } else {
          rowSums(A_spectral)
        }
        spectral_matrix <- spectral_matrix +
          regularize_tau * mean(deg_spectral) / nrow(A_spectral)
      }
      decomposition_matrix <- spectral_matrix
      spectral_order <- "magnitude"
      spectral_shift <- 0
    }

    helper_engine <- if (spectral_engine == "RSpectra") {
      "rspectra"
    } else {
      spectral_engine
    }
    protected_spectral <- c(
      "A", "d", "only_values", "nu", "nv", "use_laplacian", "engine",
      "order_by"
    )
    if (spectral_method == "eigen") {
      decomposition_arguments <- merge_options(
        defaults = list(
          A = decomposition_matrix,
          d = K,
          only_values = FALSE,
          use_laplacian = FALSE,
          engine = helper_engine,
          order_by = spectral_order
        ),
        options = spectral_options,
        protected = protected_spectral,
        option_name = "spectral_options"
      )
      decomposition <- do.call(eig_decomp, decomposition_arguments)
      U_hat <- decomposition$vectors
      spectral_values <- if (laplacian) {
        -decomposition$values
      } else {
        decomposition$values
      }
    } else {
      decomposition_arguments <- merge_options(
        defaults = list(
          A = decomposition_matrix,
          d = K,
          only_values = FALSE,
          nu = K,
          nv = 0L,
          use_laplacian = FALSE,
          engine = helper_engine,
          order_by = spectral_order
        ),
        options = spectral_options,
        protected = protected_spectral,
        option_name = "spectral_options"
      )
      decomposition <- do.call(singular_decomp, decomposition_arguments)
      U_hat <- decomposition$u
      spectral_values <- if (laplacian) {
        spectral_shift - decomposition$values
      } else {
        decomposition$values
      }
    }
    rownames(U_hat) <- rownames(A)[retained_nodes]
  }

  if (nrow(U_hat) < K) {
    stop("K cannot exceed the number of rows being clustered.",
         call. = FALSE)
  }
  if (row_normalize) {
    row_norms <- sqrt(rowSums(U_hat^2))
    positive_norm <- row_norms > 0
    U_hat[positive_norm, ] <-
      U_hat[positive_norm, , drop = FALSE] / row_norms[positive_norm]
  }

  if (ram_check && cluster_engine == "pam" && K > 1L) {
    if (!exists("report_ram_preflight", mode = "function", inherits = TRUE, envir = environment())) {
      stop(
        "Required RAM helpers are unavailable; reinstall netOP.",
        call. = FALSE
      )
    }
    clustered_n <- nrow(U_hat)
    dissimilarity_bytes <- 1.5 * 8 *
      as.double(clustered_n) * (clustered_n - 1) / 2
    report_ram_preflight(
      estimated_bytes = dissimilarity_bytes,
      operation = "Spectral clustering PAM dissimilarities",
      detail = paste0(clustered_n, " observations")
    )
    ram_formula_terms <- c(ram_formula_terms, list(list(
      estimated_bytes = dissimilarity_bytes,
      sequential_count = 1L,
      parallel_count = 1L,
      label = "PAM dissimilarities"
    )))
  }
  if (ram_check && length(ram_formula_terms) > 0L) {
    report_ram_formula(
      terms = ram_formula_terms,
      operation = "Spectral clustering conservative combined preflight",
      detail = "explicit sequential and parallel factors"
    )
  }

  if (K == 1L) {
    labels_retained <- rep.int(1L, nrow(U_hat))
    cluster_fit <- NULL
  } else if (cluster_engine == "kmeans") {
    cluster_arguments <- merge_options(
      defaults = list(
        x = U_hat,
        centers = K,
        nstart = 100L,
        iter.max = 1000L
      ),
      options = cluster_options,
      protected = c("x", "centers"),
      option_name = "cluster_options"
    )
    cluster_fit <- do.call(stats::kmeans, cluster_arguments)
    labels_retained <- cluster_fit$cluster
  } else {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop(
        "The cluster package is required for PAM or CLARA clustering.",
        call. = FALSE
      )
    }
    if (cluster_engine == "clara") {
      cluster_arguments <- merge_options(
        defaults = list(x = U_hat, k = K),
        options = cluster_options,
        protected = c("x", "k"),
        option_name = "cluster_options"
      )
      cluster_fit <- do.call(cluster::clara, cluster_arguments)
    } else {
      cluster_arguments <- merge_options(
        defaults = list(
          x = U_hat,
          k = K,
          pamonce = 6L,
          cluster.only = TRUE,
          keep.diss = FALSE,
          keep.data = FALSE
        ),
        options = cluster_options,
        protected = c("x", "k"),
        option_name = "cluster_options"
      )
      cluster_fit <- do.call(cluster::pam, cluster_arguments)
    }
    labels_retained <- if (is.atomic(cluster_fit)) {
      as.integer(cluster_fit)
    } else {
      as.integer(cluster_fit$clustering)
    }
  }

  if (supplied_U || handle_zero_degree_nodes == "none") {
    g_hat <- as.integer(labels_retained)
    names(g_hat) <- rownames(U_hat)
  } else {
    g_hat <- rep.int(NA_integer_, nrow(A))
    g_hat[retained_nodes] <- as.integer(labels_retained)
    if (handle_zero_degree_nodes == "random_label" &&
        length(zero_degree_nodes) > 0L) {
      g_hat[zero_degree_nodes] <- sample.int(
        K,
        size = length(zero_degree_nodes),
        replace = TRUE
      )
    }
    names(g_hat) <- rownames(A)
  }

  list(
    g_hat = g_hat,
    U_hat = U_hat,
    spectral_values = spectral_values,
    zero_degree_nodes = zero_degree_nodes,
    retained_nodes = retained_nodes,
    cluster_fit = cluster_fit,
    spectral_method = if (supplied_U) NULL else spectral_method,
    spectral_engine = if (supplied_U) NULL else spectral_engine,
    cluster_engine = cluster_engine
  )
}

# Fit an undirected logistic latent-space model by projected gradient ascent.
#
# Initialization uses USVT, a stable logit transform, the correct solution of
# (n I + 1 1^T) alpha = rowSums(theta), and double-centering without forming an
# n-by-n centering matrix. The diagonal is excluded from objective and gradient.
# Optional direct initial values may replace either spectral initializer.
#' @rdname embedding_estimating
#' @export
lsm_pgd <- function(
    A,
    d,
    step_size = 0.3,
    niter = 100L,
    trace = FALSE,
    Z_init = NULL,
    alpha_init = NULL,
    epsilon = 1e-6,
    use_cpp = TRUE,
    ram_check = FALSE) {
  required_helpers <- c(
    "usvt",
    "clip_probabilities",
    "eig_decomp",
    "outer_add",
    "softplus"
  )
  missing_helpers <- required_helpers[!vapply(
    required_helpers,
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE, envir = environment()
  )]
  if (length(missing_helpers) > 0L) {
    stop("Required mathematical helpers are unavailable; reinstall netOP before calling lsm_pgd().",
         call. = FALSE)
  }
  if (is.null(dim(A)) || length(dim(A)) != 2L || nrow(A) != ncol(A)) {
    stop("A must be a square matrix-like object.", call. = FALSE)
  }
  if (!(is.numeric(A) || inherits(A, "Matrix")) || any(!is.finite(A))) {
    stop("A must contain only finite numeric values.", call. = FALSE)
  }
  if (inherits(A, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE) ||
        !isTRUE(Matrix::isSymmetric(A))) {
      stop("A must be symmetric.", call. = FALSE)
    }
  } else if (!isTRUE(isSymmetric(A))) {
    stop("A must be symmetric.", call. = FALSE)
  }
  n <- nrow(A)
  if (length(d) != 1L ||
      !is.numeric(d) ||
      is.na(d) ||
      !is.finite(d) ||
      d < 1 ||
      d != floor(d) ||
      d > n) {
    stop("d must be one positive integer no larger than nrow(A).",
         call. = FALSE)
  }
  if (length(step_size) != 1L ||
      !is.numeric(step_size) ||
      is.na(step_size) ||
      !is.finite(step_size) ||
      step_size <= 0) {
    stop("step_size must be one finite positive number.", call. = FALSE)
  }
  if (length(niter) != 1L ||
      !is.numeric(niter) ||
      is.na(niter) ||
      !is.finite(niter) ||
      niter < 1 ||
      niter != floor(niter)) {
    stop("niter must be one positive integer.", call. = FALSE)
  }
  logical_inputs <- list(
    trace = trace,
    use_cpp = use_cpp,
    ram_check = ram_check
  )
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
  if (length(epsilon) != 1L ||
      !is.numeric(epsilon) ||
      is.na(epsilon) ||
      !is.finite(epsilon) ||
      epsilon <= 0 ||
      epsilon >= 0.5) {
    stop("epsilon must be one finite number strictly between 0 and 0.5.",
         call. = FALSE)
  }
  d <- as.integer(d)
  niter <- as.integer(niter)

  if (ram_check) {
    ram_helpers <- c(
      "estimate_spectral_decomp_ram",
      "estimate_matrix_product_ram",
      "report_ram_preflight",
      "report_ram_formula"
    )
    missing_ram_helpers <- ram_helpers[!vapply(
      ram_helpers,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE, envir = environment()
    )]
    if (length(missing_ram_helpers) > 0L) {
      stop(
        "Required RAM helpers are unavailable; reinstall netOP. Missing: ",
        paste(missing_ram_helpers, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    usvt_dimension <- min(ceiling(n^(1 / 3)), n)
    initialization_ram <- estimate_spectral_decomp_ram(
      n = n,
      p = n,
      K = usvt_dimension,
      method = "svd",
      engine = "irlba",
      nu = usvt_dimension,
      nv = usvt_dimension,
      dense_input = TRUE
    )$estimated_bytes
    if (is.null(Z_init)) {
      initialization_ram <- max(
        initialization_ram,
        estimate_spectral_decomp_ram(
          n = n,
          p = n,
          K = d,
          method = "eigen",
          engine = "rspectra",
          nu = d,
          nv = d,
          dense_input = TRUE
        )$estimated_bytes
      )
    }
    report_ram_preflight(
      estimated_bytes = initialization_ram,
      operation = "LSM PGD spectral initialization",
      detail = paste0(n, " x ", n)
    )
    report_ram_preflight(
      estimated_bytes = 1.25 * 8 * 6 * as.double(n)^2,
      operation = "LSM PGD dense iterative state",
      detail = "approximately six dense n x n matrices"
    )
    report_ram_preflight(
      estimated_bytes = estimate_matrix_product_ram(n, d, n),
      operation = "LSM PGD largest matrix multiplication",
      operation_count = 2 * niter + 2,
      detail = paste0(n, " x ", d, " by ", d, " x ", n)
    )
    report_ram_formula(
      terms = list(
        list(
          estimated_bytes = initialization_ram,
          sequential_count = 1L,
          parallel_count = 1L,
          label = "spectral initialization"
        ),
        list(
          estimated_bytes = 1.25 * 8 * 6 * as.double(n)^2,
          sequential_count = 1L,
          parallel_count = 1L,
          label = "dense iterative state"
        ),
        list(
          estimated_bytes = estimate_matrix_product_ram(n, d, n),
          sequential_count = 2 * niter + 2,
          parallel_count = 1L,
          label = "largest matrix multiplication"
        )
      ),
      operation = "LSM PGD conservative combined preflight",
      detail = "explicit sequential and parallel factors"
    )
  }

  # PGD forms dense theta, P_hat, and residual matrices at every iteration, so
  # converting once here avoids fragile base-method dispatch on Matrix classes.
  A <- as.matrix(A)
  if (any(!as.numeric(A) %in% c(0, 1))) {
    stop("A must contain only binary values 0 and 1.", call. = FALSE)
  }
  if (any(diag(A) != 0)) {
    stop("A must have a zero diagonal.", call. = FALSE)
  }

  P_tilde <- usvt(
    A = A,
    lower_clip = epsilon,
    upper_clip = 1 - epsilon
  )
  P_tilde <- (P_tilde + t(P_tilde)) / 2
  P_tilde <- clip_probabilities(P_tilde, eps = epsilon)
  theta_tilde <- stats::qlogis(P_tilde)

  if (is.null(alpha_init)) {
    theta_row_sums <- rowSums(theta_tilde)
    alpha_0 <- theta_row_sums / n -
      sum(theta_row_sums) / (2 * as.double(n)^2)
  } else {
    if (length(alpha_init) == 1L) alpha_init <- rep(alpha_init, n)
    if (!is.numeric(alpha_init) ||
        length(alpha_init) != n ||
        anyNA(alpha_init) ||
        any(!is.finite(alpha_init))) {
      stop("alpha_init must be one finite scalar or n finite values.",
           call. = FALSE)
    }
    alpha_0 <- as.numeric(alpha_init)
  }

  if (is.null(Z_init)) {
    G <- theta_tilde - outer_add(alpha_0, alpha_0, use_cpp = use_cpp)
    G <- sweep(G, 1L, rowMeans(G), FUN = "-")
    G <- sweep(G, 2L, colMeans(G), FUN = "-")
    G <- (G + t(G)) / 2
    decomposition <- eig_decomp(
      A = G,
      d = d,
      only_values = FALSE,
      scale_by = "none",
      use_laplacian = FALSE,
      engine = "rspectra",
      force_engine = FALSE,
      order_by = "value"
    )
    positive_values <- pmax(decomposition$values, 0)
    Z_0 <- sweep(
      decomposition$vectors,
      2L,
      sqrt(positive_values),
      FUN = "*"
    )
  } else {
    if (is.null(dim(Z_init)) ||
        length(dim(Z_init)) != 2L ||
        !is.numeric(Z_init) ||
        !identical(dim(Z_init), c(n, d)) ||
        any(!is.finite(Z_init))) {
      stop("Z_init must be a finite nrow(A)-by-d numeric matrix.",
           call. = FALSE)
    }
    Z_0 <- as.matrix(Z_init)
    Z_0 <- sweep(Z_0, 2L, colMeans(Z_0), FUN = "-")
  }

  Z_spectral_norm <- norm(Z_0, type = "2")
  Z_scale <- max(Z_spectral_norm^2, sqrt(.Machine$double.eps))
  step_size_Z <- step_size / Z_scale
  step_size_alpha <- step_size / (2 * n)

  fit_r <- function() {
    Z_hat <- Z_0
    alpha_hat <- alpha_0
    objective <- numeric(niter + 1L)
    upper_indices <- upper.tri(A, diag = FALSE)
    calculate_state <- function() {
      theta <- outer_add(
        alpha_hat,
        alpha_hat,
        use_cpp = use_cpp
      ) + tcrossprod(Z_hat)
      P_hat <- stats::plogis(theta)
      objective_value <- sum(
        A[upper_indices] * theta[upper_indices] -
          softplus(theta[upper_indices])
      )
      list(theta = theta, P_hat = P_hat, objective = objective_value)
    }
    for (iteration in seq_len(niter)) {
      state <- calculate_state()
      objective[iteration] <- state$objective
      if (trace) {
        message(
          "Iteration ", iteration,
          " objective: ", format(state$objective, digits = 10)
        )
      }
      residual <- A - state$P_hat
      diag(residual) <- 0
      Z_hat <- Z_hat + step_size_Z * crossprod(residual, Z_hat)
      alpha_hat <- alpha_hat + step_size_alpha * rowSums(residual)
      Z_hat <- sweep(Z_hat, 2L, colMeans(Z_hat), FUN = "-")
    }
    state <- calculate_state()
    objective[niter + 1L] <- state$objective
    diag(state$P_hat) <- 0
    list(
      Z_hat = Z_hat,
      alpha_hat = alpha_hat,
      P_hat = state$P_hat,
      objective = objective
    )
  }

  fit <- NULL
  if (use_cpp && exists("lsm_pgd_cpp", mode = "function", inherits = TRUE, envir = environment())) {
    compiled_lsm_pgd <- get(
      "lsm_pgd_cpp",
      mode = "function",
      inherits = TRUE, envir = environment()
    )
    compiled_result <- tryCatch(
      list(
        success = TRUE,
        value = compiled_lsm_pgd(
          A = as.matrix(A),
          Z = Z_0,
          alpha = alpha_0,
          step_size_Z = step_size_Z,
          step_size_alpha = step_size_alpha,
          niter = niter,
          trace = trace
        )
      ),
      error = function(e) {
        message(
          "lsm_pgd_cpp failed; using the R implementation: ",
          conditionMessage(e)
        )
        list(success = FALSE, value = NULL)
      }
    )
    if (compiled_result$success) fit <- compiled_result$value
  }
  if (is.null(fit)) fit <- fit_r()

  rownames(fit$Z_hat) <- rownames(A)
  dimnames(fit$P_hat) <- dimnames(A)
  fit$step_size_Z <- step_size_Z
  fit$step_size_alpha <- step_size_alpha
  fit
}
