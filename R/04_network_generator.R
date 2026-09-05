# Network generators
#
# Every network generator returns the adjacency matrix A directly. Set
# representation = "dense" for a base R matrix or representation = "sparse"
# for a Matrix::dgCMatrix. No custom network class and no igraph object is used.
# Model parameters used to generate A are stored in the lightweight
# "generator_parameters" attribute and can be retrieved with
# get_generator_parameters(). The full n-by-n probability matrix is not stored
# as an attribute because doing so would defeat the memory benefit of sparse A.

#' Network input, generation, and probability calibration
#'
#' Read and write edge lists and generate ER, SBM, DCBM, RDPG, and latent-space
#' networks. Generators return an adjacency matrix with compact truth metadata
#' available through [get_generator_parameters()].
#'
#' @name network_generators
NULL

# Validate a single logical input.
validate_generator_logical <- function(value, name) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

# Validate a positive integer-like input.
validate_generator_count <- function(value, name, allow_zero = FALSE) {
  minimum <- if (allow_zero) 0 else 1
  if (length(value) != 1L ||
      !is.numeric(value) ||
      is.na(value) ||
      !is.finite(value) ||
      value < minimum ||
      value != floor(value)) {
    qualifier <- if (allow_zero) "nonnegative" else "positive"
    stop(name, " must be one ", qualifier, " integer.", call. = FALSE)
  }
  as.integer(value)
}

# Validate and normalize an optional RNG seed.
validate_generator_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }
  if (length(seed) != 1L ||
      !is.numeric(seed) ||
      is.na(seed) ||
      !is.finite(seed) ||
      seed < 0 ||
      seed != floor(seed)) {
    stop(
      "seed must be NULL or one nonnegative integer-like value.",
      call. = FALSE
    )
  }
  as.integer(as.double(seed) %% .Machine$integer.max)
}

# Derive a valid deterministic seed without integer overflow.
offset_generator_seed <- function(seed, offset) {
  seed <- validate_generator_seed(seed)
  if (is.null(seed)) {
    return(NULL)
  }
  as.integer((as.double(seed) + as.double(offset)) %% .Machine$integer.max)
}

# Validate ncores according to CONVENTIONS.md.
validate_generator_ncores <- function(ncores) {
  validate_generator_count(ncores, "ncores")
}

# Internal apply helper. It uses 01_basic_helpers.R when available and otherwise
# falls back to lapply, keeping this file usable on its own.
generator_lapply <- function(X, FUN, ..., ncores) {
  if (ncores > 1L && exists("uni_mclapply", mode = "function", inherits = TRUE, envir = environment())) {
    return(uni_mclapply(
      X, FUN, ..., ncores = min(length(X), ncores), stop_on_error = TRUE
    ))
  }
  if (ncores > 1L) {
    message("uni_mclapply is unavailable; using sequential lapply.")
  }
  lapply(X, FUN, ...)
}

# Attach small generating parameters without changing the matrix class.
set_generator_parameters <- function(A, parameters) {
  attr(A, "generator_parameters") <- parameters
  A
}

#' Retrieve parameters attached by a network generator
#'
#' Returns the compact generating parameters stored on a network produced by a
#' netOP generator. Both dense and sparse generated matrices are supported.
#'
#' @param A A network returned by a netOP generator.
#' @return A named list of generating parameters, or `NULL` when none is attached.
#' @examples
#' A <- generate_er(n = 200, average_degree = 5, seed = 1, ncores = 1)
#' get_generator_parameters(A)
#' @seealso [generate_er()], [generate_sbm()], [generate_rdpg()], [generate_lsm()]
#' @rdname get_generator_parameters
#' @export
get_generator_parameters <- function(A) {
  attr(A, "generator_parameters", exact = TRUE)
}

# Validate an adjacency matrix used by the edge-list I/O helpers.
validate_network_for_io <- function(A) {
  if (is.null(dim(A)) || length(dim(A)) != 2L || nrow(A) != ncol(A)) {
    stop("A must be a square matrix-like object.", call. = FALSE)
  }
  if (!(is.numeric(A) || inherits(A, "Matrix"))) {
    stop("A must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(A))) {
    stop("A must contain only finite values.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check whether a network matrix is symmetric
#'
#' Tests symmetry for either a base matrix or a sparse `Matrix` object without
#' requiring users to attach Matrix.
#'
#' @param A A finite dense or sparse network matrix.
#' @param tolerance Nonnegative numerical symmetry tolerance.
#' @return A single logical value.
#' @seealso [is_symmetric_matrix()]
#' @rdname is_symmetric_network_matrix
is_symmetric_network_matrix <- function(A, tolerance = 1e-10) {
  if (inherits(A, "Matrix")) {
    return(isTRUE(Matrix::isSymmetric(A, tol = tolerance)))
  }
  isTRUE(all.equal(A, t(A), tolerance = tolerance))
}

# Resolve a delimited edge-list format from an explicit choice or extension.
resolve_network_format <- function(file, format) {
  format <- match.arg(format, c("auto", "csv", "tsv"))
  if (format != "auto") {
    return(format)
  }
  extension <- tolower(tools::file_ext(file))
  if (extension == "csv") "csv" else "tsv"
}

#' Read and write network edge lists
#'
#' `write_network()` writes a dense or sparse network as an edge list with two
#' node columns and an optional weight column. Metadata preserves network size,
#' direction, and weighting, including when isolated nodes are absent.
#'
#' @param A A finite square dense or sparse network matrix to write.
#' @param file Edge-list file path.
#' @param directed Optional direction indicator; inferred from metadata or
#'   symmetry when omitted.
#' @param weighted Optional weight indicator; inferred from edge values when
#'   omitted.
#' @param triangle For undirected output, which matrix triangle to write.
#' @param format Delimited format: `"auto"`, `"csv"`, or `"tsv"`.
#' @param include_header Whether `write_network()` writes column names.
#' @param overwrite Whether `write_network()` may replace an existing file.
#' @param tolerance Numerical tolerance used when checking symmetry and weights.
#' @param representation Whether `read_network()` returns a sparse or dense
#'   matrix.
#' @param n Optional network size used by `read_network()` when metadata cannot
#'   determine it.
#' @param has_header Whether an input edge list has a header row.
#' @return `write_network()` invisibly returns `file`; `read_network()` returns
#'   a dense matrix or sparse `dgCMatrix`, according to `representation`.
#' @examples
#' A <- generate_er(n = 200, average_degree = 5, seed = 2, ncores = 1)
#' path <- tempfile(fileext = ".csv")
#' write_network(A, path, format = "csv")
#' A_roundtrip <- read_network(path, format = "csv")
#' unlink(path)
#' dim(A_roundtrip)
#' @seealso [generate_adjacency()]
#' @rdname network_io
#' @export
write_network <- function(
    A,
    file,
    directed = NULL,
    weighted = NULL,
    triangle = c("upper", "lower"),
    format = c("auto", "csv", "tsv"),
    include_header = TRUE,
    overwrite = FALSE,
    tolerance = 1e-10) {
  validate_network_for_io(A)
  triangle <- match.arg(triangle)
  format <- resolve_network_format(file, format)
  validate_generator_logical(include_header, "include_header")
  validate_generator_logical(overwrite, "overwrite")
  if (file.exists(file) && !overwrite) {
    stop("file already exists; set overwrite = TRUE to replace it: ", file,
         call. = FALSE)
  }

  if (is.null(directed)) {
    directed <- !is_symmetric_network_matrix(A, tolerance = tolerance)
  }
  validate_generator_logical(directed, "directed")
  if (!directed && !is_symmetric_network_matrix(A, tolerance = tolerance)) {
    stop("A must be symmetric when directed = FALSE.", call. = FALSE)
  }

  if (inherits(A, "Matrix")) {
    entries <- Matrix::summary(methods::as(methods::as(A, "generalMatrix"), "dgCMatrix"))
    edge_data <- data.frame(
      node_from = entries$i,
      node_to = entries$j,
      weight = entries$x
    )
  } else {
    edge_indices <- which(A != 0, arr.ind = TRUE)
    edge_data <- data.frame(
      node_from = edge_indices[, 1L],
      node_to = edge_indices[, 2L],
      weight = if (nrow(edge_indices) == 0L) numeric() else A[edge_indices]
    )
  }
  edge_data <- edge_data[edge_data$weight != 0, , drop = FALSE]
  if (!directed) {
    triangle_mask <- if (triangle == "upper") {
      edge_data$node_from <= edge_data$node_to
    } else {
      edge_data$node_from >= edge_data$node_to
    }
    edge_data <- edge_data[triangle_mask, , drop = FALSE]
  }
  if (nrow(edge_data) > 0L) {
    edge_order <- order(edge_data$node_from, edge_data$node_to)
    edge_data <- edge_data[edge_order, , drop = FALSE]
  }

  if (is.null(weighted)) {
    weighted <- nrow(edge_data) > 0L && any(edge_data$weight != 1)
  }
  validate_generator_logical(weighted, "weighted")
  if (!weighted && nrow(edge_data) > 0L && any(edge_data$weight != 1)) {
    stop("weighted = FALSE requires every stored edge value to equal one.",
         call. = FALSE)
  }
  if (!weighted) {
    edge_data <- edge_data[, c("node_from", "node_to"), drop = FALSE]
  }

  separator <- if (format == "csv") "," else "\t"
  metadata <- c(
    paste0("# n=", nrow(A)),
    paste0("# directed=", directed),
    paste0("# weighted=", weighted),
    paste0("# triangle=", if (directed) "all" else triangle)
  )
  writeLines(metadata, con = file)
  if (include_header) {
    write(
      paste(names(edge_data), collapse = separator),
      file = file,
      append = TRUE
    )
  }
  if (nrow(edge_data) > 0L) {
    utils::write.table(
      edge_data,
      file = file,
      append = TRUE,
      sep = separator,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }
  invisible(file)
}

#' @description `read_network()` reconstructs an edge list as a dense matrix or
#'   general sparse `dgCMatrix`. For undirected input, each off-diagonal edge
#'   must occur only once.
#' @rdname network_io
#' @export
read_network <- function(
    file,
    representation = c("sparse", "dense"),
    directed = NULL,
    weighted = NULL,
    n = NULL,
    format = c("auto", "csv", "tsv"),
    has_header = TRUE) {
  if (!file.exists(file)) {
    stop("file does not exist: ", file, call. = FALSE)
  }
  representation <- match.arg(representation)
  format <- resolve_network_format(file, format)
  validate_generator_logical(has_header, "has_header")
  lines <- readLines(file, warn = FALSE)

  metadata_lines <- lines[grepl("^#[[:space:]]*[A-Za-z_]+[[:space:]]*=", lines)]
  metadata <- list()
  for (line in metadata_lines) {
    key <- sub("^#[[:space:]]*([A-Za-z_]+)[[:space:]]*=.*$", "\\1", line)
    value <- sub("^#[[:space:]]*[A-Za-z_]+[[:space:]]*=[[:space:]]*", "", line)
    metadata[[key]] <- trimws(value)
  }

  if (is.null(n) && !is.null(metadata$n)) {
    n <- suppressWarnings(as.numeric(metadata$n))
  }
  if (is.null(directed) && !is.null(metadata$directed)) {
    directed <- switch(
      toupper(metadata$directed),
      "TRUE" = TRUE,
      "FALSE" = FALSE,
      NULL
    )
  }
  if (is.null(weighted) && !is.null(metadata$weighted)) {
    weighted <- switch(
      toupper(metadata$weighted),
      "TRUE" = TRUE,
      "FALSE" = FALSE,
      NULL
    )
  }
  if (is.null(directed)) directed <- FALSE
  validate_generator_logical(directed, "directed")

  separator <- if (format == "csv") "," else "\t"
  edge_data <- tryCatch(
    utils::read.table(
      file,
      header = has_header,
      sep = separator,
      comment.char = "#",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = "numeric"
    ),
    error = function(e) {
      stop("Could not read the edge list: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  if (!(ncol(edge_data) %in% c(2L, 3L))) {
    stop("The edge list must contain exactly two or three columns.",
         call. = FALSE)
  }
  if (is.null(weighted)) weighted <- ncol(edge_data) == 3L
  validate_generator_logical(weighted, "weighted")
  expected_columns <- if (weighted) 3L else 2L
  if (ncol(edge_data) != expected_columns) {
    stop(
      if (weighted) {
        "Weighted edge lists must contain a third weight column."
      } else {
        "Unweighted edge lists must contain exactly two columns."
      },
      call. = FALSE
    )
  }

  node_from <- edge_data[[1L]]
  node_to <- edge_data[[2L]]
  weights <- if (weighted) edge_data[[3L]] else rep(1, nrow(edge_data))
  if (anyNA(node_from) ||
      anyNA(node_to) ||
      any(!is.finite(node_from)) ||
      any(!is.finite(node_to)) ||
      any(node_from < 1) ||
      any(node_to < 1) ||
      any(node_from != floor(node_from)) ||
      any(node_to != floor(node_to))) {
    stop("Node indices must be positive finite integers.", call. = FALSE)
  }
  if (anyNA(weights) || any(!is.finite(weights))) {
    stop("Edge weights must be finite and nonmissing.", call. = FALSE)
  }
  node_from <- as.integer(node_from)
  node_to <- as.integer(node_to)

  maximum_index <- if (length(node_from) == 0L) {
    0L
  } else {
    max(node_from, node_to)
  }
  if (is.null(n)) {
    n <- maximum_index
  } else {
    n <- validate_generator_count(n, "n", allow_zero = TRUE)
    if (n < maximum_index) {
      stop("n is smaller than a node index in the edge list.", call. = FALSE)
    }
  }

  duplicate_keys <- if (directed) {
    paste(node_from, node_to, sep = "_")
  } else {
    paste(pmin(node_from, node_to), pmax(node_from, node_to), sep = "_")
  }
  if (anyDuplicated(duplicate_keys)) {
    stop(
      "The edge list contains duplicate node pairs",
      if (!directed) " after ignoring orientation" else "",
      ".",
      call. = FALSE
    )
  }

  off_diagonal <- node_from != node_to
  if (representation == "dense") {
    A <- matrix(0, nrow = n, ncol = n)
    if (directed) {
      if (length(node_from) > 0L) {
        A[cbind(node_from, node_to)] <- weights
      }
      return(A)
    }
    if (any(off_diagonal)) {
      A[cbind(node_from[off_diagonal], node_to[off_diagonal])] <-
        weights[off_diagonal]
    }
    A <- A + t(A)
    if (any(!off_diagonal)) {
      diagonal_nodes <- node_from[!off_diagonal]
      A[cbind(diagonal_nodes, diagonal_nodes)] <- weights[!off_diagonal]
    }
    return(A)
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required for sparse output.", call. = FALSE)
  }
  if (directed) {
    A <- Matrix::sparseMatrix(
      i = node_from,
      j = node_to,
      x = weights,
      dims = c(n, n),
      giveCsparse = TRUE
    )
    return(methods::as(methods::as(A, "generalMatrix"), "dgCMatrix"))
  }
  A <- Matrix::sparseMatrix(
    i = node_from[off_diagonal],
    j = node_to[off_diagonal],
    x = weights[off_diagonal],
    dims = c(n, n),
    giveCsparse = TRUE
  )
  A <- A + Matrix::t(A)
  if (any(!off_diagonal)) {
    diagonal_nodes <- node_from[!off_diagonal]
    A <- A + Matrix::sparseMatrix(
      i = diagonal_nodes,
      j = diagonal_nodes,
      x = weights[!off_diagonal],
      dims = c(n, n),
      giveCsparse = TRUE
    )
  }
  methods::as(methods::as(A, "generalMatrix"), "dgCMatrix")
}

# Validate an n-by-n probability matrix.
validate_probability_matrix <- function(P, directed, tolerance = 1e-10) {
  validate_generator_logical(directed, "directed")
  if (is.null(dim(P)) || length(dim(P)) != 2L || nrow(P) != ncol(P)) {
    stop("P must be a square matrix-like object.", call. = FALSE)
  }
  if (!(is.numeric(P) || inherits(P, "Matrix"))) {
    stop("P must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(P)) || any(P < 0) || any(P > 1)) {
    stop("P must contain finite probabilities in [0, 1].", call. = FALSE)
  }
  if (!directed && !is_symmetric_network_matrix(P, tolerance = tolerance)) {
    stop("P must be symmetric when directed = FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

# Clip a generated probability matrix while remaining package-loading independent.
clip_generator_probabilities <- function(P, lower_clip, upper_clip) {
  if (length(lower_clip) != 1L ||
      length(upper_clip) != 1L ||
      !is.numeric(lower_clip) ||
      !is.numeric(upper_clip) ||
      is.na(lower_clip) ||
      is.na(upper_clip) ||
      !is.finite(lower_clip) ||
      !is.finite(upper_clip) ||
      lower_clip < 0 ||
      upper_clip > 1 ||
      lower_clip > upper_clip) {
    stop(
      "lower_clip and upper_clip must define an interval within [0, 1].",
      call. = FALSE
    )
  }
  pmin(pmax(P, lower_clip), upper_clip)
}

#' Scale probabilities to a target average degree
#'
#' These twin functions apply calibrated bisection or the historical one-step
#' scaling rule to a dense or sparse probability matrix.
#'
#' @param P Dense or sparse edge-probability matrix.
#' @param average_degree Target expected average row degree.
#' @param self_loops Whether diagonal probabilities count toward degree.
#' @param lower_clip,upper_clip Finite probability clipping bounds.
#' @param tolerance Positive bisection convergence tolerance.
#' @param max_iterations Positive maximum number of bisection iterations.
#' @return A list containing the scaled probability matrix `P` and `multiplier`.
#' @examples
#' P <- matrix(0.1, 200, 200)
#' diag(P) <- 0
#' scaled <- scale_to_average_degree(
#'   P, average_degree = 5, self_loops = FALSE,
#'   lower_clip = 0, upper_clip = 1
#' )
#' mean(rowSums(scaled$P))
#' @seealso [apply_average_degree_scaling()], [scale_lsm_to_average_degree()]
#' @rdname scale_to_average_degree
#' @export
scale_to_average_degree <- function(
    P,
    average_degree,
    self_loops,
    lower_clip,
    upper_clip,
    tolerance = 1e-8,
    max_iterations = 100L) {
  n <- nrow(P)
  maximum_degree <- n - as.integer(!self_loops)
  if (length(average_degree) != 1L ||
      !is.numeric(average_degree) ||
      is.na(average_degree) ||
      !is.finite(average_degree) ||
      average_degree < 0 ||
      average_degree > maximum_degree) {
    stop(
      "average_degree must be a finite value between 0 and ",
      maximum_degree,
      ".",
      call. = FALSE
    )
  }

  make_probability <- function(multiplier) {
    P_scaled <- clip_generator_probabilities(
      multiplier * P,
      lower_clip,
      upper_clip
    )
    if (!self_loops && n > 0L) diag(P_scaled) <- 0
    P_scaled
  }
  expected_degree <- function(multiplier) {
    sum(make_probability(multiplier)) / n
  }

  minimum_degree <- expected_degree(0)
  if (average_degree < minimum_degree - tolerance) {
    stop(
      "average_degree is below the minimum implied by lower_clip.",
      call. = FALSE
    )
  }
  if (abs(average_degree - minimum_degree) <= tolerance) {
    return(list(P = make_probability(0), multiplier = 0))
  }
  if (!any(P > 0)) {
    stop("A positive average_degree cannot be obtained from an all-zero P.",
         call. = FALSE)
  }

  lower_multiplier <- 0
  upper_multiplier <- 1
  upper_degree <- expected_degree(upper_multiplier)
  while (upper_degree < average_degree - tolerance &&
         upper_multiplier < 2^50) {
    upper_multiplier <- 2 * upper_multiplier
    upper_degree <- expected_degree(upper_multiplier)
  }
  if (upper_degree < average_degree - tolerance) {
    stop(
      "average_degree exceeds the maximum attainable under P and upper_clip.",
      call. = FALSE
    )
  }

  for (iteration in seq_len(max_iterations)) {
    midpoint <- (lower_multiplier + upper_multiplier) / 2
    midpoint_degree <- expected_degree(midpoint)
    if (abs(midpoint_degree - average_degree) <= tolerance) {
      lower_multiplier <- midpoint
      upper_multiplier <- midpoint
      break
    }
    if (midpoint_degree < average_degree) {
      lower_multiplier <- midpoint
    } else {
      upper_multiplier <- midpoint
    }
  }
  multiplier <- (lower_multiplier + upper_multiplier) / 2
  list(P = make_probability(multiplier), multiplier = multiplier)
}

#' @description `scale_to_average_degree_naive()` applies the historical
#'   one-step ratio before clipping and loop removal.
#' @rdname scale_to_average_degree
#' @export
scale_to_average_degree_naive <- function(
    P,
    average_degree,
    self_loops,
    lower_clip,
    upper_clip) {
  n <- nrow(P)
  maximum_degree <- n - as.integer(!self_loops)
  if (length(average_degree) != 1L ||
      !is.numeric(average_degree) ||
      is.na(average_degree) ||
      !is.finite(average_degree) ||
      average_degree < 0 ||
      average_degree > maximum_degree) {
    stop(
      "average_degree must be a finite value between 0 and ",
      maximum_degree,
      ".",
      call. = FALSE
    )
  }

  baseline_average_degree <- sum(P) / n
  if (average_degree == 0) {
    multiplier <- 0
  } else if (!is.finite(baseline_average_degree) ||
             baseline_average_degree <= 0) {
    stop(
      "A positive average_degree cannot be obtained when mean(rowSums(P)) ",
      "is not positive.",
      call. = FALSE
    )
  } else {
    multiplier <- average_degree / baseline_average_degree
  }
  P <- clip_generator_probabilities(
    multiplier * P,
    lower_clip,
    upper_clip
  )
  if (!self_loops && n > 0L) diag(P) <- 0
  list(P = P, multiplier = multiplier)
}

#' Apply the selected average-degree scaling method
#'
#' Dispatches to the calibrated or historical probability-scaling rule.
#'
#' @param P Dense or sparse edge-probability matrix.
#' @param average_degree Target expected average row degree.
#' @param average_degree_method Scaling method, either `"naive"` or
#'   `"calibrated"`.
#' @param self_loops Whether diagonal probabilities count toward degree.
#' @param lower_clip,upper_clip Finite probability clipping bounds.
#' @return A list containing the scaled probability matrix `P` and `multiplier`.
#' @examples
#' P <- matrix(0.1, 200, 200)
#' diag(P) <- 0
#' apply_average_degree_scaling(
#'   P, average_degree = 5, average_degree_method = "naive",
#'   self_loops = FALSE, lower_clip = 0, upper_clip = 1
#' )$multiplier
#' @seealso [scale_to_average_degree()], [scale_lsm_to_average_degree()]
#' @rdname apply_average_degree_scaling
#' @export
apply_average_degree_scaling <- function(
    P,
    average_degree,
    average_degree_method = "naive",
    self_loops,
    lower_clip,
    upper_clip) {
  if (average_degree_method == "naive") {
    return(scale_to_average_degree_naive(
      P = P,
      average_degree = average_degree,
      self_loops = self_loops,
      lower_clip = lower_clip,
      upper_clip = upper_clip
    ))
  }
  scale_to_average_degree(
    P = P,
    average_degree = average_degree,
    self_loops = self_loops,
    lower_clip = lower_clip,
    upper_clip = upper_clip
  )
}

#' Sample a network from edge probabilities
#'
#' Samples a dense or sparse adjacency matrix from `P`. Explicit row-specific
#' seeds make results reproducible across worker counts.
#'
#' @param P Dense or sparse edge-probability matrix.
#' @param n Optional network size inferred from `P` when omitted.
#' @param directed Whether to sample directed edges independently.
#' @param representation Whether to return a sparse or dense adjacency matrix.
#' @param self_loops Whether diagonal edges may be sampled.
#' @param seed Optional nonnegative reproducibility seed.
#' @param ncores Positive worker count.
#' @param validate_inputs Whether to validate `P` before sampling.
#' @return A base matrix or sparse `dgCMatrix`, according to `representation`.
#' @examples
#' P <- matrix(0.03, 200, 200)
#' diag(P) <- 0
#' A <- generate_adjacency(P, seed = 3, ncores = 1)
#' dim(A)
#' @seealso [generate_er()], [generate_rdpg()], [generate_sbm()]
#' @rdname generate_adjacency
#' @export
generate_adjacency <- function(
    P,
    representation = c("sparse", "dense"),
    directed = FALSE,
    self_loops = FALSE,
    seed = NULL,
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    validate_inputs = TRUE,
    n = NULL) {
  representation <- match.arg(representation)
  validate_generator_logical(directed, "directed")
  validate_generator_logical(self_loops, "self_loops")
  validate_generator_logical(validate_inputs, "validate_inputs")
  ncores <- validate_generator_ncores(ncores)
  seed <- validate_generator_seed(seed)
  scalar_probability <- is.null(dim(P)) && length(P) == 1L
  if (scalar_probability) {
    n <- validate_generator_count(n, "n", allow_zero = TRUE)
    if (!is.numeric(P) || is.na(P) || !is.finite(P) || P < 0 || P > 1) {
      stop("A scalar P must be one finite probability in [0, 1].",
           call. = FALSE)
    }
  } else {
    if (validate_inputs) {
      validate_probability_matrix(P, directed = directed)
    }
    n <- nrow(P)
  }
  if (n == 0L) {
    if (representation == "dense") {
      return(matrix(numeric(), nrow = 0L, ncol = 0L))
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("The Matrix package is required for sparse output.", call. = FALSE)
    }
    return(Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(0L, 0L)))
  }

  if (is.null(seed)) {
    row_seeds <- sample.int(.Machine$integer.max, n, replace = FALSE)
  } else {
    row_seeds <- as.integer(
      (as.double(seed) + seq_len(n)) %% .Machine$integer.max
    )
  }

  sampled_rows <- generator_lapply(
    seq_len(n),
    function(node_idx) {
      set.seed(row_seeds[node_idx])
      target_indices <- if (directed) {
        seq_len(n)
      } else if (self_loops) {
        seq.int(node_idx, n)
      } else if (node_idx < n) {
        seq.int(node_idx + 1L, n)
      } else {
        integer()
      }
      if (!self_loops && directed) {
        target_indices <- target_indices[target_indices != node_idx]
      }
      if (length(target_indices) == 0L) {
        return(matrix(integer(), nrow = 0L, ncol = 2L))
      }
      probabilities <- if (scalar_probability) {
        rep.int(as.numeric(P), length(target_indices))
      } else {
        as.numeric(P[node_idx, target_indices, drop = TRUE])
      }
      selected <- stats::rbinom(
        length(target_indices),
        size = 1L,
        prob = probabilities
      ) == 1L
      if (!any(selected)) {
        return(matrix(integer(), nrow = 0L, ncol = 2L))
      }
      cbind(
        rep.int(node_idx, sum(selected)),
        target_indices[selected]
      )
    },
    ncores = ncores
  )
  edge_indices <- do.call(rbind, sampled_rows)
  if (is.null(edge_indices)) {
    edge_indices <- matrix(integer(), nrow = 0L, ncol = 2L)
  }

  if (representation == "dense") {
    A <- matrix(0, nrow = n, ncol = n)
    if (nrow(edge_indices) > 0L) {
      A[edge_indices] <- 1
    }
    if (!directed) {
      diagonal_values <- diag(A)
      diag(A) <- 0
      A <- A + t(A)
      diag(A) <- diagonal_values
    }
    return(A)
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required for sparse output.", call. = FALSE)
  }
  if (directed) {
    if (nrow(edge_indices) == 0L) {
      return(Matrix::sparseMatrix(
        i = integer(),
        j = integer(),
        dims = c(n, n)
      ))
    }
    return(Matrix::sparseMatrix(
      i = edge_indices[, 1L],
      j = edge_indices[, 2L],
      x = 1,
      dims = c(n, n),
      giveCsparse = TRUE
    ))
  }

  off_diagonal <- if (nrow(edge_indices) == 0L) {
    logical()
  } else {
    edge_indices[, 1L] != edge_indices[, 2L]
  }
  A <- Matrix::sparseMatrix(
    i = edge_indices[off_diagonal, 1L],
    j = edge_indices[off_diagonal, 2L],
    x = rep(1, sum(off_diagonal)),
    dims = c(n, n),
    giveCsparse = TRUE
  )
  A <- A + Matrix::t(A)
  diagonal_indices <- if (nrow(edge_indices) == 0L) {
    integer()
  } else {
    edge_indices[!off_diagonal, 1L]
  }
  if (length(diagonal_indices) > 0L) {
    A <- A + Matrix::sparseMatrix(
      i = diagonal_indices,
      j = diagonal_indices,
      x = rep(1, length(diagonal_indices)),
      dims = c(n, n),
      giveCsparse = TRUE
    )
  }
  methods::as(methods::as(A, "generalMatrix"), "dgCMatrix")
}

#' Generate an Erdos-Renyi network
#'
#' Generates a dense or sparse Erdos-Renyi adjacency matrix and stores its
#' generating parameters as lightweight metadata.
#'
#' @param n Positive number of nodes.
#' @param p Optional edge probability.
#' @param average_degree Optional target expected average degree.
#' @param directed Whether to generate a directed network.
#' @param representation Whether to return a sparse or dense adjacency matrix.
#' @param self_loops Whether diagonal edges may be generated.
#' @param seed Optional nonnegative reproducibility seed.
#' @param ncores Positive worker count.
#' @return A generated adjacency matrix; see [get_generator_parameters()].
#' @examples
#' A <- generate_er(n = 200, average_degree = 5, seed = 1, ncores = 1)
#' get_generator_parameters(A)
#' @seealso [generate_sbm()], [generate_rdpg()]
#' @rdname generate_er
#' @export
generate_er <- function(
    n,
    p = NULL,
    average_degree = NULL,
    representation = c("sparse", "dense"),
    directed = FALSE,
    self_loops = FALSE,
    seed = NULL,
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    )) {
  n <- validate_generator_count(n, "n")
  validate_generator_logical(self_loops, "self_loops")
  if (is.null(p) == is.null(average_degree)) {
    stop("Supply exactly one of p and average_degree.", call. = FALSE)
  }
  maximum_degree <- n - as.integer(!self_loops)
  if (!is.null(average_degree)) {
    if (length(average_degree) != 1L ||
        !is.numeric(average_degree) ||
        is.na(average_degree) ||
        !is.finite(average_degree) ||
        average_degree < 0 ||
        average_degree > maximum_degree) {
      stop("average_degree must be between 0 and ", maximum_degree, ".",
           call. = FALSE)
    }
    p <- if (maximum_degree == 0L) 0 else average_degree / maximum_degree
  } else if (length(p) != 1L ||
             !is.numeric(p) ||
             is.na(p) ||
             !is.finite(p) ||
             p < 0 ||
             p > 1) {
    stop("p must be one finite probability in [0, 1].", call. = FALSE)
  }
  A <- generate_adjacency(
    P = p,
    representation = representation,
    directed = directed,
    self_loops = self_loops,
    seed = seed,
    ncores = ncores,
    validate_inputs = FALSE,
    n = n
  )
  set_generator_parameters(
    A,
    list(model = "er", n = n, p = p,
         average_degree = p * maximum_degree,
         directed = directed, self_loops = self_loops)
  )
}

#' Generate or validate latent positions
#'
#' Returns a finite `n`-by-`d` latent-position matrix, either supplied by the
#' caller or sampled from `latent_distribution`.
#'
#' @param n Positive number of nodes.
#' @param d Positive latent dimension.
#' @param Z Optional supplied latent-position matrix.
#' @param latent_distribution Function used to sample latent coordinates.
#' @param latent_args Optional named list passed to `latent_distribution`.
#' @param seed Optional nonnegative reproducibility seed.
#' @return A numeric matrix with `n` rows and `d` columns.
#' @examples
#' Z <- generate_latent_positions(n = 200, d = 3, seed = 4)
#' dim(Z)
#' @seealso [generate_rdpg()], [generate_lsm_positions()]
#' @rdname generate_latent_positions
#' @export
generate_latent_positions <- function(
    n = NULL,
    d = NULL,
    Z = NULL,
    latent_distribution = stats::runif,
    latent_args = NULL,
    seed = NULL) {
  seed <- validate_generator_seed(seed)
  if (!is.null(Z)) {
    if (is.null(dim(Z)) || length(dim(Z)) != 2L || !is.numeric(Z)) {
      stop("Z must be a numeric matrix.", call. = FALSE)
    }
    Z <- as.matrix(Z)
    if (any(!is.finite(Z))) {
      stop("Z must contain only finite values.", call. = FALSE)
    }
    if (!is.null(n) && nrow(Z) != validate_generator_count(n, "n")) {
      stop("n does not match nrow(Z).", call. = FALSE)
    }
    if (!is.null(d) && ncol(Z) != validate_generator_count(d, "d")) {
      stop("d does not match ncol(Z).", call. = FALSE)
    }
    return(Z)
  }

  if (is.null(n) || is.null(d)) {
    stop("n and d are required when Z is not supplied.", call. = FALSE)
  }
  n <- validate_generator_count(n, "n")
  d <- validate_generator_count(d, "d")
  latent_distribution <- match.fun(latent_distribution)
  if (is.null(latent_args)) {
    latent_args <- if (identical(latent_distribution, stats::runif)) {
      list(min = 0, max = 1 / sqrt(d))
    } else {
      list()
    }
  }
  if (!is.list(latent_args) ||
      (length(latent_args) > 0L && is.null(names(latent_args)))) {
    stop("latent_args must be a named list.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  values <- do.call(
    latent_distribution,
    c(list(n = n * d), latent_args)
  )
  if (!is.numeric(values) || length(values) != n * d || any(!is.finite(values))) {
    stop(
      "latent_distribution must return n * d finite numeric values.",
      call. = FALSE
    )
  }
  matrix(values, nrow = n, ncol = d)
}

#' Generate a random dot product graph
#'
#' Generates a dense or sparse RDPG adjacency matrix. Undirected models use
#' `Z`; directed models use `Z_left` and `Z_right`.
#'
#' @param n Positive number of nodes.
#' @param d Positive latent dimension.
#' @param Z Optional undirected latent-position matrix.
#' @param Z_left,Z_right Optional left and right directed latent-position
#'   matrices.
#' @param latent_distribution Function used to sample latent coordinates.
#' @param latent_args Optional named list passed to `latent_distribution`.
#' @param sparsity_multiplier Nonnegative multiplier applied to dot-product
#'   probabilities.
#' @param scale_P Whether to calibrate dot products before sampling.
#' @param average_degree Optional target expected average degree.
#' @param average_degree_method Scaling method, either `"naive"` or
#'   `"calibrated"`.
#' @param lower_clip,upper_clip Finite probability clipping bounds.
#' @param representation Whether to return a sparse or dense adjacency matrix.
#' @param directed Whether to generate a directed network.
#' @param self_loops Whether diagonal edges may be generated.
#' @param seed Optional nonnegative reproducibility seed.
#' @param ncores Positive worker count.
#' @return A generated adjacency matrix; see [get_generator_parameters()].
#' @examples
#' A <- generate_rdpg(n = 200, d = 3, seed = 2, ncores = 1)
#' get_generator_parameters(A)
#' @seealso [generate_latent_positions()], [ase()], [netcrop_rdpg()]
#' @rdname generate_rdpg
#' @export
generate_rdpg <- function(
    n = NULL,
    d = NULL,
    Z = NULL,
    Z_left = NULL,
    Z_right = NULL,
    latent_distribution = stats::runif,
    latent_args = NULL,
    sparsity_multiplier = 1,
    scale_P = TRUE,
    average_degree = NULL,
    average_degree_method = c("naive", "calibrated"),
    lower_clip = 0,
    upper_clip = 1,
    representation = c("sparse", "dense"),
    directed = FALSE,
    self_loops = FALSE,
    seed = NULL,
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    )) {
  validate_generator_logical(directed, "directed")
  validate_generator_logical(scale_P, "scale_P")
  validate_generator_logical(self_loops, "self_loops")
  average_degree_method <- match.arg(average_degree_method)
  seed <- validate_generator_seed(seed)
  if (length(sparsity_multiplier) != 1L ||
      !is.numeric(sparsity_multiplier) ||
      is.na(sparsity_multiplier) ||
      !is.finite(sparsity_multiplier) ||
      sparsity_multiplier < 0) {
    stop("sparsity_multiplier must be one nonnegative number.", call. = FALSE)
  }
  latent_distribution <- match.fun(latent_distribution)
  if (is.null(latent_args) &&
      identical(latent_distribution, stats::runif)) {
    latent_args <- list(min = 0, max = 1)
  }

  if (!directed) {
    if (!is.null(Z_left) || !is.null(Z_right)) {
      stop("Use Z, not Z_left or Z_right, for an undirected RDPG.",
           call. = FALSE)
    }
    Z <- generate_latent_positions(
      n = n,
      d = d,
      Z = Z,
      latent_distribution = latent_distribution,
      latent_args = latent_args,
      seed = seed
    )
    n <- nrow(Z)
    d <- ncol(Z)
    P <- tcrossprod(Z)
    parameters <- list(Z = Z)
  } else {
    if (!is.null(Z)) {
      stop("Use Z_left and Z_right, not Z, for a directed RDPG.",
           call. = FALSE)
    }
    supplied_Z <- if (!is.null(Z_left)) Z_left else Z_right
    if (!is.null(supplied_Z)) {
      if (is.null(n)) n <- nrow(supplied_Z)
      if (is.null(d)) d <- ncol(supplied_Z)
    }
    Z_left <- generate_latent_positions(
      n = n,
      d = d,
      Z = Z_left,
      latent_distribution = latent_distribution,
      latent_args = latent_args,
      seed = seed
    )
    Z_right <- generate_latent_positions(
      n = nrow(Z_left),
      d = ncol(Z_left),
      Z = Z_right,
      latent_distribution = latent_distribution,
      latent_args = latent_args,
      seed = offset_generator_seed(seed, 1L)
    )
    n <- nrow(Z_left)
    d <- ncol(Z_left)
    P <- tcrossprod(Z_left, Z_right)
    parameters <- list(Z_left = Z_left, Z_right = Z_right)
  }

  if (scale_P) {
    maximum_probability <- max(P)
    if (is.finite(maximum_probability) && maximum_probability > 0) {
      P <- P / maximum_probability
    }
  }
  P <- sparsity_multiplier * P
  if (is.null(average_degree)) {
    P <- clip_generator_probabilities(P, lower_clip, upper_clip)
    if (!self_loops) diag(P) <- 0
    average_degree_multiplier <- 1
  } else {
    degree_scaling <- apply_average_degree_scaling(
      P = P,
      average_degree = average_degree,
      average_degree_method = average_degree_method,
      self_loops = self_loops,
      lower_clip = lower_clip,
      upper_clip = upper_clip
    )
    P <- degree_scaling$P
    average_degree_multiplier <- degree_scaling$multiplier
  }
  if (anyNA(P)) {
    stop("Generated P contains missing or NaN probabilities.",
         call. = FALSE)
  }
  A <- generate_adjacency(
    P = P,
    representation = representation,
    directed = directed,
    self_loops = self_loops,
    seed = offset_generator_seed(seed, 2L),
    ncores = ncores,
    validate_inputs = FALSE
  )
  parameters <- c(
    list(
      model = if (directed) "generalized_rdpg" else "rdpg",
      n = n,
      d = d,
      directed = directed,
      self_loops = self_loops,
      sparsity_multiplier = sparsity_multiplier,
      scale_P = scale_P,
      average_degree_target = average_degree,
      average_degree = sum(P) / n,
      average_degree_multiplier = average_degree_multiplier,
      average_degree_method = average_degree_method,
      lower_clip = lower_clip,
      upper_clip = upper_clip
    ),
    parameters
  )
  set_generator_parameters(A, parameters)
}

#' Generate or validate community labels
#'
#' Produces community memberships for an SBM, DCBM, or latent-space model.
#'
#' @param n Positive number of nodes.
#' @param K Positive number of communities.
#' @param g_true Optional supplied community-label vector.
#' @param community_probabilities Optional length-`K` community probabilities.
#' @param seed Optional nonnegative reproducibility seed.
#' @return An integer vector of length `n` with values from `1` through `K`.
#' @examples
#' g <- generate_community_labels(n = 200, K = 3, seed = 5)
#' table(g)
#' @seealso [generate_sbm()], [generate_dcbm()], [generate_lsm_positions()]
#' @rdname generate_community_labels
#' @export
generate_community_labels <- function(
    n = NULL,
    K = NULL,
    g_true = NULL,
    community_probabilities = NULL,
    seed = NULL) {
  seed <- validate_generator_seed(seed)
  if (!is.null(g_true)) {
    if (length(g_true) == 0L || anyNA(g_true)) {
      stop("g_true must contain at least one nonmissing label.", call. = FALSE)
    }
    if (!is.null(n) && length(g_true) != validate_generator_count(n, "n")) {
      stop("n does not match length(g_true).", call. = FALSE)
    }
    label_levels <- unique(g_true)
    g_integer <- match(g_true, label_levels)
    inferred_K <- length(label_levels)
    if (!is.null(K) && validate_generator_count(K, "K") != inferred_K) {
      stop("K does not match the number of unique values in g_true.",
           call. = FALSE)
    }
    attr(g_integer, "label_levels") <- label_levels
    return(g_integer)
  }

  if (is.null(n) || is.null(K)) {
    stop("n and K are required when g_true is not supplied.", call. = FALSE)
  }
  n <- validate_generator_count(n, "n")
  K <- validate_generator_count(K, "K")
  # Match the historical SBM/DCBM default: sample labels uniformly when no
  # community proportions are supplied.
  if (is.null(community_probabilities)) {
    community_probabilities <- rep(1 / K, K)
  }
  if (!is.numeric(community_probabilities) ||
      length(community_probabilities) != K ||
      anyNA(community_probabilities) ||
      any(!is.finite(community_probabilities)) ||
      any(community_probabilities < 0) ||
      sum(community_probabilities) <= 0) {
    stop(
      "community_probabilities must be K nonnegative finite values with positive sum.",
      call. = FALSE
    )
  }
  community_probabilities <- community_probabilities /
    sum(community_probabilities)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  sample.int(n = K, size = n, replace = TRUE,
             prob = community_probabilities)
}

#' Generate inverse-beta degree parameters
#'
#' Draws the historical default degree parameters used by the DCBM generator.
#'
#' @param n Nonnegative number of degree parameters.
#' @param shape_1,shape_2 Positive beta-distribution shape parameters.
#' @return A nonnegative numeric vector of length `n`.
#' @examples
#' psi <- generate_inverse_beta_degree_parameters(200)
#' length(psi)
#' @seealso [generate_degree_parameters()], [generate_dcbm()]
#' @rdname generate_inverse_beta_degree_parameters
#' @export
generate_inverse_beta_degree_parameters <- function(
    n,
    shape_1 = 4,
    shape_2 = 1) {
  n <- validate_generator_count(n, "n", allow_zero = TRUE)
  shape_parameters <- c(shape_1 = shape_1, shape_2 = shape_2)
  if (anyNA(shape_parameters) ||
      any(!is.finite(shape_parameters)) ||
      any(shape_parameters <= 0)) {
    stop("shape_1 and shape_2 must be positive finite values.",
         call. = FALSE)
  }
  1 / stats::rbeta(n, shape1 = shape_1, shape2 = shape_2)
}

#' Generate or validate DCBM degree parameters
#'
#' Returns nonnegative node degree parameters, optionally normalized within
#' communities.
#'
#' @param n Positive number of nodes.
#' @param psi Optional supplied nonnegative degree parameters.
#' @param degree_distribution Function used to generate degree parameters.
#' @param degree_args Optional named list passed to `degree_distribution`.
#' @param degree_scale Community-wise degree-normalization method.
#' @param g_true Optional community-label vector.
#' @param seed Optional nonnegative reproducibility seed.
#' @return A nonnegative numeric vector of length `n`.
#' @examples
#' g <- generate_community_labels(n = 200, K = 3, seed = 6)
#' psi <- generate_degree_parameters(n = 200, g_true = g, seed = 7)
#' length(psi)
#' @seealso [generate_inverse_beta_degree_parameters()], [generate_dcbm()]
#' @rdname generate_degree_parameters
#' @export
generate_degree_parameters <- function(
    n = NULL,
    psi = NULL,
    degree_distribution = generate_inverse_beta_degree_parameters,
    degree_args = NULL,
    degree_scale = c("none", "max_by_community", "mean_one_by_community"),
    g_true = NULL,
    seed = NULL) {
  degree_scale <- match.arg(degree_scale)
  seed <- validate_generator_seed(seed)
  if (!is.null(psi)) {
    if (is.null(n)) n <- length(psi)
    n <- validate_generator_count(n, "n")
    if (length(psi) == 1L) psi <- rep(psi, n)
    if (!is.numeric(psi) ||
        length(psi) != n ||
        anyNA(psi) ||
        any(!is.finite(psi)) ||
        any(psi < 0)) {
      stop("psi must contain n nonnegative finite values.", call. = FALSE)
    }
  } else {
    if (is.null(n)) {
      stop("n is required when psi is not supplied.", call. = FALSE)
    }
    n <- validate_generator_count(n, "n")
    degree_distribution <- match.fun(degree_distribution)
    if (is.null(degree_args)) {
      degree_args <- if (identical(degree_distribution, stats::runif)) {
        list(min = 0.5, max = 1)
      } else {
        list()
      }
    }
    if (!is.list(degree_args) ||
        (length(degree_args) > 0L && is.null(names(degree_args)))) {
      stop("degree_args must be a named list.", call. = FALSE)
    }
    if (!is.null(seed)) set.seed(seed)
    psi <- do.call(degree_distribution, c(list(n = n), degree_args))
    if (!is.numeric(psi) ||
        length(psi) != n ||
        anyNA(psi) ||
        any(!is.finite(psi)) ||
        any(psi < 0)) {
      stop(
        "degree_distribution must return n nonnegative finite values.",
        call. = FALSE
      )
    }
  }

  if (degree_scale != "none") {
    if (is.null(g_true) || length(g_true) != n) {
      stop("g_true with length n is required for community-wise scaling.",
           call. = FALSE)
    }
    for (community in unique(g_true)) {
      community_mask <- g_true == community
      denominator <- if (degree_scale == "max_by_community") {
        max(psi[community_mask])
      } else {
        mean(psi[community_mask])
      }
      if (!is.finite(denominator) || denominator <= 0) {
        stop("Cannot scale a community whose degree parameters are all zero.",
             call. = FALSE)
      }
      psi[community_mask] <- psi[community_mask] / denominator
    }
  }
  unname(psi)
}

#' Generate block-model networks
#'
#' `generate_dcbm()` generates a degree-corrected stochastic block model;
#' `generate_sbm()` is its ordinary-block-model twin with unit degree effects.
#' Both return dense or sparse adjacency matrices.
#'
#' @param n Positive number of nodes.
#' @param K Positive number of communities.
#' @param g_true Optional supplied community-label vector.
#' @param community_probabilities Optional length-`K` community probabilities.
#' @param P_block Optional `K`-by-`K` block-probability matrix.
#' @param alpha,beta Within- and between-community probabilities used when
#'   `P_block` is omitted.
#' @param psi Optional supplied DCBM degree parameters.
#' @param degree_distribution Function used to generate DCBM degree parameters.
#' @param degree_args Optional named list passed to `degree_distribution`.
#' @param degree_scale Community-wise degree-normalization method.
#' @param sparsity_multiplier Nonnegative multiplier applied to probabilities.
#' @param average_degree Optional target expected average degree.
#' @param average_degree_method Scaling method, either `"naive"` or
#'   `"calibrated"`.
#' @param lower_clip,upper_clip Finite probability clipping bounds.
#' @param representation Whether to return a sparse or dense adjacency matrix.
#' @param directed Whether to generate a directed network.
#' @param self_loops Whether diagonal edges may be generated.
#' @param seed Optional nonnegative reproducibility seed.
#' @param ncores Positive worker count.
#' @return A generated adjacency matrix; see [get_generator_parameters()].
#' @examples
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.4, beta = 0.1,
#'                   seed = 3, ncores = 1)
#' get_generator_parameters(A)
#' @seealso [generate_community_labels()], [estimate_sbm()], [estimate_dcbm()]
#' @rdname generate_blockmodel
#' @export
generate_dcbm <- function(
    n = NULL,
    K = NULL,
    g_true = NULL,
    community_probabilities = NULL,
    P_block = NULL,
    alpha = 0.2,
    beta = 0.2,
    psi = NULL,
    degree_distribution = generate_inverse_beta_degree_parameters,
    degree_args = NULL,
    degree_scale = c("none", "max_by_community", "mean_one_by_community"),
    sparsity_multiplier = 1,
    average_degree = NULL,
    average_degree_method = c("naive", "calibrated"),
    lower_clip = 0,
    upper_clip = 1,
    representation = c("sparse", "dense"),
    directed = FALSE,
    self_loops = FALSE,
    seed = NULL,
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    )) {
  validate_generator_logical(directed, "directed")
  validate_generator_logical(self_loops, "self_loops")
  degree_scale <- match.arg(degree_scale)
  average_degree_method <- match.arg(average_degree_method)
  seed <- validate_generator_seed(seed)
  if (is.null(n) && !is.null(g_true)) n <- length(g_true)
  g_true <- generate_community_labels(
    n = n,
    K = K,
    g_true = g_true,
    community_probabilities = community_probabilities,
    seed = seed
  )
  n <- length(g_true)
  K <- if (is.null(K)) {
    max(g_true)
  } else {
    validate_generator_count(K, "K")
  }

  if (is.null(P_block)) {
    scalar_parameters <- c(alpha = alpha, beta = beta)
    if (anyNA(scalar_parameters) ||
        any(!is.finite(scalar_parameters)) ||
        alpha < 0 ||
        alpha > 1 ||
        beta < 0) {
      stop("alpha must be in [0, 1] and beta must be nonnegative.",
           call. = FALSE)
    }
    P_block <- matrix(alpha * beta, nrow = K, ncol = K)
    diag(P_block) <- alpha
  } else {
    if (is.null(dim(P_block)) ||
        !is.numeric(P_block) ||
        !identical(dim(P_block), c(K, K)) ||
        any(!is.finite(P_block)) ||
        any(P_block < 0) ||
        any(P_block > 1)) {
      stop("P_block must be a finite K-by-K probability matrix.",
           call. = FALSE)
    }
    P_block <- as.matrix(P_block)
  }
  if (!directed &&
      !isTRUE(all.equal(P_block, t(P_block), tolerance = 1e-10))) {
    stop("P_block must be symmetric when directed = FALSE.", call. = FALSE)
  }
  if (length(sparsity_multiplier) != 1L ||
      !is.numeric(sparsity_multiplier) ||
      is.na(sparsity_multiplier) ||
      !is.finite(sparsity_multiplier) ||
      sparsity_multiplier < 0) {
    stop("sparsity_multiplier must be one nonnegative number.", call. = FALSE)
  }

  psi <- generate_degree_parameters(
    n = n,
    psi = psi,
    degree_distribution = degree_distribution,
    degree_args = degree_args,
    degree_scale = degree_scale,
    g_true = g_true,
    seed = offset_generator_seed(seed, 1L)
  )
  P <- P_block[g_true, g_true, drop = FALSE]
  P <- P * psi
  P <- sweep(P, 2L, psi, `*`)
  P <- sparsity_multiplier * P
  if (is.null(average_degree)) {
    P <- clip_generator_probabilities(P, lower_clip, upper_clip)
    if (!self_loops) diag(P) <- 0
    average_degree_multiplier <- 1
  } else {
    degree_scaling <- apply_average_degree_scaling(
      P = P,
      average_degree = average_degree,
      average_degree_method = average_degree_method,
      self_loops = self_loops,
      lower_clip = lower_clip,
      upper_clip = upper_clip
    )
    P <- degree_scaling$P
    average_degree_multiplier <- degree_scaling$multiplier
  }

  if (anyNA(P)) {
    stop("Generated P contains missing or NaN probabilities.",
         call. = FALSE)
  }

  A <- generate_adjacency(
    P = P,
    representation = representation,
    directed = directed,
    self_loops = self_loops,
    seed = offset_generator_seed(seed, 2L),
    ncores = ncores,
    validate_inputs = FALSE
  )
  set_generator_parameters(
    A,
    list(
      model = "dcbm",
      n = n,
      K = K,
      g_true = g_true,
      psi = psi,
      P_block = P_block,
      directed = directed,
      self_loops = self_loops,
      sparsity_multiplier = sparsity_multiplier,
      average_degree_target = average_degree,
      average_degree = sum(P) / n,
      average_degree_multiplier = average_degree_multiplier,
      average_degree_method = average_degree_method,
      degree_scale = degree_scale
    )
  )
}

#' @description `generate_sbm()` fixes every degree parameter to one.
#' @rdname generate_blockmodel
#' @export
generate_sbm <- function(
    n = NULL,
    K = NULL,
    g_true = NULL,
    community_probabilities = NULL,
    P_block = NULL,
    alpha = 0.2,
    beta = 0.2,
    sparsity_multiplier = 1,
    average_degree = NULL,
    average_degree_method = c("naive", "calibrated"),
    lower_clip = 0,
    upper_clip = 1,
    representation = c("sparse", "dense"),
    directed = FALSE,
    self_loops = FALSE,
    seed = NULL,
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    )) {
  if (is.null(n) && !is.null(g_true)) n <- length(g_true)
  if (is.null(n)) {
    stop("n is required when g_true is not supplied.", call. = FALSE)
  }
  A <- generate_dcbm(
    n = n,
    K = K,
    g_true = g_true,
    community_probabilities = community_probabilities,
    P_block = P_block,
    alpha = alpha,
    beta = beta,
    psi = rep(1, n),
    degree_scale = "none",
    sparsity_multiplier = sparsity_multiplier,
    average_degree = average_degree,
    average_degree_method = average_degree_method,
    lower_clip = lower_clip,
    upper_clip = upper_clip,
    representation = representation,
    directed = directed,
    self_loops = self_loops,
    seed = seed,
    ncores = ncores
  )
  parameters <- get_generator_parameters(A)
  parameters$model <- "sbm"
  set_generator_parameters(A, parameters)
}

#' Generate truncated standard-normal values
#'
#' Draws standard-normal noise on a finite interval by inverse-CDF sampling.
#'
#' @param n Positive number of values to generate.
#' @param lower_bound,upper_bound Finite truncation bounds with
#'   `lower_bound < upper_bound`.
#' @return A numeric vector of length `n`.
#' @examples
#' x <- generate_truncated_normal(200, lower_bound = -2, upper_bound = 2)
#' range(x)
#' @seealso [generate_lsm_positions()]
#' @rdname generate_truncated_normal
#' @export
generate_truncated_normal <- function(
    n,
    lower_bound = -2,
    upper_bound = 2) {
  n <- validate_generator_count(n, "n")
  if (length(lower_bound) != 1L ||
      length(upper_bound) != 1L ||
      !is.numeric(lower_bound) ||
      !is.numeric(upper_bound) ||
      is.na(lower_bound) ||
      is.na(upper_bound) ||
      !is.finite(lower_bound) ||
      !is.finite(upper_bound) ||
      lower_bound >= upper_bound) {
    stop("lower_bound and upper_bound must be finite with lower_bound < upper_bound.",
         call. = FALSE)
  }
  lower_probability <- stats::pnorm(lower_bound)
  upper_probability <- stats::pnorm(upper_bound)
  stats::qnorm(stats::runif(n, lower_probability, upper_probability))
}

#' Normalize latent-space positions
#'
#' Centers and scales a latent-position matrix using the LSM convention.
#'
#' @param Z A finite numeric latent-position matrix.
#' @return A centered and normalized matrix with the same dimensions as `Z`.
#' @examples
#' Z <- matrix(stats::rnorm(200 * 3), 200, 3)
#' dim(normalize_lsm_positions(Z))
#' @seealso [generate_lsm_positions()], [generate_lsm()]
#' @rdname normalize_lsm_positions
#' @export
normalize_lsm_positions <- function(Z) {
  Z <- sweep(Z, 2L, colMeans(Z), FUN = "-")
  G_small <- crossprod(Z)
  normalization <- sqrt(sqrt(sum(G_small^2)) / nrow(Z))
  if (!is.finite(normalization) || normalization <= 0) {
    stop("Latent positions cannot be normalized because they have zero scale.",
         call. = FALSE)
  }
  Z / normalization
}

#' Generate latent-space-model positions
#'
#' Generates an `n`-by-`d` latent-position matrix from a finite Gaussian
#' mixture with `K` communities.
#'
#' @param n Positive number of nodes.
#' @param d Positive latent dimension.
#' @param K Positive number of communities.
#' @param g_true Community-label vector of length `n`.
#' @param mean_lower,mean_upper Finite bounds for community means.
#' @param noise_lower,noise_upper Finite truncation bounds for latent noise.
#' @param seed Optional nonnegative reproducibility seed.
#' @return A numeric matrix with `n` rows and `d` columns.
#' @examples
#' g <- generate_community_labels(n = 200, K = 3, seed = 4)
#' Z <- generate_lsm_positions(n = 200, d = 3, K = 3, g_true = g, seed = 5)
#' dim(Z)
#' @seealso [generate_lsm()], [generate_lsm_alpha()], [normalize_lsm_positions()]
#' @rdname generate_lsm_positions
#' @export
generate_lsm_positions <- function(
    n,
    d,
    K,
    g_true,
    mean_lower = -1,
    mean_upper = 1,
    noise_lower = -2,
    noise_upper = 2,
    seed = NULL) {
  n <- validate_generator_count(n, "n")
  d <- validate_generator_count(d, "d")
  K <- validate_generator_count(K, "K")
  if (length(g_true) != n || anyNA(g_true) ||
      any(g_true < 1L) || any(g_true > K)) {
    stop("g_true must contain n labels in 1, ..., K.", call. = FALSE)
  }
  if (length(mean_lower) != 1L ||
      length(mean_upper) != 1L ||
      !is.numeric(mean_lower) ||
      !is.numeric(mean_upper) ||
      !is.finite(mean_lower) ||
      !is.finite(mean_upper) ||
      mean_lower >= mean_upper) {
    stop("mean_lower and mean_upper must be finite with mean_lower < mean_upper.",
         call. = FALSE)
  }
  seed <- validate_generator_seed(seed)
  if (!is.null(seed)) set.seed(seed)
  community_means <- matrix(
    stats::runif(K * d, mean_lower, mean_upper),
    nrow = K,
    ncol = d
  )
  noise <- matrix(
    generate_truncated_normal(n * d, noise_lower, noise_upper),
    nrow = n,
    ncol = d
  )
  community_means[g_true, , drop = FALSE] + noise
}

#' Generate latent-space-model intercepts
#'
#' Validates or generates the node intercept used by the latent-space model.
#'
#' @param n Positive number of nodes.
#' @param alpha Optional scalar or length-`n` intercept vector.
#' @param seed Optional nonnegative reproducibility seed.
#' @return A finite numeric vector of length `n`.
#' @examples
#' alpha <- generate_lsm_alpha(n = 200, seed = 6)
#' length(alpha)
#' @seealso [generate_lsm()], [generate_lsm_positions()]
#' @rdname generate_lsm_alpha
#' @export
generate_lsm_alpha <- function(n, alpha = NULL, seed = NULL) {
  n <- validate_generator_count(n, "n")
  seed <- validate_generator_seed(seed)
  if (is.null(alpha)) {
    if (!is.null(seed)) set.seed(seed)
    return(-stats::runif(n, 1, 3) / 2)
  }
  if (length(alpha) == 1L) alpha <- rep(alpha, n)
  if (!is.numeric(alpha) ||
      length(alpha) != n ||
      anyNA(alpha) ||
      any(!is.finite(alpha))) {
    stop("alpha must be NULL, one finite scalar, or n finite values.",
         call. = FALSE)
  }
  unname(alpha)
}

#' Scale latent-space logits to a target average degree
#'
#' Shifts an LSM logit matrix using the historical iterative method or calibrated
#' bisection.
#'
#' @param theta Finite square latent-space logit matrix.
#' @param average_degree Target expected average degree.
#' @param average_degree_method Scaling method, either `"naive"` or
#'   `"calibrated"`.
#' @param self_loops Whether diagonal probabilities count toward degree.
#' @param lower_clip,upper_clip Finite probability clipping bounds.
#' @param naive_iterations Positive iteration count for historical scaling.
#' @param tolerance Positive bisection convergence tolerance.
#' @param max_iterations Positive maximum number of bisection iterations.
#' @return A list containing the probability matrix `P` and `intercept_shift`.
#' @examples
#' Z <- matrix(stats::rnorm(200 * 3), 200, 3)
#' theta <- -tcrossprod(Z)
#' scaled <- scale_lsm_to_average_degree(theta, average_degree = 5)
#' mean(rowSums(scaled$P))
#' @seealso [generate_lsm()], [scale_to_average_degree()]
#' @rdname scale_lsm_to_average_degree
#' @export
scale_lsm_to_average_degree <- function(
    theta,
    average_degree,
    average_degree_method = c("naive", "calibrated"),
    self_loops = FALSE,
    lower_clip = 0,
    upper_clip = 1,
    naive_iterations = 10L,
    tolerance = 1e-8,
    max_iterations = 100L) {
  average_degree_method <- match.arg(average_degree_method)
  n <- nrow(theta)
  maximum_degree <- n - as.integer(!self_loops)
  if (length(average_degree) != 1L ||
      !is.numeric(average_degree) ||
      is.na(average_degree) ||
      !is.finite(average_degree) ||
      average_degree < 0 ||
      average_degree > maximum_degree) {
    stop("average_degree must be a finite value between 0 and ",
         maximum_degree, ".", call. = FALSE)
  }
  naive_iterations <- validate_generator_count(
    naive_iterations,
    "naive_iterations"
  )

  make_probability <- function(intercept_shift, remove_diagonal = TRUE) {
    P <- stats::plogis(theta + intercept_shift)
    P <- clip_generator_probabilities(P, lower_clip, upper_clip)
    if (remove_diagonal && !self_loops && n > 0L) diag(P) <- 0
    P
  }

  if (average_degree_method == "naive") {
    if (average_degree == 0) {
      P <- matrix(0, nrow = n, ncol = n)
      return(list(P = P, intercept_shift = -Inf))
    }
    intercept_shift <- 0
    for (iteration in seq_len(naive_iterations)) {
      P_iteration <- make_probability(
        intercept_shift,
        remove_diagonal = FALSE
      )
      current_average_degree <- sum(P_iteration) / n
      if (!is.finite(current_average_degree) || current_average_degree <= 0) {
        stop("Naive LSM degree scaling encountered a nonpositive degree.",
             call. = FALSE)
      }
      intercept_shift <- intercept_shift -
        log(current_average_degree / average_degree)
    }
    return(list(
      P = make_probability(intercept_shift, remove_diagonal = TRUE),
      intercept_shift = intercept_shift
    ))
  }

  expected_degree <- function(intercept_shift) {
    sum(make_probability(intercept_shift, remove_diagonal = TRUE)) / n
  }
  minimum_degree <- expected_degree(-Inf)
  maximum_attainable_degree <- expected_degree(Inf)
  if (average_degree < minimum_degree - tolerance ||
      average_degree > maximum_attainable_degree + tolerance) {
    stop(
      "average_degree is unattainable under the requested probability clips.",
      call. = FALSE
    )
  }
  if (abs(average_degree - minimum_degree) <= tolerance) {
    return(list(P = make_probability(-Inf), intercept_shift = -Inf))
  }
  if (abs(average_degree - maximum_attainable_degree) <= tolerance) {
    return(list(P = make_probability(Inf), intercept_shift = Inf))
  }

  lower_shift <- -1
  upper_shift <- 1
  while (expected_degree(lower_shift) > average_degree) {
    lower_shift <- 2 * lower_shift
  }
  while (expected_degree(upper_shift) < average_degree) {
    upper_shift <- 2 * upper_shift
  }
  for (iteration in seq_len(max_iterations)) {
    midpoint <- (lower_shift + upper_shift) / 2
    midpoint_degree <- expected_degree(midpoint)
    if (abs(midpoint_degree - average_degree) <= tolerance) {
      lower_shift <- midpoint
      upper_shift <- midpoint
      break
    }
    if (midpoint_degree < average_degree) {
      lower_shift <- midpoint
    } else {
      upper_shift <- midpoint
    }
  }
  intercept_shift <- (lower_shift + upper_shift) / 2
  list(
    P = make_probability(intercept_shift),
    intercept_shift = intercept_shift
  )
}

#' Generate a latent-space-model network
#'
#' Generates a dense or sparse Hoff/Ma-style latent-space-model adjacency
#' matrix from supplied or generated positions, memberships, and intercepts.
#'
#' @param n Positive number of nodes.
#' @param d Positive latent dimension.
#' @param K Positive number of communities.
#' @param g_true Optional supplied community-label vector.
#' @param community_probabilities Optional length-`K` community probabilities.
#' @param alpha Optional scalar or node-specific intercepts.
#' @param Z Optional undirected latent-position matrix.
#' @param Z_left,Z_right Optional left and right directed latent-position
#'   matrices.
#' @param normalize_Z Whether to center and scale latent positions.
#' @param distance_adjustment Whether to apply the historical distance scaling.
#' @param mean_lower,mean_upper Finite bounds for generated community means.
#' @param noise_lower,noise_upper Finite truncation bounds for generated latent
#'   noise.
#' @param average_degree Optional target expected average degree.
#' @param average_degree_method Scaling method, either `"naive"` or
#'   `"calibrated"`.
#' @param naive_iterations Positive iteration count for historical scaling.
#' @param lower_clip,upper_clip Finite probability clipping bounds.
#' @param representation Whether to return a sparse or dense adjacency matrix.
#' @param directed Whether to generate a directed network.
#' @param self_loops Whether diagonal edges may be generated.
#' @param seed Optional nonnegative reproducibility seed.
#' @param ncores Positive worker count.
#' @return A generated adjacency matrix; see [get_generator_parameters()].
#' @examples
#' A <- generate_lsm(n = 200, d = 3, K = 3, seed = 7, ncores = 1)
#' get_generator_parameters(A)
#' @seealso [generate_lsm_positions()], [generate_lsm_alpha()], [lsm_pgd()],
#'   [netcrop_lsm()]
#' @rdname generate_lsm
#' @export
generate_lsm <- function(
    n = NULL,
    d = NULL,
    K = 1L,
    g_true = NULL,
    community_probabilities = NULL,
    alpha = NULL,
    Z = NULL,
    Z_left = NULL,
    Z_right = NULL,
    normalize_Z = TRUE,
    distance_adjustment = TRUE,
    mean_lower = -1,
    mean_upper = 1,
    noise_lower = -2,
    noise_upper = 2,
    average_degree = NULL,
    average_degree_method = c("naive", "calibrated"),
    naive_iterations = 10L,
    lower_clip = 0,
    upper_clip = 1,
    representation = c("sparse", "dense"),
    directed = FALSE,
    self_loops = FALSE,
    seed = NULL,
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    )) {
  validate_generator_logical(directed, "directed")
  validate_generator_logical(normalize_Z, "normalize_Z")
  validate_generator_logical(distance_adjustment, "distance_adjustment")
  validate_generator_logical(self_loops, "self_loops")
  average_degree_method <- match.arg(average_degree_method)
  seed <- validate_generator_seed(seed)
  if (!is.null(g_true) && missing(K)) {
    K <- length(unique(g_true))
  }

  supplied_Z <- if (!is.null(Z)) {
    Z
  } else if (!is.null(Z_left)) {
    Z_left
  } else {
    Z_right
  }
  if (!is.null(supplied_Z)) {
    if (is.null(n)) n <- nrow(supplied_Z)
    if (is.null(d)) d <- ncol(supplied_Z)
  }
  if (is.null(n) || is.null(d)) {
    stop("n and d are required when no latent-position matrix is supplied.",
         call. = FALSE)
  }
  n <- validate_generator_count(n, "n")
  d <- validate_generator_count(d, "d")
  K <- validate_generator_count(K, "K")

  needs_generated_positions <- if (directed) {
    is.null(Z_left) || is.null(Z_right)
  } else {
    is.null(Z)
  }
  if (needs_generated_positions || !is.null(g_true)) {
    g_true <- generate_community_labels(
      n = n,
      K = K,
      g_true = g_true,
      community_probabilities = community_probabilities,
      seed = seed
    )
  }

  if (!directed) {
    if (!is.null(Z_left) || !is.null(Z_right)) {
      stop("Use Z, not Z_left or Z_right, for an undirected LSM.",
           call. = FALSE)
    }
    if (is.null(Z)) {
      Z <- generate_lsm_positions(
        n = n,
        d = d,
        K = K,
        g_true = g_true,
        mean_lower = mean_lower,
        mean_upper = mean_upper,
        noise_lower = noise_lower,
        noise_upper = noise_upper,
        seed = offset_generator_seed(seed, 1L)
      )
    } else {
      Z <- generate_latent_positions(n = n, d = d, Z = Z)
    }
    if (normalize_Z) Z <- normalize_lsm_positions(Z)
    Z_left <- Z
    Z_right <- Z
  } else {
    if (!is.null(Z)) {
      stop("Use Z_left and Z_right, not Z, for a directed LSM.",
           call. = FALSE)
    }
    if (is.null(Z_left)) {
      Z_left <- generate_lsm_positions(
        n = n,
        d = d,
        K = K,
        g_true = g_true,
        mean_lower = mean_lower,
        mean_upper = mean_upper,
        noise_lower = noise_lower,
        noise_upper = noise_upper,
        seed = offset_generator_seed(seed, 1L)
      )
    } else {
      Z_left <- generate_latent_positions(n = n, d = d, Z = Z_left)
    }
    if (is.null(Z_right)) {
      Z_right <- generate_lsm_positions(
        n = n,
        d = d,
        K = K,
        g_true = g_true,
        mean_lower = mean_lower,
        mean_upper = mean_upper,
        noise_lower = noise_lower,
        noise_upper = noise_upper,
        seed = offset_generator_seed(seed, 2L)
      )
    } else {
      Z_right <- generate_latent_positions(n = n, d = d, Z = Z_right)
    }
    if (normalize_Z) {
      Z_left <- normalize_lsm_positions(Z_left)
      Z_right <- normalize_lsm_positions(Z_right)
    }
  }

  alpha <- generate_lsm_alpha(
    n = n,
    alpha = alpha,
    seed = offset_generator_seed(seed, 3L)
  )
  alpha_left <- alpha / 2
  alpha_right <- alpha / 2
  if (distance_adjustment) {
    # alpha_left <- alpha_left - sqrt(rowSums(Z_left^2)) / 2
    # alpha_right <- alpha_right - sqrt(rowSums(Z_right^2)) / 2
    alpha_left <- alpha_left - rowSums(Z_left^2)
    alpha_right <- alpha_right - rowSums(Z_right^2)
  }
  # theta <- tcrossprod(Z_left, Z_right)
  theta <- 2 * tcrossprod(Z_left, Z_right)
  theta <- sweep(theta, 1L, alpha_left, `+`)
  theta <- sweep(theta, 2L, alpha_right, `+`)

  if (is.null(average_degree)) {
    P <- clip_generator_probabilities(
      stats::plogis(theta),
      lower_clip,
      upper_clip
    )
    if (!self_loops) diag(P) <- 0
    intercept_shift <- 0
  } else {
    degree_scaling <- scale_lsm_to_average_degree(
      theta = theta,
      average_degree = average_degree,
      average_degree_method = average_degree_method,
      self_loops = self_loops,
      lower_clip = lower_clip,
      upper_clip = upper_clip,
      naive_iterations = naive_iterations
    )
    P <- degree_scaling$P
    intercept_shift <- degree_scaling$intercept_shift
  }
  if (anyNA(P)) {
    stop("Generated P contains missing or NaN probabilities.",
         call. = FALSE)
  }
  A <- generate_adjacency(
    P = P,
    representation = representation,
    directed = directed,
    self_loops = self_loops,
    seed = offset_generator_seed(seed, 4L),
    ncores = ncores,
    validate_inputs = FALSE
  )

  parameters <- list(
    model = "lsm",
    n = n,
    d = d,
    K = K,
    g_true = g_true,
    alpha = alpha,
    alpha_left = alpha_left + intercept_shift / 2,
    alpha_right = alpha_right + intercept_shift / 2,
    average_degree_target = average_degree,
    average_degree = sum(P) / n,
    average_degree_method = average_degree_method,
    average_degree_intercept_shift = intercept_shift,
    normalize_Z = normalize_Z,
    distance_adjustment = distance_adjustment,
    directed = directed,
    self_loops = self_loops
  )
  if (directed) {
    parameters$Z_left <- Z_left
    parameters$Z_right <- Z_right
  } else {
    parameters$Z <- Z
  }
  set_generator_parameters(A, parameters)
}
