# Network utility functions
#
# Compiled accelerators are registered when the package is installed.
#
# All R-facing functions have dependency-light pure-R fallbacks and do not use
# igraph. A finite nonzero adjacency entry represents an edge. Connected
# components use the undirected union of A and t(A). Shortest paths may preserve
# direction and may use either unit edge lengths or nonnegative matrix weights.

# Validate and normalize a dense or Matrix-package adjacency matrix.
validate_adjacency <- function(A) {
  if (is.null(dim(A)) || length(dim(A)) != 2L || nrow(A) != ncol(A)) {
    stop("A must be a square matrix-like object.", call. = FALSE)
  }
  if (!(is.numeric(A) || inherits(A, "Matrix"))) {
    stop("A must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(A))) {
    stop("A must contain only finite values.", call. = FALSE)
  }
  if (inherits(A, "Matrix") &&
      !requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required for Matrix inputs.", call. = FALSE)
  }

  invisible(TRUE)
}

# Build adjacency-list neighbors for the pure-R graph algorithms.
#' @rdname adjacency_neighbors
#' @export
adjacency_neighbors <- function(
    A,
    directed = FALSE,
    self_loops = c("ignore", "include")) {
  validate_adjacency(A)
  if (length(directed) != 1L || !is.logical(directed) || is.na(directed)) {
    stop("directed must be TRUE or FALSE.", call. = FALSE)
  }
  self_loops <- match.arg(self_loops)

  n <- nrow(A)
  neighbors <- vector("list", n)
  if (n == 0L) {
    return(neighbors)
  }

  if (inherits(A, "Matrix")) {
    A_c <- if (inherits(A, "dgCMatrix")) {
      A
    } else {
      methods::as(methods::as(A, "generalMatrix"), "dgCMatrix")
    }
    entries <- Matrix::summary(A_c)
    entries <- entries[entries$x != 0, , drop = FALSE]
    if (self_loops == "ignore") {
      entries <- entries[entries$i != entries$j, , drop = FALSE]
    }
    from <- entries$i
    to <- entries$j
    if (!directed) {
      from_original <- from
      from <- c(from, to)
      to <- c(to, from_original)
    }
    if (length(from) > 0L) {
      split_neighbors <- split(to, factor(from, levels = seq_len(n)))
      neighbors <- lapply(split_neighbors, function(nodes) sort(unique(nodes)))
    }
  } else {
    A <- as.matrix(A)
    for (node in seq_len(n)) {
      connected <- A[node, ] != 0
      if (!directed) {
        connected <- connected | A[, node] != 0
      }
      if (self_loops == "ignore") {
        connected[node] <- FALSE
      }
      neighbors[[node]] <- which(connected)
    }
  }
  neighbors
}

# Compute connected components as a list of one-based node-index vectors.
#
# Components are undirected/weak: an edge in either direction connects nodes.
# Edge weights do not affect membership. Components are returned in discovery
# order, and nodes within each component follow breadth-first search order.
#' @rdname connected_components
#' @export
connected_components <- function(
    A,
    self_loops = c("ignore", "include"),
    use_cpp = TRUE) {
  validate_adjacency(A)
  logical_inputs <- list(use_cpp = use_cpp)
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
  self_loops <- match.arg(self_loops)

  components <- NULL
  is_sparse <- inherits(A, "Matrix")
  cpp_name <- if (is_sparse) {
    "connected_components_sparse_cpp"
  } else {
    "connected_components_dense_cpp"
  }
  if (use_cpp && exists(cpp_name, mode = "function", inherits = TRUE, envir = environment())) {
    cpp_function <- get(cpp_name, mode = "function", inherits = TRUE, envir = environment())
    A_cpp <- if (is_sparse && !inherits(A, "dgCMatrix")) {
      methods::as(methods::as(A, "generalMatrix"), "dgCMatrix")
    } else if (is_sparse) {
      A
    } else {
      as.matrix(A)
    }
    compiled_result <- tryCatch(
      list(
        success = TRUE,
        value = cpp_function(
          A_cpp,
          self_loops == "include"
        )
      ),
      error = function(e) {
        message(cpp_name, " failed; using the R implementation: ",
                conditionMessage(e))
        list(success = FALSE, value = NULL)
      }
    )
    if (compiled_result$success) {
      components <- compiled_result$value
    }
  }

  if (is.null(components)) {
    neighbors <- adjacency_neighbors(
      A,
      directed = FALSE,
      self_loops = self_loops
    )
    n <- length(neighbors)
    visited <- rep(FALSE, n)
    components <- list()
    for (start_node in seq_len(n)) {
      if (visited[start_node]) {
        next
      }
      queue <- integer(n)
      head <- 1L
      tail <- 1L
      queue[tail] <- start_node
      visited[start_node] <- TRUE
      component <- integer()
      while (head <= tail) {
        node <- queue[head]
        head <- head + 1L
        component <- c(component, node)
        for (neighbor in neighbors[[node]]) {
          if (!visited[neighbor]) {
            visited[neighbor] <- TRUE
            tail <- tail + 1L
            queue[tail] <- neighbor
          }
        }
      }
      components[[length(components) + 1L]] <- component
    }
  }

  if (length(components) > 0L) {
    names(components) <- paste0("component_", seq_along(components))
  }
  components
}

# Return the first largest connected component and its induced adjacency matrix.
#
# When components tie in size, discovery order determines the winner, matching
# which.max behavior from the pasted implementation.
#' @rdname connected_components
#' @export
largest_connected_component <- function(
    A,
    self_loops = c("ignore", "include"),
    use_cpp = TRUE,
    sort_nodes = TRUE) {
  validate_adjacency(A)
  if (length(sort_nodes) != 1L ||
      !is.logical(sort_nodes) ||
      is.na(sort_nodes)) {
    stop("sort_nodes must be TRUE or FALSE.", call. = FALSE)
  }

  self_loops <- match.arg(self_loops)
  components <- connected_components(
    A,
    self_loops = self_loops,
    use_cpp = use_cpp
  )
  if (length(components) == 0L) {
    nodes <- integer()
  } else {
    sizes <- lengths(components)
    nodes <- components[[which.max(sizes)]]
  }
  if (sort_nodes) {
    nodes <- sort(nodes)
  }
  submatrix <- A[nodes, nodes, drop = FALSE]
  if (is.null(rownames(submatrix))) {
    rownames(submatrix) <- as.character(nodes)
  }
  if (is.null(colnames(submatrix))) {
    colnames(submatrix) <- as.character(nodes)
  }

  list(
    nodes = nodes,
    size = length(nodes),
    submatrix = submatrix
  )
}

# Build weighted adjacency lists without densifying sparse matrices.
#' @rdname adjacency_neighbors
#' @export
adjacency_weighted_neighbors <- function(
    A,
    directed = FALSE,
    self_loops = c("ignore", "include")) {
  validate_adjacency(A)
  self_loops <- match.arg(self_loops)
  n <- nrow(A)
  nodes <- vector("list", n)
  weights <- vector("list", n)

  if (inherits(A, "Matrix")) {
    A_c <- if (inherits(A, "dgCMatrix")) {
      A
    } else {
      methods::as(methods::as(A, "generalMatrix"), "dgCMatrix")
    }
    entries <- Matrix::summary(A_c)
    entries <- entries[entries$x != 0, , drop = FALSE]
    if (self_loops == "ignore") {
      entries <- entries[entries$i != entries$j, , drop = FALSE]
    }
    from <- entries$i
    to <- entries$j
    edge_weights <- entries$x
    if (!directed) {
      from_original <- from
      from <- c(from, to)
      to <- c(to, from_original)
      edge_weights <- rep(edge_weights, 2L)
    }
    if (length(from) > 0L) {
      edge_groups <- split(
        seq_along(from),
        factor(from, levels = seq_len(n))
      )
      for (node in seq_len(n)) {
        indices <- edge_groups[[node]]
        if (length(indices) == 0L) {
          next
        }
        node_targets <- to[indices]
        target_groups <- split(edge_weights[indices], node_targets)
        nodes[[node]] <- as.integer(names(target_groups))
        weights[[node]] <- vapply(target_groups, min, numeric(1))
      }
    }
  } else {
    A <- as.matrix(A)
    for (node in seq_len(n)) {
      if (directed) {
        node_weights <- A[node, ]
        connected <- node_weights != 0
      } else {
        outgoing <- A[node, ]
        incoming <- A[, node]
        connected <- outgoing != 0 | incoming != 0
        node_weights <- ifelse(
          outgoing == 0,
          incoming,
          ifelse(incoming == 0, outgoing, pmin(outgoing, incoming))
        )
      }
      if (self_loops == "ignore") {
        connected[node] <- FALSE
      }
      nodes[[node]] <- which(connected)
      weights[[node]] <- unname(node_weights[connected])
    }
  }
  list(nodes = nodes, weights = weights)
}

# Compute the n-by-n matrix of shortest-path distances.
#
# Rows are source nodes and columns are destination nodes. Diagonal values are
# zero and unreachable pairs are Inf. With directed = FALSE, an edge in either
# direction is traversable in both directions. If reciprocal edge weights
# differ, the smaller nonzero weight is used. Weighted paths require
# nonnegative weights; zero denotes absence of an edge.
#' @rdname shortest_path_distances
#' @export
shortest_path_distances <- function(
    A,
    directed = FALSE,
    weighted = FALSE,
    self_loops = c("ignore", "include"),
    use_cpp = TRUE) {
  validate_adjacency(A)
  logical_inputs <- list(
    directed = directed,
    weighted = weighted,
    use_cpp = use_cpp
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
  self_loops <- match.arg(self_loops)
  if (weighted && any(A < 0)) {
    stop(
      "weighted shortest paths require nonnegative edge weights.",
      call. = FALSE
    )
  }

  distances <- NULL
  is_sparse <- inherits(A, "Matrix")
  cpp_name <- if (is_sparse) {
    "shortest_path_distances_sparse_cpp"
  } else {
    "shortest_path_distances_dense_cpp"
  }
  if (use_cpp && exists(cpp_name, mode = "function", inherits = TRUE, envir = environment())) {
    cpp_function <- get(cpp_name, mode = "function", inherits = TRUE, envir = environment())
    A_cpp <- if (is_sparse && !inherits(A, "dgCMatrix")) {
      methods::as(methods::as(A, "generalMatrix"), "dgCMatrix")
    } else if (is_sparse) {
      A
    } else {
      as.matrix(A)
    }
    compiled_result <- tryCatch(
      list(
        success = TRUE,
        value = cpp_function(
          A_cpp,
          directed,
          weighted,
          self_loops == "include"
        )
      ),
      error = function(e) {
        message(cpp_name, " failed; using the R implementation: ",
                conditionMessage(e))
        list(success = FALSE, value = NULL)
      }
    )
    if (compiled_result$success) {
      distances <- compiled_result$value
    }
  }

  if (is.null(distances)) {
    if (weighted) {
      adjacency <- adjacency_weighted_neighbors(
        A,
        directed = directed,
        self_loops = self_loops
      )
    } else {
      neighbors <- adjacency_neighbors(
        A,
        directed = directed,
        self_loops = self_loops
      )
    }
    n <- nrow(A)
    distances <- matrix(Inf, nrow = n, ncol = n)
    for (source in seq_len(n)) {
      source_distances <- rep(Inf, n)
      source_distances[source] <- 0
      if (weighted) {
        visited <- rep(FALSE, n)
        for (iteration in seq_len(n)) {
          candidates <- which(!visited)
          if (length(candidates) == 0L) {
            break
          }
          node <- candidates[which.min(source_distances[candidates])]
          if (!is.finite(source_distances[node])) {
            break
          }
          visited[node] <- TRUE
          node_neighbors <- adjacency$nodes[[node]]
          node_weights <- adjacency$weights[[node]]
          for (index in seq_along(node_neighbors)) {
            neighbor <- node_neighbors[index]
            candidate_distance <- source_distances[node] + node_weights[index]
            if (candidate_distance < source_distances[neighbor]) {
              source_distances[neighbor] <- candidate_distance
            }
          }
        }
      } else {
        queue <- integer(n)
        head <- 1L
        tail <- 1L
        queue[tail] <- source
        while (head <= tail) {
          node <- queue[head]
          head <- head + 1L
          for (neighbor in neighbors[[node]]) {
            if (is.infinite(source_distances[neighbor])) {
              source_distances[neighbor] <- source_distances[node] + 1
              tail <- tail + 1L
              queue[tail] <- neighbor
            }
          }
        }
      }
      distances[source, ] <- source_distances
    }
  }

  dimnames(distances) <- dimnames(A)
  distances
}
