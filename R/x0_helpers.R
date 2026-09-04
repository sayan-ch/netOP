# Shared SONNET and NETCROP helpers
#
# SONNET and both NETCROP variants share the core op_splitter() partitioning
# rule defined in this file.

# Generate repeated overlap-plus-partition node subsets.
#
# For each repetition, sample common overlap nodes and partition the rest. When
# s * m is smaller than n - o, the default augments the overlap with all
# remainder nodes. Alternatively, extras can be balanced randomly across the
# non-overlap pieces or, with a warning, ignored.
op_splitter <- function(
    n,
    s,
    o,
    n_repetitions = 1L,
    m = floor((n - o) / s),
    remainder_handling = c(
      "augment_overlap",
      "distribute_randomly",
      "ignore"
    ),
    seed = NULL) {
  remainder_handling <- match.arg(remainder_handling)
  count_inputs <- list(
    n = n,
    s = s,
    o = o,
    n_repetitions = n_repetitions,
    m = m
  )
  invalid_count <- vapply(
    count_inputs,
    function(value) {
      length(value) != 1L ||
        !is.numeric(value) ||
        is.na(value) ||
        !is.finite(value) ||
        value < 0 ||
        value != floor(value)
    },
    logical(1)
  )
  invalid_count[c("n", "s", "n_repetitions")] <-
    invalid_count[c("n", "s", "n_repetitions")] |
    unlist(count_inputs[c("n", "s", "n_repetitions")]) < 1
  if (any(invalid_count)) {
    stop(
      paste(names(count_inputs)[invalid_count], collapse = ", "),
      " must be valid integer counts; n, s, and n_repetitions must be positive.",
      call. = FALSE
    )
  }
  n <- as.integer(n)
  s <- as.integer(s)
  o <- as.integer(o)
  n_repetitions <- as.integer(n_repetitions)
  m <- as.integer(m)
  if (o > n) {
    stop("o cannot exceed n.", call. = FALSE)
  }
  if (as.double(s) * m > n - o) {
    stop("s * m cannot exceed n - o.", call. = FALSE)
  }
  remainder_count <- n - o - s * m
  o_updated <- if (remainder_handling == "augment_overlap") {
    o + remainder_count
  } else {
    o
  }
  if (remainder_handling == "ignore" && remainder_count > 0L) {
    warning(
      remainder_count,
      " node(s) will be ignored in each repetition; this is not recommended.",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (!is.null(seed)) {
    if (length(seed) != 1L ||
        !is.numeric(seed) ||
        is.na(seed) ||
        !is.finite(seed) ||
        seed < 0 ||
        seed != floor(seed) ||
        seed > .Machine$integer.max) {
      stop(
        "seed must be NULL or an integer from 0 to .Machine$integer.max.",
        call. = FALSE
      )
    }
    seed <- as.double(seed)
  }

  node_ids <- seq_len(n)
  repetitions <- lapply(seq_len(n_repetitions), function(repetition) {
    if (!is.null(seed)) {
      repetition_seed <- as.integer(
        (seed + repetition - 1) %% .Machine$integer.max
      )
      set.seed(repetition_seed)
    }
    overlap_nodes <- if (o == 0L) {
      integer(0)
    } else {
      sample(node_ids, size = o, replace = FALSE)
    }
    nonoverlap_pool <- setdiff(node_ids, overlap_nodes)
    ordered_nonoverlap <- if (length(nonoverlap_pool) == 0L) {
      integer(0)
    } else {
      sample(nonoverlap_pool, size = length(nonoverlap_pool), replace = FALSE)
    }
    selected_count <- s * m
    split_extras <- replicate(s, integer(0), simplify = FALSE)

    if (remainder_handling == "augment_overlap") {
      additional_overlap <- if (remainder_count == 0L) {
        integer(0)
      } else {
        ordered_nonoverlap[seq_len(remainder_count)]
      }
      overlap_nodes <- c(overlap_nodes, additional_overlap)
      selected_nonoverlap <- if (selected_count == 0L) {
        integer(0)
      } else {
        start <- remainder_count + 1L
        ordered_nonoverlap[start + seq_len(selected_count) - 1L]
      }
    } else {
      selected_nonoverlap <- if (selected_count == 0L) {
        integer(0)
      } else {
        ordered_nonoverlap[seq_len(selected_count)]
      }
    }

    if (remainder_handling == "distribute_randomly" &&
        remainder_count > 0L) {
      remaining_nodes <- ordered_nonoverlap[
        selected_count + seq_len(remainder_count)
      ]
      extra_counts <- rep.int(remainder_count %/% s, s)
      bonus_count <- remainder_count %% s
      if (bonus_count > 0L) {
        bonus_splits <- sample.int(s, size = bonus_count, replace = FALSE)
        extra_counts[bonus_splits] <- extra_counts[bonus_splits] + 1L
      }
      cursor <- 1L
      for (split_id in seq_len(s)) {
        count <- extra_counts[split_id]
        if (count > 0L) {
          split_extras[[split_id]] <- remaining_nodes[
            cursor + seq_len(count) - 1L
          ]
          cursor <- cursor + count
        }
      }
    }

    splits <- lapply(seq_len(s), function(split_id) {
      split_nonoverlap <- if (m == 0L) {
        integer(0)
      } else {
        first <- (split_id - 1L) * m + 1L
        last <- split_id * m
        selected_nonoverlap[first:last]
      }
      c(overlap_nodes, split_nonoverlap, split_extras[[split_id]])
    })
    names(splits) <- paste0("split_", seq_len(s))
    splits
  })
  names(repetitions) <- paste0("repetition_", seq_len(n_repetitions))
  list(
    splits = repetitions,
    o = o_updated,
    remainder_count = remainder_count,
    remainder_handling = remainder_handling
  )
}

# Select SONNET subnetwork and overlap counts by constrained grid search.
#
# Computational cost is evaluated without forming n^theta explicitly:
#   s * ((o + (n - o) / s) / n)^theta.
# Feasibility requires this cost to be no larger than q * ncores. Candidate
# counts are integers and (n - o) must be divisible by s. By default the
# overlap starts at min(K^3, n - s_lower), while its upper bound n - s is
# applied separately for each candidate s.
#' @rdname model_selection
#' @export
sonnet_param_select <- function(
    n,
    K = 5L,
    theta = 3,
    q = 0.005,
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    s_lower = 2,
    s_upper = 50,
    o_lower = NULL,
    o_upper = NULL) {
  validate_finite_scalar <- function(value, name) {
    if (length(value) != 1L ||
        !is.numeric(value) ||
        is.na(value) ||
        !is.finite(value)) {
      stop(name, " must be one finite numeric value.", call. = FALSE)
    }
    invisible(TRUE)
  }

  validate_finite_scalar(n, "n")
  if (n < 3 || n != floor(n) || n > .Machine$integer.max) {
    stop(
      "n must be one integer from 3 through .Machine$integer.max.",
      call. = FALSE
    )
  }
  n <- as.integer(n)

  validate_finite_scalar(K, "K")
  if (K < 1 || K != floor(K) || K > n) {
    stop("K must be one positive integer no larger than n.",
         call. = FALSE)
  }
  K <- as.integer(K)

  validate_finite_scalar(theta, "theta")
  if (theta <= 0) {
    stop("theta must be positive.", call. = FALSE)
  }
  validate_finite_scalar(q, "q")
  if (q <= 0 || q > 1) {
    stop("q must be in (0, 1].", call. = FALSE)
  }
  validate_finite_scalar(ncores, "ncores")
  if (ncores < 1 || ncores != floor(ncores) ||
      ncores > .Machine$integer.max) {
    stop("ncores must be one positive integer.", call. = FALSE)
  }
  ncores <- as.integer(ncores)
  validate_finite_scalar(s_lower, "s_lower")
  validate_finite_scalar(s_upper, "s_upper")
  if (!is.null(o_lower)) {
    validate_finite_scalar(o_lower, "o_lower")
  }
  if (!is.null(o_upper)) {
    validate_finite_scalar(o_upper, "o_upper")
  }

  s_lower_integer <- max(ceiling(s_lower), 2L)
  s_upper_integer <- min(floor(s_upper), n - 1L)
  o_lower_integer <- if (is.null(o_lower)) {
    min(K^3, n - s_lower_integer)
  } else {
    max(ceiling(o_lower), 1L)
  }
  o_upper_integer <- if (is.null(o_upper)) {
    n - 2L
  } else {
    min(floor(o_upper), n - 2L)
  }
  if (s_lower_integer > s_upper_integer) {
    stop("The normalized s bounds contain no admissible integer values.",
         call. = FALSE)
  }
  if (o_lower_integer > o_upper_integer) {
    stop("The normalized o bounds contain no admissible integer values.",
         call. = FALSE)
  }

  best_objective <- -Inf
  best_num_subnetworks <- NA_integer_
  best_overlap_size <- NA_integer_
  best_piece_size <- NA_real_
  best_computational_fraction <- NA_real_
  fallback_objective <- -Inf
  fallback_num_subnetworks <- NA_integer_
  fallback_overlap_size <- NA_integer_
  fallback_piece_size <- NA_real_
  fallback_computational_fraction <- Inf
  feasible_count <- 0
  evaluated_count <- 0

  for (num_subnetworks_candidate in
       seq.int(s_lower_integer, s_upper_integer)) {
    overlap_upper_candidate <- min(
      o_upper_integer,
      n - num_subnetworks_candidate
    )
    if (o_lower_integer > overlap_upper_candidate) next

    overlap_candidates <- seq.int(
      o_lower_integer,
      overlap_upper_candidate
    )
    evaluated_count <- evaluated_count + length(overlap_candidates)
    divisible_mask <-
      (n - overlap_candidates) %% num_subnetworks_candidate == 0
    overlap_candidates <- overlap_candidates[divisible_mask]
    if (length(overlap_candidates) == 0L) next

    piece_sizes <- (n - overlap_candidates) / num_subnetworks_candidate
    subnetwork_fractions <-
      (overlap_candidates + piece_sizes) / n
    computational_fractions <- num_subnetworks_candidate *
      subnetwork_fractions^theta
    objectives_all <- 1 - ((n - overlap_candidates) / n)^2 *
      ((num_subnetworks_candidate - 1) / num_subnetworks_candidate)
    fallback_idx <- which.min(computational_fractions)
    fallback_fraction_candidate <-
      computational_fractions[fallback_idx]
    fallback_objective_candidate <- objectives_all[fallback_idx]
    if (fallback_fraction_candidate < fallback_computational_fraction ||
        (fallback_fraction_candidate == fallback_computational_fraction &&
         fallback_objective_candidate > fallback_objective)) {
      fallback_objective <- fallback_objective_candidate
      fallback_num_subnetworks <- as.integer(
        ceiling(num_subnetworks_candidate)
      )
      fallback_overlap_size <- as.integer(overlap_candidates[fallback_idx])
      fallback_piece_size <- piece_sizes[fallback_idx]
      fallback_computational_fraction <- fallback_fraction_candidate
    }
    feasible_mask <- computational_fractions <= q * ncores
    feasible_count <- feasible_count + sum(feasible_mask)
    if (!any(feasible_mask)) next

    feasible_overlaps <- overlap_candidates[feasible_mask]
    feasible_piece_sizes <- piece_sizes[feasible_mask]
    feasible_computational_fractions <-
      computational_fractions[feasible_mask]
    objectives <- objectives_all[feasible_mask]
    candidate_idx <- which.max(objectives)
    candidate_objective <- objectives[candidate_idx]

    if (candidate_objective > best_objective) {
      best_objective <- candidate_objective
      best_num_subnetworks <- as.integer(
        ceiling(num_subnetworks_candidate)
      )
      best_overlap_size <- as.integer(feasible_overlaps[candidate_idx])
      best_piece_size <- feasible_piece_sizes[candidate_idx]
      best_computational_fraction <-
        feasible_computational_fractions[candidate_idx]
    }
  }

  q_requested <- q
  if (!is.finite(best_objective) &&
      is.finite(fallback_computational_fraction)) {
    best_objective <- fallback_objective
    best_num_subnetworks <- fallback_num_subnetworks
    best_overlap_size <- fallback_overlap_size
    best_piece_size <- fallback_piece_size
    best_computational_fraction <- fallback_computational_fraction
    q <- fallback_computational_fraction / ncores
    feasible_count <- 1
  }
  if (!is.finite(best_objective)) {
    stop(
      "No (num_subnetworks, overlap_size) pair satisfies the constraint ",
      "within the requested bounds.",
      call. = FALSE
    )
  }

  list(
    num_subnetworks = best_num_subnetworks,
    overlap_size = best_overlap_size,
    m = best_piece_size,
    objective = best_objective,
    data_use = 100 * best_objective,
    computational_fraction = best_computational_fraction,
    parallel_runtime_fraction = best_computational_fraction / ncores,
    K = K,
    theta = theta,
    q = q,
    q_requested = q_requested,
    ncores = ncores,
    feasible_count = feasible_count,
    evaluated_count = evaluated_count,
    bounds = list(
      s_lower = as.integer(s_lower_integer),
      s_upper = as.integer(s_upper_integer),
      o_lower = as.integer(o_lower_integer),
      o_upper = if (is.null(o_upper)) NULL else as.integer(o_upper_integer)
    )
  )
}

# User-provided dplyr prototype, retained as a commented reference.
#
# library(dplyr)
#
# sonnet_param_select <- function(
#     n,
#     theta = 3,
#     q = 0.005,
#     K = 5,
#     s_lower = 2,
#     s_upper = 50,
#     o_lower = min(K^3, n - s_lower),
#     o_upper = n - s_lower,
#     ncore = 20
# ) {
#   tab <- expand.grid(
#     s = seq.int(s_lower, s_upper),
#     o = seq.int(o_lower, o_upper)
#   )
#
#   tab2 <- tab |>
#     filter(o <= n - s) |>
#     filter((n-o) %% s == 0) |>
#     mutate(
#       m = (n - o) / s,
#       q = q
#     )
#
#   tab3 <- tab2 |>
#     mutate(
#       constraint = s * ((o+m) / n)^theta
#     ) |>
#     filter(
#       constraint <= ncore * q
#     )
#
#   tab3 |>
#     mutate(
#       obj = 100 * (1 - (1 - o/n)^2 * (s-1) / s)
#     ) |>
#     slice_max(
#       obj
#     )
# }

# Select NETCROP subnetwork and overlap parameters over a requested range.
#
# The overlap-range positions interpolate between the lower and upper overlap
# proportions implied by test_prop. As in the supplied method, an o_range value
# of exactly one is replaced by 0.8 to avoid the singular upper endpoint. When
# n is supplied, remainder nodes augment the overlap so the remaining nodes
# divide evenly among the selected subnetworks.
#' @rdname model_selection
#' @export
netcrop_param_select <- function(
    test_prop = 0.02,
    n = NULL,
    o_range = c(0, 0.8)) {
  if (length(test_prop) != 1L ||
      !is.numeric(test_prop) ||
      is.na(test_prop) ||
      !is.finite(test_prop) ||
      test_prop <= 0 ||
      test_prop > 0.5) {
    stop("test_prop must be one finite number in (0, 0.5].",
         call. = FALSE)
  }
  if (!is.numeric(o_range) ||
      length(o_range) == 0L ||
      anyNA(o_range) ||
      any(!is.finite(o_range)) ||
      any(o_range < 0 | o_range > 1)) {
    stop("o_range must contain finite values in [0, 1].",
         call. = FALSE)
  }
  if (!is.null(n)) {
    if (length(n) != 1L ||
        !is.numeric(n) ||
        is.na(n) ||
        !is.finite(n) ||
        n < 3 ||
        n != floor(n) ||
        n > .Machine$integer.max) {
      stop(
        "n must be NULL or one integer from 3 through .Machine$integer.max.",
        call. = FALSE
      )
    }
    n <- as.integer(n)
  }

  o_range_adjusted <- o_range
  endpoint_mask <- o_range_adjusted == 1
  if (any(endpoint_mask)) {
    warning(
      "o_range values equal to 1 were replaced by 0.8 to avoid the ",
      "singular upper endpoint.",
      call. = FALSE,
      immediate. = TRUE
    )
    o_range_adjusted[endpoint_mask] <- 0.8
  }
  overlap_upper <- 1 - sqrt(test_prop)
  overlap_lower <- 1 - sqrt(2 * test_prop)
  overlap_proportion <- overlap_lower +
    (overlap_upper - overlap_lower) * o_range_adjusted
  nonoverlap_proportion_squared <- (1 - overlap_proportion)^2
  num_subnetworks_unrounded <-
    nonoverlap_proportion_squared /
      (nonoverlap_proportion_squared - test_prop)
  # rounding_tolerance <- sqrt(.Machine$double.eps) *
  #   pmax(1, abs(num_subnetworks_unrounded))
  rounding_tolerance <- 0
  num_subnetworks <- ceiling(
    num_subnetworks_unrounded - rounding_tolerance
  )
  if (any(!is.finite(num_subnetworks)) ||
      any(num_subnetworks < 2) ||
      any(num_subnetworks > .Machine$integer.max)) {
    stop(
      "The requested parameter range does not yield valid subnetwork counts.",
      call. = FALSE
    )
  }
  num_subnetworks <- as.integer(num_subnetworks)

  if (is.null(n)) {
    return(data.frame(
      test_prop = rep(test_prop, length(overlap_proportion)),
      overlap_proportion = overlap_proportion,
      num_subnetworks = num_subnetworks
    ))
  }

  overlap_size <- ceiling(n * overlap_proportion)
  piece_size <- floor((n - overlap_size) / num_subnetworks)
  if (any(piece_size < 1L)) {
    stop(
      "n is too small for at least one selected parameter pair to retain ",
      "one non-overlap node per subnetwork.",
      call. = FALSE
    )
  }
  overlap_extra <- n - overlap_size - num_subnetworks * piece_size
  overlap_size <- as.integer(overlap_size + overlap_extra)

  data.frame(
    test_prop = rep(test_prop, length(overlap_proportion)),
    overlap_proportion = overlap_proportion,
    n = rep(n, length(overlap_proportion)),
    num_subnetworks = num_subnetworks,
    overlap_size = overlap_size
  )
}

# Generate SONNET partitions with one overlap shared across all repetitions.
#
# Remainder nodes always augment the shared overlap. The first repetition's
# subnetworks include that overlap; extra repetitions contain only freshly
# permuted non-overlap pieces, matching the original SONNET flow.
sonnet_splitter <- function(
    n,
    num_subnetworks,
    overlap_size,
    extra_nrep = 0L,
    m = floor((n - overlap_size) / num_subnetworks),
    seed = NULL) {
  if (length(extra_nrep) != 1L ||
      !is.numeric(extra_nrep) ||
      is.na(extra_nrep) ||
      !is.finite(extra_nrep) ||
      extra_nrep < 0 ||
      extra_nrep != floor(extra_nrep)) {
    stop("extra_nrep must be one non-negative integer.", call. = FALSE)
  }
  extra_nrep <- as.integer(extra_nrep)
  first_split <- op_splitter(
    n = n,
    s = num_subnetworks,
    o = overlap_size,
    n_repetitions = 1L,
    m = m,
    remainder_handling = "augment_overlap",
    seed = seed
  )
  first_subnetworks <- first_split$splits[[1L]]
  overlap_nodes <- Reduce(intersect, first_subnetworks)
  nonoverlap_nodes <- setdiff(seq_len(n), overlap_nodes)
  subnetworks <- vector("list", extra_nrep + 1L)
  subnetworks[[1L]] <- first_subnetworks

  if (extra_nrep > 0L) {
    for (repetition in seq_len(extra_nrep)) {
      if (!is.null(seed)) {
        repetition_seed <- as.integer(
          (as.double(seed) + repetition) %% .Machine$integer.max
        )
        set.seed(repetition_seed)
      }
      ordered_nodes <- if (length(nonoverlap_nodes) == 0L) {
        integer(0)
      } else {
        sample(
          nonoverlap_nodes,
          size = length(nonoverlap_nodes),
          replace = FALSE
        )
      }
      pieces <- lapply(seq_len(num_subnetworks), function(split_id) {
        if (m == 0L) {
          return(integer(0))
        }
        first <- (split_id - 1L) * m + 1L
        last <- split_id * m
        ordered_nodes[first:last]
      })
      names(pieces) <- paste0("split_", seq_len(num_subnetworks))
      subnetworks[[repetition + 1L]] <- pieces
    }
  }
  names(subnetworks) <- paste0(
    "repetition_",
    seq_len(extra_nrep + 1L)
  )
  list(
    subnetworks = subnetworks,
    overlap_nodes = overlap_nodes,
    nonoverlap_nodes = nonoverlap_nodes,
    overlap_size = length(overlap_nodes),
    requested_overlap_size = as.integer(overlap_size),
    remainder_count = first_split$remainder_count,
    m = as.integer(m),
    remainder_handling = "augment_overlap"
  )
}

# Generate SONNET partitions with independently sampled overlap per repetition.
#
# Remainder nodes augment the overlap in every repetition, so every node is
# represented and every repetition produces a complete labeling.
sonnet_splitter_independent <- function(
    n,
    num_subnetworks,
    overlap_size,
    extra_nrep = 0L,
    m = floor((n - overlap_size) / num_subnetworks),
    seed = NULL) {
  split <- op_splitter(
    n = n,
    s = num_subnetworks,
    o = overlap_size,
    n_repetitions = extra_nrep + 1L,
    m = m,
    remainder_handling = "augment_overlap",
    seed = seed
  )
  overlap_nodes <- lapply(split$splits, function(subnetworks) {
    Reduce(intersect, subnetworks)
  })
  nonoverlap_nodes <- lapply(overlap_nodes, function(overlap) {
    setdiff(seq_len(n), overlap)
  })
  list(
    subnetworks = split$splits,
    overlap_nodes = overlap_nodes,
    nonoverlap_nodes = nonoverlap_nodes,
    overlap_size = split$o,
    requested_overlap_size = as.integer(overlap_size),
    remainder_count = split$remainder_count,
    m = as.integer(m),
    remainder_handling = "augment_overlap"
  )
}

# Generate independent overlap partitions for NETCROP repetitions.
#
# Remainder nodes augment the overlap in every repetition. Consequently every
# node is retained, every non-overlap piece has the same size, and the returned
# effective overlap may be larger than overlap_size.
netcrop_splitter <- function(
    n,
    num_subnetworks,
    overlap_size,
    nrep = 1L,
    seed = NULL) {
  split <- op_splitter(
    n = n,
    s = num_subnetworks,
    o = overlap_size,
    n_repetitions = nrep,
    m = floor((n - overlap_size) / num_subnetworks),
    remainder_handling = "augment_overlap",
    seed = seed
  )
  overlap_nodes <- lapply(split$splits, function(subnetworks) {
    Reduce(intersect, subnetworks)
  })
  nonoverlap_pieces <- lapply(seq_along(split$splits), function(repetition) {
    lapply(split$splits[[repetition]], function(nodes) {
      setdiff(nodes, overlap_nodes[[repetition]])
    })
  })
  names(nonoverlap_pieces) <- names(split$splits)
  list(
    subnetworks = split$splits,
    overlap_nodes = overlap_nodes,
    nonoverlap_pieces = nonoverlap_pieces,
    overlap_size = split$o,
    requested_overlap_size = as.integer(overlap_size),
    remainder_count = split$remainder_count,
    piece_size = as.integer(
      floor((n - overlap_size) / num_subnetworks)
    ),
    nrep = as.integer(nrep),
    remainder_handling = "augment_overlap"
  )
}
