# NETCROP dimension selection for symmetric RDPGs
#
# This file preserves the original four-stage computational flow: decompose
# overlapping subnetworks, construct every candidate embedding, align each
# non-reference embedding through the overlap, and evaluate every unordered
# pair of non-overlap pieces. Package loading resolves all internal helpers.

#' @rdname netcrop_rdpg
#' @export
netcrop_rdpg <- function(
    A,
    d_candidates,
    num_subnetworks = NULL,
    overlap_size = NULL,
    nrep = 1L,
    losses = "sse",
    eig_options = list(),
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    seed = NULL,
    verbose = TRUE,
    force_windows = FALSE,
    ram_check = FALSE,
    parameter_select_options = list(),
    retain_intermediates = c("all", "minimal")) {
  call <- match.call()
  retain_intermediates <- match.arg(retain_intermediates)
  required_helpers <- c(
    "uni_mclapply", "is_symmetric_matrix", "eig_decomp", "procrustes",
    "netcrop_param_select", "netcrop_splitter", "modal",
    "estimate_spectral_decomp_ram",
    "estimate_matrix_product_ram", "report_ram_formula"
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
  validate_named_list <- function(value, name) {
    if (!is.list(value) ||
        (length(value) > 0L &&
         (is.null(names(value)) || any(!nzchar(names(value)))))) {
      stop(name, " must be a named list.", call. = FALSE)
    }
    invisible(TRUE)
  }

  if (is.null(dim(A)) || length(dim(A)) != 2L ||
      nrow(A) != ncol(A) || nrow(A) < 2L) {
    stop("A must be a square matrix-like object with at least two rows.",
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
    stop("A must be symmetric for symmetric RDPG NETCROP.", call. = FALSE)
  }
  A_diagonal <- if (inherits(A, "Matrix")) Matrix::diag(A) else diag(A)
  if (any(A_diagonal != 0)) {
    stop("A must have a zero diagonal (no self-loops).", call. = FALSE)
  }
  if (any(A_values < 0)) {
    stop("A must be non-negative.", call. = FALSE)
  }
  if (sum(A_values) == 0) {
    stop("A is an empty graph; RDPG NETCROP is undefined.", call. = FALSE)
  }
  if (inherits(A, "symmetricMatrix")) {
    A <- methods::as(A, "generalMatrix")
  }

  n <- nrow(A)
  nrep <- validate_count(nrep, "nrep", minimum = 1L)
  ncores <- validate_count(ncores, "ncores", minimum = 1L)
  if (is.null(num_subnetworks) || is.null(overlap_size)) {
    validate_named_list(parameter_select_options, "parameter_select_options")
    if ("n" %in% names(parameter_select_options)) {
      stop("parameter_select_options cannot override n.", call. = FALSE)
    }
    unknown_parameter_options <- setdiff(
      names(parameter_select_options),
      names(formals(netcrop_param_select))
    )
    if (length(unknown_parameter_options) > 0L) {
      stop(
        "parameter_select_options contains unsupported component(s): ",
        paste(unknown_parameter_options, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    selected_parameters <- do.call(
      netcrop_param_select,
      utils::modifyList(
        list(test_prop = 0.02, n = n, o_range = 0),
        parameter_select_options,
        keep.null = TRUE
      )
    )
    if (nrow(selected_parameters) != 1L) {
      stop(
        "parameter_select_options must select exactly one NETCROP pair.",
        call. = FALSE
      )
    }
    if (is.null(num_subnetworks)) {
      num_subnetworks <- selected_parameters$num_subnetworks[[1L]]
    }
    if (is.null(overlap_size)) {
      overlap_size <- selected_parameters$overlap_size[[1L]]
    }
  }
  num_subnetworks <- validate_count(
    num_subnetworks,
    "num_subnetworks",
    minimum = 2L
  )
  overlap_size <- validate_count(overlap_size, "overlap_size", minimum = 1L)
  verbose <- validate_flag(verbose, "verbose")
  force_windows <- validate_flag(force_windows, "force_windows")
  ram_check <- validate_flag(ram_check, "ram_check")
  if (overlap_size >= n) {
    stop("overlap_size must be smaller than nrow(A).", call. = FALSE)
  }
  if (floor((n - overlap_size) / num_subnetworks) < 1L) {
    stop("Every non-overlap piece must contain at least one node.",
         call. = FALSE)
  }
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
        !is.finite(seed) || seed < 0 || seed != floor(seed) ||
        seed > .Machine$integer.max) {
      stop(
        "seed must be NULL or an integer from 0 to .Machine$integer.max.",
        call. = FALSE
      )
    }
    seed <- as.double(seed)
  }
  if (!is.numeric(d_candidates) || length(d_candidates) < 1L ||
      anyNA(d_candidates) || any(!is.finite(d_candidates)) ||
      any(d_candidates < 1) || any(d_candidates != floor(d_candidates))) {
    stop("d_candidates must contain positive integers.", call. = FALSE)
  }
  d_candidates <- unique(as.integer(d_candidates))
  d_max <- max(d_candidates)
  if (!is.character(losses) || length(losses) < 1L ||
      anyNA(losses) || any(!nzchar(losses))) {
    stop("losses must contain non-empty function names.", call. = FALSE)
  }
  losses <- unique(losses)
  validate_named_list(eig_options, "eig_options")
  protected_eig_options <- c(
    "A", "d", "only_values", "order_by", "validate_inputs"
  )
  conflicts <- intersect(names(eig_options), protected_eig_options)
  if (length(conflicts) > 0L) {
    stop(
      "eig_options cannot override ", paste(conflicts, collapse = ", "), ".",
      call. = FALSE
    )
  }
  unsupported_eig_options <- setdiff(
    names(eig_options),
    names(formals(eig_decomp))
  )
  if (length(unsupported_eig_options) > 0L) {
    stop(
      "eig_options contains unsupported argument(s): ",
      paste(unsupported_eig_options, collapse = ", "), ".",
      call. = FALSE
    )
  }
  resolved_eig_options <- utils::modifyList(
    list(
      scale_by = "none",
      use_laplacian = FALSE,
      engine = "rspectra",
      force_engine = TRUE,
      safe_d_multiplier = 1
    ),
    eig_options,
    keep.null = TRUE
  )

  selected_loss_functions <- lapply(losses, function(loss_name) {
    if (!exists(loss_name, mode = "function", inherits = TRUE, envir = environment())) {
      stop("Loss function not found: ", loss_name, ".", call. = FALSE)
    }
    get(loss_name, mode = "function", inherits = TRUE, envir = environment())
  })
  names(selected_loss_functions) <- losses
  bundled_loss_names <- intersect(
    c(
      "validate_error_inputs", "clip_values", "clip_probabilities",
      "sse", "sae", "mse", "mae", "bin_dev", "bin_dev_mean",
      "auc", "auc_as_loss"
    ),
    ls(envir = .GlobalEnv, all.names = TRUE)
  )
  loss_environment <- new.env(parent = baseenv())
  for (function_name in bundled_loss_names) {
    loss_function <- get(
      function_name,
      envir = .GlobalEnv,
      mode = "function",
      inherits = FALSE
    )
    environment(loss_function) <- loss_environment
    assign(function_name, loss_function, envir = loss_environment)
  }
  loss_functions <- lapply(seq_along(losses), function(loss_id) {
    loss_name <- losses[loss_id]
    if (exists(
      loss_name,
      envir = loss_environment,
      mode = "function",
      inherits = FALSE
    )) {
      get(
        loss_name,
        envir = loss_environment,
        mode = "function",
        inherits = FALSE
      )
    } else {
      selected_loss_functions[[loss_id]]
    }
  })
  names(loss_functions) <- losses

  splitter <- netcrop_splitter(
    n = n,
    num_subnetworks = num_subnetworks,
    overlap_size = overlap_size,
    nrep = nrep,
    seed = seed
  )
  effective_overlap_size <- splitter$overlap_size
  piece_size <- splitter$piece_size
  subgraph_size <- effective_overlap_size + piece_size
  if (d_max >= subgraph_size) {
    stop(
      "max(d_candidates) must be smaller than the effective subgraph size (",
      subgraph_size, ").",
      call. = FALSE
    )
  }

  degree_values <- if (inherits(A, "Matrix")) {
    Matrix::rowSums(A)
  } else {
    rowSums(A)
  }
  rho_hat <- mean(degree_values) / (n - 1)
  A_spectral <- A
  raw_grid <- data.frame(
    repetition = rep(seq_len(nrep), each = num_subnetworks),
    split = rep(seq_len(num_subnetworks), times = nrep)
  )
  candidate_grid <- expand.grid(
    candidate_id = seq_along(d_candidates),
    split = seq_len(num_subnetworks),
    repetition = seq_len(nrep),
    KEEP.OUT.ATTRS = FALSE
  )
  pair_grid <- t(utils::combn(num_subnetworks, 2L))
  loss_grid <- expand.grid(
    pair_id = seq_len(nrow(pair_grid)),
    candidate_id = seq_along(d_candidates),
    repetition = seq_len(nrep),
    KEEP.OUT.ATTRS = FALSE
  )
  raw_ncores <- min(ncores, nrow(raw_grid))
  embedding_ncores <- min(ncores, nrow(candidate_grid))
  alignment_ncores <- min(ncores, nrow(candidate_grid))
  loss_ncores <- min(ncores, nrow(loss_grid))
  raw_task_lookup <- matrix(
    seq_len(nrow(raw_grid)),
    nrow = nrep,
    ncol = num_subnetworks,
    byrow = TRUE
  )
  candidate_task_lookup <- array(
    seq_len(nrow(candidate_grid)),
    dim = c(length(d_candidates), num_subnetworks, nrep)
  )
  loss_accepts_prevalidation <- vapply(
    loss_functions,
    function(fun) "validate_inputs" %in% names(formals(fun)),
    logical(1)
  )

  manage_future_plan <- TRUE
  os_type <- if (force_windows) "windows" else .Platform$OS.type
  worker_future_packages <- if (inherits(A, "Matrix")) "Matrix" else NULL
  maximum_stage_ncores <- max(
    raw_ncores, embedding_ncores, alignment_ncores, loss_ncores
  )
  if (os_type != "unix" && maximum_stage_ncores > 1L) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop(
        "future and future.apply are required for Windows parallel execution.",
        call. = FALSE
      )
    }
    previous_future_plan <- future::plan()
    previous_renv_sync_check <- Sys.getenv(
      "RENV_CONFIG_SYNCHRONIZED_CHECK",
      unset = NA_character_
    )
    Sys.setenv(RENV_CONFIG_SYNCHRONIZED_CHECK = "FALSE")
    on.exit(
      try(future::plan(previous_future_plan), silent = TRUE),
      add = TRUE
    )
    on.exit({
      if (is.na(previous_renv_sync_check)) {
        Sys.unsetenv("RENV_CONFIG_SYNCHRONIZED_CHECK")
      } else {
        Sys.setenv(
          RENV_CONFIG_SYNCHRONIZED_CHECK = previous_renv_sync_check
        )
      }
    }, add = TRUE)
    future::plan(future::multisession, workers = maximum_stage_ncores)
    manage_future_plan <- FALSE
  }

  ram_report <- NULL
  if (ram_check) {
    decomposition_ram <- do.call(
      estimate_spectral_decomp_ram,
      c(
        list(
          n = subgraph_size,
          p = subgraph_size,
          K = d_max,
          method = "eigen",
          dense_input = !inherits(A, "sparseMatrix")
        ),
        resolved_eig_options[c("engine", "force_engine", "safe_d_multiplier")]
      )
    )
    embedding_ram <- 8 * as.double(subgraph_size) * d_max
    alignment_ram <- estimate_matrix_product_ram(
      nrow_left = piece_size,
      shared_dimension = d_max,
      ncol_right = d_max
    )
    prediction_ram <- estimate_matrix_product_ram(
      nrow_left = piece_size,
      shared_dimension = d_max,
      ncol_right = piece_size
    )
    ram_report <- report_ram_formula(
      terms = list(
        list(
          estimated_bytes = decomposition_ram$estimated_bytes,
          sequential_count = 1L,
          parallel_count = raw_ncores,
          label = "subnetwork eigendecomposition"
        ),
        list(
          estimated_bytes = embedding_ram,
          sequential_count = 1L,
          parallel_count = embedding_ncores,
          label = "candidate embedding"
        ),
        list(
          estimated_bytes = alignment_ram,
          sequential_count = 1L,
          parallel_count = alignment_ncores,
          label = "non-overlap rotation"
        ),
        list(
          estimated_bytes = prediction_ram,
          sequential_count = 1L,
          parallel_count = loss_ncores,
          label = "held-out probability product"
        )
      ),
      operation = "RDPG NETCROP conservative combined preflight",
      detail = "four sequential parallel stages; peak is overestimated additively"
    )
  }

  if (verbose) {
    message(
      "Stage 1/4: decomposing ", nrow(raw_grid),
      " subnetworks with ", raw_ncores, " worker(s)."
    )
  }
  raw_time <- system.time({
    raw_output <- uni_mclapply(
      seq_len(nrow(raw_grid)),
      function(task_id) {
        repetition <- raw_grid$repetition[task_id]
        split_id <- raw_grid$split[task_id]
        nodes <- splitter$subnetworks[[repetition]][[split_id]]
        A_subnetwork <- A_spectral[nodes, nodes, drop = FALSE]
        decomposition <- do.call(
          eig_decomp,
          c(
            list(
              A = A_subnetwork,
              d = d_max,
              only_values = FALSE,
              order_by = "magnitude",
              validate_inputs = FALSE
            ),
            resolved_eig_options
          )
        )
        list(
          U = decomposition$vectors,
          values = decomposition$values,
          negative_eigenvalues = sum(decomposition$values < 0)
        )
      },
      ncores = raw_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      manage_future_plan = manage_future_plan,
      future_packages = worker_future_packages
    )
  })
  eigen_diagnostics <- do.call(rbind, lapply(seq_along(raw_output), function(id) {
    data.frame(
      repetition = raw_grid$repetition[id],
      split = raw_grid$split[id],
      negative_eigenvalues = raw_output[[id]]$negative_eigenvalues,
      stringsAsFactors = FALSE
    )
  }))
  if (verbose) {
    message(
      "Stage 2/4: constructing ", nrow(candidate_grid),
      " candidate embeddings with ", embedding_ncores, " worker(s)."
    )
  }
  embedding_time <- system.time({
    embedding_output <- uni_mclapply(
      seq_len(nrow(candidate_grid)),
      function(task_id) {
        repetition <- candidate_grid$repetition[task_id]
        split_id <- candidate_grid$split[task_id]
        candidate_id <- candidate_grid$candidate_id[task_id]
        d <- d_candidates[candidate_id]
        raw_id <- raw_task_lookup[repetition, split_id]
        U <- raw_output[[raw_id]]$U[, seq_len(d), drop = FALSE]
        scales <- sqrt(abs(raw_output[[raw_id]]$values[seq_len(d)]))
        sweep(U, 2L, scales, `*`)
      },
      ncores = embedding_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      manage_future_plan = manage_future_plan,
      future_packages = worker_future_packages
    )
  })
  if (retain_intermediates == "minimal") {
    raw_output <- NULL
  }
  if (verbose) {
    message(
      "Stage 3/4: aligning ", nrow(candidate_grid),
      " non-overlap embeddings with ", alignment_ncores, " worker(s)."
    )
  }
  alignment_time <- system.time({
    aligned_output <- uni_mclapply(
      seq_len(nrow(candidate_grid)),
      function(task_id) {
        repetition <- candidate_grid$repetition[task_id]
        split_id <- candidate_grid$split[task_id]
        candidate_id <- candidate_grid$candidate_id[task_id]
        X <- embedding_output[[task_id]]
        overlap_count <- length(splitter$overlap_nodes[[repetition]])
        if (split_id == 1L) {
          return(X[-seq_len(overlap_count), , drop = FALSE])
        }
        standard_id <- candidate_task_lookup[candidate_id, 1L, repetition]
        rotation <- procrustes(
          X = X[seq_len(overlap_count), , drop = FALSE],
          X_star = embedding_output[[standard_id]][
            seq_len(overlap_count), , drop = FALSE
          ],
          validate_inputs = FALSE
        )$rotation
        tcrossprod(
          X[-seq_len(overlap_count), , drop = FALSE],
          t(rotation)
        )
      },
      ncores = alignment_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      manage_future_plan = manage_future_plan,
      future_packages = worker_future_packages
    )
  })
  if (retain_intermediates == "minimal") {
    embedding_output <- NULL
  }

  pair_count <- nrow(pair_grid)
  if (verbose) {
    message(
      "Stage 4/4: evaluating ", nrow(loss_grid),
      " held-out subnetwork pairs with ", loss_ncores, " worker(s)."
    )
  }
  loss_time <- system.time({
    loss_output <- uni_mclapply(
      seq_len(nrow(loss_grid)),
      function(task_id) {
        repetition <- loss_grid$repetition[task_id]
        candidate_id <- loss_grid$candidate_id[task_id]
        pair_id <- loss_grid$pair_id[task_id]
        split_left <- pair_grid[pair_id, 1L]
        split_right <- pair_grid[pair_id, 2L]
        nodes_left <- splitter$nonoverlap_pieces[[repetition]][[split_left]]
        nodes_right <- splitter$nonoverlap_pieces[[repetition]][[split_right]]
        A_test <- A[nodes_left, nodes_right, drop = FALSE]
        left_id <- candidate_task_lookup[
          candidate_id, split_left, repetition
        ]
        right_id <- candidate_task_lookup[
          candidate_id, split_right, repetition
        ]
        P_hat <- tcrossprod(
          aligned_output[[left_id]],
          aligned_output[[right_id]]
        )
        P_hat <- pmax(P_hat, 1e-6)
        P_hat <- pmin(P_hat, 1 - 1e-6)
        A_test_numeric <- as.numeric(A_test)
        P_hat_numeric <- as.numeric(P_hat)
        records <- matrix(NA_real_, nrow = length(losses), ncol = 6L)
        for (loss_id in seq_along(losses)) {
          loss_name <- losses[loss_id]
          loss_arguments <- list(A_test_numeric, P_hat_numeric)
          if (loss_accepts_prevalidation[loss_id]) {
            loss_arguments$validate_inputs <- FALSE
          }
          value <- do.call(
            loss_functions[[loss_name]], loss_arguments
          ) / pair_count
          if (length(value) != 1L || !is.numeric(value) ||
              is.na(value) || !is.finite(value)) {
            stop(
              "Non-finite loss for repetition ", repetition,
              ", d = ", d_candidates[candidate_id],
              ", loss = ", loss_name,
              ", split pair = ", split_left, "-", split_right, ".",
              call. = FALSE
            )
          }
          records[loss_id, ] <- c(
            repetition,
            d_candidates[candidate_id],
            loss_id,
            split_left,
            split_right,
            as.numeric(value)
          )
        }
        records
      },
      ncores = loss_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      manage_future_plan = manage_future_plan,
      future_packages = worker_future_packages
    )
  })
  if (retain_intermediates == "minimal") {
    aligned_output <- NULL
  }

  loss_matrix <- do.call(rbind, loss_output)
  cv_all_loss <- data.frame(
    repetition = as.integer(loss_matrix[, 1L]),
    d = as.integer(loss_matrix[, 2L]),
    loss = losses[as.integer(loss_matrix[, 3L])],
    split_left = as.integer(loss_matrix[, 4L]),
    split_right = as.integer(loss_matrix[, 5L]),
    loss_value = loss_matrix[, 6L],
    stringsAsFactors = FALSE
  )
  grouping <- interaction(
    cv_all_loss$repetition,
    cv_all_loss$d,
    cv_all_loss$loss,
    drop = TRUE,
    lex.order = TRUE
  )
  cv_loss <- do.call(rbind, lapply(split(cv_all_loss, grouping), function(x) {
    data.frame(
      repetition = x$repetition[1L],
      d = x$d[1L],
      loss = x$loss[1L],
      average_loss = sum(x$loss_value),
      stringsAsFactors = FALSE
    )
  }))
  rownames(cv_loss) <- NULL
  cv_loss <- cv_loss[order(
    cv_loss$repetition,
    match(cv_loss$loss, losses),
    match(cv_loss$d, d_candidates)
  ), ]
  best_dimension_cv <- do.call(rbind, lapply(seq_len(nrep), function(repetition) {
    do.call(rbind, lapply(losses, function(loss_name) {
      candidates <- cv_loss[
        cv_loss$repetition == repetition & cv_loss$loss == loss_name,
        ,
        drop = FALSE
      ]
      best_id <- which.min(candidates$average_loss)
      data.frame(
        repetition = repetition,
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
      best_dimension_cv$loss == loss_name,
      ,
      drop = FALSE
    ]
    data.frame(
      loss = loss_name,
      d_hat = modal(selected$d_hat),
      mean_average_loss = mean(selected$average_loss),
      stringsAsFactors = FALSE
    )
  }))
  rownames(overall_best) <- NULL
  unordered_test_proportion <-
    choose(num_subnetworks, 2) * piece_size^2 / choose(n, 2)

  out <- list(
    algorithm = "NETCROP",
    cv_loss = cv_loss,
    cv_all_loss = cv_all_loss,
    best_dimension_cv = best_dimension_cv,
    overall_best = overall_best,
    d_candidates = d_candidates,
    num_subnetworks = num_subnetworks,
    requested_overlap_size = overlap_size,
    effective_overlap_size = effective_overlap_size,
    piece_size = piece_size,
    unordered_test_proportion = unordered_test_proportion,
    nrep = nrep,
    losses = losses,
    rho_hat = rho_hat,
    eig_options = resolved_eig_options,
    splitter = splitter,
    raw_grid = raw_grid,
    candidate_grid = candidate_grid,
    loss_grid = loss_grid,
    raw_output = if (retain_intermediates == "all") raw_output else NULL,
    embedding_output = if (retain_intermediates == "all") {
      embedding_output
    } else {
      NULL
    },
    aligned_output = if (retain_intermediates == "all") {
      aligned_output
    } else {
      NULL
    },
    retain_intermediates = retain_intermediates,
    eigen_diagnostics = eigen_diagnostics,
    ncores = list(
      requested = ncores,
      decomposition = raw_ncores,
      embedding = embedding_ncores,
      alignment = alignment_ncores,
      loss = loss_ncores
    ),
    timing = c(
      decomposition = unname(raw_time["elapsed"]),
      embedding = unname(embedding_time["elapsed"]),
      alignment = unname(alignment_time["elapsed"]),
      loss = unname(loss_time["elapsed"]),
      total = unname(
        raw_time["elapsed"] + embedding_time["elapsed"] +
          alignment_time["elapsed"] + loss_time["elapsed"]
      )
    ),
    ram_preflight = ram_report,
    seed = seed,
    call = call
  )
  class(out) <- "netcrop_rdpg"
  out
}

# Print the selected RDPG dimension for every requested loss.
#' @rdname netcrop_rdpg
#' @export
print.netcrop_rdpg <- function(x, ...) {
  algorithm <- if (is.null(x$algorithm)) "NETCROP" else toupper(x$algorithm)
  cat(algorithm, "results for symmetric RDPG\n")
  cat(strrep("-", nchar(algorithm) + 27L), "\n", sep = "")
  print(x$overall_best, row.names = FALSE)
  invisible(x)
}

# Summarize a symmetric-RDPG NETCROP fit.
#' @rdname netcrop_rdpg
#' @export
summary.netcrop_rdpg <- function(object, ...) {
  algorithm <- if (is.null(object$algorithm)) {
    "NETCROP"
  } else {
    toupper(object$algorithm)
  }
  if (algorithm == "ECV") {
    result <- list(
      algorithm = algorithm,
      call = object$call,
      d_candidates = object$d_candidates,
      nrep = object$nrep,
      completed_repetitions = object$completed_repetitions,
      valid_repetitions = object$valid_repetitions,
      cv = object$cv,
      train_proportion = object$train_proportion,
      holdout_proportion = object$holdout_proportion,
      holdout_count = object$holdout_count,
      best_dimension_cv = object$best_dimension_cv,
      overall_best = object$overall_best,
      failure_diagnostics = object$failure_diagnostics,
      ncores = object$ncores,
      timing = object$timing
    )
    class(result) <- "summary.netcrop_rdpg"
    return(result)
  }
  result <- list(
    algorithm = algorithm,
    call = object$call,
    d_candidates = object$d_candidates,
    nrep = object$nrep,
    num_subnetworks = object$num_subnetworks,
    requested_overlap_size = object$requested_overlap_size,
    effective_overlap_size = object$effective_overlap_size,
    piece_size = object$piece_size,
    unordered_test_proportion = object$unordered_test_proportion,
    rho_hat = object$rho_hat,
    best_dimension_cv = object$best_dimension_cv,
    overall_best = object$overall_best,
    negative_eigenvalues = sum(
      object$eigen_diagnostics$negative_eigenvalues
    ),
    ncores = object$ncores,
    timing = object$timing
  )
  class(result) <- "summary.netcrop_rdpg"
  result
}

# Print a symmetric-RDPG NETCROP summary.
#' @rdname netcrop_rdpg
#' @export
print.summary.netcrop_rdpg <- function(x, ...) {
  algorithm <- if (is.null(x$algorithm)) "NETCROP" else toupper(x$algorithm)
  if (algorithm == "ECV") {
    cat("Summary of ECV symmetric-RDPG dimension selection\n")
    cat("-------------------------------------------------\n")
    cat("Candidate d:", paste(x$d_candidates, collapse = ", "), "\n")
    cat(
      "Completed repetitions:", x$completed_repetitions,
      "of", x$nrep, "requested\n"
    )
    cat("CV folds per repetition:", x$cv, "\n")
    cat(sprintf("Training proportion: %.2f%%\n", 100 * x$train_proportion))
    cat(sprintf("Holdout proportion: %.2f%%\n", 100 * x$holdout_proportion))
    cat("Held-out unordered pairs per fold:", x$holdout_count, "\n")
    cat("Best dimensions per repetition:\n")
    print(x$best_dimension_cv, row.names = FALSE)
    cat("Overall best dimensions:\n")
    print(x$overall_best, row.names = FALSE)
    if (!is.null(x$failure_diagnostics) &&
        nrow(x$failure_diagnostics) > 0L) {
      cat("Failed repetitions:\n")
      print(x$failure_diagnostics, row.names = FALSE)
    }
    cat("Timing (seconds):\n")
    print(x$timing)
    return(invisible(x))
  }
  cat("Summary of NETCROP symmetric-RDPG dimension selection\n")
  cat("----------------------------------------------------\n")
  cat("Candidate d:", paste(x$d_candidates, collapse = ", "), "\n")
  cat("Repetitions:", x$nrep, "\n")
  cat("Subnetworks per repetition:", x$num_subnetworks, "\n")
  cat(
    "Overlap size:", x$effective_overlap_size,
    if (x$effective_overlap_size != x$requested_overlap_size) {
      paste0(" (requested ", x$requested_overlap_size, ")")
    } else {
      ""
    },
    "\n"
  )
  cat("Non-overlap piece size:", x$piece_size, "\n")
  cat(sprintf("Estimated sparsity rho_hat: %.6g\n", x$rho_hat))
  cat(sprintf(
    "Held-out unordered-pair proportion: %.2f%%\n",
    100 * x$unordered_test_proportion
  ))
  cat(
    "Negative retained eigenvalues across decompositions:",
    x$negative_eigenvalues, "\n"
  )
  cat("Best dimensions per repetition:\n")
  print(x$best_dimension_cv, row.names = FALSE)
  cat("Overall best dimensions:\n")
  print(x$overall_best, row.names = FALSE)
  cat("Timing (seconds):\n")
  print(x$timing)
  invisible(x)
}

# Plot RDPG CV loss curves, optionally aggregating across repetitions.
#' @rdname netcrop_rdpg
#' @export
plot.netcrop_rdpg <- function(x, aggregate = TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting.", call. = FALSE)
  }
  dots <- list(...)
  dot_is_result <- vapply(dots, inherits, logical(1), what = "netcrop_rdpg")
  if (inherits(aggregate, "netcrop_rdpg") || any(dot_is_result)) {
    if (!exists("plot_rdpg_comparison", mode = "function", inherits = TRUE, envir = environment())) {
      stop(
        "Required comparison helpers are unavailable; reinstall netOP.",
        call. = FALSE
      )
    }
    comparison_results <- c(
      list(x),
      if (inherits(aggregate, "netcrop_rdpg")) list(aggregate) else NULL,
      dots[dot_is_result]
    )
    comparison_options <- dots[!dot_is_result]
    loss_scale <- if ("loss_scale" %in% names(comparison_options)) {
      comparison_options$loss_scale
    } else {
      "relative"
    }
    return(do.call(
      plot_rdpg_comparison,
      c(comparison_results, list(loss_scale = loss_scale))
    ))
  }
  if (length(aggregate) != 1L || !is.logical(aggregate) || is.na(aggregate)) {
    stop("aggregate must be TRUE or FALSE.", call. = FALSE)
  }
  d_breaks <- sort(unique(x$cv_loss$d))
  algorithm <- if (is.null(x$algorithm)) "NETCROP" else toupper(x$algorithm)
  valid_cv_loss <- x$cv_loss[is.finite(x$cv_loss$average_loss), , drop = FALSE]
  if (nrow(valid_cv_loss) == 0L) {
    stop("No finite CV losses are available for plotting.", call. = FALSE)
  }
  if (!aggregate) {
    plot_data <- valid_cv_loss
    plot_data$repetition <- factor(plot_data$repetition)
    return(
      ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
          x = d,
          y = average_loss,
          group = repetition,
          color = repetition
        )
      ) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::scale_x_continuous(breaks = d_breaks) +
        ggplot2::facet_wrap(~loss, scales = "free_y") +
        ggplot2::labs(
          title = paste(algorithm, "CV loss by RDPG dimension"),
          x = "Latent dimension (d)",
          y = "Average CV loss",
          color = "Repetition"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")
    )
  }
  grouping <- interaction(
    valid_cv_loss$d,
    valid_cv_loss$loss,
    drop = TRUE,
    lex.order = TRUE
  )
  plot_data <- do.call(rbind, lapply(split(valid_cv_loss, grouping), function(z) {
    standard_deviation <- stats::sd(z$average_loss)
    if (is.na(standard_deviation)) {
      standard_deviation <- 0
    }
    data.frame(
      d = z$d[1L],
      loss = z$loss[1L],
      mean_loss = mean(z$average_loss),
      sd_loss = standard_deviation,
      stringsAsFactors = FALSE
    )
  }))
  plot_data$lower <- plot_data$mean_loss - plot_data$sd_loss
  plot_data$upper <- plot_data$mean_loss + plot_data$sd_loss
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = d,
      y = mean_loss,
      ymin = lower,
      ymax = upper
    )
  ) +
    ggplot2::geom_ribbon(alpha = 0.2) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = d_breaks) +
    ggplot2::facet_wrap(~loss, scales = "free_y") +
    ggplot2::labs(
      title = paste(algorithm, "CV loss by RDPG dimension"),
      x = "Latent dimension (d)",
      y = "Mean CV loss (plus or minus one SD)"
    ) +
    ggplot2::theme_minimal()
}

# system.time(
#   net <- generate_rdpg(
#     n = 10000, d = 10, ncores = 5
#   )
# )
#
# # system.time(
# #   net <- RDPG.gen(n = 10000, d = 10, X = NULL, rho = 0.65,
# #                   ncore = 5)$A
# # )
#
# system.time(
#   rdpg_out <- netcrop_rdpg(
#     A = net,
#     d_candidates = 1:20,
#     num_subnetworks = 3L,
#     overlap_size = 8002L,
#     nrep = 1L,
#     losses = "bin_dev",
#     ncores = 5L,
#     verbose = TRUE,
#     ram_check = FALSE
#   )
# )
# rdpg_out
# summary(rdpg_out)
# plot(rdpg_out)
