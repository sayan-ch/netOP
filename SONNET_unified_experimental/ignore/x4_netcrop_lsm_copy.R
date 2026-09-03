# Copy of NETCROP dimension selection for reparameterized latent-space models
#
# This implementation preserves the pasted three-stage computational flow:
# fit every candidate on every overlapping subnetwork, exactly reparameterize
# the fitted inner-product LSM as a squared-distance LSM, rigidly align each fit
# to the first split through the overlap, and evaluate every unordered pair of
# non-overlap pieces. Source the numbered helpers, x0_helpers.R, and
# x1_sonnet.R before this file.

# Select an LSM latent dimension by overlapping-subnetwork cross-validation.
netcrop_lsm_reparameterized <- function(
    A,
    d_candidates,
    num_subnetworks = NULL,
    overlap_size = NULL,
    nrep = 1L,
    losses = c("sse", "bin_dev", "auc_as_loss"),
    lsm_options = list(),
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
    "uni_mclapply", "is_symmetric_matrix", "lsm_pgd", "procrustes",
    "netcrop_param_select", "netcrop_splitter", "modal",
    "estimate_spectral_decomp_ram", "estimate_matrix_product_ram",
    "report_ram_formula"
  )
  missing_helpers <- required_helpers[!vapply(
    required_helpers,
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE
  )]
  if (length(missing_helpers) > 0L) {
    stop(
      "Source the numbered helper files, x0_helpers.R, and x1_sonnet.R ",
      "before x4_netcrop_lsm_copy.R. Missing: ",
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
  if (any(!A_values %in% c(0, 1))) {
    stop("A must contain only binary values 0 and 1.", call. = FALSE)
  }
  if (!is_symmetric_matrix(A)) {
    stop("A must be symmetric for LSM NETCROP.", call. = FALSE)
  }
  A_diagonal <- if (inherits(A, "Matrix")) Matrix::diag(A) else diag(A)
  if (any(A_diagonal != 0)) {
    stop("A must have a zero diagonal (no self-loops).", call. = FALSE)
  }
  if (sum(A_values) == 0) {
    stop("A is an empty graph; LSM NETCROP is undefined.", call. = FALSE)
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
        paste(unknown_parameter_options, collapse = ", "), ".",
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

  validate_named_list(lsm_options, "lsm_options")
  protected_lsm_options <- c("A", "d", "Z_init", "alpha_init", "ram_check")
  lsm_conflicts <- intersect(names(lsm_options), protected_lsm_options)
  if (length(lsm_conflicts) > 0L) {
    stop(
      "lsm_options cannot override ",
      paste(lsm_conflicts, collapse = ", "), ".",
      call. = FALSE
    )
  }
  allowed_lsm_options <- c("step_size", "niter", "trace", "epsilon", "use_cpp")
  unsupported_lsm_options <- setdiff(names(lsm_options), allowed_lsm_options)
  if (length(unsupported_lsm_options) > 0L) {
    stop(
      "lsm_options contains unsupported component(s): ",
      paste(unsupported_lsm_options, collapse = ", "), ".",
      call. = FALSE
    )
  }
  resolved_lsm_options <- utils::modifyList(
    list(
      step_size = 0.3,
      niter = 100L,
      trace = FALSE,
      epsilon = 1e-6,
      use_cpp = TRUE
    ),
    lsm_options,
    keep.null = TRUE
  )

  selected_loss_functions <- lapply(losses, function(loss_name) {
    if (!exists(loss_name, mode = "function", inherits = TRUE)) {
      stop("Loss function not found: ", loss_name, ".", call. = FALSE)
    }
    get(loss_name, mode = "function", inherits = TRUE)
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
  if (effective_overlap_size <= d_max) {
    warning(
      "The effective overlap size is no larger than max(d_candidates); ",
      "translated Procrustes alignment may be rank-deficient or unstable.",
      call. = FALSE,
      immediate. = TRUE
    )
  }

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
  candidate_task_lookup <- array(
    seq_len(nrow(candidate_grid)),
    dim = c(length(d_candidates), num_subnetworks, nrep)
  )
  fit_ncores <- min(ncores, nrow(candidate_grid))
  alignment_ncores <- min(ncores, nrow(candidate_grid))
  loss_ncores <- min(ncores, nrow(loss_grid))
  loss_accepts_prevalidation <- vapply(
    loss_functions,
    function(fun) "validate_inputs" %in% names(formals(fun)),
    logical(1)
  )

  manage_future_plan <- TRUE
  os_type <- if (force_windows) "windows" else .Platform$OS.type
  worker_future_packages <- if (inherits(A, "Matrix")) "Matrix" else NULL
  maximum_stage_ncores <- max(fit_ncores, alignment_ncores, loss_ncores)
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
    on.exit(try(future::plan(previous_future_plan), silent = TRUE), add = TRUE)
    on.exit({
      if (is.na(previous_renv_sync_check)) {
        Sys.unsetenv("RENV_CONFIG_SYNCHRONIZED_CHECK")
      } else {
        Sys.setenv(RENV_CONFIG_SYNCHRONIZED_CHECK = previous_renv_sync_check)
      }
    }, add = TRUE)
    future::plan(future::multisession, workers = maximum_stage_ncores)
    manage_future_plan <- FALSE
  }

  ram_report <- NULL
  if (ram_check) {
    usvt_dimension <- min(ceiling(subgraph_size^(1 / 3)), subgraph_size)
    initialization_ram <- max(
      estimate_spectral_decomp_ram(
        n = subgraph_size,
        p = subgraph_size,
        K = usvt_dimension,
        method = "svd",
        engine = "irlba",
        nu = usvt_dimension,
        nv = usvt_dimension,
        dense_input = TRUE
      )$estimated_bytes,
      estimate_spectral_decomp_ram(
        n = subgraph_size,
        p = subgraph_size,
        K = d_max,
        method = "eigen",
        engine = "rspectra",
        nu = d_max,
        nv = d_max,
        dense_input = TRUE
      )$estimated_bytes
    )
    dense_fit_state_ram <- 1.25 * 8 * 6 * as.double(subgraph_size)^2
    fit_product_ram <- estimate_matrix_product_ram(
      subgraph_size,
      d_max,
      subgraph_size
    )
    alignment_ram <- estimate_matrix_product_ram(
      piece_size,
      d_max,
      d_max
    )
    prediction_ram <- estimate_matrix_product_ram(
      piece_size,
      d_max,
      piece_size
    )
    ram_report <- report_ram_formula(
      terms = list(
        list(
          estimated_bytes = initialization_ram + dense_fit_state_ram +
            fit_product_ram,
          sequential_count = 1L,
          parallel_count = fit_ncores,
          label = "dense subnetwork LSM fit"
        ),
        list(
          estimated_bytes = alignment_ram,
          sequential_count = 1L,
          parallel_count = alignment_ncores,
          label = "non-overlap translated alignment"
        ),
        list(
          estimated_bytes = prediction_ram,
          sequential_count = 1L,
          parallel_count = loss_ncores,
          label = "held-out probability construction"
        )
      ),
      operation = "LSM NETCROP conservative combined preflight",
      detail = "three sequential parallel stages; peak is overestimated additively"
    )
  }

  if (verbose) {
    message(
      "Stage 1/3: fitting ", nrow(candidate_grid),
      " candidate subnetworks with ", fit_ncores, " worker(s)."
    )
  }
  fit_time <- system.time({
    fit_output <- uni_mclapply(
      seq_len(nrow(candidate_grid)),
      function(task_id) {
        repetition <- candidate_grid$repetition[task_id]
        split_id <- candidate_grid$split[task_id]
        candidate_id <- candidate_grid$candidate_id[task_id]
        d <- d_candidates[candidate_id]
        if (!is.null(seed)) {
          task_seed <- as.integer(
            (seed + 100000 + task_id) %% .Machine$integer.max
          )
          set.seed(task_seed)
        }
        nodes <- splitter$subnetworks[[repetition]][[split_id]]
        A_subnetwork <- A[nodes, nodes, drop = FALSE]
        fit <- do.call(
          lsm_pgd,
          c(
            list(A = A_subnetwork, d = d, ram_check = FALSE),
            resolved_lsm_options
          )
        )
        X_hat <- fit$Z_hat
        Z_hat <- X_hat / sqrt(2)
        gamma_hat <- as.numeric(fit$alpha_hat) + rowSums(X_hat^2) / 2
        list(
          X_hat = X_hat,
          beta_hat = as.numeric(fit$alpha_hat),
          Z_hat = Z_hat,
          gamma_hat = gamma_hat,
          objective = fit$objective,
          step_size_Z = fit$step_size_Z,
          step_size_alpha = fit$step_size_alpha
        )
      },
      ncores = fit_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      manage_future_plan = manage_future_plan,
      future_packages = worker_future_packages
    )
  })

  if (verbose) {
    message(
      "Stage 2/3: aligning ", nrow(candidate_grid),
      " fitted subnetworks with ", alignment_ncores, " worker(s)."
    )
  }
  alignment_time <- system.time({
    aligned_output <- uni_mclapply(
      seq_len(nrow(candidate_grid)),
      function(task_id) {
        repetition <- candidate_grid$repetition[task_id]
        split_id <- candidate_grid$split[task_id]
        candidate_id <- candidate_grid$candidate_id[task_id]
        fit <- fit_output[[task_id]]
        overlap_indices <- seq_len(effective_overlap_size)
        nonoverlap_indices <- effective_overlap_size + seq_len(piece_size)
        if (split_id == 1L) {
          return(list(
            Z_aligned = fit$Z_hat[nonoverlap_indices, , drop = FALSE],
            gamma_hat = fit$gamma_hat[nonoverlap_indices],
            rotation = diag(d_candidates[candidate_id]),
            translation = rep(0, d_candidates[candidate_id])
          ))
        }
        standard_id <- candidate_task_lookup[candidate_id, 1L, repetition]
        alignment <- procrustes(
          X = fit$Z_hat[overlap_indices, , drop = FALSE],
          X_star = fit_output[[standard_id]]$Z_hat[
            overlap_indices, , drop = FALSE
          ],
          translate = TRUE,
          dilate = FALSE,
          validate_inputs = FALSE
        )
        Z_rotated <- tcrossprod(
          fit$Z_hat[nonoverlap_indices, , drop = FALSE],
          t(alignment$rotation)
        )
        Z_aligned <- sweep(
          Z_rotated,
          2L,
          alignment$translation,
          `+`
        )

        # The original inner-product alignment adjusted beta after translation:
        # alpha_aligned <- fit$beta_hat[nonoverlap_indices] +
        #   sum(alignment$translation^2) / 2 +
        #   drop(Z_rotated %*% alignment$translation)

        # After the exact distance reparameterization, gamma_i is invariant
        # under rotation and translation, and rigid alignment preserves all
        # pairwise distances. Therefore gamma must not be adjusted here.
        gamma_aligned <- fit$gamma_hat[nonoverlap_indices]

        list(
          Z_aligned = Z_aligned,
          gamma_hat = gamma_aligned,
          rotation = alignment$rotation,
          translation = alignment$translation
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
    fit_output <- NULL
  }

  pair_count <- nrow(pair_grid)
  if (verbose) {
    message(
      "Stage 3/3: evaluating ", nrow(loss_grid),
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
        d <- d_candidates[candidate_id]
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
        left_fit <- aligned_output[[left_id]]
        right_fit <- aligned_output[[right_id]]
        left_norms <- rowSums(left_fit$Z_aligned^2)
        right_norms <- rowSums(right_fit$Z_aligned^2)
        distance_squared <- outer(left_norms, right_norms, FUN = "+") -
          2 * tcrossprod(left_fit$Z_aligned, right_fit$Z_aligned)
        distance_squared <- pmax(distance_squared, 0)
        theta_hat <- outer(
          left_fit$gamma_hat,
          right_fit$gamma_hat,
          FUN = "+"
        ) - distance_squared
        P_hat <- stats::plogis(theta_hat)
        A_test_numeric <- as.numeric(A_test)
        P_hat_numeric <- as.numeric(P_hat)

        # Preserve the pasted penalty but use d itself rather than the
        # candidate's position in d_candidates.
        penalty <- 2 * (effective_overlap_size + piece_size + d)
        records <- matrix(NA_real_, nrow = length(losses), ncol = 8L)
        for (loss_id in seq_along(losses)) {
          loss_name <- losses[loss_id]
          loss_arguments <- list(A_test_numeric, P_hat_numeric)
          if (loss_accepts_prevalidation[loss_id]) {
            loss_arguments$validate_inputs <- FALSE
          }
          raw_loss <- do.call(loss_functions[[loss_name]], loss_arguments)
          value <- (raw_loss + penalty) / pair_count
          if (length(value) != 1L || !is.numeric(value) ||
              is.na(value) || !is.finite(value)) {
            stop(
              "Non-finite loss for repetition ", repetition,
              ", d = ", d,
              ", loss = ", loss_name,
              ", split pair = ", split_left, "-", split_right,
              ". For AUC-based loss, a held-out block may contain only one ",
              "edge class.",
              call. = FALSE
            )
          }
          records[loss_id, ] <- c(
            repetition,
            d,
            loss_id,
            split_left,
            split_right,
            as.numeric(raw_loss),
            penalty,
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
    raw_loss = loss_matrix[, 6L],
    penalty = loss_matrix[, 7L],
    loss_value = loss_matrix[, 8L],
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
    lsm_options = resolved_lsm_options,
    splitter = splitter,
    candidate_grid = candidate_grid,
    loss_grid = loss_grid,
    fit_output = if (retain_intermediates == "all") fit_output else NULL,
    aligned_output = if (retain_intermediates == "all") {
      aligned_output
    } else {
      NULL
    },
    retain_intermediates = retain_intermediates,
    ncores = list(
      requested = ncores,
      fit = fit_ncores,
      alignment = alignment_ncores,
      loss = loss_ncores
    ),
    timing = c(
      fit = unname(fit_time["elapsed"]),
      alignment = unname(alignment_time["elapsed"]),
      loss = unname(loss_time["elapsed"]),
      total = unname(
        fit_time["elapsed"] + alignment_time["elapsed"] +
          loss_time["elapsed"]
      )
    ),
    ram_preflight = ram_report,
    seed = seed,
    call = call
  )
  class(out) <- "netcrop_lsm_reparameterized"
  out
}

# Print the selected LSM dimension for every requested loss.
print.netcrop_lsm_reparameterized <- function(x, ...) {
  cat("NETCROP results for latent-space models\n")
  cat("---------------------------------------\n")
  print(x$overall_best, row.names = FALSE)
  invisible(x)
}

# Summarize an LSM NETCROP fit.
summary.netcrop_lsm_reparameterized <- function(object, ...) {
  result <- list(
    call = object$call,
    d_candidates = object$d_candidates,
    nrep = object$nrep,
    num_subnetworks = object$num_subnetworks,
    requested_overlap_size = object$requested_overlap_size,
    effective_overlap_size = object$effective_overlap_size,
    piece_size = object$piece_size,
    unordered_test_proportion = object$unordered_test_proportion,
    best_dimension_cv = object$best_dimension_cv,
    overall_best = object$overall_best,
    lsm_options = object$lsm_options,
    ncores = object$ncores,
    timing = object$timing
  )
  class(result) <- "summary.netcrop_lsm_reparameterized"
  result
}

# Print an LSM NETCROP summary.
print.summary.netcrop_lsm_reparameterized <- function(x, ...) {
  cat("Summary of NETCROP latent-space-model dimension selection\n")
  cat("---------------------------------------------------------\n")
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
  cat(sprintf(
    "Held-out unordered-pair proportion: %.2f%%\n",
    100 * x$unordered_test_proportion
  ))
  cat("Best dimensions per repetition:\n")
  print(x$best_dimension_cv, row.names = FALSE)
  cat("Overall best dimensions:\n")
  print(x$overall_best, row.names = FALSE)
  cat("Timing (seconds):\n")
  print(x$timing)
  invisible(x)
}

# Plot LSM CV loss curves, optionally aggregating across repetitions.
plot.netcrop_lsm_reparameterized <- function(x, aggregate = TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting.", call. = FALSE)
  }
  if (length(aggregate) != 1L || !is.logical(aggregate) || is.na(aggregate)) {
    stop("aggregate must be TRUE or FALSE.", call. = FALSE)
  }
  d_breaks <- sort(unique(x$cv_loss$d))
  if (!aggregate) {
    plot_data <- x$cv_loss
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
          title = "NETCROP CV loss by LSM dimension",
          x = "Latent dimension (d)",
          y = "Penalized average CV loss",
          color = "Repetition"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")
    )
  }
  grouping <- interaction(
    x$cv_loss$d,
    x$cv_loss$loss,
    drop = TRUE,
    lex.order = TRUE
  )
  plot_data <- do.call(rbind, lapply(split(x$cv_loss, grouping), function(z) {
    standard_deviation <- stats::sd(z$average_loss)
    if (is.na(standard_deviation)) standard_deviation <- 0
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
      title = "NETCROP CV loss by LSM dimension",
      x = "Latent dimension (d)",
      y = "Mean penalized CV loss (plus or minus one SD)"
    ) +
    ggplot2::theme_minimal()
}



# Benchmark code (intentionally not run when this file is sourced).
# library(Matrix)
# system.time(
#   net <- generate_lsm(
#     n = 1e3, d = 5, alpha = 0
#   )
# )
#
# system.time(
#   nc_out <- netcrop_lsm_reparameterized(
#     A = net,
#     d_candidates = 1:10,
#     num_subnetworks = NULL,
#     overlap_size = NULL,
#     nrep = 1L,
#     lsm_options = list(niter = 100L),
#     losses = "sse",
#     ncores = 8L,
#     force_windows = FALSE
#   )
# )
# nc_out
# summary(nc_out)
# plot(nc_out)
