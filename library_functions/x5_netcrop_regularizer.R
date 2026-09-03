# NETCROP regularization-parameter selection for SBM and DCBM
#
# Source 01_basic_helpers.R through 06_estimators.R and x0_helpers.R first.
# This implementation uses observed-network cross-validation only.

# Compare two pairs of community-label vectors using negative NMI.
pair_nmi_loss <- function(g_1, g_2, g_3, g_4) {
  vectors <- list(g_1 = g_1, g_2 = g_2, g_3 = g_3, g_4 = g_4)
  invalid <- vapply(vectors, function(g) {
    !is.numeric(g) || length(g) < 1L || anyNA(g) ||
      any(!is.finite(g)) || any(g < 1) || any(g != floor(g))
  }, logical(1))
  if (any(invalid) || length(g_1) * length(g_2) !=
      length(g_3) * length(g_4)) {
    stop("The two label-pair Cartesian products must have matching sizes.",
         call. = FALSE)
  }
  labels_1 <- (rep(g_1, times = length(g_2)) - 1L) * max(g_2) +
    rep(g_2, each = length(g_1))
  labels_2 <- (rep(g_3, times = length(g_4)) - 1L) * max(g_4) +
    rep(g_4, each = length(g_3))
  if (length(unique(labels_1)) == 1L && length(unique(labels_2)) == 1L) {
    return(-1)
  }
  value <- nmi(labels_1, labels_2)
  if (is.finite(value)) -value else NA_real_
}

# Compare two pairs of community-label vectors using matched Hamming loss.
pair_hamming_loss <- function(g_1, g_2, g_3, g_4) {
  vectors <- list(g_1 = g_1, g_2 = g_2, g_3 = g_3, g_4 = g_4)
  invalid <- vapply(vectors, function(g) {
    !is.numeric(g) || length(g) < 1L || anyNA(g) ||
      any(!is.finite(g)) || any(g < 1) || any(g != floor(g))
  }, logical(1))
  if (any(invalid) || length(g_1) * length(g_2) !=
      length(g_3) * length(g_4)) {
    stop("The two label-pair Cartesian products must have matching sizes.",
         call. = FALSE)
  }
  labels_1 <- (rep(g_1, times = length(g_2)) - 1L) * max(g_2) +
    rep(g_2, each = length(g_1))
  labels_2 <- (rep(g_3, times = length(g_4)) - 1L) * max(g_4) +
    rep(g_4, each = length(g_3))
  label_match_greedy(
    labels_1,
    labels_2,
    K = max(c(labels_1, labels_2)),
    algorithm = "hungarian"
  )$mismatch_rate
}

# Select a spectral regularization parameter by overlapping-subnetwork CV.
netcrop_tune_regularizer <- function(
    A,
    K,
    tau_candidates,
    use_dcbm = FALSE,
    num_subnetworks = NULL,
    overlap_size = NULL,
    nrep = 1L,
    use_laplacian = FALSE,
    dcbm_est_method = c("plugin", "spectral"),
    losses = "sse",
    loss_types = NULL,
    label_reference = c("full_network", "leave_pair_out"),
    spectral_options = list(),
    cluster_options = list(),
    estimator_options = list(),
    matching_method = c("greedy", "hungarian", "brute_force"),
    confirm_large = NULL,
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
  dcbm_est_method <- match.arg(dcbm_est_method)
  matching_method <- match.arg(matching_method)
  label_reference <- match.arg(label_reference)
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
  validate_named_list <- function(value, name) {
    if (!is.list(value) || (length(value) > 0L &&
                            (is.null(names(value)) || any(!nzchar(names(value)))))) {
      stop(name, " must be a named list.", call. = FALSE)
    }
    invisible(TRUE)
  }

  K <- validate_count(K, "K")
  if (K == 1L) {
    message("K = 1; cross-validation is not necessary.")
    return(invisible(NULL))
  }
  required_helpers <- c(
    "uni_mclapply", "netcrop_param_select", "netcrop_splitter",
    "is_symmetric_matrix", "spectral_cluster", "estimate_sbm",
    "estimate_dcbm", "label_match_greedy", "label_match_brute_force",
    "sse", "bin_dev", "auc_as_loss", "nmi"
  )
  missing_helpers <- required_helpers[!vapply(
    required_helpers, exists, logical(1), mode = "function", inherits = TRUE
  )]
  if (length(missing_helpers) > 0L) {
    stop(
      "Source the numbered helpers and x0_helpers.R first. Missing: ",
      paste(missing_helpers, collapse = ", "), ".",
      call. = FALSE
    )
  }
  required_packages <- c("RSpectra", "tibble")
  if (isTRUE(use_dcbm)) required_packages <- c(required_packages, "cluster")
  if (inherits(A, "Matrix")) required_packages <- c(required_packages, "Matrix")
  missing_packages <- unique(required_packages)[!vapply(
    unique(required_packages), requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing_packages) > 0L) {
    stop("Required package(s) are not installed: ",
         paste(missing_packages, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(dim(A)) || length(dim(A)) != 2L || nrow(A) != ncol(A) ||
      nrow(A) < 3L || !(is.numeric(A) || inherits(A, "Matrix"))) {
    stop("A must be a numeric square matrix with at least three nodes.",
         call. = FALSE)
  }
  A_values <- if (inherits(A, "sparseMatrix") &&
                  "x" %in% methods::slotNames(A)) {
    methods::slot(A, "x")
  } else {
    as.numeric(A)
  }
  if (anyNA(A_values) || any(!is.finite(A_values)) || any(A_values < 0)) {
    stop("A must contain only finite non-negative values.", call. = FALSE)
  }
  n <- nrow(A)
  nrep <- validate_count(nrep, "nrep")
  ncores <- validate_count(ncores, "ncores")
  use_dcbm <- validate_flag(use_dcbm, "use_dcbm")
  use_laplacian <- validate_flag(use_laplacian, "use_laplacian")
  verbose <- validate_flag(verbose, "verbose")
  force_windows <- validate_flag(force_windows, "force_windows")
  ram_check <- validate_flag(ram_check, "ram_check")
  if (!is.null(confirm_large)) {
    confirm_large <- validate_flag(confirm_large, "confirm_large")
  }
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
        !is.finite(seed) || seed < 0 || seed != floor(seed) ||
        seed > .Machine$integer.max) {
      stop("seed must be NULL or a non-negative integer within R's range.",
           call. = FALSE)
    }
    seed <- as.double(seed)
  }
  if (!is.numeric(tau_candidates) || length(tau_candidates) < 1L ||
      anyNA(tau_candidates) || any(!is.finite(tau_candidates)) ||
      any(tau_candidates < 0) || anyDuplicated(tau_candidates)) {
    stop("tau_candidates must contain unique finite non-negative numbers.",
         call. = FALSE)
  }
  tau_candidates <- as.numeric(tau_candidates)
  if (!is.character(losses) || length(losses) < 1L || anyNA(losses) ||
      any(!nzchar(losses)) || anyDuplicated(losses)) {
    stop("losses must contain unique, non-empty function names.",
         call. = FALSE)
  }
  validate_named_list(spectral_options, "spectral_options")
  validate_named_list(cluster_options, "cluster_options")
  validate_named_list(estimator_options, "estimator_options")
  validate_named_list(parameter_select_options, "parameter_select_options")

  if (!is_symmetric_matrix(A)) {
    warning("A is asymmetric; using the undirected binarized union A + t(A).",
            call. = FALSE)
    A_nonzero <- A != 0
    A <- (A_nonzero + if (inherits(A_nonzero, "Matrix")) {
      Matrix::t(A_nonzero)
    } else {
      t(A_nonzero)
    }) > 0
  } else if (!all(A_values %in% c(0, 1))) {
    warning("A is weighted; binarizing its nonzero entries.", call. = FALSE)
    A <- (A != 0) * 1
  } else if (is.logical(A)) {
    A <- A * 1
  }
  if (!is.numeric(A)) A <- A * 1
  if (inherits(A, "symmetricMatrix")) {
    A <- methods::as(A, "generalMatrix")
  }
  A_diagonal <- if (inherits(A, "Matrix")) Matrix::diag(A) else diag(A)
  if (any(A_diagonal != 0)) {
    warning("Diagonal entries of A were set to zero.", call. = FALSE)
    if (inherits(A, "Matrix")) {
      A <- A - Matrix::Diagonal(x = A_diagonal)
    } else {
      diag(A) <- 0
    }
  }
  input_was_sparse <- inherits(A, "sparseMatrix")

  if (is.null(num_subnetworks) || is.null(overlap_size)) {
    if ("n" %in% names(parameter_select_options)) {
      stop("parameter_select_options cannot override n.", call. = FALSE)
    }
    unknown <- setdiff(names(parameter_select_options),
                       names(formals(netcrop_param_select)))
    if (length(unknown) > 0L) {
      stop("Unsupported parameter_select_options: ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    selected <- do.call(
      netcrop_param_select,
      utils::modifyList(
        list(test_prop = 0.02, n = n, o_range = 0),
        parameter_select_options,
        keep.null = TRUE
      )
    )
    if (nrow(selected) != 1L) {
      stop("Parameter selection must return exactly one split pair.",
           call. = FALSE)
    }
    if (is.null(num_subnetworks)) {
      num_subnetworks <- selected$num_subnetworks[[1L]]
    }
    if (is.null(overlap_size)) overlap_size <- selected$overlap_size[[1L]]
  }
  num_subnetworks <- validate_count(num_subnetworks, "num_subnetworks", 2L)
  overlap_size <- validate_count(overlap_size, "overlap_size")
  if (overlap_size >= n ||
      floor((n - overlap_size) / num_subnetworks) < 1L) {
    stop("The split must retain at least one node per non-overlap piece.",
         call. = FALSE)
  }
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
  if (K >= subgraph_size) {
    stop("K must be smaller than the effective subgraph size.", call. = FALSE)
  }

  probability_loss_names <- c("sse", "bin_dev", "auc_as_loss")
  label_loss_names <- c("pair_nmi_loss", "pair_hamming_loss")
  if (is.null(loss_types)) {
    unknown_losses <- setdiff(losses,
                              c(probability_loss_names, label_loss_names))
    if (length(unknown_losses) > 0L) {
      stop("loss_types must classify custom loss(es): ",
           paste(unknown_losses, collapse = ", "), ".", call. = FALSE)
    }
    loss_types <- ifelse(losses %in% probability_loss_names,
                         "probability", "label")
    names(loss_types) <- losses
  } else {
    if (!is.character(loss_types) || is.null(names(loss_types)) ||
        anyDuplicated(names(loss_types)) || !setequal(names(loss_types), losses) ||
        anyNA(loss_types) ||
        any(!loss_types %in% c("probability", "label"))) {
      stop("loss_types must classify every loss as 'probability' or 'label'.",
           call. = FALSE)
    }
    loss_types <- loss_types[losses]
  }
  loss_functions <- lapply(losses, function(loss_name) {
    if (!exists(loss_name, mode = "function", inherits = TRUE)) {
      stop("Loss function not found: ", loss_name, ".", call. = FALSE)
    }
    get(loss_name, mode = "function", inherits = TRUE)
  })
  names(loss_functions) <- losses
  has_probability_loss <- any(loss_types == "probability")
  has_label_loss <- any(loss_types == "label")
  use_full_label_reference <- has_label_loss &&
    label_reference == "full_network"

  protected_spectral <- intersect(
    names(spectral_options),
    c("A", "K", "regularize_tau", "row_normalize", "cluster_engine",
      "cluster_options", "ram_check", "validate_inputs")
  )
  if (length(protected_spectral) > 0L) {
    stop("spectral_options cannot override: ",
         paste(protected_spectral, collapse = ", "), ".", call. = FALSE)
  }
  unsupported_spectral <- setdiff(
    names(spectral_options), names(formals(spectral_cluster))
  )
  cluster_function <- if (use_dcbm) cluster::clara else stats::kmeans
  unsupported_cluster <- setdiff(
    names(cluster_options),
    setdiff(names(formals(cluster_function)), c("x", "k", "centers"))
  )
  if (length(unsupported_spectral) > 0L ||
      length(unsupported_cluster) > 0L) {
    stop(
      "Unsupported spectral or clustering option(s): ",
      paste(c(unsupported_spectral, unsupported_cluster), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  estimator_function <- if (use_dcbm) estimate_dcbm else estimate_sbm
  protected_estimator <- intersect(
    names(estimator_options),
    c("A", "g", "K", "method", "psi_omit", "validate_inputs")
  )
  unsupported_estimator <- setdiff(names(estimator_options),
                                   names(formals(estimator_function)))
  if (length(protected_estimator) > 0L ||
      length(unsupported_estimator) > 0L) {
    stop("Invalid estimator_options: ",
         paste(unique(c(protected_estimator, unsupported_estimator)),
               collapse = ", "), ".", call. = FALSE)
  }
  if (matching_method == "brute_force" && K > 8L &&
      !isTRUE(confirm_large)) {
    stop("Brute-force matching above K = 8 requires confirm_large = TRUE.",
         call. = FALSE)
  }

  match_labels <- function(match_this, standard) {
    if (matching_method == "brute_force") {
      return(label_match_brute_force(
        match_this, standard, K = K, return_mapping = TRUE,
        confirm_large = confirm_large
      ))
    }
    label_match_greedy(
      match_this, standard, K = K, algorithm = matching_method,
      return_mapping = TRUE
    )
  }
  set_task_seed <- function(stage, task_idx) {
    if (is.null(seed)) return(invisible(NULL))
    task_seed <- as.integer((seed + as.double(stage) * 1000003 + task_idx) %%
                              .Machine$integer.max)
    set.seed(task_seed)
    invisible(task_seed)
  }
  fit_labels <- function(A_work, tau, task_seed = NULL) {
    if (!is.null(task_seed)) set.seed(task_seed)
    tryCatch({
      resolved_spectral_options <- spectral_options
      if (inherits(A_work, "sparseMatrix") && tau == 0) {
        resolved_spectral_options$force_engine <- TRUE
      }
      cluster_defaults <- if (use_dcbm) {
        list(metric = "manhattan", cluster.only = TRUE, samples = 5L)
      } else {
        list(nstart = 100L, iter.max = 10^7)
      }
      fit_arguments <- c(
        list(
          A = A_work, K = K, laplacian = use_laplacian,
          normalize_laplacian = TRUE, regularize_tau = tau,
          handle_zero_degree_nodes = if (tau == 0) "random_label" else "none",
          row_normalize = use_dcbm, spectral_engine = "RSpectra",
          cluster_engine = if (use_dcbm) "clara" else "kmeans",
          cluster_options = utils::modifyList(
            cluster_defaults, cluster_options, keep.null = TRUE
          ),
          ram_check = FALSE, validate_inputs = FALSE
        ),
        resolved_spectral_options
      )
      fit <- if (verbose) {
        do.call(spectral_cluster, fit_arguments)
      } else {
        suppressMessages(do.call(spectral_cluster, fit_arguments))
      }
      list(success = TRUE, labels = as.integer(fit$g_hat), error = NULL)
    }, error = function(error) {
      list(success = FALSE, labels = NULL, error = conditionMessage(error))
    })
  }

  pair_count <- choose(num_subnetworks, 2L)
  fit_grid <- expand.grid(
    repetition = seq_len(nrep), split = seq_len(num_subnetworks),
    tau_idx = seq_along(tau_candidates), KEEP.OUT.ATTRS = FALSE
  )
  fit_grid <- fit_grid[order(fit_grid$repetition, fit_grid$tau_idx,
                             fit_grid$split), ]
  pair_grid <- do.call(rbind, lapply(seq_len(nrep), function(rep_idx) {
    do.call(rbind, lapply(seq_along(tau_candidates), function(tau_idx) {
      pairs <- utils::combn(num_subnetworks, 2L)
      data.frame(repetition = rep_idx, tau_idx = tau_idx,
                 piece_1 = pairs[1L, ], piece_2 = pairs[2L, ])
    }))
  }))
  fit_workers <- min(nrow(fit_grid), ncores)
  loss_workers <- min(nrow(pair_grid), ncores)
  reference_workers <- if (use_full_label_reference) {
    min(length(tau_candidates), ncores)
  } else {
    1L
  }

  manage_future_plan <- TRUE
  os_type <- if (force_windows) "windows" else .Platform$OS.type
  worker_future_packages <- if (input_was_sparse) "Matrix" else NULL
  if (os_type != "unix" &&
      max(fit_workers, reference_workers, loss_workers) > 1L) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("future and future.apply are required for Windows parallel execution.",
           call. = FALSE)
    }
    previous_plan <- future::plan()
    previous_renv_check <- Sys.getenv("RENV_CONFIG_SYNCHRONIZED_CHECK",
                                      unset = NA_character_)
    Sys.setenv(RENV_CONFIG_SYNCHRONIZED_CHECK = "FALSE")
    on.exit(try(future::plan(previous_plan), silent = TRUE), add = TRUE)
    on.exit({
      if (is.na(previous_renv_check)) {
        Sys.unsetenv("RENV_CONFIG_SYNCHRONIZED_CHECK")
      } else {
        Sys.setenv(RENV_CONFIG_SYNCHRONIZED_CHECK = previous_renv_check)
      }
    }, add = TRUE)
    future::plan(
      future::multisession,
      workers = max(fit_workers, reference_workers, loss_workers)
    )
    manage_future_plan <- FALSE
  }

  ram_report <- NULL
  if (ram_check) {
    ram_helpers <- c("estimate_spectral_decomp_ram",
                     "estimate_matrix_product_ram", "report_ram_formula")
    missing_ram <- ram_helpers[!vapply(
      ram_helpers, exists, logical(1), mode = "function", inherits = TRUE
    )]
    if (length(missing_ram) > 0L) {
      stop("Missing RAM helper(s): ", paste(missing_ram, collapse = ", "),
           ".", call. = FALSE)
    }
    subgraph_ram <- estimate_spectral_decomp_ram(
      n = subgraph_size, p = subgraph_size, K = K, method = "eigen",
      engine = "rspectra", dense_input = any(tau_candidates > 0)
    )$estimated_bytes
    heldout_ram <- estimate_matrix_product_ram(piece_size, 1L, piece_size)
    ram_terms <- list(
      list(estimated_bytes = subgraph_ram, sequential_count = 1L,
           parallel_count = fit_workers, label = "subgraph spectral fit"),
      list(estimated_bytes = heldout_ram,
           sequential_count = if (has_probability_loss) 3L else 1L,
           parallel_count = loss_workers,
           label = "held-out adjacency and probability matrices")
    )
    if (has_label_loss) {
      reference_ram <- estimate_spectral_decomp_ram(
        n = n, p = n, K = K, method = "eigen", engine = "rspectra",
        dense_input = any(tau_candidates > 0)
      )$estimated_bytes + 8 * as.double(n)^2
      ram_terms <- c(ram_terms, list(list(
        estimated_bytes = reference_ram,
        sequential_count = 1L,
        parallel_count = if (use_full_label_reference) {
          reference_workers
        } else {
          loss_workers
        },
        label = if (use_full_label_reference) {
          "full-network label reference fit"
        } else {
          "leave-pair-out label reference fit"
        }
      )))
    }
    ram_report <- report_ram_formula(
      terms = ram_terms,
      operation = "NETCROP regularizer conservative preflight",
      detail = "stage-specific concurrent matrix allocations"
    )
  }

  if (verbose) {
    message("Stage 1/3: fitting ", nrow(fit_grid),
            " subnetwork/tau combinations with ", fit_workers, " worker(s).")
  }
  fit_time <- system.time({
    raw_fits <- uni_mclapply(
      seq_len(nrow(fit_grid)),
      function(task_idx) {
        set_task_seed(1L, task_idx)
        rep_idx <- fit_grid$repetition[task_idx]
        split_idx <- fit_grid$split[task_idx]
        tau_idx <- fit_grid$tau_idx[task_idx]
        nodes <- splitter$subnetworks[[rep_idx]][[split_idx]]
        fit <- fit_labels(A[nodes, nodes, drop = FALSE],
                          tau_candidates[tau_idx])
        c(fit, list(nodes = nodes))
      },
      ncores = fit_workers,
      force_windows = force_windows,
      stop_on_error = TRUE,
      manage_future_plan = manage_future_plan,
      future_packages = worker_future_packages
    )
  })
  fit_lookup <- array(
    NA_integer_,
    dim = c(nrep, length(tau_candidates), num_subnetworks)
  )
  for (task_idx in seq_len(nrow(fit_grid))) {
    fit_lookup[fit_grid$repetition[task_idx], fit_grid$tau_idx[task_idx],
               fit_grid$split[task_idx]] <- task_idx
  }

  if (verbose) {
    message("Stage 2/3: aligning labels and estimating model parameters.")
  }
  estimate_time <- system.time({
    estimates <- vector("list", nrow(fit_grid))
    matching_records <- vector("list", nrow(fit_grid))
    for (task_idx in seq_len(nrow(fit_grid))) {
      rep_idx <- fit_grid$repetition[task_idx]
      tau_idx <- fit_grid$tau_idx[task_idx]
      split_idx <- fit_grid$split[task_idx]
      fit <- raw_fits[[task_idx]]
      if (!fit$success) {
        estimates[[task_idx]] <- list(success = FALSE, error = fit$error)
        next
      }
      labels <- fit$labels
      underidentified <- FALSE
      mismatch_rate <- 0
      if (split_idx != 1L) {
        reference_fit <- raw_fits[[fit_lookup[rep_idx, tau_idx, 1L]]]
        if (!reference_fit$success) {
          estimates[[task_idx]] <- list(
            success = FALSE,
            error = paste("Reference split failed:", reference_fit$error)
          )
          next
        }
        overlap <- splitter$overlap_nodes[[rep_idx]]
        local_positions <- match(overlap, fit$nodes)
        reference_positions <- match(overlap, reference_fit$nodes)
        local_labels <- labels[local_positions]
        reference_labels <- reference_fit$labels[reference_positions]
        underidentified <- length(unique(local_labels)) < K ||
          length(unique(reference_labels)) < K
        matched <- match_labels(local_labels, reference_labels)
        labels <- unname(matched$mapping[labels])
        mismatch_rate <- matched$mismatch_rate
      }
      piece_nodes <- splitter$nonoverlap_pieces[[rep_idx]][[split_idx]]
      piece_positions <- match(piece_nodes, fit$nodes)
      B_hat <- NULL
      psi_hat <- NULL
      if (has_probability_loss) {
        estimation_error <- NULL
        model_fit <- tryCatch({
          if (use_dcbm) {
            do.call(
              estimate_dcbm,
              c(list(
                A = A[fit$nodes, fit$nodes, drop = FALSE],
                g = labels,
                K = K,
                method = dcbm_est_method,
                psi_omit = effective_overlap_size,
                validate_inputs = FALSE
              ), estimator_options)
            )
          } else {
            list(
              B_hat = do.call(
                estimate_sbm,
                c(list(
                  A = A[fit$nodes, fit$nodes, drop = FALSE],
                  g = labels,
                  K = K,
                  validate_inputs = FALSE
                ), estimator_options)
              ),
              psi_hat = rep(1, piece_size)
            )
          }
        }, error = function(error) {
          estimation_error <<- conditionMessage(error)
          NULL
        })
        if (is.null(model_fit)) {
          estimates[[task_idx]] <- list(
            success = FALSE,
            error = paste("Estimation failed:", estimation_error)
          )
          next
        }
        B_hat <- model_fit$B_hat
        psi_hat <- model_fit$psi_hat
      }
      estimates[[task_idx]] <- list(
        success = TRUE,
        labels = labels[piece_positions],
        B_hat = B_hat,
        psi_hat = psi_hat,
        nodes = piece_nodes,
        error = NULL
      )
      matching_records[[task_idx]] <- data.frame(
        repetition = rep_idx,
        tau = tau_candidates[tau_idx],
        split = split_idx,
        mismatch_rate = mismatch_rate,
        underidentified = underidentified
      )
    }
  })

  average_B_hat <- array(
    vector("list", nrep * length(tau_candidates)),
    dim = c(nrep, length(tau_candidates))
  )
  if (has_probability_loss) {
    for (rep_idx in seq_len(nrep)) {
      for (tau_idx in seq_along(tau_candidates)) {
        selected_estimates <- estimates[fit_lookup[rep_idx, tau_idx, ]]
        if (all(vapply(selected_estimates, `[[`, logical(1), "success"))) {
          matrices <- lapply(selected_estimates, `[[`, "B_hat")
          if (all(vapply(matrices, function(B_hat) {
            is.matrix(B_hat) && all(is.finite(B_hat))
          }, logical(1)))) {
            average_B_hat[[rep_idx, tau_idx]] <-
              Reduce(`+`, matrices) / num_subnetworks
          }
        }
      }
    }
  }

  reference_time <- c(elapsed = 0)
  complete_reference_fits <- NULL
  if (use_full_label_reference) {
    if (verbose) {
      message("Reference stage: fitting the complete network at each tau.")
    }
    reference_time <- system.time({
      complete_reference_fits <- uni_mclapply(
        seq_along(tau_candidates),
        function(tau_idx) {
          set_task_seed(2L, tau_idx)
          fit_labels(A, tau_candidates[tau_idx])
        },
        ncores = reference_workers,
        force_windows = force_windows,
        stop_on_error = TRUE,
        manage_future_plan = manage_future_plan,
        future_packages = worker_future_packages
      )
    })
  }

  evaluate_probability_loss <- function(loss_name, A_test, P_hat) {
    A_vec <- as.numeric(A_test)
    P_vec <- as.numeric(P_hat)
    if (loss_name == "sse") return(sum((A_vec - P_vec)^2))
    if (loss_name == "bin_dev") {
      P_vec <- pmin(pmax(P_vec, 1e-10), 1 - 1e-10)
      return(-sum(A_vec * log(P_vec) + (1 - A_vec) * log1p(-P_vec)))
    }
    if (loss_name == "auc_as_loss") {
      if (length(unique(A_vec)) < 2L) return(NA_real_)
      ranks <- rank(P_vec, ties.method = "average")
      n_positive <- sum(A_vec == 1)
      n_negative <- sum(A_vec == 0)
      auc_value <- (sum(ranks[A_vec == 1]) -
                      n_positive * (n_positive + 1) / 2) /
        (n_positive * n_negative)
      return(1 - auc_value)
    }
    loss_functions[[loss_name]](A_vec, P_vec)
  }
  evaluate_label_loss <- function(loss_name, g_1, g_2, g_3, g_4) {
    if (loss_name == "pair_hamming_loss") {
      labels_1 <- (rep(g_1, times = length(g_2)) - 1L) * max(g_2) +
        rep(g_2, each = length(g_1))
      labels_2 <- (rep(g_3, times = length(g_4)) - 1L) * max(g_4) +
        rep(g_4, each = length(g_3))
      return(label_match_greedy(
        labels_1, labels_2, K = max(c(labels_1, labels_2)),
        algorithm = "hungarian"
      )$mismatch_rate)
    }
    if (loss_name == "pair_nmi_loss") {
      labels_1 <- (rep(g_1, times = length(g_2)) - 1L) * max(g_2) +
        rep(g_2, each = length(g_1))
      labels_2 <- (rep(g_3, times = length(g_4)) - 1L) * max(g_4) +
        rep(g_4, each = length(g_3))
      if (length(unique(labels_1)) == 1L &&
          length(unique(labels_2)) == 1L) return(-1)
      contingency <- table(labels_1, labels_2)
      probabilities <- contingency / sum(contingency)
      row_probabilities <- rowSums(probabilities)
      col_probabilities <- colSums(probabilities)
      observed <- probabilities > 0
      mutual_information <- sum(probabilities[observed] * log(
        probabilities[observed] /
          (row_probabilities[row(probabilities)][observed] *
             col_probabilities[col(probabilities)][observed])
      ))
      joint_probabilities <- as.numeric(probabilities)
      joint_probabilities <- joint_probabilities[joint_probabilities > 0]
      joint_entropy <- -sum(joint_probabilities * log(joint_probabilities))
      if (joint_entropy == 0) return(NA_real_)
      return(-mutual_information / joint_entropy)
    }
    loss_functions[[loss_name]](g_1, g_2, g_3, g_4)
  }

  if (verbose) {
    message("Stage 3/3: evaluating ", nrow(pair_grid),
            " held-out piece pairs with ", loss_workers, " worker(s).")
  }
  loss_time <- system.time({
    pair_results <- uni_mclapply(
      seq_len(nrow(pair_grid)),
      function(task_idx) {
        set_task_seed(3L, task_idx)
        rep_idx <- pair_grid$repetition[task_idx]
        tau_idx <- pair_grid$tau_idx[task_idx]
        piece_1 <- pair_grid$piece_1[task_idx]
        piece_2 <- pair_grid$piece_2[task_idx]
        estimate_1 <- estimates[[fit_lookup[rep_idx, tau_idx, piece_1]]]
        estimate_2 <- estimates[[fit_lookup[rep_idx, tau_idx, piece_2]]]
        values <- stats::setNames(rep(NA_real_, length(losses)), losses)
        errors <- character()
        if (!estimate_1$success || !estimate_2$success) {
          errors <- c(if (!estimate_1$success) estimate_1$error,
                      if (!estimate_2$success) estimate_2$error)
          return(list(values = values, errors = unique(errors)))
        }
        nodes_1 <- estimate_1$nodes
        nodes_2 <- estimate_2$nodes
        A_test <- A[nodes_1, nodes_2, drop = FALSE]
        P_hat <- NULL
        if (has_probability_loss) {
          B_hat <- average_B_hat[[rep_idx, tau_idx]]
          if (is.null(B_hat)) {
            errors <- c(errors, "No finite averaged block estimate.")
          } else {
            P_hat <- B_hat[estimate_1$labels, estimate_2$labels, drop = FALSE] *
              tcrossprod(estimate_1$psi_hat, estimate_2$psi_hat)
            P_hat <- pmin(pmax(P_hat, 1e-6), 1 - 1e-6)
            if (any(!is.finite(P_hat))) {
              errors <- c(errors, "Predicted probabilities are non-finite.")
              P_hat <- NULL
            }
          }
        }
        reference_labels <- NULL
        if (has_label_loss) {
          reference_fit <- if (use_full_label_reference) {
            complete_reference_fits[[tau_idx]]
          } else {
            A_reference <- A
            A_reference[nodes_1, nodes_2] <- 0
            A_reference[nodes_2, nodes_1] <- 0
            fit_labels(A_reference, tau_candidates[tau_idx])
          }
          if (reference_fit$success) {
            reference_labels <- reference_fit$labels
          } else {
            errors <- c(errors,
                        paste("Label reference failed:",
                              reference_fit$error))
          }
        }
        for (loss_idx in seq_along(losses)) {
          loss_name <- losses[loss_idx]
          value <- tryCatch({
            if (loss_types[loss_name] == "probability") {
              if (is.null(P_hat)) NA_real_ else
                evaluate_probability_loss(loss_name, A_test, P_hat)
            } else {
              if (is.null(reference_labels)) NA_real_ else
                evaluate_label_loss(
                  loss_name, estimate_1$labels, estimate_2$labels,
                  reference_labels[nodes_1], reference_labels[nodes_2]
                )
            }
          }, error = function(error) {
            errors <<- c(errors, paste(loss_name, conditionMessage(error)))
            NA_real_
          })
          if (length(value) != 1L || !is.numeric(value) ||
              !is.finite(value)) {
            errors <- c(errors,
                        paste(loss_name, "returned a non-finite scalar."))
            value <- NA_real_
          }
          values[loss_name] <- value
        }
        list(values = values, errors = unique(errors))
      },
      ncores = loss_workers,
      force_windows = force_windows,
      stop_on_error = TRUE,
      manage_future_plan = manage_future_plan,
      future_packages = worker_future_packages
    )
  })

  cv_records <- vector("list", nrep * length(tau_candidates) * length(losses))
  record_idx <- 1L
  for (rep_idx in seq_len(nrep)) {
    for (tau_idx in seq_along(tau_candidates)) {
      pair_ids <- which(pair_grid$repetition == rep_idx &
                          pair_grid$tau_idx == tau_idx)
      for (loss_name in losses) {
        values <- vapply(
          pair_results[pair_ids],
          function(result) result$values[[loss_name]],
          numeric(1)
        )
        valid <- all(is.finite(values))
        cv_records[[record_idx]] <- data.frame(
          repetition = rep_idx,
          tau = tau_candidates[tau_idx],
          loss = loss_name,
          average_loss = if (valid) mean(values) else NA_real_,
          valid_pairs = sum(is.finite(values)),
          total_pairs = length(values),
          valid = valid
        )
        record_idx <- record_idx + 1L
      }
    }
  }
  cv_loss <- do.call(rbind, cv_records)
  rownames(cv_loss) <- NULL

  select_finite_minimum <- function(tau, values) {
    finite <- is.finite(values)
    if (!any(finite)) return(NA_real_)
    minimum <- min(values[finite])
    min(tau[finite & values == minimum])
  }
  best_records <- vector("list", nrep * length(losses))
  best_idx <- 1L
  for (rep_idx in seq_len(nrep)) {
    for (loss_name in losses) {
      rows <- cv_loss$repetition == rep_idx & cv_loss$loss == loss_name
      best_records[[best_idx]] <- data.frame(
        repetition = rep_idx,
        loss = loss_name,
        tau_hat = select_finite_minimum(
          cv_loss$tau[rows], cv_loss$average_loss[rows]
        )
      )
      best_idx <- best_idx + 1L
    }
  }
  best_tau_each_rep <- do.call(rbind, best_records)
  overall_records <- lapply(losses, function(loss_name) {
    average_by_tau <- vapply(tau_candidates, function(tau) {
      values <- cv_loss$average_loss[
        cv_loss$loss == loss_name & cv_loss$tau == tau
      ]
      if (length(values) != nrep || any(!is.finite(values))) NA_real_ else
        mean(values)
    }, numeric(1))
    data.frame(
      loss = loss_name,
      tau_hat = select_finite_minimum(tau_candidates, average_by_tau),
      average_loss = if (any(is.finite(average_by_tau))) {
        min(average_by_tau, na.rm = TRUE)
      } else {
        NA_real_
      }
    )
  })
  overall_best <- do.call(rbind, overall_records)

  loss_table <- data.frame(candidate_tau = tau_candidates)
  for (loss_name in losses) {
    for (rep_idx in seq_len(nrep)) {
      rows <- cv_loss$loss == loss_name & cv_loss$repetition == rep_idx
      loss_table[[paste0(loss_name, "_rep_", rep_idx)]] <-
        cv_loss$average_loss[rows][match(tau_candidates, cv_loss$tau[rows])]
    }
  }
  loss_table <- tibble::as_tibble(loss_table)
  selection <- list(candidate_tau = tau_candidates)
  for (loss_name in losses) {
    selection[[paste0("tau_hat_each_rep_", loss_name)]] <-
      best_tau_each_rep$tau_hat[best_tau_each_rep$loss == loss_name]
    selection[[paste0("tau_hat_", loss_name)]] <-
      overall_best$tau_hat[overall_best$loss == loss_name]
  }

  diagnostic_rows <- lapply(seq_along(pair_results), function(task_idx) {
    errors <- pair_results[[task_idx]]$errors
    if (length(errors) == 0L) return(NULL)
    data.frame(
      repetition = pair_grid$repetition[task_idx],
      tau = tau_candidates[pair_grid$tau_idx[task_idx]],
      piece_1 = pair_grid$piece_1[task_idx],
      piece_2 = pair_grid$piece_2[task_idx],
      message = errors
    )
  })
  diagnostics <- do.call(rbind, diagnostic_rows)
  matching_diagnostics <- do.call(
    rbind,
    matching_records[!vapply(matching_records, is.null, logical(1))]
  )
  if (is.null(matching_diagnostics)) matching_diagnostics <- data.frame()

  if (verbose) {
    message("Regularizer CV complete; ", sum(!cv_loss$valid),
            " candidate/repetition/loss combination(s) were non-finite.")
  }
  out <- list(
    algorithm = "NETCROP",
    loss = loss_table,
    selection = selection,
    cv_loss = cv_loss,
    best_tau_each_rep = best_tau_each_rep,
    overall_best = overall_best,
    tau_candidates = tau_candidates,
    K = K,
    model = if (use_dcbm) "DCBM" else "SBM",
    losses = data.frame(loss = losses, type = unname(loss_types)),
    nrep = nrep,
    num_subnetworks = num_subnetworks,
    requested_overlap_size = overlap_size,
    effective_overlap_size = effective_overlap_size,
    piece_size = piece_size,
    unordered_test_proportion = pair_count * piece_size^2 / choose(n, 2L),
    matching_diagnostics = if (retain_intermediates == "all") {
      matching_diagnostics
    } else NULL,
    diagnostics = diagnostics,
    partitions = if (retain_intermediates == "all") splitter else NULL,
    raw_fits = if (retain_intermediates == "all") raw_fits else NULL,
    estimates = if (retain_intermediates == "all") estimates else NULL,
    pair_results = if (retain_intermediates == "all") pair_results else NULL,
    retain_intermediates = retain_intermediates,
    options = list(
      use_laplacian = use_laplacian,
      dcbm_est_method = dcbm_est_method,
      spectral_options = spectral_options,
      cluster_options = cluster_options,
      estimator_options = estimator_options,
      matching_method = matching_method,
      label_reference = label_reference,
      parameter_select_options = parameter_select_options
    ),
    ncores = list(requested = ncores, fitting = fit_workers,
                  reference = reference_workers, loss = loss_workers),
    timing = c(
      fitting = unname(fit_time["elapsed"]),
      estimation = unname(estimate_time["elapsed"]),
      reference = unname(reference_time["elapsed"]),
      loss = unname(loss_time["elapsed"]),
      total = unname(fit_time["elapsed"] + estimate_time["elapsed"] +
                       reference_time["elapsed"] + loss_time["elapsed"])
    ),
    ram_preflight = ram_report,
    seed = seed,
    call = call
  )
  class(out) <- "netcrop_regularizer"
  out
}

# Print selected regularization parameters.
print.netcrop_regularizer <- function(x, ...) {
  algorithm <- if (is.null(x$algorithm)) "NETCROP" else toupper(x$algorithm)
  cat(algorithm, "spectral-regularizer selection\n")
  cat(strrep("-", nchar(algorithm) + 31L), "\n", sep = "")
  cat("Model:", x$model, "  K:", x$K, "\n")
  print(x$overall_best, row.names = FALSE)
  invisible(x)
}

# Summarize a NETCROP regularizer fit.
summary.netcrop_regularizer <- function(object, ...) {
  algorithm <- if (is.null(object$algorithm)) {
    "NETCROP"
  } else {
    toupper(object$algorithm)
  }
  if (algorithm == "DKEST") {
    result <- list(
      algorithm = algorithm,
      call = object$call,
      model = object$model,
      K = object$K,
      tau_candidates = object$tau_candidates,
      use_laplacian = object$options$use_laplacian,
      dcbm_est_method = object$options$dcbm_est_method,
      tau_hat = object$tau_hat,
      selected_dk_stat = object$selected_dk_stat,
      overall_best = object$overall_best,
      valid_candidates = sum(object$dk_statistic$valid),
      invalid_candidates = sum(!object$dk_statistic$valid),
      diagnostics = object$diagnostics,
      ncores = object$ncores,
      timing = object$timing
    )
    class(result) <- "summary.netcrop_regularizer"
    return(result)
  }
  result <- list(
    algorithm = algorithm,
    call = object$call,
    model = object$model,
    K = object$K,
    tau_candidates = object$tau_candidates,
    nrep = object$nrep,
    num_subnetworks = object$num_subnetworks,
    requested_overlap_size = object$requested_overlap_size,
    effective_overlap_size = object$effective_overlap_size,
    piece_size = object$piece_size,
    unordered_test_proportion = object$unordered_test_proportion,
    best_tau_each_rep = object$best_tau_each_rep,
    overall_best = object$overall_best,
    invalid_combinations = sum(!object$cv_loss$valid),
    diagnostic_count = if (is.null(object$diagnostics)) 0L else
      nrow(object$diagnostics),
    ncores = object$ncores,
    timing = object$timing
  )
  class(result) <- "summary.netcrop_regularizer"
  result
}

# Print a NETCROP regularizer summary.
print.summary.netcrop_regularizer <- function(x, ...) {
  algorithm <- if (is.null(x$algorithm)) "NETCROP" else toupper(x$algorithm)
  if (algorithm == "DKEST") {
    cat("Summary of DKEST spectral-regularizer selection\n")
    cat("-----------------------------------------------\n")
    cat("Model:", x$model, "  K:", x$K, "\n")
    cat("Candidates:", paste(x$tau_candidates, collapse = ", "), "\n")
    cat("Laplacian scaling:", x$use_laplacian, "\n")
    if (x$model == "DCBM") {
      cat("DCBM estimator:", x$dcbm_est_method, "\n")
    }
    cat("Valid candidates:", x$valid_candidates, "\n")
    cat("Invalid candidates:", x$invalid_candidates, "\n")
    cat("Overall selection:\n")
    print(x$overall_best, row.names = FALSE)
    if (!is.null(x$diagnostics) && nrow(x$diagnostics) > 0L) {
      cat("Diagnostics:\n")
      print(x$diagnostics, row.names = FALSE)
    }
    cat("Timing (seconds):\n")
    print(x$timing)
    return(invisible(x))
  }
  cat("Summary of NETCROP spectral-regularizer selection\n")
  cat("-------------------------------------------------\n")
  cat("Model:", x$model, "  K:", x$K, "\n")
  cat("Candidates:", paste(x$tau_candidates, collapse = ", "), "\n")
  cat("Repetitions:", x$nrep, "\n")
  cat("Subnetworks per repetition:", x$num_subnetworks, "\n")
  cat("Effective overlap size:", x$effective_overlap_size, "\n")
  cat("Non-overlap piece size:", x$piece_size, "\n")
  cat(sprintf("Held-out unordered-pair proportion: %.2f%%\n",
              100 * x$unordered_test_proportion))
  cat("Invalid CV combinations:", x$invalid_combinations, "\n")
  cat("Diagnostics:", x$diagnostic_count, "\n")
  cat("Overall selections:\n")
  print(x$overall_best, row.names = FALSE)
  cat("Timing (seconds):\n")
  print(x$timing)
  invisible(x)
}

# Plot regularization-parameter CV loss curves.
plot.netcrop_regularizer <- function(x, aggregate = TRUE, ...) {
  algorithm <- if (is.null(x$algorithm)) "NETCROP" else toupper(x$algorithm)
  if (algorithm == "DKEST") {
    stop(
      "DKEST directly estimates tau_hat and has no CV loss to plot.",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting.", call. = FALSE)
  }
  if (length(aggregate) != 1L || !is.logical(aggregate) || is.na(aggregate)) {
    stop("aggregate must be TRUE or FALSE.", call. = FALSE)
  }
  plot_data <- x$cv_loss[
    is.finite(x$cv_loss$average_loss), , drop = FALSE
  ]
  if (nrow(plot_data) == 0L) {
    stop("No finite regularizer criteria are available for plotting.",
         call. = FALSE)
  }
  if (aggregate) {
    plot_data <- stats::aggregate(
      average_loss ~ tau + loss,
      data = plot_data,
      FUN = function(values) {
        if (any(!is.finite(values))) NA_real_ else mean(values)
      },
      na.action = stats::na.pass
    )
    mapping <- ggplot2::aes(x = tau, y = average_loss)
  } else {
    plot_data$repetition <- factor(plot_data$repetition)
    mapping <- ggplot2::aes(
      x = tau, y = average_loss, group = repetition, color = repetition
    )
  }
  ggplot2::ggplot(plot_data, mapping) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~loss, scales = "free_y") +
    ggplot2::labs(
      title = paste(algorithm, "criterion by spectral regularization"),
      x = "Regularization parameter (tau)",
      y = "Average CV loss",
      color = "Repetition"
    ) +
    ggplot2::theme_minimal()
}




# Benchmark code (intentionally not run when this file is sourced).
# library(Matrix)
# system.time(
#   net1 <- generate_dcbm(
#     n = 1000,
#     K = 5,
#     P_block = 0.3 * diag(2 / 3, 5, 5) + 0.3 * matrix(1 / 3, 5, 5),
#     psi = sample(rbeta(300, 1, 4), 1000, replace = TRUE),
#     ncores = 8L
#   )
# )
# system.time(
#   net2 <- generate_dcbm(
#     n = 1000,
#     K = 5,
#     P_block = 0.3 * diag(2 / 3, 5, 5) + 0.3 * matrix(1 / 3, 5, 5),
#     psi = sample(rbeta(300, 1, 4), 1000, replace = TRUE),
#     ncores = 8L
#   )
# )
# 
# 
# system.time(
#   nc_out1 <- netcrop_tune_regularizer(
#     A = net1,
#     K = 5,
#     tau_candidates = seq(0, 2, by = 0.05),
#     use_dcbm = TRUE,
#     nrep = 5L,
#     use_laplacian = TRUE,
#     losses = "sse",
#     # label_reference = "leave_pair_out",
#     ncores = 8L,
#     force_windows = FALSE
#   )
# )
# 
# system.time(
#   nc_out2 <- netcrop_tune_regularizer(
#     A = net2,
#     K = 5,
#     tau_candidates = seq(0, 2, by = 0.05),
#     use_dcbm = TRUE,
#     nrep = 5L,
#     use_laplacian = TRUE,
#     losses = "sse",
#     # label_reference = "leave_pair_out",
#     ncores = 8L,
#     force_windows = FALSE
#   )
# )
# 
# # nc_out
# # summary(nc_out)
# # plot(nc_out)
# 
# 
# system.time(
#   dk.out1 <- dkest_tune_regularizer(
#     A = net1,
#     K = 5,
#     tau_candidates = seq(0, 2, by = 0.05),
#     use_laplacian = TRUE,
#     use_dcbm = TRUE,
#     ncores = 8L
#   )
# )
# 
# system.time(
#   dk.out2 <- dkest_tune_regularizer(
#     A = net2,
#     K = 5,
#     tau_candidates = seq(0, 2, by = 0.05),
#     use_laplacian = TRUE,
#     use_dcbm = TRUE,
#     ncores = 8L
#   )
# )
# 
# # dk.out
# # summary(dk.out)
# # plot(dk.out)
# 
# system.time(
#   or.out <- oracle_plotter(
#     A = list(net1, net2),
#     g_true = list(attributes(net1)$generator_parameters$g_true,
#                   attributes(net2)$generator_parameters$g_true),
#     tau_candidates = seq(0, 2, by = 0.05),
#     netcrop_outcomes = list(nc_out1, nc_out2),
#     dkest_outcomes = list(dk.out1, dk.out2),
#     include_netcrop_mean = TRUE,
#     include_netcrop_mode = TRUE,
#     losses = NULL,
#     engines = c("sonnet", "spectral_cluster"),
#     matching_method = "greedy",
#     ncores = 8L,
#     force_windows = FALSE
#   )
# )
# 
# or.out
# attributes(or.out)
