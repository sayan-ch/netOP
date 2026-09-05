# NETCROP block-model selection
#
# This file preserves the original three-stage algorithm: obtain subnetwork
# spectral labels, estimate block-model parameters, and evaluate every pair of
# non-overlap pieces. Package loading resolves all internal helpers.

# Select an SBM or DCBM and number of communities by overlapping-subnetwork CV.
#' Select a block model with NETCROP
#'
#' Selects the number of communities and either the stochastic block model
#' (SBM), degree-corrected block model (DCBM), or both by overlapping-subnetwork
#' cross-validation. The procedure obtains spectral labels on each subnetwork,
#' estimates block-model parameters, aligns labels through the common overlap,
#' and evaluates predictions between every pair of non-overlap pieces.
#'
#' @param A Finite square symmetric binary adjacency matrix with zero diagonal,
#'   including a supported sparse `Matrix` object.
#' @param K_candidates Non-empty vector of positive integer community counts;
#'   duplicates are discarded. Conventionally this is `1:5`.
#' @param num_subnetworks Optional integer number of subnetworks, at least 2. If
#'   it or `overlap_size` is `NULL`, missing values are selected with
#'   [netcrop_param_select()].
#' @param overlap_size Optional positive integer overlap size. It must leave
#'   enough non-overlap nodes for the requested subnetworks and candidates.
#' @param nrep Positive integer number of independent NETCROP repetitions.
#' @param losses Character vector naming prediction-loss functions available in
#'   the calling environment, such as `"sse"`.
#' @param model_candidates Character vector selecting `"SBM"`, `"DCBM"`, or
#'   both candidate model families.
#' @param sbm_est_options,dcbm_est_options Named lists of additional arguments
#'   for the corresponding block-model estimators.
#' @param matching_method Label-alignment method: `"greedy"`, `"hungarian"`, or
#'   `"brute_force"`.
#' @param confirm_large Optional logical passed to label matching to confirm an
#'   unusually large brute-force matching problem.
#' @param ncores Positive integer worker count used across decomposition,
#'   estimation, alignment, and loss tasks.
#' @param seed Optional nonnegative integer-like reproducibility seed.
#' @param verbose Whether to print progress messages.
#' @param force_windows Whether to force the Windows-compatible parallel
#'   backend, including for backend testing on other platforms.
#' @param ram_check Whether to run RAM preflight checks before major operations.
#' @param parameter_select_options Named list of additional options passed to
#'   [netcrop_param_select()] when partition parameters are selected
#'   automatically.
#' @param retain_intermediates Either `"all"` or `"minimal"`, controlling how
#'   many fitted models and intermediate arrays are retained.
#' @param laplacian Whether to convert subnetwork adjacency matrices to graph
#'   Laplacians before spectral clustering.
#' @param regularize_tau Nonnegative Laplacian regularization strength shared by
#'   SBM and DCBM candidates.
#' @param x A `netcrop_blockmodel` or `summary.netcrop_blockmodel` object. For
#'   plotting, it is the first result to display or compare.
#' @param object A fitted `netcrop_blockmodel` object to summarize.
#' @param aggregate Whether a plot should aggregate loss curves across
#'   repetitions. A `netcrop_blockmodel` object in this position instead starts
#'   a comparison plot.
#' @param ... Additional plotting arguments or further `netcrop_blockmodel`
#'   objects for comparison; ignored by print and summary methods.
#'
#' @return `netcrop_blockmodel()` returns an object of class
#'   `netcrop_blockmodel` containing candidate-wise cross-validation losses,
#'   repetition and overall selections, split and alignment diagnostics,
#'   retained intermediates, resource information, and timing.
#'   `summary.netcrop_blockmodel()` returns a compact summary object; print
#'   methods return their input invisibly, and the plot method returns a
#'   `ggplot` object.
#' @seealso [ecv_stability_blockmodel()], [ncv_stability_blockmodel()]
#' @examples
#' \donttest{
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.5, beta = 0.1,
#'                   seed = 5, ncores = 1)
#' fit <- netcrop_blockmodel(
#'   A, K_candidates = 1:5, num_subnetworks = 2, overlap_size = 30,
#'   model_candidates = "SBM", ncores = 1, seed = 6, verbose = FALSE
#' )
#' fit
#' summary(fit)
#' }
#' @rdname netcrop_blockmodel
#' @export
netcrop_blockmodel <- function(
    A,
    K_candidates,
    num_subnetworks = NULL,
    overlap_size = NULL,
    nrep = 1L,
    losses = "sse",
    model_candidates = c("SBM", "DCBM"),
    sbm_est_options = list(),
    dcbm_est_options = list(),
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
    retain_intermediates = c("all", "minimal"),
    laplacian = FALSE,
    regularize_tau = 0) {
  call <- match.call()
  laplacian_missing <- missing(laplacian)
  regularize_tau_missing <- missing(regularize_tau)
  matching_method <- match.arg(matching_method)
  retain_intermediates <- match.arg(retain_intermediates)
  required_helpers <- c(
    "uni_mclapply", "netcrop_param_select", "netcrop_splitter",
    "is_symmetric_matrix", "eig_decomp",
    "singular_decomp", "graph_laplacian", "spectral_cluster",
    "estimate_sbm", "estimate_dcbm", "label_match_greedy",
    "label_match_brute_force", "clip_probabilities", "modal"
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
      paste(missing_helpers, collapse = ", "),
      ".",
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
  merge_options <- function(defaults, supplied, name) {
    validate_named_list(supplied, name)
    unknown <- setdiff(names(supplied), names(defaults))
    if (length(unknown) > 0L) {
      stop(
        name, " contains unsupported component(s): ",
        paste(unknown, collapse = ", "), ".",
        call. = FALSE
      )
    }
    utils::modifyList(defaults, supplied, keep.null = TRUE)
  }
  compatible_cluster_options <- function(options, engine) {
    cluster_function <- switch(
      engine,
      kmeans = stats::kmeans,
      clara = cluster::clara,
      pam = cluster::pam
    )
    accepted <- setdiff(
      names(formals(cluster_function)),
      c("x", "centers", "k")
    )
    options[intersect(names(options), accepted)]
  }

  if (is.null(dim(A)) || length(dim(A)) != 2L ||
      nrow(A) != ncol(A) || nrow(A) < 1L) {
    stop("A must be a non-empty square matrix-like object.", call. = FALSE)
  }
  is_numeric_matrix <- is.numeric(A) || inherits(A, "Matrix")
  stored_values_are_finite <- if (inherits(A, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(A)) {
      all(is.finite(methods::slot(A, "x")))
    } else {
      TRUE
    }
  } else {
    all(is.finite(A))
  }
  if (!is_numeric_matrix || !stored_values_are_finite) {
    stop("A must contain only finite numeric values.", call. = FALSE)
  }
  input_was_sparse <- inherits(A, "sparseMatrix")
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
    stop(
      "Every non-overlap piece must contain at least one node.",
      call. = FALSE
    )
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
  if (!is.null(confirm_large) &&
      (length(confirm_large) != 1L || !is.logical(confirm_large) ||
       is.na(confirm_large))) {
    stop("confirm_large must be NULL, TRUE, or FALSE.", call. = FALSE)
  }

  if (!is.numeric(K_candidates) || length(K_candidates) < 1L ||
      anyNA(K_candidates) || any(!is.finite(K_candidates)) ||
      any(K_candidates < 1) || any(K_candidates != floor(K_candidates))) {
    stop("K_candidates must contain positive integers.", call. = FALSE)
  }
  K_candidates <- unique(as.integer(K_candidates))
  K_max <- max(K_candidates)
  model_candidates <- unique(toupper(as.character(model_candidates)))
  if (length(model_candidates) < 1L ||
      any(!model_candidates %in% c("SBM", "DCBM"))) {
    stop(
      "model_candidates must contain one or both of 'SBM' and 'DCBM'.",
      call. = FALSE
    )
  }
  if (!is.character(losses) || length(losses) < 1L ||
      anyNA(losses) || any(!nzchar(losses))) {
    stop("losses must contain non-empty function names.", call. = FALSE)
  }
  losses <- unique(losses)
  selected_loss_functions <- lapply(losses, function(loss_name) {
    if (!exists(loss_name, mode = "function", inherits = TRUE, envir = environment())) {
      stop("Loss function not found: ", loss_name, ".", call. = FALSE)
    }
    get(loss_name, mode = "function", inherits = TRUE, envir = environment())
  })
  names(selected_loss_functions) <- losses
  # future/multisession cannot discover helpers reached dynamically through a
  # list of function objects. Bundle the standard loss functions and their
  # dependencies in one serializable lexical environment so Windows workers
  # receive exactly the same implementations as the main R process.
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

  A_was_symmetric <- is_symmetric_matrix(A)
  A_values <- if (inherits(A, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(A)) {
      methods::slot(A, "x")
    } else {
      rep.int(1, length(methods::slot(A, "i")))
    }
  } else {
    as.numeric(A)
  }
  A_was_binary <- all(A_values %in% c(0, 1))
  A_binary <- NULL
  if (!A_was_symmetric) {
    warning(
      "A is asymmetric; using the undirected binarized union A + t(A).",
      call. = FALSE
    )
    A_binary <- A != 0
    A <- (A_binary + if (inherits(A_binary, "Matrix")) {
      Matrix::t(A_binary)
    } else {
      t(A_binary)
    }) > 0
  } else {
    if (!A_was_binary) {
      warning("A is weighted; binarizing its nonzero entries.", call. = FALSE)
      A <- (A != 0) * 1
    } else if (is.logical(A)) {
      A <- A * 1
    }
  }
  if (!is.numeric(A)) A <- A * 1
  if (inherits(A, "symmetricMatrix")) {
    A <- methods::as(A, "generalMatrix")
  }

  spectral_defaults <- list(
    laplacian = FALSE,
    normalize_laplacian = TRUE,
    regularize_tau = 0,
    handle_zero_degree_nodes = "none",
    row_normalize = FALSE,
    spectral_method = "eigen",
    spectral_engine = "RSpectra",
    spectral_options = list(),
    cluster_engine = "clara",
    cluster_options = list(samples = 5L)
  )
  sbm_defaults <- list(
    spectral_cluster = utils::modifyList(
      spectral_defaults,
      list(
        row_normalize = FALSE,
        cluster_options = list(
          metric = "euclidean",
          cluster.only = TRUE,
          samples = 5L
        )
      )
    ),
    estimate = list()
  )
  dcbm_defaults <- list(
    spectral_cluster = utils::modifyList(
      spectral_defaults,
      list(
        row_normalize = TRUE,
        cluster_options = list(
          metric = "manhattan",
          cluster.only = TRUE,
          samples = 5L
        )
      )
    ),
    estimate = list(method = "plugin")
  )
  legacy_spectral_values <- function(option_name) {
    values <- list()
    if ("SBM" %in% model_candidates &&
        !is.null(sbm_est_options$spectral_cluster[[option_name]])) {
      values$SBM <- sbm_est_options$spectral_cluster[[option_name]]
    }
    if ("DCBM" %in% model_candidates &&
        !is.null(dcbm_est_options$spectral_cluster[[option_name]])) {
      values$DCBM <- dcbm_est_options$spectral_cluster[[option_name]]
    }
    values
  }
  resolve_shared_direct_option <- function(value, was_missing, option_name) {
    legacy <- legacy_spectral_values(option_name)
    if (length(legacy) == 0L) {
      return(value)
    }
    if (length(legacy) > 1L &&
        !all(vapply(legacy[-1L], identical, logical(1), legacy[[1L]]))) {
      stop(
        "SBM and DCBM ", option_name, " options must agree.",
        call. = FALSE
      )
    }
    legacy_value <- legacy[[1L]]
    if (!was_missing && !identical(value, legacy_value)) {
      stop(
        "Specify ", option_name, " directly or in estimator options, not both ",
        "with different values.",
        call. = FALSE
      )
    }
    if (was_missing) legacy_value else value
  }
  laplacian <- resolve_shared_direct_option(
    laplacian, laplacian_missing, "laplacian"
  )
  regularize_tau <- resolve_shared_direct_option(
    regularize_tau, regularize_tau_missing, "regularize_tau"
  )
  laplacian <- validate_flag(laplacian, "laplacian")
  if (length(regularize_tau) != 1L || !is.numeric(regularize_tau) ||
      is.na(regularize_tau) || !is.finite(regularize_tau) ||
      regularize_tau < 0) {
    stop("regularize_tau must be one finite non-negative number.",
         call. = FALSE)
  }
  sbm_options <- merge_options(
    sbm_defaults,
    sbm_est_options,
    "sbm_est_options"
  )
  dcbm_options <- merge_options(
    dcbm_defaults,
    dcbm_est_options,
    "dcbm_est_options"
  )
  sbm_options$spectral_cluster <- merge_options(
    sbm_defaults$spectral_cluster,
    sbm_options$spectral_cluster,
    "sbm_est_options$spectral_cluster"
  )
  dcbm_options$spectral_cluster <- merge_options(
    dcbm_defaults$spectral_cluster,
    dcbm_options$spectral_cluster,
    "dcbm_est_options$spectral_cluster"
  )
  sbm_options$spectral_cluster$laplacian <- laplacian
  dcbm_options$spectral_cluster$laplacian <- laplacian
  sbm_options$spectral_cluster$regularize_tau <- regularize_tau
  dcbm_options$spectral_cluster$regularize_tau <- regularize_tau
  sbm_options$spectral_cluster$cluster_options <- compatible_cluster_options(
    sbm_options$spectral_cluster$cluster_options,
    sbm_options$spectral_cluster$cluster_engine
  )
  dcbm_options$spectral_cluster$cluster_options <- compatible_cluster_options(
    dcbm_options$spectral_cluster$cluster_options,
    dcbm_options$spectral_cluster$cluster_engine
  )
  validate_named_list(sbm_options$estimate, "sbm_est_options$estimate")
  validate_named_list(dcbm_options$estimate, "dcbm_est_options$estimate")
  dcbm_estimation_method <- if (is.null(dcbm_options$estimate$method)) {
    "plugin"
  } else {
    match.arg(dcbm_options$estimate$method, c("plugin", "spectral"))
  }
  sbm_estimate_conflicts <- intersect(
    names(sbm_options$estimate),
    c("A", "g", "K", "directed", "self_loops", "validate_inputs")
  )
  dcbm_estimate_conflicts <- intersect(
    names(dcbm_options$estimate),
    c("A", "g", "K", "row_norm", "psi_omit", "validate_inputs")
  )
  if (length(sbm_estimate_conflicts) > 0L ||
      length(dcbm_estimate_conflicts) > 0L) {
    stop(
      "Estimator options cannot override protected arguments: ",
      paste(c(sbm_estimate_conflicts, dcbm_estimate_conflicts),
            collapse = ", "), ".",
      call. = FALSE
    )
  }
  unsupported_sbm_estimate <- setdiff(
    names(sbm_options$estimate), names(formals(estimate_sbm))
  )
  unsupported_dcbm_estimate <- setdiff(
    names(dcbm_options$estimate), names(formals(estimate_dcbm))
  )
  if (length(unsupported_sbm_estimate) > 0L ||
      length(unsupported_dcbm_estimate) > 0L) {
    stop(
      "Unsupported estimator option(s): ",
      paste(
        c(unsupported_sbm_estimate, unsupported_dcbm_estimate),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  if (!is.null(dcbm_options$estimate$stabilizer)) {
    stabilizer <- dcbm_options$estimate$stabilizer
    if (length(stabilizer) != 1L || !is.numeric(stabilizer) ||
        is.na(stabilizer) || !is.finite(stabilizer) || stabilizer < 0) {
      stop("stabilizer must be one finite non-negative number.",
           call. = FALSE)
    }
  }
  if (!is.null(dcbm_options$estimate$spectral_options)) {
    validate_named_list(
      dcbm_options$estimate$spectral_options,
      "dcbm_est_options$estimate$spectral_options"
    )
  }
  for (option_name in c("spectral_options", "cluster_options")) {
    validate_named_list(
      sbm_options$spectral_cluster[[option_name]],
      paste0("sbm_est_options$spectral_cluster$", option_name)
    )
    validate_named_list(
      dcbm_options$spectral_cluster[[option_name]],
      paste0("dcbm_est_options$spectral_cluster$", option_name)
    )
  }

  shared_names <- c(
    "laplacian", "normalize_laplacian", "regularize_tau",
    "spectral_method", "spectral_engine", "spectral_options"
  )
  active_spectral_options <- if (identical(model_candidates, "DCBM")) {
    dcbm_options$spectral_cluster
  } else {
    sbm_options$spectral_cluster
  }
  if (all(c("SBM", "DCBM") %in% model_candidates)) {
    conflicts <- shared_names[!vapply(shared_names, function(name) {
      identical(
        sbm_options$spectral_cluster[[name]],
        dcbm_options$spectral_cluster[[name]]
      )
    }, logical(1))]
    if (length(conflicts) > 0L) {
      stop(
        paste0(
          "SBM and DCBM decomposition options must agree to reuse one ",
          "decomposition. Conflicts: ", paste(conflicts, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
  }
  active_spectral_options$laplacian <- match.arg(
    as.character(active_spectral_options$laplacian),
    c("TRUE", "FALSE")
  ) == "TRUE"
  active_spectral_options$normalize_laplacian <- match.arg(
    as.character(active_spectral_options$normalize_laplacian),
    c("TRUE", "FALSE")
  ) == "TRUE"
  active_spectral_options$spectral_method <- match.arg(
    tolower(active_spectral_options$spectral_method),
    c("eigen", "svd")
  )
  active_spectral_options$spectral_engine <- match.arg(
    active_spectral_options$spectral_engine,
    c("RSpectra", "irlba", "base")
  )
  if (length(active_spectral_options$regularize_tau) != 1L ||
      !is.numeric(active_spectral_options$regularize_tau) ||
      is.na(active_spectral_options$regularize_tau) ||
      !is.finite(active_spectral_options$regularize_tau) ||
      active_spectral_options$regularize_tau < 0) {
    stop("regularize_tau must be one finite non-negative number.",
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
  if (K_max >= subgraph_size) {
    stop(
      "max(K_candidates) must be smaller than the effective subgraph size (",
      subgraph_size, ").",
      call. = FALSE
    )
  }

  if (matching_method == "brute_force" && K_max > 8L) {
    warning_text <- paste0(
      "Brute-force matching will check as many as ",
      format(factorial(K_max), scientific = FALSE, trim = TRUE),
      " permutations per match and can be very slow. ",
      "Greedy matching is preferred."
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
        c("No, cancel", "Yes, continue"),
        title = paste0(warning_text, " Continue?")
      )
      if (choice != 2L) {
        stop("Brute-force matching was cancelled.", call. = FALSE)
      }
    } else if (!confirm_large) {
      stop("Brute-force matching was not confirmed.", call. = FALSE)
    }
  }

  raw_task_count <- nrep * num_subnetworks
  estimation_combination_count <- raw_task_count * length(K_candidates)
  estimation_task_count <- raw_task_count
  pair_count <- choose(num_subnetworks, 2L)
  loss_task_count <- nrep * length(K_candidates) * pair_count
  raw_ncores <- min(raw_task_count, ncores)
  estimation_ncores <- min(estimation_task_count, ncores)
  loss_ncores <- min(loss_task_count, ncores)
  manage_future_plan <- TRUE
  os_type <- if (force_windows) "windows" else .Platform$OS.type
  worker_future_packages <- if (inherits(A, "Matrix")) {
    "Matrix"
  } else {
    NULL
  }
  maximum_stage_ncores <- max(raw_ncores, estimation_ncores, loss_ncores)
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
    future::plan(
      future::multisession,
      workers = maximum_stage_ncores
    )
    manage_future_plan <- FALSE
  }

  ram_report <- NULL
  if (ram_check) {
    ram_helpers <- c(
      "estimate_spectral_decomp_ram", "estimate_matrix_product_ram",
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
        paste(missing_ram_helpers, collapse = ", "), ".",
        call. = FALSE
      )
    }
    dense_spectral_input <- isTRUE(active_spectral_options$laplacian) ||
      active_spectral_options$regularize_tau > 0
    spectral_ram <- estimate_spectral_decomp_ram(
      n = subgraph_size,
      p = subgraph_size,
      K = K_max,
      method = active_spectral_options$spectral_method,
      engine = tolower(active_spectral_options$spectral_engine),
      force_engine = isTRUE(
        active_spectral_options$spectral_options$force_engine
      ),
      safe_d_multiplier = if (is.null(
        active_spectral_options$spectral_options$safe_d_multiplier
      )) 1 else active_spectral_options$spectral_options$safe_d_multiplier,
      nu = K_max,
      nv = K_max,
      dense_input = dense_spectral_input
    )$estimated_bytes
    ram_terms <- list(list(
      estimated_bytes = spectral_ram,
      sequential_count = 1L,
      parallel_count = raw_ncores,
      label = "shared subgraph spectral decomposition"
    ))
    if (isTRUE(active_spectral_options$laplacian)) {
      ram_terms <- c(ram_terms, list(list(
        estimated_bytes = estimate_matrix_product_ram(
          subgraph_size,
          1L,
          subgraph_size
        ),
        sequential_count = 1L,
        parallel_count = raw_ncores,
        label = "dense Laplacian or normalized-adjacency construction"
      )))
    }
    loss_matrix_count <- if (all(c("SBM", "DCBM") %in%
                                 model_candidates)) 5L else if (
                                   "DCBM" %in% model_candidates
                                 ) 4L else 2L
    retained_storage <-
      raw_task_count * (
        8 * as.double(subgraph_size) * K_max +
          4 * as.double(subgraph_size) * length(K_candidates) *
            length(model_candidates)
      ) +
      estimation_combination_count * length(model_candidates) * (
        4 * as.double(piece_size) +
          8 * as.double(K_max)^2 +
          if ("DCBM" %in% model_candidates) {
            8 * as.double(piece_size)
          } else {
            0
          }
      )
    ram_terms <- c(ram_terms, list(
      list(
        estimated_bytes = 8 * as.double(subgraph_size)^2,
        sequential_count = 1L,
        parallel_count = estimation_ncores,
        label = "estimation-stage subgraph copy"
      ),
      list(
        estimated_bytes = estimate_matrix_product_ram(
          piece_size,
          1L,
          piece_size
        ),
        sequential_count = loss_matrix_count,
        parallel_count = loss_ncores,
        label = "largest held-out matrix multiplication or allocation"
      ),
      list(
        estimated_bytes = as.numeric(utils::object.size(A)),
        sequential_count = 1L,
        parallel_count = max(raw_ncores, estimation_ncores, loss_ncores),
        label = "full adjacency copy per worker"
      ),
      list(
        estimated_bytes = retained_storage,
        sequential_count = 1L,
        parallel_count = 1L,
        label = "retained labels, estimates, and diagnostics"
      )
    ))
    ram_report <- report_ram_formula(
      terms = ram_terms,
      operation = "NETCROP block-model conservative preflight",
      detail = "explicit stage-specific sequential and parallel factors"
    )
  }

  set_task_seed <- function(stage, task_id) {
    if (is.null(seed)) {
      return(invisible(NULL))
    }
    task_seed <- as.integer(
      (seed + as.double(stage) * 1000003 + task_id) %%
        .Machine$integer.max
    )
    set.seed(task_seed)
    invisible(task_seed)
  }
  task_grid <- expand.grid(
    repetition = seq_len(nrep),
    split = seq_len(num_subnetworks),
    KEEP.OUT.ATTRS = FALSE
  )
  task_grid <- task_grid[order(task_grid$repetition, task_grid$split), ]

  build_spectral_representation <- function(A_subgraph) {
    tau <- active_spectral_options$regularize_tau
    laplacian <- active_spectral_options$laplacian
    normalized <- active_spectral_options$normalize_laplacian
    method <- active_spectral_options$spectral_method
    engine <- tolower(active_spectral_options$spectral_engine)
    spectral_options <- active_spectral_options$spectral_options
    strict_sparse_path <- input_was_sparse && tau == 0
    if (strict_sparse_path && engine == "base") {
      stop(
        paste0(
          "engine = 'base' would densify a sparse adjacency matrix. ",
          "Use RSpectra for eigen decompositions or RSpectra/irlba for SVD, ",
          "or provide a dense A explicitly."
        ),
        call. = FALSE
      )
    }
    if (strict_sparse_path) {
      # Low-level decomposition helpers normally fall back to base R for small
      # matrices or unavailable partial engines. That fallback is forbidden
      # here because it would silently densify a sparse NETCROP subgraph.
      spectral_options$force_engine <- TRUE
    }
    protected <- c(
      "A", "d", "only_values", "nu", "nv", "engine",
      "use_laplacian", "order_by", "validate_inputs"
    )
    conflicts <- intersect(names(spectral_options), protected)
    if (length(conflicts) > 0L) {
      stop(
        "spectral_options cannot override ",
        paste(conflicts, collapse = ", "), ".",
        call. = FALSE
      )
    }
    deg <- if (inherits(A_subgraph, "Matrix")) {
      Matrix::rowSums(A_subgraph)
    } else {
      rowSums(A_subgraph)
    }
    average_degree <- mean(deg)
    if (laplacian && normalized) {
      spectral_matrix <- A_subgraph
      if (tau > 0) {
        spectral_matrix <- spectral_matrix +
          tau * average_degree / nrow(A_subgraph)
      }
      inverse_sqrt_deg <- numeric(length(deg))
      regularized_deg <- deg + tau * average_degree
      positive <- regularized_deg > 0
      inverse_sqrt_deg[positive] <- 1 / sqrt(regularized_deg[positive])
      if (inherits(spectral_matrix, "sparseMatrix")) {
        scaling_matrix <- Matrix::Diagonal(x = inverse_sqrt_deg)
        spectral_matrix <- crossprod(
          scaling_matrix,
          tcrossprod(spectral_matrix, scaling_matrix)
        )
      } else {
        spectral_matrix <- spectral_matrix * tcrossprod(inverse_sqrt_deg)
      }
      spectral_matrix[!is.finite(spectral_matrix)] <- 0
      order_by <- "magnitude"
    } else if (laplacian) {
      spectral_matrix <- -graph_laplacian(
        A_subgraph,
        normalized = FALSE,
        tau = tau
      )
      order_by <- "value"
    } else {
      spectral_matrix <- A_subgraph
      if (tau > 0) {
        spectral_matrix <- spectral_matrix +
          tau * average_degree / nrow(A_subgraph)
      }
      order_by <- "magnitude"
    }
    if (method == "eigen") {
      arguments <- utils::modifyList(
        list(
          A = spectral_matrix,
          d = K_max,
          only_values = FALSE,
          use_laplacian = FALSE,
          engine = engine,
          order_by = order_by,
          validate_inputs = FALSE
        ),
        spectral_options,
        keep.null = TRUE
      )
      decomposition <- do.call(eig_decomp, arguments)
      return(list(U = decomposition$vectors, values = decomposition$values))
    }
    arguments <- utils::modifyList(
      list(
        A = spectral_matrix,
        d = K_max,
        only_values = FALSE,
        nu = K_max,
        nv = 0L,
        use_laplacian = FALSE,
        engine = engine,
        order_by = "value",
        validate_inputs = FALSE
      ),
      spectral_options,
      keep.null = TRUE
    )
    decomposition <- do.call(singular_decomp, arguments)
    list(U = decomposition$u, values = decomposition$values)
  }

  cluster_from_U <- function(U, K, options) {
    do.call(
      spectral_cluster,
      list(
        U = U[, seq_len(K), drop = FALSE],
        K = K,
        row_normalize = options$row_normalize,
        cluster_engine = options$cluster_engine,
        cluster_options = options$cluster_options,
        ram_check = FALSE,
        validate_inputs = FALSE
      )
    )
  }

  if (verbose) {
    message(
      "Stage 1/3: clustering ", raw_task_count,
      " subnetworks with ", raw_ncores, " worker(s)."
    )
  }
  raw_time <- system.time({
    raw_output <- uni_mclapply(
      seq_len(raw_task_count),
      function(task_id) {
        if (input_was_sparse &&
            !requireNamespace("Matrix", quietly = TRUE)) {
          stop("The Matrix package is required for sparse A.", call. = FALSE)
        }
        set_task_seed(1L, task_id)
        repetition <- task_grid$repetition[task_id]
        split_id <- task_grid$split[task_id]
        nodes <- splitter$subnetworks[[repetition]][[split_id]]
        A_subgraph <- A[nodes, nodes, drop = FALSE]
        representation <- build_spectral_representation(A_subgraph)
        labels_by_model <- stats::setNames(vector("list", length(model_candidates)),
                                    model_candidates)
        retain_row_norms <- "DCBM" %in% model_candidates &&
          dcbm_estimation_method == "spectral"
        row_norms <- if (retain_row_norms) {
          vector("list", length(K_candidates))
        } else {
          NULL
        }
        for (candidate_id in seq_along(K_candidates)) {
          K <- K_candidates[candidate_id]
          U_K <- representation$U[, seq_len(K), drop = FALSE]
          if (K == 1L) {
            constant_labels <- rep.int(1L, nrow(U_K))
            if ("SBM" %in% model_candidates) {
              labels_by_model$SBM[[candidate_id]] <- constant_labels
            }
            if ("DCBM" %in% model_candidates) {
              labels_by_model$DCBM[[candidate_id]] <- constant_labels
              if (retain_row_norms) {
                row_norms[[candidate_id]] <- pmax(abs(U_K[, 1L]), 1e-6)
              }
            }
            next
          }
          if ("SBM" %in% model_candidates) {
            labels_by_model$SBM[[candidate_id]] <-
              cluster_from_U(U_K, K, sbm_options$spectral_cluster)$g_hat
          }
          if ("DCBM" %in% model_candidates) {
            row_norm <- pmax(sqrt(rowSums(U_K^2)), 1e-6)
            U_normalized <- U_K / row_norm
            dcbm_cluster_options <- dcbm_options$spectral_cluster
            dcbm_cluster_options$row_normalize <- FALSE
            labels_by_model$DCBM[[candidate_id]] <-
              cluster_from_U(
                U_normalized, K, dcbm_cluster_options
              )$g_hat
            if (retain_row_norms) row_norms[[candidate_id]] <- row_norm
          }
        }
        list(
          labels = labels_by_model,
          row_norms = row_norms,
          spectral_values = representation$values,
          nodes = nodes
        )
      },
      ncores = raw_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      manage_future_plan = manage_future_plan,
      future_packages = worker_future_packages
    )
  })

  raw_task_lookup <- matrix(
    seq_len(raw_task_count),
    nrow = nrep,
    ncol = num_subnetworks,
    byrow = TRUE
  )
  match_labels <- function(match_this, standard, K) {
    if (matching_method == "brute_force") {
      label_match_brute_force(
        match_this = match_this,
        standard = standard,
        K = K,
        return_mapping = TRUE,
        confirm_large = TRUE
      )
    } else {
      label_match_greedy(
        match_this = match_this,
        standard = standard,
        K = K,
        algorithm = matching_method,
        return_mapping = TRUE
      )
    }
  }
  if (verbose) {
    message(
      "Stage 2/3: estimating ", estimation_combination_count,
      " subnetwork/candidate combinations in ", estimation_task_count,
      " grouped task(s) with ", estimation_ncores, " worker(s)."
    )
  }
  estimation_time <- system.time({
    estimation_output <- uni_mclapply(
      seq_len(estimation_task_count),
      function(raw_task_id) {
        if (input_was_sparse &&
            !requireNamespace("Matrix", quietly = TRUE)) {
          stop("The Matrix package is required for sparse A.", call. = FALSE)
        }
        repetition <- task_grid$repetition[raw_task_id]
        split_id <- task_grid$split[raw_task_id]
        nodes <- splitter$subnetworks[[repetition]][[split_id]]
        A_subgraph <- A[nodes, nodes, drop = FALSE]
        overlap_count <- length(splitter$overlap_nodes[[repetition]])
        standard_task_id <- raw_task_lookup[repetition, 1L]
        candidate_results <- vector("list", length(K_candidates))
        for (candidate_id in seq_along(K_candidates)) {
          set_task_seed(
            2L,
            (raw_task_id - 1L) * length(K_candidates) + candidate_id
          )
          K <- K_candidates[candidate_id]
          model_results <- stats::setNames(
            vector("list", length(model_candidates)), model_candidates
          )
          diagnostics <- vector("list", length(model_candidates))
          for (model_id in seq_along(model_candidates)) {
            model <- model_candidates[model_id]
            labels <- raw_output[[raw_task_id]]$labels[[model]][[candidate_id]]
            standard <- raw_output[[standard_task_id]]$labels[[model]][[candidate_id]]
            overlap_labels <- labels[seq_len(overlap_count)]
            overlap_standard <- standard[seq_len(overlap_count)]
            if (split_id == 1L || K == 1L) {
              matched <- list(
                matched_labels = as.integer(labels),
                mismatch_rate = 0,
                mapping = seq_len(K),
                overlap = unclass(table(
                  factor(overlap_labels, levels = seq_len(K)),
                  factor(overlap_standard, levels = seq_len(K))
                ))
              )
            } else {
              matched <- match_labels(
                match_this = overlap_labels,
                standard = overlap_standard,
                K = K
              )
              matched$matched_labels <- unname(matched$mapping[labels])
            }
            matched_labels <- as.integer(matched$matched_labels)
            missing_match <- which(rowSums(matched$overlap) == 0L)
            missing_standard <- which(colSums(matched$overlap) == 0L)
            diagnostics[[model_id]] <- list(
              repetition = repetition,
              split = split_id,
              K = K,
              model = model,
              matching_method = matching_method,
              mismatch_rate = matched$mismatch_rate,
              missing_match_communities = paste(missing_match, collapse = ","),
              missing_standard_communities = paste(
                missing_standard, collapse = ","
              ),
              underidentified = length(missing_match) > 0L ||
                length(missing_standard) > 0L,
              effective_overlap_size = overlap_count
            )
            if (model == "SBM") {
              fit <- if (K == 1L) {
                denominator <- nrow(A_subgraph) * (nrow(A_subgraph) - 1L)
                matrix(
                  sum(A_subgraph) / denominator,
                  nrow = 1L,
                  dimnames = list("community_1", "community_1")
                )
              } else {
                do.call(
                  estimate_sbm,
                  c(
                    list(
                      A = A_subgraph,
                      g = matched_labels,
                      K = K,
                      directed = FALSE,
                      self_loops = FALSE,
                      validate_inputs = FALSE
                    ),
                    sbm_options$estimate
                  )
                )
              }
            # The historical estimator returned the observed diagonal value
            # for a singleton community, whereas estimate_sbm() correctly has
            # no loop-free within-community dyad and returns NA. Retain the
            # historical NETCROP behavior without changing estimate_sbm().
            missing_blocks <- which(is.na(fit), arr.ind = TRUE)
            if (nrow(missing_blocks) > 0L) {
              for (missing_id in seq_len(nrow(missing_blocks))) {
                block_row <- missing_blocks[missing_id, 1L]
                block_col <- missing_blocks[missing_id, 2L]
                if (block_row != block_col) {
                  stop(
                    "SBM estimation produced an unobserved off-diagonal block.",
                    call. = FALSE
                  )
                }
                block_nodes <- which(matched_labels == block_row)
                if (length(block_nodes) != 1L) {
                  stop(
                    "SBM estimation produced an unexpected missing block.",
                    call. = FALSE
                  )
                }
                fit[block_row, block_col] <-
                  as.numeric(A_subgraph[block_nodes, block_nodes])
              }
            }
              model_results[[model]] <- list(
                g = matched_labels[-seq_len(overlap_count)],
                B_hat = fit,
                psi_hat = NULL
              )
            } else {
              estimate_options <- dcbm_options$estimate
              row_norm <- if (dcbm_estimation_method == "spectral") {
                raw_output[[raw_task_id]]$row_norms[[candidate_id]]
              } else {
                NULL
              }
              fit <- if (K == 1L && dcbm_estimation_method == "plugin") {
                stabilizer <- if (is.null(estimate_options$stabilizer)) {
                  0.01
                } else {
                  estimate_options$stabilizer
                }
                B_hat <- matrix(
                  sum(A_subgraph) + stabilizer,
                  nrow = 1L,
                  dimnames = list("community_1", "community_1")
                )
                degree <- if (inherits(A_subgraph, "Matrix")) {
                  Matrix::rowSums(A_subgraph)
                } else {
                  rowSums(A_subgraph)
                }
                list(
                  B_hat = B_hat,
                  psi_hat = as.numeric(degree[-seq_len(overlap_count)] /
                    B_hat[1L, 1L])
                )
              } else if (K == 1L) {
                denominator <- sum(row_norm)^2
                B_hat <- matrix(
                  sum(A_subgraph) / denominator,
                  nrow = 1L,
                  dimnames = list("community_1", "community_1")
                )
                list(
                  B_hat = B_hat,
                  psi_hat = as.numeric(row_norm[-seq_len(overlap_count)])
                )
              } else {
                do.call(
                  estimate_dcbm,
                  c(
                    list(
                      A = A_subgraph,
                      g = matched_labels,
                      K = K,
                      row_norm = row_norm,
                      psi_omit = overlap_count,
                      validate_inputs = FALSE
                    ),
                    estimate_options
                  )
                )
              }
              model_results[[model]] <- list(
                g = matched_labels[-seq_len(overlap_count)],
                B_hat = fit$B_hat,
                psi_hat = fit$psi_hat
              )
            }
          }
          candidate_results[[candidate_id]] <- list(
            models = model_results,
            diagnostics = diagnostics
          )
        }
        list(candidates = candidate_results)
      },
      ncores = estimation_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE,
      manage_future_plan = manage_future_plan,
      future_packages = worker_future_packages
    )
  })
  if (retain_intermediates == "minimal") {
    raw_output <- NULL
  }

  estimate_lookup <- function(repetition, split_id, candidate_id, model) {
    raw_task_id <- raw_task_lookup[repetition, split_id]
    estimation_output[[raw_task_id]]$candidates[[candidate_id]]$models[[model]]
  }
  averaged_estimates <- vector("list", nrep)
  for (repetition in seq_len(nrep)) {
    averaged_estimates[[repetition]] <- vector(
      "list",
      length(K_candidates)
    )
    for (candidate_id in seq_along(K_candidates)) {
      averaged_estimates[[repetition]][[candidate_id]] <- stats::setNames(
        vector("list", length(model_candidates)),
        model_candidates
      )
      for (model in model_candidates) {
        B_sum <- 0
        for (split_id in seq_len(num_subnetworks)) {
          B_sum <- B_sum + estimate_lookup(
            repetition,
            split_id,
            candidate_id,
            model
          )$B_hat / num_subnetworks
        }
        averaged_estimates[[repetition]][[candidate_id]][[model]] <- B_sum
      }
    }
  }

  pair_grid <- t(utils::combn(num_subnetworks, 2L))
  loss_grid <- expand.grid(
    pair_id = seq_len(nrow(pair_grid)),
    candidate_id = seq_along(K_candidates),
    repetition = seq_len(nrep),
    KEEP.OUT.ATTRS = FALSE
  )
  loss_grid <- loss_grid[
    order(loss_grid$repetition, loss_grid$candidate_id, loss_grid$pair_id),
  ]
  loss_accepts_prevalidation <- vapply(
    loss_functions,
    function(loss_function) "validate_inputs" %in% names(formals(loss_function)),
    logical(1)
  )
  if (verbose) {
    message(
      "Stage 3/3: evaluating ", loss_task_count,
      " held-out subnetwork pairs with ", loss_ncores, " worker(s)."
    )
  }
  loss_time <- system.time({
    loss_output <- uni_mclapply(
      seq_len(loss_task_count),
      function(task_id) {
        if (input_was_sparse &&
            !requireNamespace("Matrix", quietly = TRUE)) {
          stop("The Matrix package is required for sparse A.", call. = FALSE)
        }
        set_task_seed(3L, task_id)
        repetition <- loss_grid$repetition[task_id]
        candidate_id <- loss_grid$candidate_id[task_id]
        pair_id <- loss_grid$pair_id[task_id]
        split_left <- pair_grid[pair_id, 1L]
        split_right <- pair_grid[pair_id, 2L]
        nodes_left <- splitter$nonoverlap_pieces[[repetition]][[split_left]]
        nodes_right <- splitter$nonoverlap_pieces[[repetition]][[split_right]]
        A_test <- A[nodes_left, nodes_right, drop = FALSE]
        A_test_numeric <- as.numeric(A_test)
        records <- matrix(
          NA_real_,
          nrow = length(model_candidates) * length(losses),
          ncol = 7L
        )
        record_id <- 1L
        for (model_id in seq_along(model_candidates)) {
          model <- model_candidates[model_id]
          left_raw_id <- raw_task_lookup[repetition, split_left]
          right_raw_id <- raw_task_lookup[repetition, split_right]
          left <- estimation_output[[left_raw_id]]$candidates[[candidate_id]]$
            models[[model]]
          right <- estimation_output[[right_raw_id]]$candidates[[candidate_id]]$
            models[[model]]
          B_hat <- averaged_estimates[[repetition]][[candidate_id]][[model]]
          P_hat <- B_hat[left$g, right$g, drop = FALSE]
          if (model == "DCBM") {
            P_hat <- P_hat * tcrossprod(left$psi_hat, right$psi_hat)
            P_hat <- clip_probabilities(P_hat, eps = 1e-6)
          }
          P_hat_numeric <- as.numeric(P_hat)
          for (loss_id in seq_along(losses)) {
            loss_name <- losses[loss_id]
            loss_arguments <- list(A_test_numeric, P_hat_numeric)
            if (loss_accepts_prevalidation[loss_id]) {
              loss_arguments$validate_inputs <- FALSE
            }
            loss_value <- do.call(
              loss_functions[[loss_name]], loss_arguments
            ) / pair_count
            if (length(loss_value) != 1L || !is.numeric(loss_value) ||
                is.na(loss_value) || !is.finite(loss_value)) {
              stop(
                "Non-finite loss for repetition ", repetition,
                ", K = ", K_candidates[candidate_id],
                ", model = ", model,
                ", loss = ", loss_name,
                ", split pair = ", split_left, "-", split_right,
                ".",
                call. = FALSE
              )
            }
            records[record_id, ] <- c(
              repetition,
              K_candidates[candidate_id],
              model_id,
              loss_id,
              split_left,
              split_right,
              as.numeric(loss_value)
            )
            record_id <- record_id + 1L
          }
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
  loss_matrix <- do.call(rbind, loss_output)
  cv_all_loss <- data.frame(
    repetition = as.integer(loss_matrix[, 1L]),
    K = as.integer(loss_matrix[, 2L]),
    model = model_candidates[as.integer(loss_matrix[, 3L])],
    loss = losses[as.integer(loss_matrix[, 4L])],
    split_left = as.integer(loss_matrix[, 5L]),
    split_right = as.integer(loss_matrix[, 6L]),
    loss_value = loss_matrix[, 7L],
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
    data.frame(
      repetition = x$repetition[1L],
      K = x$K[1L],
      model = x$model[1L],
      loss = x$loss[1L],
      average_loss = sum(x$loss_value),
      stringsAsFactors = FALSE
    )
  }))
  rownames(cv_loss) <- NULL
  cv_loss <- cv_loss[order(
    cv_loss$repetition,
    match(cv_loss$loss, losses),
    match(cv_loss$model, model_candidates),
    match(cv_loss$K, K_candidates)
  ), ]
  best_model_cv <- do.call(rbind, lapply(seq_len(nrep), function(repetition) {
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
        model = candidates$model[best_id],
        K = candidates$K[best_id],
        average_loss = candidates$average_loss[best_id],
        best_model = paste0(
          candidates$model[best_id],
          "-",
          candidates$K[best_id]
        ),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(best_model_cv) <- NULL
  overall_best <- do.call(rbind, lapply(losses, function(loss_name) {
    selected <- best_model_cv[best_model_cv$loss == loss_name, , drop = FALSE]
    best <- modal(selected$best_model)
    data.frame(
      loss = loss_name,
      best_model = best,
      mean_average_loss = mean(selected$average_loss),
      stringsAsFactors = FALSE
    )
  }))
  rownames(overall_best) <- NULL
  diagnostic_records <- unlist(lapply(estimation_output, function(raw_result) {
    unlist(lapply(
      raw_result$candidates,
      `[[`,
      "diagnostics"
    ), recursive = FALSE)
  }), recursive = FALSE)
  matching_diagnostics <- do.call(rbind, lapply(diagnostic_records, function(x) {
    as.data.frame(x, stringsAsFactors = FALSE)
  }))
  rownames(matching_diagnostics) <- NULL
  matching_diagnostics$repetition <- as.integer(
    matching_diagnostics$repetition
  )
  matching_diagnostics$split <- as.integer(matching_diagnostics$split)
  matching_diagnostics$K <- as.integer(matching_diagnostics$K)
  matching_diagnostics$mismatch_rate <- as.numeric(
    matching_diagnostics$mismatch_rate
  )
  matching_diagnostics$underidentified <- as.logical(
    matching_diagnostics$underidentified
  )
  matching_diagnostics$effective_overlap_size <- as.integer(
    matching_diagnostics$effective_overlap_size
  )
  matching_summary <- list(
    underidentified_matches = sum(matching_diagnostics$underidentified),
    total_matches = nrow(matching_diagnostics)
  )
  if (retain_intermediates == "minimal") {
    estimation_output <- NULL
    averaged_estimates <- NULL
    matching_diagnostics <- NULL
  }
  candidate_models <- expand.grid(
    model = model_candidates,
    K = K_candidates,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  candidate_models <- candidate_models[order(
    match(candidate_models$model, model_candidates),
    match(candidate_models$K, K_candidates)
  ), ]
  historical_train_proportion <-
    num_subnetworks * subgraph_size * (subgraph_size - 1) /
    (n * (n - 1))
  historical_test_proportion <-
    num_subnetworks * (num_subnetworks - 1) * piece_size^2 /
    (2 * n * (n - 1))
  unordered_test_proportion <-
    choose(num_subnetworks, 2) * piece_size^2 / choose(n, 2)

  out <- list(
    algorithm = "NETCROP",
    cv_loss = cv_loss,
    cv_all_loss = cv_all_loss,
    best_model_cv = best_model_cv,
    overall_best = overall_best,
    candidate_models = candidate_models,
    num_subnetworks = num_subnetworks,
    requested_overlap_size = overlap_size,
    effective_overlap_size = effective_overlap_size,
    piece_size = piece_size,
    historical_train_proportion = historical_train_proportion,
    historical_test_proportion = historical_test_proportion,
    unordered_test_proportion = unordered_test_proportion,
    nrep = nrep,
    losses = losses,
    model_candidates = model_candidates,
    sbm_est_options = sbm_options,
    dcbm_est_options = dcbm_options,
    resolved_spectral_options = active_spectral_options,
    splitter = splitter,
    raw_output = if (retain_intermediates == "all") raw_output else NULL,
    estimation_output = if (retain_intermediates == "all") {
      estimation_output
    } else {
      NULL
    },
    averaged_estimates = if (retain_intermediates == "all") {
      averaged_estimates
    } else {
      NULL
    },
    matching_diagnostics = if (retain_intermediates == "all") {
      matching_diagnostics
    } else {
      NULL
    },
    matching_summary = matching_summary,
    retain_intermediates = retain_intermediates,
    matching_method = matching_method,
    ncores = list(
      requested = ncores,
      clustering = raw_ncores,
      estimation = estimation_ncores,
      loss = loss_ncores
    ),
    timing = c(
      clustering = unname(raw_time["elapsed"]),
      estimation = unname(estimation_time["elapsed"]),
      loss = unname(loss_time["elapsed"]),
      total = unname(raw_time["elapsed"] + estimation_time["elapsed"] +
                       loss_time["elapsed"])
    ),
    ram_preflight = ram_report,
    seed = seed,
    call = call
  )
  class(out) <- "netcrop_blockmodel"
  out
}

# Print the selected block model for every requested loss.
#' Print the selected block model for every requested loss.
#' @rdname netcrop_blockmodel
#' @export
print.netcrop_blockmodel <- function(x, ...) {
  algorithm <- if (is.null(x$algorithm)) "NETCROP" else toupper(x$algorithm)
  cat(algorithm, "results for block models\n")
  cat(strrep("-", nchar(algorithm) + 25L), "\n", sep = "")
  print(x$overall_best, row.names = FALSE)
  invisible(x)
}

# Summarize a NETCROP block-model fit.
#' Summarize selections, split diagnostics, and timing.
#' @rdname netcrop_blockmodel
#' @export
summary.netcrop_blockmodel <- function(object, ...) {
  algorithm <- if (is.null(object$algorithm)) {
    "NETCROP"
  } else {
    toupper(object$algorithm)
  }
  if (algorithm == "ECV") {
    result <- list(
      algorithm = algorithm,
      call = object$call,
      candidate_models = object$candidate_models,
      nrep = object$nrep,
      completed_repetitions = object$completed_repetitions,
      valid_repetitions = object$valid_repetitions,
      cv = object$cv,
      train_proportion = object$train_proportion,
      holdout_proportion = object$holdout_proportion,
      holdout_count = object$holdout_count,
      best_model_cv = object$best_model_cv,
      overall_best = object$overall_best,
      selection_diagnostics = object$selection_diagnostics,
      failure_diagnostics = object$failure_diagnostics,
      ncores = object$ncores,
      timing = object$timing
    )
    class(result) <- "summary.netcrop_blockmodel"
    return(result)
  }
  if (algorithm == "NCV") {
    result <- list(
      algorithm = algorithm,
      call = object$call,
      candidate_models = object$candidate_models,
      nrep = object$nrep,
      completed_repetitions = object$completed_repetitions,
      valid_repetitions = object$valid_repetitions,
      cv = object$cv,
      fold_sizes = object$fold_sizes,
      dc_est = object$dc_est,
      tau = object$tau,
      use_laplacian = object$use_laplacian,
      best_model_cv = object$best_model_cv,
      overall_best = object$overall_best,
      failure_diagnostics = object$failure_diagnostics,
      ncores = object$ncores,
      timing = object$timing
    )
    class(result) <- "summary.netcrop_blockmodel"
    return(result)
  }
  matching_summary <- object$matching_summary
  if (is.null(matching_summary) && !is.null(object$matching_diagnostics)) {
    matching_summary <- list(
      underidentified_matches = sum(
        object$matching_diagnostics$underidentified
      ),
      total_matches = nrow(object$matching_diagnostics)
    )
  }
  result <- list(
    algorithm = algorithm,
    call = object$call,
    candidate_models = object$candidate_models,
    nrep = object$nrep,
    num_subnetworks = object$num_subnetworks,
    requested_overlap_size = object$requested_overlap_size,
    effective_overlap_size = object$effective_overlap_size,
    piece_size = object$piece_size,
    unordered_test_proportion = object$unordered_test_proportion,
    best_model_cv = object$best_model_cv,
    overall_best = object$overall_best,
    underidentified_matches = matching_summary$underidentified_matches,
    total_matches = matching_summary$total_matches,
    matching_method = object$matching_method,
    ncores = object$ncores,
    timing = object$timing
  )
  class(result) <- "summary.netcrop_blockmodel"
  result
}

# Print a NETCROP block-model summary.
#' Print a `summary.netcrop_blockmodel` object.
#' @rdname netcrop_blockmodel
#' @export
print.summary.netcrop_blockmodel <- function(x, ...) {
  algorithm <- if (is.null(x$algorithm)) "NETCROP" else toupper(x$algorithm)
  if (algorithm == "ECV") {
    cat("Summary of ECV block-model selection\n")
    cat("------------------------------------\n")
    cat("Candidate models:", paste(unique(x$candidate_models$model),
                                    collapse = ", "), "\n")
    cat("Candidate K:", paste(unique(x$candidate_models$K),
                               collapse = ", "), "\n")
    cat(
      "Completed repetitions:", x$completed_repetitions,
      "of", x$nrep, "requested\n"
    )
    cat("CV folds per repetition:", x$cv, "\n")
    cat(sprintf("Training proportion: %.2f%%\n", 100 * x$train_proportion))
    cat(sprintf("Holdout proportion: %.2f%%\n", 100 * x$holdout_proportion))
    cat("Held-out unordered pairs per fold:", x$holdout_count, "\n")
    cat("Best models per repetition:\n")
    print(x$best_model_cv, row.names = FALSE)
    cat("Overall best models:\n")
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
  if (algorithm == "NCV") {
    cat("Summary of NCV block-model selection\n")
    cat("------------------------------------\n")
    cat("Candidate models:", paste(unique(x$candidate_models$model),
                                    collapse = ", "), "\n")
    cat("Candidate K:", paste(unique(x$candidate_models$K),
                               collapse = ", "), "\n")
    cat(
      "Completed repetitions:", x$completed_repetitions,
      "of", x$nrep, "requested\n"
    )
    cat("CV folds per repetition:", x$cv, "\n")
    cat("Fold sizes:", paste(x$fold_sizes, collapse = ", "), "\n")
    cat("DCBM estimator:", x$dc_est, "\n")
    cat("Laplacian scaling:", x$use_laplacian, "\n")
    cat("Best models per repetition:\n")
    print(x$best_model_cv, row.names = FALSE)
    cat("Overall best models:\n")
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
  cat("Summary of NETCROP block-model selection\n")
  cat("------------------------------------------\n")
  cat("Candidate models:", paste(unique(x$candidate_models$model),
                                  collapse = ", "), "\n")
  cat("Candidate K:", paste(unique(x$candidate_models$K),
                             collapse = ", "), "\n")
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
  cat(
    "Underidentified overlap matches:", x$underidentified_matches,
    "of", x$total_matches, "\n"
  )
  cat("Best models per repetition:\n")
  print(x$best_model_cv, row.names = FALSE)
  cat("Overall best models:\n")
  print(x$overall_best, row.names = FALSE)
  cat("Timing (seconds):\n")
  print(x$timing)
  invisible(x)
}

# Plot CV loss curves, optionally aggregating across repetitions.
#' Plot loss curves, optionally aggregating repetitions or comparing fits.
#' @rdname netcrop_blockmodel
#' @export
plot.netcrop_blockmodel <- function(x, aggregate = TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for plotting.", call. = FALSE)
  }
  dots <- list(...)
  dot_is_result <- vapply(
    dots, inherits, logical(1), what = "netcrop_blockmodel"
  )
  if (inherits(aggregate, "netcrop_blockmodel") || any(dot_is_result)) {
    if (!exists(
      "plot_blockmodel_comparison", mode = "function", inherits = TRUE, envir = environment()
    )) {
      stop(
        "Required comparison helpers are unavailable; reinstall netOP.",
        call. = FALSE
      )
    }
    comparison_results <- c(
      list(x),
      if (inherits(aggregate, "netcrop_blockmodel")) list(aggregate) else NULL,
      dots[dot_is_result]
    )
    comparison_options <- dots[!dot_is_result]
    loss_scale <- if ("loss_scale" %in% names(comparison_options)) {
      comparison_options$loss_scale
    } else {
      "relative"
    }
    return(do.call(
      plot_blockmodel_comparison,
      c(comparison_results, list(loss_scale = loss_scale))
    ))
  }
  if (length(aggregate) != 1L || !is.logical(aggregate) || is.na(aggregate)) {
    stop("aggregate must be TRUE or FALSE.", call. = FALSE)
  }
  K_breaks <- sort(unique(x$cv_loss$K))
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
          x = K,
          y = average_loss,
          group = interaction(repetition, model),
          color = model,
          shape = repetition
        )
      ) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::scale_x_continuous(breaks = K_breaks) +
        ggplot2::facet_wrap(~loss, scales = "free_y") +
        ggplot2::labs(
          title = paste(algorithm, "CV loss by number of communities"),
          x = "Number of communities (K)",
          y = "Average CV loss",
          color = "Model",
          shape = "Repetition"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")
    )
  }
  grouping <- interaction(
    valid_cv_loss$K,
    valid_cv_loss$model,
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
      K = z$K[1L],
      model = z$model[1L],
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
      x = K,
      y = mean_loss,
      ymin = lower,
      ymax = upper,
      group = model,
      color = model,
      fill = model
    )
  ) +
    ggplot2::geom_ribbon(alpha = 0.2, color = NA) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = K_breaks) +
    ggplot2::facet_wrap(~loss, scales = "free_y") +
    ggplot2::labs(
      title = paste(algorithm, "CV loss by number of communities"),
      x = "Number of communities (K)",
      y = "Mean CV loss (plus or minus one SD)",
      color = "Model",
      fill = "Model"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

# library(Matrix)
# system.time(
#   net <- generate_dcbm(n = 10000, K = 5, ncores = 5, average_degree = 80)
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
#     K = 20,
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
#     K_candidates = 1:20,
#     num_subnetworks = 3,
#     overlap_size = 8002,
#     nrep = 1,
#     losses = "mse",
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
