# SONNET fitting and S3 methods

# Fit SONNET using overlapping subnetwork spectral clusterings.
#
# This implementation preserves the original partition, alignment, and modal
# aggregation flow. Additional outputs expose every intermediate result without
# changing how the final labels are computed. Arguments in ... are forwarded
# to spectral_cluster(); A, U, and K are reserved by sonnet().
.sonnet_fit <- function(
    A,
    K,
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
    share_overlap = TRUE,
    parameter_select_options = list(),
    ...) {
  call <- match.call()
  matching_method <- match.arg(matching_method)
  required_helpers <- c(
    "uni_mclapply",
    "spectral_cluster",
    "label_match_greedy",
    "label_match_brute_force",
    "modal",
    "sonnet_param_select",
    "sonnet_splitter",
    "sonnet_splitter_independent"
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
  if (is.null(dim(A)) ||
      length(dim(A)) != 2L ||
      nrow(A) != ncol(A) ||
      nrow(A) < 1L ||
      !(is.numeric(A) || inherits(A, "Matrix")) ||
      any(!is.finite(A))) {
    stop("A must be a non-empty finite square numeric matrix.", call. = FALSE)
  }
  n <- nrow(A)
  if (is.null(num_subnetworks) || is.null(overlap_size)) {
    if (!is.list(parameter_select_options) ||
        (length(parameter_select_options) > 0L &&
         (is.null(names(parameter_select_options)) ||
          any(!nzchar(names(parameter_select_options)))))) {
      stop("parameter_select_options must be a named list.", call. = FALSE)
    }
    protected_parameter_options <- intersect(
      names(parameter_select_options),
      c("n", "K", "ncores")
    )
    if (length(protected_parameter_options) > 0L) {
      stop(
        "parameter_select_options cannot override ",
        paste(protected_parameter_options, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    unknown_parameter_options <- setdiff(
      names(parameter_select_options),
      names(formals(sonnet_param_select))
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
      sonnet_param_select,
      utils::modifyList(
        list(n = n, K = K, theta = 3, q = 0.005, ncores = ncores),
        parameter_select_options,
        keep.null = TRUE
      )
    )
    if (is.null(num_subnetworks)) {
      num_subnetworks <- selected_parameters$num_subnetworks
    }
    if (is.null(overlap_size)) {
      overlap_size <- selected_parameters$overlap_size
    }
  }
  integer_inputs <- list(
    K = K,
    num_subnetworks = num_subnetworks,
    overlap_size = overlap_size,
    extra_nrep = extra_nrep,
    ncores = ncores
  )
  invalid_integer <- vapply(
    integer_inputs,
    function(value) {
      length(value) != 1L ||
        !is.numeric(value) ||
        is.na(value) ||
        !is.finite(value) ||
        value != floor(value)
    },
    logical(1)
  )
  if (any(invalid_integer)) {
    stop(
      paste(names(integer_inputs)[invalid_integer], collapse = ", "),
      " must be integer-valued scalars.",
      call. = FALSE
    )
  }
  K <- as.integer(K)
  num_subnetworks <- as.integer(num_subnetworks)
  overlap_size <- as.integer(overlap_size)
  extra_nrep <- as.integer(extra_nrep)
  ncores <- as.integer(ncores)
  if (K < 1L || K > n) {
    stop("K must be between 1 and nrow(A).", call. = FALSE)
  }
  if (num_subnetworks < 2L) {
    stop("num_subnetworks must be at least 2.", call. = FALSE)
  }
  if (overlap_size < 1L || overlap_size > n) {
    stop("overlap_size must be between 1 and nrow(A).", call. = FALSE)
  }
  if (extra_nrep < 0L) {
    stop("extra_nrep must be non-negative.", call. = FALSE)
  }
  if (ncores < 1L) {
    stop("ncores must be positive.", call. = FALSE)
  }
  logical_inputs <- list(
    confirm_large = confirm_large,
    verbose = verbose,
    force_windows = force_windows,
    ram_check = ram_check,
    share_overlap = share_overlap
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
  spectral_arguments <- list(...)
  if (length(spectral_arguments) > 0L &&
      (is.null(names(spectral_arguments)) ||
       any(!nzchar(names(spectral_arguments))))) {
    stop("All arguments in ... must be named.", call. = FALSE)
  }
  reserved_spectral <- intersect(
    names(spectral_arguments),
    c("A", "U", "K")
  )
  if (length(reserved_spectral) > 0L) {
    stop(
      "... cannot override ",
      paste(reserved_spectral, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  allowed_spectral <- setdiff(
    names(formals(spectral_cluster)),
    c("A", "U", "K")
  )
  unknown_spectral <- setdiff(names(spectral_arguments), allowed_spectral)
  if (length(unknown_spectral) > 0L) {
    stop(
      "Unknown spectral_cluster() argument(s) in ...: ",
      paste(unknown_spectral, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  resolved_spectral_arguments <- utils::modifyList(
    list(
      laplacian = FALSE,
      normalize_laplacian = TRUE,
      regularize_tau = 0,
      handle_zero_degree_nodes = "none",
      row_normalize = FALSE,
      spectral_method = "eigen",
      spectral_engine = "RSpectra",
      spectral_options = list(),
      cluster_engine = "clara",
      cluster_options = list()
    ),
    spectral_arguments,
    keep.null = TRUE
  )
  
  m <- floor((n - overlap_size) / num_subnetworks)
  remainder_count <- n - overlap_size - num_subnetworks * m
  effective_overlap_size <- overlap_size + remainder_count
  first_subgraph_size <- effective_overlap_size + m
  if (first_subgraph_size < K) {
    stop(
      "Each first-repetition subgraph must contain at least K nodes.",
      call. = FALSE
    )
  }
  if (share_overlap && extra_nrep > 0L && m < K) {
    stop(
      "Each extra-repetition subgraph must contain at least K nodes.",
      call. = FALSE
    )
  }
  
  if (matching_method == "brute_force" && K > 8L && confirm_large) {
    warning_text <- paste0(
      "Brute-force matching will be used repeatedly across ",
      factorial(K),
      " permutations per match and may be extremely slow and memory ",
      "consuming. label_match_greedy() is preferred."
    )
    if (!interactive()) {
      stop(
        warning_text,
        " Set confirm_large = FALSE to bypass this prompt.",
        call. = FALSE
      )
    }
    choice <- utils::menu(
      c("No, cancel", "Yes, continue"),
      title = paste0(warning_text, " Do you want to continue?")
    )
    if (choice != 2L) {
      stop("SONNET was cancelled before computation.", call. = FALSE)
    }
  }
  
  total_start <- proc.time()[["elapsed"]]
  if (verbose) {
    message(
      "Starting SONNET with n = ", n,
      ", K = ", K,
      ", subnetworks = ", num_subnetworks,
      ", extra repetitions = ", extra_nrep, "."
    )
  }
  splitter_time <- system.time({
    splitter <- if (share_overlap) {
      sonnet_splitter(
        n = n,
        num_subnetworks = num_subnetworks,
        overlap_size = overlap_size,
        extra_nrep = extra_nrep,
        m = m,
        seed = seed
      )
    } else {
      sonnet_splitter_independent(
        n = n,
        num_subnetworks = num_subnetworks,
        overlap_size = overlap_size,
        extra_nrep = extra_nrep,
        m = m,
        seed = seed
      )
    }
  })[["elapsed"]]
  overlap_nodes <- if (share_overlap) {
    splitter$overlap_nodes
  } else {
    splitter$overlap_nodes[[1L]]
  }
  nonoverlap_nodes <- if (share_overlap) {
    splitter$nonoverlap_nodes
  } else {
    splitter$nonoverlap_nodes[[1L]]
  }
  if (verbose && splitter$remainder_count > 0L) {
    message(
      "Augmented overlap_size from ",
      overlap_size,
      " to ",
      splitter$overlap_size,
      " using ",
      splitter$remainder_count,
      " remainder node(s)."
    )
  }
  
  set_task_seed <- function(offset) {
    if (!is.null(seed)) {
      task_seed <- as.integer(
        (seed + as.double(offset)) %% .Machine$integer.max
      )
      set.seed(task_seed)
    }
    invisible(NULL)
  }
  match_labels <- function(match_this, standard, return_mapping = FALSE) {
    if (matching_method == "greedy") {
      return(label_match_greedy(
        match_this = match_this,
        standard = standard,
        K = K,
        algorithm = "greedy",
        return_mapping = return_mapping
      ))
    }
    if (matching_method == "hungarian") {
      return(label_match_greedy(
        match_this = match_this,
        standard = standard,
        K = K,
        algorithm = "hungarian",
        return_mapping = return_mapping
      ))
    }
    label_match_brute_force(
      match_this = match_this,
      standard = standard,
      K = K,
      return_mapping = return_mapping,
      confirm_large = TRUE
    )
  }
  membership_from_labels <- function(label_matrix, fallback = NULL) {
    prediction_count <- rowSums(label_matrix > 0L)
    membership <- matrix(0, nrow = n, ncol = K)
    for (community in seq_len(K)) {
      membership[, community] <- rowSums(label_matrix == community)
    }
    predicted <- prediction_count > 0L
    membership[predicted, ] <-
      membership[predicted, , drop = FALSE] / prediction_count[predicted]
    if (any(!predicted)) {
      if (is.null(fallback) ||
          any(fallback[!predicted] < 1L | fallback[!predicted] > K)) {
        stop(
          "Some nodes have no valid predicted label in a repetition.",
          call. = FALSE
        )
      }
      missing_rows <- which(!predicted)
      membership[cbind(missing_rows, fallback[missing_rows])] <- 1
    }
    rownames(membership) <- rownames(A)
    colnames(membership) <- paste0("community_", seq_len(K))
    membership
  }

  check_parallel_spectral_ram <- function(
      matrix_dimension,
      worker_count) {
    if (!ram_check) {
      return(invisible(NULL))
    }
    required_ram_helpers <- c(
      "available_ram",
      "format_bytes",
      "estimate_rspectra_ram",
      "estimate_dense_eigen_ram",
      "estimate_partial_svd_ram",
      "estimate_dense_svd_ram",
      "estimate_matrix_product_ram",
      "report_ram_preflight",
      "report_ram_formula"
    )
    missing_ram_helpers <- required_ram_helpers[!vapply(
      required_ram_helpers,
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

    method <- resolved_spectral_arguments$spectral_method
    engine <- resolved_spectral_arguments$spectral_engine
    options <- resolved_spectral_arguments$spectral_options
    dense_laplacian <- isTRUE(resolved_spectral_arguments$laplacian)
    dense_input <- dense_laplacian ||
      resolved_spectral_arguments$regularize_tau > 0
    force_engine <- isTRUE(options$force_engine)
    safe_d_multiplier <- if (is.null(options$safe_d_multiplier)) {
      1
    } else {
      options$safe_d_multiplier
    }
    d_compute <- min(
      ceiling(safe_d_multiplier * K),
      matrix_dimension
    )
    engine_lower <- tolower(engine)
    if (method == "eigen" && engine_lower == "irlba" && !force_engine) {
      engine_lower <- if (requireNamespace("RSpectra", quietly = TRUE)) {
        "rspectra"
      } else {
        "base"
      }
    }
    if (matrix_dimension <= 200L && engine_lower != "base" &&
        !force_engine) {
      engine_lower <- "base"
    }
    package_name <- if (engine_lower == "rspectra") "RSpectra" else engine
    if (engine_lower != "base" && !force_engine &&
        !requireNamespace(package_name, quietly = TRUE)) {
      engine_lower <- "base"
    }
    if (engine_lower != "base" && d_compute >= matrix_dimension) {
      if (K >= matrix_dimension) {
        engine_lower <- "base"
      } else {
        d_compute <- matrix_dimension - 1L
      }
    }

    if (method == "eigen") {
      if (engine_lower == "base") {
        expected_per_worker <- estimate_dense_eigen_ram(
          n = matrix_dimension,
          input_already_counted = !dense_input
        )
        operation <- "SONNET dense eigendecomposition"
      } else {
        ncv_estimate <- min(
          matrix_dimension,
          max(2L * d_compute + 1L, 20L)
        )
        expected_per_worker <- estimate_rspectra_ram(
          n = matrix_dimension,
          K = d_compute,
          ncv = ncv_estimate
        )
        if (dense_input) {
          expected_per_worker <- expected_per_worker +
            8 * as.double(matrix_dimension)^2
        }
        operation <- "SONNET RSpectra eigendecomposition"
      }
    } else if (engine_lower == "base") {
      expected_per_worker <- estimate_dense_svd_ram(
        n = matrix_dimension,
        p = matrix_dimension,
        nu = K,
        nv = 0L,
        input_already_counted = !dense_input
      )
      operation <- "SONNET dense singular decomposition"
    } else {
      ncv_estimate <- min(
        matrix_dimension,
        max(2L * d_compute + 1L, 20L)
      )
      expected_per_worker <- estimate_partial_svd_ram(
        n = matrix_dimension,
        p = matrix_dimension,
        K = d_compute,
        ncv = ncv_estimate
      )
      if (dense_input) {
        expected_per_worker <- expected_per_worker +
          8 * as.double(matrix_dimension)^2
      }
      operation <- paste0(
        "SONNET ",
        if (engine_lower == "rspectra") "RSpectra" else "irlba",
        " singular decomposition"
      )
    }

    expected_total <- expected_per_worker * worker_count
    available_bytes <- available_ram()
    report <- paste0(
      operation,
      " for matrices up to ", matrix_dimension, " x ", matrix_dimension,
      if (dense_laplacian) " with a dense Laplacian input" else "",
      ": expected RAM use ", format_bytes(expected_per_worker),
      " x ", worker_count, " core(s) = ", format_bytes(expected_total),
      " out of ", format_bytes(available_bytes), " available RAM."
    )
    if (is.na(available_bytes) || expected_total > available_bytes) {
      warning(report, call. = FALSE)
    } else {
      message(report)
    }
    ram_formula_terms <- list(list(
      estimated_bytes = expected_per_worker,
      sequential_count = 1L,
      parallel_count = worker_count,
      label = "spectral decomposition"
    ))
    if (dense_laplacian) {
      multiplication_count <- if (isTRUE(
        resolved_spectral_arguments$normalize_laplacian
      )) 2L else 1L
      multiplication_ram <- estimate_matrix_product_ram(
        matrix_dimension,
        matrix_dimension,
        matrix_dimension
      )
      report_ram_preflight(
        estimated_bytes = multiplication_ram,
        operation = "SONNET largest matrix operation",
        operation_count = multiplication_count,
        ncores = worker_count,
        detail = paste0(
          matrix_dimension, " x ", matrix_dimension,
          " dense Laplacian construction"
        )
      )
      ram_formula_terms <- c(ram_formula_terms, list(list(
        estimated_bytes = multiplication_ram,
        sequential_count = multiplication_count,
        parallel_count = worker_count,
        label = "largest dense Laplacian matrix operation"
      )))
    }
    full_adjacency_copy <- as.numeric(utils::object.size(A))
    ram_formula_terms <- c(ram_formula_terms, list(list(
      estimated_bytes = full_adjacency_copy,
      sequential_count = 1L,
      parallel_count = worker_count,
      label = "full adjacency copy"
    )))
    repetition_count <- extra_nrep + 1L
    task_count <- num_subnetworks * repetition_count
    label_storage <-
      3 * 4 * as.double(n) * task_count +
      8 * as.double(n) * K * repetition_count +
      4 * as.double(n) * repetition_count
    ram_formula_terms <- c(ram_formula_terms, list(list(
      estimated_bytes = label_storage,
      sequential_count = 1L,
      parallel_count = 1L,
      label = "label and membership storage"
    )))
    report_ram_formula(
      terms = ram_formula_terms,
      operation = "SONNET conservative combined preflight",
      detail = paste0(
        "explicit sequential and parallel factors; retained result storage is ",
        "not multiplied by worker count"
      )
    )
    invisible(list(
      expected_per_worker = expected_per_worker,
      worker_count = worker_count,
      expected_total = expected_total,
      available = available_bytes
    ))
  }
  
  total_tasks <- num_subnetworks * (extra_nrep + 1L)
  clustering_ncores <- min(total_tasks, ncores)
  task_dimensions <- vapply(seq_len(total_tasks), function(task_id) {
    repetition <- ceiling(task_id / num_subnetworks)
    position <- task_id %% num_subnetworks + 1L
    length(splitter$subnetworks[[repetition]][[position]])
  }, integer(1))
  check_parallel_spectral_ram(
    matrix_dimension = max(task_dimensions),
    worker_count = clustering_ncores
  )
  if (verbose) {
    message(
      "Clustering ", total_tasks,
      " subgraphs with ", clustering_ncores,
      " worker(s)."
    )
  }
  worker_spectral_arguments <- spectral_arguments
  if ("spectral_options" %in% names(worker_spectral_arguments) &&
      is.list(worker_spectral_arguments$spectral_options)) {
    worker_spectral_arguments$spectral_options$ram_check <- NULL
  }
  clustering_time <- system.time({
    clustered <- uni_mclapply(
      seq_len(total_tasks),
      function(task_id) {
        repetition <- ceiling(task_id / num_subnetworks)
        position <- task_id %% num_subnetworks + 1L
        index <- splitter$subnetworks[[repetition]][[position]]
        set_task_seed(1000L + task_id)
        fit <- do.call(
          spectral_cluster,
          c(
            list(
              A = A[index, index, drop = FALSE],
              K = K,
              ram_check = FALSE
            ),
            worker_spectral_arguments
          )
        )
        g_hat <- fit$g_hat
        if (length(g_hat) != length(index) ||
            anyNA(g_hat) ||
            any(!is.finite(g_hat)) ||
            any(g_hat < 1L | g_hat > K) ||
            any(g_hat != floor(g_hat))) {
          stop(
            "spectral_cluster() returned invalid g_hat for subgraph task ",
            task_id,
            ".",
            call. = FALSE
          )
        }
        labels <- integer(n)
        labels[index] <- as.integer(g_hat)
        labels
      },
      ncores = clustering_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE
    )
    raw_subnetwork_labels <- do.call(cbind, clustered)
  })[["elapsed"]]
  
  if (share_overlap) {
  first_alignment_time <- system.time({
    first_alignment <- lapply(2:num_subnetworks, function(column_id) {
      match <- match_labels(
        match_this = raw_subnetwork_labels[overlap_nodes, column_id],
        standard = raw_subnetwork_labels[overlap_nodes, 1L],
        return_mapping = TRUE
      )
      aligned <- integer(n)
      observed <- raw_subnetwork_labels[, column_id] != 0L
      aligned[observed] <- unname(
        match$mapping[raw_subnetwork_labels[observed, column_id]]
      )
      list(labels = aligned, match = match)
    })
    aligned_first <- do.call(
      cbind,
      c(
        list(raw_subnetwork_labels[, 1L]),
        lapply(first_alignment, `[[`, "labels")
      )
    )
    standard_initial <- integer(n)
    standard_initial[overlap_nodes] <- apply(
      aligned_first[overlap_nodes, , drop = FALSE],
      1L,
      modal
    )
    standard_initial[nonoverlap_nodes] <- rowSums(
      aligned_first[nonoverlap_nodes, , drop = FALSE]
    )
    first_membership <- membership_from_labels(aligned_first)
  })[["elapsed"]]
  if (verbose) {
    message("Constructed the first-repetition standard labeling.")
  }
  
  extra_alignment <- list()
  aligned_extra_by_repetition <- list()
  extra_partition_labels <- list()
  extra_memberships <- list()
  repetition_matching_ncores <- 0L
  repetition_alignment_time <- 0
  if (extra_nrep > 0L) {
    extra_task_ids <- seq.int(
      num_subnetworks + 1L,
      total_tasks
    )
    repetition_matching_ncores <- min(length(extra_task_ids), ncores)
    if (verbose) {
      message(
        "Aligning ", length(extra_task_ids),
        " extra subgraphs with ", repetition_matching_ncores,
        " worker(s)."
      )
    }
    repetition_alignment_time <- system.time({
      extra_alignment <- uni_mclapply(
        extra_task_ids,
        function(column_id) {
          observed <- raw_subnetwork_labels[, column_id] != 0L
          match <- match_labels(
            match_this = raw_subnetwork_labels[observed, column_id],
            standard = standard_initial[observed],
            return_mapping = TRUE
          )
          aligned <- integer(n)
          aligned[observed] <- match$matched_labels
          list(labels = aligned, match = match)
        },
        ncores = repetition_matching_ncores,
        force_windows = force_windows,
        stop_on_error = TRUE
      )
      aligned_extra <- do.call(
        cbind,
        lapply(extra_alignment, `[[`, "labels")
      )
      aligned_extra_by_repetition <- lapply(
        seq_len(extra_nrep),
        function(repetition) {
          first <- (repetition - 1L) * num_subnetworks + 1L
          last <- repetition * num_subnetworks
          aligned_extra[, first:last, drop = FALSE]
        }
      )
      extra_partition_labels <- lapply(
        aligned_extra_by_repetition,
        rowSums
      )
      extra_memberships <- lapply(
        aligned_extra_by_repetition,
        membership_from_labels,
        fallback = standard_initial
      )
    })[["elapsed"]]
  }
  
  aggregation_time <- system.time({
    labels <- standard_initial
    if (extra_nrep > 0L) {
      extra_label_matrix <- do.call(cbind, extra_partition_labels)
      labels[nonoverlap_nodes] <- apply(
        cbind(
          standard_initial[nonoverlap_nodes],
          extra_label_matrix[nonoverlap_nodes, , drop = FALSE]
        ),
        1L,
        modal
      )
    }
    repetition_labels <- c(
      list(standard_initial),
      lapply(extra_partition_labels, function(x) {
        x[overlap_nodes] <- standard_initial[overlap_nodes]
        x
      })
    )
    membership_matrices <- c(list(first_membership), extra_memberships)
  })[["elapsed"]]
  } else {
    align_repetition <- function(repetition) {
      column_ids <- (repetition - 1L) * num_subnetworks +
        seq_len(num_subnetworks)
      overlap <- splitter$overlap_nodes[[repetition]]
      nonoverlap <- splitter$nonoverlap_nodes[[repetition]]
      repetition_raw <- raw_subnetwork_labels[, column_ids, drop = FALSE]
      intra_alignment <- lapply(2:num_subnetworks, function(column_id) {
        match <- match_labels(
          match_this = repetition_raw[overlap, column_id],
          standard = repetition_raw[overlap, 1L],
          return_mapping = TRUE
        )
        aligned <- integer(n)
        observed <- repetition_raw[, column_id] != 0L
        aligned[observed] <- unname(
          match$mapping[repetition_raw[observed, column_id]]
        )
        list(labels = aligned, match = match)
      })
      aligned <- do.call(
        cbind,
        c(
          list(repetition_raw[, 1L]),
          lapply(intra_alignment, `[[`, "labels")
        )
      )
      repetition_labels <- integer(n)
      repetition_labels[overlap] <- apply(
        aligned[overlap, , drop = FALSE],
        1L,
        modal
      )
      repetition_labels[nonoverlap] <- rowSums(
        aligned[nonoverlap, , drop = FALSE]
      )
      list(
        labels = repetition_labels,
        aligned = aligned,
        intra_alignment = intra_alignment
      )
    }

    first_alignment_time <- system.time({
      first_fit <- align_repetition(1L)
      standard_initial <- first_fit$labels
      aligned_first <- first_fit$aligned
      first_alignment <- first_fit$intra_alignment
      first_membership <- membership_from_labels(aligned_first)
    })[["elapsed"]]
    if (verbose) {
      message("Constructed the zeroth-repetition standard labeling.")
    }

    extra_alignment <- list()
    aligned_extra_by_repetition <- list()
    extra_partition_labels <- list()
    extra_memberships <- list()
    repetition_matching_ncores <- 0L
    repetition_alignment_time <- 0
    if (extra_nrep > 0L) {
      repetition_ids <- seq.int(2L, extra_nrep + 1L)
      repetition_matching_ncores <- min(length(repetition_ids), ncores)
      if (verbose) {
        message(
          "Aligning ", length(repetition_ids),
          " complete extra-repetition labelings with ",
          repetition_matching_ncores, " worker(s)."
        )
      }
      repetition_alignment_time <- system.time({
        extra_fits <- uni_mclapply(
          repetition_ids,
          function(repetition) {
            fit <- align_repetition(repetition)
            cross_match <- match_labels(
              match_this = fit$labels,
              standard = standard_initial,
              return_mapping = TRUE
            )
            observed <- fit$aligned != 0L
            aligned_to_standard <- fit$aligned
            aligned_to_standard[observed] <- unname(
              cross_match$mapping[fit$aligned[observed]]
            )
            list(
              labels = as.integer(cross_match$matched_labels),
              aligned = aligned_to_standard,
              intra_alignment = fit$intra_alignment,
              cross_repetition_match = cross_match
            )
          },
          ncores = repetition_matching_ncores,
          force_windows = force_windows,
          stop_on_error = TRUE
        )
        extra_partition_labels <- lapply(extra_fits, `[[`, "labels")
        aligned_extra_by_repetition <- lapply(extra_fits, `[[`, "aligned")
        extra_memberships <- lapply(
          aligned_extra_by_repetition,
          membership_from_labels
        )
        extra_alignment <- lapply(extra_fits, function(fit) {
          list(
            within_repetition = fit$intra_alignment,
            against_zeroth_repetition = fit$cross_repetition_match
          )
        })
      })[["elapsed"]]
    }

    aggregation_time <- system.time({
      repetition_labels <- c(list(standard_initial), extra_partition_labels)
      label_votes <- do.call(cbind, repetition_labels)
      labels <- apply(label_votes, 1L, modal)
      membership_matrices <- c(list(first_membership), extra_memberships)
    })[["elapsed"]]
  }
  total_time <- proc.time()[["elapsed"]] - total_start
  
  result <- list(
    labels = as.integer(labels),
    membership_matrices = membership_matrices,
    repetition_labels = repetition_labels,
    raw_subnetwork_labels = raw_subnetwork_labels,
    aligned_subnetwork_labels = c(
      list(aligned_first),
      aligned_extra_by_repetition
    ),
    alignment = list(
      first_repetition = first_alignment,
      extra_repetitions = extra_alignment
    ),
    splitter = splitter,
    overlap_nodes = splitter$overlap_nodes,
    subnetwork_nodes = splitter$subnetworks,
    parameters = list(
      n = n,
      K = K,
      num_subnetworks = num_subnetworks,
      overlap_size = splitter$overlap_size,
      requested_overlap_size = overlap_size,
      extra_nrep = extra_nrep,
      m = m,
      seed = seed,
      matching_method = matching_method,
      force_windows = force_windows,
      ram_check = ram_check,
      share_overlap = share_overlap,
      spectral_arguments = spectral_arguments,
      resolved_spectral_arguments = resolved_spectral_arguments
    ),
    ncores = list(
      requested = ncores,
      clustering = clustering_ncores,
      repetition_matching = repetition_matching_ncores
    ),
    timing = c(
      splitter = splitter_time,
      clustering = clustering_time,
      first_alignment = first_alignment_time,
      repetition_alignment = repetition_alignment_time,
      aggregation = aggregation_time,
      total = total_time
    ),
    call = call
  )
  class(result) <- "sonnet"
  if (verbose) {
    message("SONNET finished in ", format(total_time, digits = 4), " seconds.")
  }
  result
}

# Fit SONNET with one overlap shared by all repetitions.
#' @rdname sonnet_overlap_variants
#' @export
sonnet_shared_overlap <- function(...) {
  result <- .sonnet_fit(..., share_overlap = TRUE)
  result$call <- match.call()
  result
}

# Fit SONNET with an independently sampled overlap in every repetition.
#' @rdname sonnet_overlap_variants
#' @export
sonnet_independent_overlap <- function(...) {
  result <- .sonnet_fit(..., share_overlap = FALSE)
  result$call <- match.call()
  result
}

# Fit either SONNET overlap variant.
#' @rdname sonnet
#' @export
sonnet <- function(
    A,
    K,
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
    ...) {
  result <- .sonnet_fit(
    A = A,
    K = K,
    num_subnetworks = num_subnetworks,
    overlap_size = overlap_size,
    extra_nrep = extra_nrep,
    ncores = ncores,
    seed = seed,
    matching_method = matching_method,
    confirm_large = confirm_large,
    verbose = verbose,
    force_windows = force_windows,
    ram_check = ram_check,
    share_overlap = share_overlap,
    parameter_select_options = parameter_select_options,
    ...
  )
  result$call <- match.call()
  result
}

# Print a compact SONNET fit overview.
#' @rdname sonnet
#' @export
print.sonnet <- function(x, ...) {
  cat("SONNET fit\n")
  cat("  Nodes:", x$parameters$n, "\n")
  cat("  Communities:", x$parameters$K, "\n")
  cat("  Subnetworks per repetition:", x$parameters$num_subnetworks, "\n")
  cat("  Extra repetitions:", x$parameters$extra_nrep, "\n")
  cat(
    "  Overlap across repetitions:",
    if (isTRUE(x$parameters$share_overlap)) "shared" else "independent",
    "\n"
  )
  cat(
    "  Overlap size:",
    x$parameters$overlap_size,
    if (x$parameters$overlap_size != x$parameters$requested_overlap_size) {
      paste0(" (requested ", x$parameters$requested_overlap_size, ")")
    } else {
      ""
    },
    "\n"
  )
  cat("  Matching:", x$parameters$matching_method, "\n")
  cat(
    "  Spectral clustering:",
    x$parameters$resolved_spectral_arguments$spectral_method,
    "/",
    x$parameters$resolved_spectral_arguments$spectral_engine,
    "/",
    x$parameters$resolved_spectral_arguments$cluster_engine,
    "\n"
  )
  cat("  Final community sizes:\n")
  print(table(
    factor(x$labels, levels = seq_len(x$parameters$K))
  ))
  invisible(x)
}

# Summarize a SONNET fit.
#' @rdname sonnet
#' @export
summary.sonnet <- function(object, ...) {
  maximum_membership <- vapply(
    object$membership_matrices,
    function(membership) mean(apply(membership, 1L, max)),
    numeric(1)
  )
  result <- list(
    call = object$call,
    n = object$parameters$n,
    K = object$parameters$K,
    community_sizes = table(factor(
      object$labels,
      levels = seq_len(object$parameters$K)
    )),
    num_subnetworks = object$parameters$num_subnetworks,
    repetitions = object$parameters$extra_nrep + 1L,
    share_overlap = object$parameters$share_overlap,
    overlap_size = object$parameters$overlap_size,
    requested_overlap_size = object$parameters$requested_overlap_size,
    remainder_count = object$splitter$remainder_count,
    matching_method = object$parameters$matching_method,
    spectral_method = object$parameters$resolved_spectral_arguments$spectral_method,
    spectral_engine = object$parameters$resolved_spectral_arguments$spectral_engine,
    cluster_engine = object$parameters$resolved_spectral_arguments$cluster_engine,
    mean_maximum_membership = maximum_membership,
    ncores = object$ncores,
    timing = object$timing
  )
  class(result) <- "summary.sonnet"
  result
}

# Print a SONNET summary.
#' @rdname sonnet
#' @export
print.summary.sonnet <- function(x, ...) {
  cat("Summary of SONNET fit\n")
  cat("  Nodes:", x$n, "\n")
  cat("  Communities:", x$K, "\n")
  cat("  Repetitions:", x$repetitions, "\n")
  cat("  Subnetworks per repetition:", x$num_subnetworks, "\n")
  cat(
    "  Overlap across repetitions:",
    if (isTRUE(x$share_overlap)) "shared" else "independent",
    "\n"
  )
  cat("  Effective overlap:", x$overlap_size, "\n")
  cat("  Matching:", x$matching_method, "\n")
  cat(
    "  Spectral clustering:",
    x$spectral_method,
    "/",
    x$spectral_engine,
    "/",
    x$cluster_engine,
    "\n"
  )
  cat("  Community sizes:\n")
  print(x$community_sizes)
  cat("  Mean maximum membership by repetition:\n")
  print(x$mean_maximum_membership)
  cat("  Timing (seconds):\n")
  print(x$timing)
  invisible(x)
}


# system.time(
#   net <- generate_dcbm(n = 10000, K = 20, ncores = 5, average_degree = 500)
# )
# 
# g_true <- attributes(net)$generator_parameters$g_true
# 
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
#     cluster_engine = "kmeans"
#     # cluster_options = list("metric" = "manhattan")
#   )})
# )
# 
# pram$Peak_RAM_Used_MiB / 1024^2
# label_match_greedy(sc.out$g_hat, g_true)$mismatch_rate
# 
# system.time(
#   pram2 <- peakRAM::peakRAM({son.out <- sonnet(
#     A = net,
#     K = 20,
#     num_subnetworks = 5,
#     overlap_size = 1000,
#     extra_nrep = 0,
#     ncores = 5,
#     seed = 52,
#     matching_method = "greedy",
#     verbose = TRUE,
#     share_overlap = TRUE,
#     laplacian = TRUE,
#     regularize_tau = 0,
#     row_normalize = TRUE,
#     spectral_method = "eigen",
#     cluster_engine = "clara",
#     cluster_options = list("metric" = "manhattan", sampsize = 200, samples = 50),
#     force_windows = FALSE
#   )})
# )
# 
# pram2$Peak_RAM_Used_MiB / 1024^2
# label_match_greedy(son.out$labels, g_true)$mismatch_rate
# 
# 
# system.time(
#   ase.out <- ase(A = net, d = 20)
# )
# 
# ase.out
# 
# 
# 
# 
# 
