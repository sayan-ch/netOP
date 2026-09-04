# Basic cross-platform helpers
#
# This file follows NAMING_CONVENTION.md. The original proposed name
# `uni.mclapply` is therefore exposed as `uni_mclapply`, and `mc.cores` is
# exposed as `ncores`. The dotted `mc.cores`, `.Platform$OS.type`, and
# `future.seed` names below are external R/package interfaces and cannot be
# renamed here.

# Apply a function over X using a platform-appropriate parallel backend.
#
# On Unix-like systems, this helper first attempts parallel::mclapply. On
# Windows, or when force_windows is TRUE, it uses
# future.apply::future_lapply. By default it creates a temporary multisession
# plan and restores the caller's plan; high-level multi-stage algorithms may
# set manage_future_plan = FALSE after creating one reusable pool. Any backend
# error causes a safe fallback to base::lapply unless strict mode is requested.
#
# Arguments:
# - X: A vector, list, or other object accepted by lapply-style functions.
# - FUN: The function to apply to each element of X.
# - ...: Additional arguments forwarded to FUN.
# - ncores: Positive integer number of worker processes. By default, use half
#   of the detected logical cores, rounded down, with a minimum of one.
# - force_windows: Logical flag for exercising the Windows/future backend on
#   another operating system.
# - stop_on_error: If TRUE, stop with failed task indices and worker messages
#   instead of falling back or returning worker-level error objects.
# - manage_future_plan: Create and restore a local multisession plan. Set FALSE
#   only when the caller owns an already configured future plan.
# - future_packages: Optional packages loaded before future workers evaluate
#   tasks; this is important for S4 method registration (for example Matrix).
# - suppress_worker_renv_sync_check: Temporarily disable renv's synchronized
#   project check while multisession workers start, avoiding one identical
#   project-status message per worker. The parent environment is restored.
#
# Returns:
# - A list with one result per element of X.
#' @rdname uni_mclapply
#' @export
uni_mclapply <- function(
    X,
    FUN,
    ...,
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    force_windows = FALSE,
    stop_on_error = FALSE,
    manage_future_plan = TRUE,
    future_packages = NULL,
    suppress_worker_renv_sync_check = TRUE) {
  if (!is.function(FUN)) {
    stop("FUN must be a function.", call. = FALSE)
  }

  if (length(ncores) != 1L ||
      !is.numeric(ncores) ||
      is.na(ncores) ||
      !is.finite(ncores) ||
      ncores < 1 ||
      ncores != floor(ncores)) {
    stop("ncores must be one positive integer.", call. = FALSE)
  }
  ncores <- as.integer(ncores)

  if (length(force_windows) != 1L ||
      !is.logical(force_windows) ||
      is.na(force_windows)) {
    stop("force_windows must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(stop_on_error) != 1L ||
      !is.logical(stop_on_error) ||
      is.na(stop_on_error)) {
    stop("stop_on_error must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(manage_future_plan) != 1L ||
      !is.logical(manage_future_plan) || is.na(manage_future_plan)) {
    stop("manage_future_plan must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(future_packages) &&
      (!is.character(future_packages) || anyNA(future_packages) ||
       any(!nzchar(future_packages)))) {
    stop("future_packages must be NULL or package names.", call. = FALSE)
  }
  if (length(suppress_worker_renv_sync_check) != 1L ||
      !is.logical(suppress_worker_renv_sync_check) ||
      is.na(suppress_worker_renv_sync_check)) {
    stop("suppress_worker_renv_sync_check must be TRUE or FALSE.",
         call. = FALSE)
  }

  loop_length <- length(X)
  if (loop_length == 0L) {
    return(vector("list", 0L))
  }
  ncores <- min(ncores, loop_length)

  worker_function <- if (!stop_on_error) {
    FUN
  } else {
    function(value, ...) {
      tryCatch(
        structure(
          list(success = TRUE, value = FUN(value, ...)),
          class = "uni_mclapply_worker_result"
        ),
        error = function(e) {
          structure(
            list(
              success = FALSE,
              message = conditionMessage(e),
              call = conditionCall(e),
              classes = class(e)
            ),
            class = "uni_mclapply_worker_result"
          )
        }
      )
    }
  }

  unwrap_results <- function(results) {
    if (!stop_on_error) {
      return(results)
    }
    backend_failures <- which(vapply(
      results,
      inherits,
      logical(1),
      what = "try-error"
    ))
    if (length(backend_failures) > 0L) {
      stop(
        "Parallel worker failure at task(s) ",
        paste(backend_failures, collapse = ", "),
        ": ",
        paste(vapply(
          results[backend_failures],
          as.character,
          character(1)
        ), collapse = " | "),
        call. = FALSE
      )
    }
    malformed <- which(!vapply(
      results,
      inherits,
      logical(1),
      what = "uni_mclapply_worker_result"
    ))
    if (length(malformed) > 0L) {
      stop(
        "Parallel backend returned malformed results for task(s) ",
        paste(malformed, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    failed <- which(!vapply(results, `[[`, logical(1), "success"))
    if (length(failed) > 0L) {
      messages <- vapply(results[failed], `[[`, character(1), "message")
      stop(
        "Worker error in task(s) ",
        paste(failed, collapse = ", "),
        ": ",
        paste(messages, collapse = " | "),
        call. = FALSE
      )
    }
    lapply(results, `[[`, "value")
  }

  os_type <- if (force_windows) "windows" else .Platform$OS.type

  if (ncores == 1L) {
    return(unwrap_results(lapply(X, worker_function, ...)))
  }

  if (os_type == "unix") {
    results <- tryCatch(
      parallel::mclapply(X, worker_function, ..., mc.cores = ncores),
      error = function(e) {
        if (stop_on_error) {
          stop(
            "Parallel backend failed: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
        message("mclapply failed; falling back to lapply: ", e$message)
        lapply(X, FUN, ...)
      }
    )
    return(unwrap_results(results))
  }

  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    if (stop_on_error) {
      stop(
        paste0(
          "future and future.apply are required for strict parallel ",
          "Windows execution."
        ),
        call. = FALSE
      )
    }
    message(
      "future and future.apply are required for parallel Windows execution; ",
      "falling back to lapply."
    )
    return(lapply(X, FUN, ...))
  }

  if (suppress_worker_renv_sync_check) {
    previous_renv_sync_check <- Sys.getenv(
      "RENV_CONFIG_SYNCHRONIZED_CHECK",
      unset = NA_character_
    )
    Sys.setenv(RENV_CONFIG_SYNCHRONIZED_CHECK = "FALSE")
    on.exit({
      if (is.na(previous_renv_sync_check)) {
        Sys.unsetenv("RENV_CONFIG_SYNCHRONIZED_CHECK")
      } else {
        Sys.setenv(
          RENV_CONFIG_SYNCHRONIZED_CHECK = previous_renv_sync_check
        )
      }
    }, add = TRUE)
  }

  results <- tryCatch({
    if (manage_future_plan) {
      previous_plan <- future::plan()
      on.exit(
        try(future::plan(previous_plan), silent = TRUE),
        add = TRUE
      )
      future::plan(future::multisession, workers = ncores)
    }
    future.apply::future_lapply(
      X,
      worker_function,
      ...,
      future.seed = TRUE,
      future.packages = future_packages
    )
  }, error = function(e) {
    if (stop_on_error) {
      stop(
        "Parallel backend failed: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
    message("future_lapply failed; falling back to lapply: ", e$message)
    lapply(X, FUN, ...)
  })
  unwrap_results(results)
}

# Return currently available system RAM in bytes.
#
# ps is deliberately optional. Different ps versions and operating systems use
# slightly different field names, so check the known alternatives.
#' @rdname ram_reporting
#' @export
available_ram <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) {
    return(NA_real_)
  }
  info <- tryCatch(ps::ps_system_memory(), error = function(e) NULL)
  if (is.null(info)) {
    return(NA_real_)
  }
  for (name in c("available", "avail", "free")) {
    if (name %in% names(info)) {
      value <- suppressWarnings(as.numeric(info[[name]]))
      if (length(value) == 1L && is.finite(value) && value >= 0) {
        return(value)
      }
    }
  }
  NA_real_
}

# Format a non-negative byte count using binary units.
#' @rdname ram_reporting
#' @export
format_bytes <- function(x) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x)) {
    return("unknown")
  }
  if (x < 0) {
    stop("x must be non-negative.", call. = FALSE)
  }
  if (is.infinite(x)) {
    return("Inf")
  }
  units <- c("B", "KiB", "MiB", "GiB", "TiB", "PiB")
  unit_id <- if (x == 0) {
    0L
  } else {
    min(floor(log(x, base = 1024)), length(units) - 1L)
  }
  sprintf("%.2f %s", x / 1024^unit_id, units[unit_id + 1L])
}

# Report estimated versus available RAM and warn above a conservative limit.
#' @rdname ram_reporting
#' @export
warn_if_insufficient_ram <- function(
    estimated_bytes,
    max_fraction = 0.60,
    reserve_bytes = 2 * 1024^3,
    operation = "operation") {
  if (length(estimated_bytes) != 1L ||
      !is.numeric(estimated_bytes) ||
      is.na(estimated_bytes) ||
      !is.finite(estimated_bytes) ||
      estimated_bytes < 0) {
    stop("estimated_bytes must be one finite non-negative number.",
         call. = FALSE)
  }
  if (length(max_fraction) != 1L ||
      !is.numeric(max_fraction) ||
      is.na(max_fraction) ||
      !is.finite(max_fraction) ||
      max_fraction <= 0 ||
      max_fraction > 1) {
    stop("max_fraction must be in (0, 1].", call. = FALSE)
  }
  if (length(reserve_bytes) != 1L ||
      !is.numeric(reserve_bytes) ||
      is.na(reserve_bytes) ||
      !is.finite(reserve_bytes) ||
      reserve_bytes < 0) {
    stop("reserve_bytes must be one finite non-negative number.",
         call. = FALSE)
  }
  if (length(operation) != 1L || !is.character(operation) ||
      is.na(operation) || !nzchar(operation)) {
    stop("operation must be one non-empty string.", call. = FALSE)
  }

  available_bytes <- available_ram()
  estimate_text <- format_bytes(estimated_bytes)
  if (is.na(available_bytes)) {
    warning(
      operation, ": expected RAM use is ", estimate_text,
      " out of unknown available RAM because package 'ps' is unavailable ",
      "or the platform did not report available memory.",
      call. = FALSE
    )
    return(invisible(NA))
  }

  usable_bytes <- max(
    0,
    min(max_fraction * available_bytes, available_bytes - reserve_bytes)
  )
  report <- paste0(
    operation, ": expected RAM use is ", estimate_text,
    " out of ", format_bytes(available_bytes), " available RAM."
  )
  if (estimated_bytes > usable_bytes) {
    warning(
      report,
      " Conservative usable-memory threshold: ",
      format_bytes(usable_bytes), ".",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  message(report)
  invisible(TRUE)
}

# Estimate RAM for an RSpectra Lanczos/Arnoldi decomposition.
#' @rdname ram_estimators
#' @export
estimate_rspectra_ram <- function(
    n,
    K,
    ncv = min(n, max(2 * K + 1, 20)),
    safety_factor = 2) {
  values <- c(n = n, K = K, ncv = ncv)
  if (any(!is.finite(values)) || any(values != floor(values)) ||
      n < 1 || K < 1 || K >= n || ncv <= K || ncv > n) {
    stop("Require integer n >= ncv > K >= 1.", call. = FALSE)
  }
  if (length(safety_factor) != 1L || !is.numeric(safety_factor) ||
      is.na(safety_factor) || !is.finite(safety_factor) ||
      safety_factor <= 0) {
    stop("safety_factor must be one positive finite number.", call. = FALSE)
  }
  doubles <- as.double(n) * ncv + as.double(n) * K +
    6 * as.double(n) + as.double(ncv) * (ncv + 8)
  safety_factor * 8 * doubles
}

# Estimate RAM for a dense symmetric eigendecomposition.
#' @rdname ram_estimators
#' @export
estimate_dense_eigen_ram <- function(
    n,
    input_already_counted = TRUE,
    safety_factor = 1.5) {
  if (length(n) != 1L || !is.numeric(n) || is.na(n) ||
      !is.finite(n) || n < 1 || n != floor(n)) {
    stop("n must be one positive integer.", call. = FALSE)
  }
  if (length(input_already_counted) != 1L ||
      !is.logical(input_already_counted) || is.na(input_already_counted)) {
    stop("input_already_counted must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(safety_factor) != 1L || !is.numeric(safety_factor) ||
      is.na(safety_factor) || !is.finite(safety_factor) ||
      safety_factor <= 0) {
    stop("safety_factor must be one positive finite number.", call. = FALSE)
  }
  dense_matrices <- if (input_already_counted) 2 else 3
  safety_factor * 8 * as.double(n)^2 * dense_matrices
}

# Estimate RAM for a partial rectangular singular decomposition.
#' @rdname ram_estimators
#' @export
estimate_partial_svd_ram <- function(
    n,
    p,
    K,
    ncv = min(min(n, p), max(2 * K + 1, 20)),
    safety_factor = 2) {
  values <- c(n = n, p = p, K = K, ncv = ncv)
  min_dim <- min(n, p)
  if (any(!is.finite(values)) || any(values != floor(values)) ||
      n < 1 || p < 1 || K < 1 || K >= min_dim ||
      ncv <= K || ncv > min_dim) {
    stop("Require positive integer dimensions and min(n, p) >= ncv > K.",
         call. = FALSE)
  }
  if (length(safety_factor) != 1L || !is.numeric(safety_factor) ||
      is.na(safety_factor) || !is.finite(safety_factor) ||
      safety_factor <= 0) {
    stop("safety_factor must be one positive finite number.", call. = FALSE)
  }
  doubles <- as.double(n + p) * ncv + as.double(n + p) * K +
    6 * as.double(n + p) + as.double(ncv) * (ncv + 8)
  safety_factor * 8 * doubles
}

# Estimate RAM for a dense rectangular singular decomposition.
#' @rdname ram_estimators
#' @export
estimate_dense_svd_ram <- function(
    n,
    p,
    nu = min(n, p),
    nv = min(n, p),
    input_already_counted = TRUE,
    safety_factor = 1.5) {
  values <- c(n = n, p = p, nu = nu, nv = nv)
  if (any(!is.finite(values)) || any(values != floor(values)) ||
      n < 1 || p < 1 || nu < 0 || nv < 0 || nu > n || nv > p) {
    stop("n and p must be positive integer dimensions; nu and nv must fit.",
         call. = FALSE)
  }
  if (length(input_already_counted) != 1L ||
      !is.logical(input_already_counted) || is.na(input_already_counted)) {
    stop("input_already_counted must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(safety_factor) != 1L || !is.numeric(safety_factor) ||
      is.na(safety_factor) || !is.finite(safety_factor) ||
      safety_factor <= 0) {
    stop("safety_factor must be one positive finite number.", call. = FALSE)
  }
  input_doubles <- as.double(n) * p *
    if (input_already_counted) 1 else 2
  output_and_work <- as.double(n) * nu + as.double(p) * nv +
    as.double(min(n, p))^2
  safety_factor * 8 * (input_doubles + output_and_work)
}

# Estimate the output and packing workspace for one dense matrix product.
#' @rdname ram_estimators
#' @export
estimate_matrix_product_ram <- function(
    nrow_left,
    shared_dimension,
    ncol_right,
    safety_factor = 1.25) {
  dimensions <- c(nrow_left, shared_dimension, ncol_right)
  if (any(!is.numeric(dimensions)) || anyNA(dimensions) ||
      any(!is.finite(dimensions)) || any(dimensions < 1) ||
      any(dimensions != floor(dimensions))) {
    stop("Matrix-product dimensions must be positive integers.",
         call. = FALSE)
  }
  if (length(safety_factor) != 1L || !is.numeric(safety_factor) ||
      is.na(safety_factor) || !is.finite(safety_factor) ||
      safety_factor <= 0) {
    stop("safety_factor must be one positive finite number.", call. = FALSE)
  }
  output_doubles <- as.double(nrow_left) * ncol_right
  packing_doubles <- min(
    as.double(nrow_left) * shared_dimension,
    as.double(shared_dimension) * ncol_right
  )
  safety_factor * 8 * (output_doubles + packing_doubles)
}

# Estimate the decomposition selected by the audited fallback rules.
#' @rdname ram_estimators
#' @export
estimate_spectral_decomp_ram <- function(
    n,
    p = n,
    K,
    method = c("eigen", "svd"),
    engine = c("rspectra", "irlba", "base"),
    force_engine = FALSE,
    safe_d_multiplier = 1,
    nu = K,
    nv = if (method[1L] == "eigen") K else 0L,
    dense_input = FALSE) {
  method <- match.arg(method)
  engine <- match.arg(tolower(engine), c("rspectra", "irlba", "base"))
  counts <- c(n = n, p = p, K = K, nu = nu, nv = nv)
  if (any(!is.numeric(counts)) || anyNA(counts) ||
      any(!is.finite(counts)) || any(counts != floor(counts)) ||
      n < 1 || p < 1 || K < 1 || K > min(n, p) || nu < 0 || nv < 0) {
    stop("Invalid decomposition dimensions.", call. = FALSE)
  }
  logical_values <- c(force_engine = force_engine, dense_input = dense_input)
  if (anyNA(logical_values) || !all(logical_values %in% c(TRUE, FALSE))) {
    stop("force_engine and dense_input must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (length(safe_d_multiplier) != 1L ||
      !is.numeric(safe_d_multiplier) || is.na(safe_d_multiplier) ||
      !is.finite(safe_d_multiplier) || safe_d_multiplier < 1) {
    stop("safe_d_multiplier must be at least one.", call. = FALSE)
  }

  min_dim <- min(n, p)
  d_compute <- min(ceiling(safe_d_multiplier * K), min_dim)
  resolved_engine <- engine
  if (method == "eigen" && resolved_engine == "irlba" && !force_engine) {
    resolved_engine <- if (requireNamespace("RSpectra", quietly = TRUE)) {
      "rspectra"
    } else {
      "base"
    }
  }
  package_name <- if (resolved_engine == "rspectra") {
    "RSpectra"
  } else {
    resolved_engine
  }
  if (min_dim <= 200L && resolved_engine != "base" && !force_engine) {
    resolved_engine <- "base"
  } else if (resolved_engine != "base" && !force_engine &&
             !requireNamespace(package_name, quietly = TRUE)) {
    resolved_engine <- "base"
  }
  if (resolved_engine != "base" && d_compute >= min_dim) {
    if (K >= min_dim) {
      resolved_engine <- "base"
    } else {
      d_compute <- min_dim - 1L
    }
  }

  if (resolved_engine == "base") {
    estimated_bytes <- if (method == "eigen") {
      estimate_dense_eigen_ram(
        n = n,
        input_already_counted = !dense_input
      )
    } else {
      estimate_dense_svd_ram(
        n = n,
        p = p,
        nu = min(nu, n, K),
        nv = min(nv, p, K),
        input_already_counted = !dense_input
      )
    }
  } else {
    ncv <- min(min_dim, max(2 * d_compute + 1, 20))
    estimated_bytes <- if (method == "eigen") {
      estimate_rspectra_ram(n = n, K = d_compute, ncv = ncv)
    } else {
      estimate_partial_svd_ram(
        n = n,
        p = p,
        K = d_compute,
        ncv = ncv
      )
    }
    if (dense_input) {
      estimated_bytes <- estimated_bytes + 8 * as.double(n) * p
    }
  }
  list(
    estimated_bytes = estimated_bytes,
    resolved_engine = resolved_engine,
    computed_dimension = d_compute,
    dense_input = dense_input
  )
}

# Report an algorithm-level RAM estimate, optionally across parallel workers.
#' @rdname ram_reporting
#' @export
report_ram_preflight <- function(
    estimated_bytes,
    operation,
    operation_count = 1L,
    ncores = 1L,
    detail = NULL) {
  if (length(estimated_bytes) != 1L || !is.numeric(estimated_bytes) ||
      is.na(estimated_bytes) || !is.finite(estimated_bytes) ||
      estimated_bytes < 0) {
    stop("estimated_bytes must be one finite non-negative number.",
         call. = FALSE)
  }
  multipliers <- c(operation_count = operation_count, ncores = ncores)
  if (any(!is.numeric(multipliers)) || anyNA(multipliers) ||
      any(!is.finite(multipliers)) || any(multipliers < 1) ||
      any(multipliers != floor(multipliers))) {
    stop("operation_count and ncores must be positive integers.",
         call. = FALSE)
  }
  if (length(operation) != 1L || !is.character(operation) ||
      is.na(operation) || !nzchar(operation)) {
    stop("operation must be one non-empty string.", call. = FALSE)
  }
  total_bytes <- estimated_bytes * operation_count * ncores
  available_bytes <- available_ram()
  report <- paste0(
    operation,
    if (!is.null(detail)) paste0(" (", detail, ")") else "",
    ": expected RAM use ", format_bytes(estimated_bytes),
    " x ", operation_count, " operation(s)",
    " x ", ncores, " core(s) = ", format_bytes(total_bytes),
    " out of ", format_bytes(available_bytes), " available RAM."
  )
  if (is.na(available_bytes) || total_bytes > available_bytes) {
    warning(report, call. = FALSE)
  } else {
    message(report)
  }
  invisible(list(
    expected_per_operation = estimated_bytes,
    operation_count = as.integer(operation_count),
    ncores = as.integer(ncores),
    expected_total = total_bytes,
    available = available_bytes
  ))
}

# Report an additive high-level RAM formula with explicit operation factors.
#' @rdname ram_reporting
#' @export
report_ram_formula <- function(terms, operation, detail = NULL) {
  if (!is.list(terms) || length(terms) < 1L) {
    stop("terms must be a non-empty list.", call. = FALSE)
  }
  if (length(operation) != 1L || !is.character(operation) ||
      is.na(operation) || !nzchar(operation)) {
    stop("operation must be one non-empty string.", call. = FALSE)
  }
  normalized_terms <- lapply(seq_along(terms), function(term_id) {
    term <- terms[[term_id]]
    if (!is.list(term)) {
      stop("Every RAM formula term must be a list.", call. = FALSE)
    }
    bytes <- term$estimated_bytes
    sequential_count <- if (is.null(term$sequential_count)) {
      1L
    } else {
      term$sequential_count
    }
    parallel_count <- if (is.null(term$parallel_count)) {
      1L
    } else {
      term$parallel_count
    }
    label <- if (is.null(term$label)) {
      paste0("term ", term_id)
    } else {
      term$label
    }
    counts <- c(sequential_count, parallel_count)
    if (length(bytes) != 1L || !is.numeric(bytes) || is.na(bytes) ||
        !is.finite(bytes) || bytes < 0 ||
        any(!is.numeric(counts)) || anyNA(counts) ||
        any(!is.finite(counts)) || any(counts < 1) ||
        any(counts != floor(counts)) ||
        length(label) != 1L || !is.character(label) ||
        is.na(label) || !nzchar(label)) {
      stop("Invalid RAM formula term.", call. = FALSE)
    }
    list(
      estimated_bytes = bytes,
      sequential_count = sequential_count,
      parallel_count = parallel_count,
      label = label,
      subtotal = bytes * sequential_count * parallel_count
    )
  })
  total_bytes <- sum(vapply(
    normalized_terms,
    `[[`,
    numeric(1),
    "subtotal"
  ))
  formula <- paste(vapply(normalized_terms, function(term) {
    paste0(
      "(", format_bytes(term$estimated_bytes),
      " [", term$label, "] x ", term$sequential_count,
      " sequential operation(s)) x ", term$parallel_count,
      " parallel operation(s)"
    )
  }, character(1)), collapse = " + ")
  available_bytes <- available_ram()
  report <- paste0(
    operation,
    if (!is.null(detail)) paste0(" (", detail, ")") else "",
    ": ", formula,
    " = ", format_bytes(total_bytes),
    " out of ", format_bytes(available_bytes), " available RAM."
  )
  if (is.na(available_bytes) || total_bytes > available_bytes) {
    warning(report, call. = FALSE)
  } else {
    message(report)
  }
  invisible(list(
    terms = normalized_terms,
    expected_total = total_bytes,
    available = available_bytes
  ))
}

# Run an expression or function through peakRAM::peakRAM.
#
# peakRAM returns only a metrics data frame and discards the evaluated result.
# This wrapper passes it a closure that captures the result, allowing existing
# callers to receive both the value and the peakRAM measurements without
# evaluating the process twice.
#
# Usage:
#
# measure_peak_ram({
#   mean(1:5000)
# })
#
# measure_peak_ram(mean, x = 1:5000)
#
# Arguments:
# - process: Unevaluated R expression or a function to execute and measure.
# - ...: Additional arguments forwarded when process evaluates to a function.
#
# Returns:
# - A list containing the process result and a one-row metrics data frame with
#   elapsed_seconds, total_ram_used_mib, and peak_ram_used_mib.
#' @rdname measure_peak_ram
#' @export
measure_peak_ram <- function(process, ...) {
  if (missing(process)) {
    stop("process must be an expression or function.", call. = FALSE)
  }
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    stop(
      "The peakRAM package is required by measure_peak_ram(). ",
      "Install it with install.packages(\"peakRAM\").",
      call. = FALSE
    )
  }

  process_expression <- substitute(process)
  process_environment <- parent.frame()
  process_args <- list(...)
  process_result <- NULL

  run_process <- function() {
    process_value <- eval(process_expression, envir = process_environment)
    if (is.function(process_value)) {
      return(do.call(process_value, process_args))
    }
    if (length(process_args) > 0L) {
      stop(
        "Additional arguments can only be supplied when process is a function.",
        call. = FALSE
      )
    }
    process_value
  }

  measured_process <- function() {
    process_result <<- run_process()
    invisible(NULL)
  }

  peak_ram_metrics <- peakRAM::peakRAM(measured_process)
  required_columns <- c(
    "Elapsed_Time_sec",
    "Total_RAM_Used_MiB",
    "Peak_RAM_Used_MiB"
  )
  missing_columns <- setdiff(required_columns, names(peak_ram_metrics))
  if (length(missing_columns) > 0L || nrow(peak_ram_metrics) != 1L) {
    stop(
      "Unexpected output from peakRAM::peakRAM().",
      call. = FALSE
    )
  }

  list(
    result = process_result,
    metrics = data.frame(
      elapsed_seconds = peak_ram_metrics$Elapsed_Time_sec[[1L]],
      total_ram_used_mib = peak_ram_metrics$Total_RAM_Used_MiB[[1L]],
      peak_ram_used_mib = peak_ram_metrics$Peak_RAM_Used_MiB[[1L]],
      check.names = FALSE
    )
  )
}

# Run a simulation function repeatedly with optional parallelism and resuming.
#
# one_simulation must accept a named `simulation` argument containing the
# one-based simulation number. Additional inputs supplied through `...` are
# forwarded unchanged to one_simulation.
#
# When results_file is supplied, results are stored as an RDS list with one
# record per simulation. Each record contains the simulation number, success
# flag, elapsed time, optional peak RAM, returned value, and error message. When
# resource measurement is enabled, elapsed_seconds and peak_ram_used_mib are
# also appended to data-frame or list results before the records are saved.
# The replace, resume, and archive actions follow the workflow used by the
# simulation scripts from which this helper was generalized.
#
# Arguments:
# - one_simulation: Function that performs one simulation.
# - nsim: Positive integer number of simulations.
# - ...: Additional arguments forwarded to one_simulation.
# - use_parallel_simulations: Whether to parallelize the outer loop over nsim.
#   The default is FALSE because one_simulation commonly uses parallel methods
#   internally, and nested parallelism can oversubscribe the machine.
# - ncores_outer: Positive integer number of outer-loop workers, used only when
#   use_parallel_simulations is TRUE. A `ncores` argument in `...` remains
#   available to one_simulation and its internal methods.
# - seed: Optional non-negative integer-like base seed. Simulation i uses
#   seed + i - 1, reduced safely to R's supported integer range.
# - seeds: Optional vector of nsim explicit non-negative integer-like seeds.
#   Supply either seed or seeds, not both. The seed used is stored in each
#   simulation record.
# - results_file: Optional path to an RDS file used for durable results.
# - action: How to handle an existing results_file: replace it, resume its
#   incomplete simulations, or archive it and start again.
# - retry_failed: When resuming, rerun records whose success field is FALSE.
# - continue_on_error: Return all records when any simulation fails. If FALSE,
#   save available records first and then stop with a summary error.
# - measure_resources: Measure per-simulation elapsed time and peak additional
#   R-managed RAM, and append those values to data-frame or list results.
# - show_progress: Print a compact status line for every newly run simulation.
#   Sequential runs print immediately after each simulation. Parallel runs with
#   more than one outer worker may omit these messages.
# - force_windows: Forwarded to uni_mclapply for backend testing/selection.
#
# Returns:
# - A list of nsim simulation records, invisibly when no work remains during a
#   resume and normally otherwise.
#' @rdname run_simulations
#' @export
run_simulations <- function(
    one_simulation,
    nsim,
    ...,
    use_parallel_simulations = FALSE,
    ncores_outer = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    seed = NULL,
    seeds = NULL,
    results_file = NULL,
    action = c("replace", "resume", "archive"),
    retry_failed = TRUE,
    continue_on_error = TRUE,
    measure_resources = FALSE,
    show_progress = TRUE,
    force_windows = FALSE) {
  if (!is.function(one_simulation)) {
    stop("one_simulation must be a function.", call. = FALSE)
  }

  if (length(nsim) != 1L ||
      !is.numeric(nsim) ||
      is.na(nsim) ||
      !is.finite(nsim) ||
      nsim < 1 ||
      nsim != floor(nsim)) {
    stop("nsim must be one positive integer.", call. = FALSE)
  }
  nsim <- as.integer(nsim)

  if (length(ncores_outer) != 1L ||
      !is.numeric(ncores_outer) ||
      is.na(ncores_outer) ||
      !is.finite(ncores_outer) ||
      ncores_outer < 1 ||
      ncores_outer != floor(ncores_outer)) {
    stop("ncores_outer must be one positive integer.", call. = FALSE)
  }
  ncores_outer <- as.integer(ncores_outer)

  if (!is.null(seed) && !is.null(seeds)) {
    stop("Supply either seed or seeds, not both.", call. = FALSE)
  }
  if (!is.null(seed) &&
      (length(seed) != 1L ||
       !is.numeric(seed) ||
       is.na(seed) ||
       !is.finite(seed) ||
       seed < 0 ||
       seed != floor(seed))) {
    stop(
      "seed must be NULL or one non-negative finite integer-like value.",
      call. = FALSE
    )
  }
  if (!is.null(seeds) &&
      (length(seeds) != nsim ||
       !is.numeric(seeds) ||
       anyNA(seeds) ||
       any(!is.finite(seeds)) ||
       any(seeds < 0) ||
       any(seeds != floor(seeds)))) {
    stop(
      "seeds must be NULL or nsim non-negative finite integer-like values.",
      call. = FALSE
    )
  }

  simulation_seeds <- rep(NA_integer_, nsim)
  if (!is.null(seed)) {
    simulation_seeds <- as.integer(
      (as.double(seed) + seq_len(nsim) - 1) %% .Machine$integer.max
    )
  } else if (!is.null(seeds)) {
    simulation_seeds <- as.integer(
      as.double(seeds) %% .Machine$integer.max
    )
  }

  logical_inputs <- list(
    use_parallel_simulations = use_parallel_simulations,
    retry_failed = retry_failed,
    continue_on_error = continue_on_error,
    measure_resources = measure_resources,
    show_progress = show_progress,
    force_windows = force_windows
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
      paste0(
        paste(names(logical_inputs)[invalid_logical], collapse = ", "),
        " must be TRUE or FALSE."
      ),
      call. = FALSE
    )
  }

  action <- match.arg(action)
  if (!is.null(results_file) &&
      (length(results_file) != 1L ||
       !is.character(results_file) ||
       is.na(results_file) ||
       !nzchar(results_file))) {
    stop("results_file must be NULL or one non-empty path.", call. = FALSE)
  }

  archive_file <- function(file) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    archived_file <- paste0(file, "_", timestamp)
    suffix <- 1L
    while (file.exists(archived_file)) {
      archived_file <- paste0(file, "_", timestamp, "_", suffix)
      suffix <- suffix + 1L
    }
    if (!file.rename(file, archived_file)) {
      stop("Could not archive results file: ", file, call. = FALSE)
    }
    message("Archived ", file, " to ", archived_file)
    invisible(archived_file)
  }

  results <- vector("list", nsim)
  if (!is.null(results_file) && file.exists(results_file)) {
    if (action == "archive") {
      archive_file(results_file)
    } else if (action == "resume") {
      results <- readRDS(results_file)
      if (!is.list(results)) {
        stop(
          "Cannot resume because results_file does not contain a list.",
          call. = FALSE
        )
      }
      length(results) <- nsim
    }
  }

  if (action != "resume") {
    results <- vector("list", nsim)
  } else if (retry_failed) {
    failed <- vapply(
      results,
      function(record) {
        is.list(record) &&
          !is.null(record$success) &&
          identical(record$success, FALSE)
      },
      logical(1)
    )
    results[failed] <- vector("list", sum(failed))
  }

  simulations <- which(vapply(results, is.null, logical(1)))
  if (length(simulations) == 0L) {
    message("All ", nsim, " simulations are already complete.")
    return(invisible(results))
  }

  run_one_simulation <- function(simulation) {
    simulation_seed <- simulation_seeds[simulation]
    if (!is.na(simulation_seed)) {
      set.seed(simulation_seed)
    }

    execute_simulation <- function() {
      tryCatch(
        {
          value <- one_simulation(simulation = simulation, ...)
          list(success = TRUE, result = value, error = NULL)
        },
        error = function(e) {
          list(success = FALSE, result = NULL, error = conditionMessage(e))
        }
      )
    }

    if (measure_resources) {
      resource_measurement <- measure_peak_ram({
        execute_simulation()
      })
      simulation_result <- resource_measurement$result
      elapsed_seconds <- resource_measurement$metrics$elapsed_seconds
      peak_ram_used_mib <- resource_measurement$metrics$peak_ram_used_mib
    } else {
      start_time <- proc.time()[["elapsed"]]
      simulation_result <- execute_simulation()
      elapsed_seconds <- unname(proc.time()[["elapsed"]] - start_time)
      peak_ram_used_mib <- NA_real_
    }

    if (measure_resources && simulation_result$success) {
      if (is.data.frame(simulation_result$result)) {
        result_rows <- nrow(simulation_result$result)
        simulation_result$result$elapsed_seconds <- rep(
          elapsed_seconds,
          result_rows
        )
        simulation_result$result$peak_ram_used_mib <- rep(
          peak_ram_used_mib,
          result_rows
        )
      } else if (is.list(simulation_result$result)) {
        simulation_result$result$elapsed_seconds <- elapsed_seconds
        simulation_result$result$peak_ram_used_mib <- peak_ram_used_mib
      }
    }

    list(
      simulation = simulation,
      seed = if (is.na(simulation_seed)) NULL else simulation_seed,
      success = simulation_result$success,
      elapsed_seconds = elapsed_seconds,
      peak_ram_used_mib = peak_ram_used_mib,
      result = simulation_result$result,
      error = simulation_result$error
    )
  }

  report_simulation <- function(record) {
    status <- if (record$success) "success" else "failed"
    message(
      sprintf(
        "Simulation %d/%d finished | status=%s | elapsed=%.3f s",
        record$simulation,
        nsim,
        status,
        record$elapsed_seconds
      )
    )
    if (!record$success) {
      message("  error=", record$error)
    }
    invisible(NULL)
  }

  use_parallel_outer_loop <- use_parallel_simulations && ncores_outer > 1L
  if (use_parallel_outer_loop) {
    new_results <- uni_mclapply(
      simulations,
      run_one_simulation,
      ncores = ncores_outer,
      force_windows = force_windows
    )
  } else {
    new_results <- vector("list", length(simulations))
    for (simulation_index in seq_along(simulations)) {
      new_results[[simulation_index]] <- run_one_simulation(
        simulations[simulation_index]
      )
      if (show_progress) {
        report_simulation(new_results[[simulation_index]])
      }
    }
  }
  results[simulations] <- new_results

  if (!is.null(results_file)) {
    results_dir <- dirname(results_file)
    if (!dir.exists(results_dir) &&
        !dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create results directory: ", results_dir,
           call. = FALSE)
    }
    saveRDS(results, results_file)
  }

  failed <- vapply(
    results,
    function(record) {
      is.list(record) &&
        !is.null(record$success) &&
        identical(record$success, FALSE)
    },
    logical(1)
  )
  if (any(failed) && !continue_on_error) {
    stop(
      sum(failed), " of ", nsim, " simulations failed. ",
      "Inspect the returned/saved records for error messages.",
      call. = FALSE
    )
  }

  results
}
