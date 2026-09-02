# NETCROP dimension selection for symmetric RDPGs
#
# This file preserves the original four-stage computational flow: decompose
# overlapping subnetworks, construct every candidate embedding, align each
# non-reference embedding through the overlap, and evaluate every unordered
# pair of non-overlap pieces. Source the numbered helpers, x1_sonnet.R, and
# xx_helpers.R before this file.

netcrop_rdpg <- function(
    A,
    d_candidates,
    num_subnetworks,
    overlap_size,
    nrep = 1L,
    losses = c("sse", "bin_dev", "auc_as_loss"),
    eig_options = list(),
    ncores = max(
      floor(parallel::detectCores() / 2),
      1L,
      na.rm = TRUE
    ),
    seed = NULL,
    verbose = FALSE,
    force_windows = FALSE,
    ram_check = TRUE) {
  call <- match.call()
  required_helpers <- c(
    "uni_mclapply", "is_symmetric_matrix", "eig_decomp", "procrustes",
    "netcrop_splitter", "modal", "estimate_spectral_decomp_ram",
    "estimate_matrix_product_ram", "report_ram_formula"
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
      "Source the numbered helper files, x1_sonnet.R, and xx_helpers.R ",
      "before x3_netcrop_rdpg.R. Missing: ",
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
  if (!(is.numeric(A) || inherits(A, "Matrix")) || any(!is.finite(A))) {
    stop("A must contain only finite numeric values.", call. = FALSE)
  }
  if (!is_symmetric_matrix(A)) {
    stop("A must be symmetric for symmetric RDPG NETCROP.", call. = FALSE)
  }
  if (any(diag(A) != 0)) {
    stop("A must have a zero diagonal (no self-loops).", call. = FALSE)
  }
  if (any(A < 0)) {
    stop("A must be non-negative.", call. = FALSE)
  }
  if (sum(A) == 0) {
    stop("A is an empty graph; RDPG NETCROP is undefined.", call. = FALSE)
  }
  
  n <- nrow(A)
  num_subnetworks <- validate_count(
    num_subnetworks,
    "num_subnetworks",
    minimum = 2L
  )
  overlap_size <- validate_count(overlap_size, "overlap_size", minimum = 1L)
  nrep <- validate_count(nrep, "nrep", minimum = 1L)
  ncores <- validate_count(ncores, "ncores", minimum = 1L)
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
  
  rho_hat <- mean(rowSums(A)) / (n - 1)
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
      stop_on_error = TRUE
    )
  })
  
  raw_lookup <- function(repetition, split_id) {
    which(raw_grid$repetition == repetition & raw_grid$split == split_id)
  }
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
        raw_id <- raw_lookup(repetition, split_id)
        U <- raw_output[[raw_id]]$U[, seq_len(d), drop = FALSE]
        scales <- sqrt(abs(raw_output[[raw_id]]$values[seq_len(d)]))
        sweep(U, 2L, scales, `*`)
      },
      ncores = embedding_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE
    )
  })
  
  embedding_lookup <- function(repetition, split_id, candidate_id) {
    which(
      candidate_grid$repetition == repetition &
        candidate_grid$split == split_id &
        candidate_grid$candidate_id == candidate_id
    )
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
        standard_id <- embedding_lookup(repetition, 1L, candidate_id)
        rotation <- procrustes(
          X = X[seq_len(overlap_count), , drop = FALSE],
          X_star = embedding_output[[standard_id]][
            seq_len(overlap_count), , drop = FALSE
          ]
        )$rotation
        tcrossprod(
          X[-seq_len(overlap_count), , drop = FALSE],
          t(rotation)
        )
      },
      ncores = alignment_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE
    )
  })
  
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
        left_id <- embedding_lookup(repetition, split_left, candidate_id)
        right_id <- embedding_lookup(repetition, split_right, candidate_id)
        P_hat <- tcrossprod(
          aligned_output[[left_id]],
          aligned_output[[right_id]]
        ) * rho_hat
        P_hat <- pmax(P_hat, 1e-6)
        P_hat <- pmin(P_hat, 1 - 1e-6)
        records <- lapply(losses, function(loss_name) {
          value <- loss_functions[[loss_name]](
            as.numeric(A_test),
            as.numeric(P_hat)
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
          data.frame(
            repetition = repetition,
            d = d_candidates[candidate_id],
            loss = loss_name,
            split_left = split_left,
            split_right = split_right,
            loss_value = as.numeric(value),
            stringsAsFactors = FALSE
          )
        })
        do.call(rbind, records)
      },
      ncores = loss_ncores,
      force_windows = force_windows,
      stop_on_error = TRUE
    )
  })
  
  cv_all_loss <- do.call(rbind, loss_output)
  rownames(cv_all_loss) <- NULL
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
  eigen_diagnostics <- do.call(rbind, lapply(seq_along(raw_output), function(id) {
    data.frame(
      repetition = raw_grid$repetition[id],
      split = raw_grid$split[id],
      negative_eigenvalues = raw_output[[id]]$negative_eigenvalues,
      stringsAsFactors = FALSE
    )
  }))
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
    rho_hat = rho_hat,
    eig_options = resolved_eig_options,
    splitter = splitter,
    raw_grid = raw_grid,
    candidate_grid = candidate_grid,
    loss_grid = loss_grid,
    raw_output = raw_output,
    embedding_output = embedding_output,
    aligned_output = aligned_output,
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
print.netcrop_rdpg <- function(x, ...) {
  cat("NETCROP results for symmetric RDPG\n")
  cat("-----------------------------------\n")
  print(x$overall_best, row.names = FALSE)
  invisible(x)
}

# Summarize a symmetric-RDPG NETCROP fit.
summary.netcrop_rdpg <- function(object, ...) {
  result <- list(
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
print.summary.netcrop_rdpg <- function(x, ...) {
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
plot.netcrop_rdpg <- function(x, aggregate = TRUE, ...) {
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
          title = "NETCROP CV loss by RDPG dimension",
          x = "Latent dimension (d)",
          y = "Average CV loss",
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
      title = "NETCROP CV loss by RDPG dimension",
      x = "Latent dimension (d)",
      y = "Mean CV loss (plus or minus one SD)"
    ) +
    ggplot2::theme_minimal()
}

system.time(
  net <- generate_rdpg(
    n = 10000, d = 10, ncores = 5
  )
)

system.time(
  net <- RDPG.gen(n = 10000, d = 10, X = NULL, rho = 0.65,
                  ncore = 5)$A
)

system.time(
  rdpg_out <- netcrop_rdpg(
    A = net,
    d_candidates = 1:20,
    num_subnetworks = 3L,
    overlap_size = 8002L,
    nrep = 1L,
    losses = "bin_dev",
    ncores = 5L,
    verbose = TRUE,
    ram_check = FALSE
  )
)
rdpg_out
summary(rdpg_out)
plot(rdpg_out)










