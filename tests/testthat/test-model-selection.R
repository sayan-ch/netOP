small_block_network <- function(seed = 20) {
  generate_sbm(n = 24, K = 2, alpha = 0.55, beta = 0.08,
               representation = "dense", seed = seed, ncores = 1)
}

test_that("canonical model-selection openings support positional use", {
  expect_identical(names(formals(netcrop_blockmodel))[1:2],
                   c("A", "K_candidates"))
  expect_identical(names(formals(netcrop_rdpg))[1:2],
                   c("A", "d_candidates"))
  expect_identical(names(formals(netcrop_lsm))[1:2],
                   c("A", "d_candidates"))
  expect_identical(names(formals(netcrop_tune_regularizer))[1:3],
                   c("A", "K", "tau_candidates"))
  expect_false("tau_candidates" %in% names(formals(oracle_plotter)))
  expect_identical(formals(oracle_plotter)$include_sonnet_tau_zero, TRUE)
  expect_identical(formals(oracle_plotter)$include_spectral_tau_zero, TRUE)
  expect_identical(eval(formals(netcrop_tune_regularizer)$cluster_engine)[1L],
                   "clara")
  expect_identical(names(formals(ecv_stability_blockmodel))[1:2],
                   c("A", "max_K"))
  expect_identical(names(formals(ecv_stability_rdpg))[1:2],
                   c("A", "max_d"))
  expect_identical(names(formals(ncv_stability_blockmodel))[1:2],
                   c("A", "max_K"))
  expect_identical(formals(netcrop_blockmodel)$losses, "sse")
  expect_identical(formals(netcrop_rdpg)$losses, "sse")
  expect_identical(formals(netcrop_lsm)$losses, "sse")
  expect_identical(formals(ecv_stability_blockmodel)$losses, "sse")
  expect_identical(formals(ncv_stability_blockmodel)$losses, "sse")
  expect_identical(formals(sonnet)$laplacian, FALSE)
  expect_identical(formals(sonnet)$regularize_tau, 0)
  expect_identical(formals(sonnet)$regularize_subnetworks, TRUE)
})

regularizer_outcome <- function(
    algorithm,
    candidates,
    selected,
    K = 2L,
    model = "DCBM") {
  output <- list(
    algorithm = algorithm,
    K = K,
    model = model,
    tau_candidates = candidates,
    tau_hat = selected,
    seed = 101L,
    options = list(
      use_laplacian = FALSE,
      matching_method = "greedy",
      spectral_options = list(),
      cluster_engine = "clara",
      cluster_options = list()
    ),
    num_subnetworks = 2L,
    effective_overlap_size = 8L
  )
  if (algorithm == "NETCROP") {
    output$best_tau_each_rep <- data.frame(
      repetition = 1L,
      loss = "sse",
      tau_hat = selected
    )
  }
  class(output) <- "netcrop_regularizer"
  output
}

test_that("oracle derives a common tuner grid and keeps zero engine-specific", {
  skip_if_not_installed("ggplot2")
  A <- small_block_network(51)
  truth <- get_generator_parameters(A)$g_true
  netcrop <- regularizer_outcome("NETCROP", c(0.1, 0.2), 0.2)
  dkest <- regularizer_outcome("DKEST", c(0.2, 0.3), 0.2)

  expect_warning(
    plot <- oracle_plotter(
      A,
      g_true = truth,
      netcrop_outcomes = netcrop,
      dkest_outcomes = dkest,
      engines = "spectral_cluster",
      include_sonnet_tau_zero = FALSE,
      include_spectral_tau_zero = TRUE,
      spectral_cluster_options = list(
        spectral_engine = "base",
        cluster_engine = "clara"
      ),
      ncores = 1,
      seed = 52,
      verbose = FALSE
    ),
    "using their intersection"
  )
  metadata <- attr(plot, "metadata")
  expect_equal(metadata$tau_candidates, 0.2)
  expect_false(metadata$candidate_grids_match)
  expect_true(metadata$include_spectral_tau_zero)
  expect_true(any(attr(plot, "accuracy_data")$method == "tau = 0"))
  fit_engines <- vapply(
    metadata$fit_parameters,
    function(fit) fit$parameters$cluster_engine,
    character(1)
  )
  expect_true(all(fit_engines == "clara"))
})

test_that("oracle rejects missing and disjoint tuner candidate grids", {
  skip_if_not_installed("ggplot2")
  A <- small_block_network(53)
  truth <- get_generator_parameters(A)$g_true
  expect_error(
    oracle_plotter(A, g_true = truth, engines = "spectral_cluster"),
    "At least one"
  )
  expect_error(
    oracle_plotter(
      A,
      g_true = truth,
      netcrop_outcomes = regularizer_outcome("NETCROP", 0.1, 0.1),
      dkest_outcomes = regularizer_outcome("DKEST", 0.2, 0.2),
      engines = "spectral_cluster",
      verbose = FALSE
    ),
    "no common tau candidates"
  )
})

test_that("NETCROP regularizer accepts each clustering backend", {
  A <- small_block_network(54)
  configurations <- list(
    clara = list(samples = 2L),
    kmeans = list(nstart = 2L, iter.max = 20L),
    pam = list(do.swap = FALSE, pamonce = 6L)
  )
  for (engine in names(configurations)) {
    fit <- netcrop_tune_regularizer(
      A,
      K = 2,
      tau_candidates = 0.1,
      num_subnetworks = 2,
      overlap_size = 8,
      nrep = 1,
      spectral_options = list(spectral_engine = "base"),
      cluster_engine = engine,
      cluster_options = configurations[[engine]],
      ncores = 1,
      seed = 55,
      verbose = FALSE,
      retain_intermediates = "minimal"
    )
    expect_identical(fit$options$cluster_engine, engine)
  }
})

test_that("SONNET runs deterministically on a benchmark-derived small case", {
  A <- small_block_network()
  fit_1 <- sonnet(A, 2, num_subnetworks = 2, overlap_size = 8,
                  ncores = 1, seed = 21, verbose = FALSE,
                  spectral_engine = "base")
  fit_2 <- sonnet(A, 2, num_subnetworks = 2, overlap_size = 8,
                  ncores = 1, seed = 21, verbose = FALSE,
                  spectral_engine = "base")
  expect_s3_class(fit_1, "sonnet")
  expect_identical(fit_1$g_hat, fit_2$g_hat)
})

test_that("NETCROP block-model selection returns a classed result", {
  A <- small_block_network(22)
  fit <- netcrop_blockmodel(
    A, 1:2, num_subnetworks = 2, overlap_size = 8, nrep = 1,
    losses = "sse", ncores = 1, seed = 23, verbose = FALSE,
    sbm_est_options = list(spectral_cluster = list(spectral_engine = "base")),
    dcbm_est_options = list(spectral_cluster = list(spectral_engine = "base"))
  )
  expect_s3_class(fit, "netcrop_blockmodel")
  expect_identical(fit$algorithm, "NETCROP")
})

test_that("NETCROP block-model regularization controls are shared", {
  A <- small_block_network(27)
  direct <- netcrop_blockmodel(
    A, 1:2, num_subnetworks = 2, overlap_size = 8, nrep = 1,
    laplacian = TRUE, regularize_tau = 0.2,
    ncores = 1, seed = 28, verbose = FALSE,
    sbm_est_options = list(spectral_cluster = list(spectral_engine = "base")),
    dcbm_est_options = list(spectral_cluster = list(spectral_engine = "base"))
  )
  expect_true(direct$sbm_est_options$spectral_cluster$laplacian)
  expect_equal(
    direct$sbm_est_options$spectral_cluster$regularize_tau, 0.2
  )
  expect_identical(
    direct$sbm_est_options$spectral_cluster$laplacian,
    direct$dcbm_est_options$spectral_cluster$laplacian
  )

  legacy <- netcrop_blockmodel(
    A, 1:2, num_subnetworks = 2, overlap_size = 8, nrep = 1,
    ncores = 1, seed = 28, verbose = FALSE,
    sbm_est_options = list(spectral_cluster = list(
      laplacian = TRUE, regularize_tau = 0.2, spectral_engine = "base"
    )),
    dcbm_est_options = list(spectral_cluster = list(
      laplacian = TRUE, regularize_tau = 0.2, spectral_engine = "base"
    ))
  )
  expect_equal(legacy$cv_loss, direct$cv_loss)
})

test_that("SONNET supports subnetwork and whole-network regularization", {
  A <- small_block_network(29)
  subnet <- sonnet(
    A, 2, num_subnetworks = 2, overlap_size = 8, ncores = 1,
    seed = 30, verbose = FALSE, spectral_engine = "base",
    laplacian = TRUE, regularize_tau = 0.1
  )
  whole <- sonnet(
    A, 2, num_subnetworks = 2, overlap_size = 8, ncores = 1,
    seed = 30, verbose = FALSE, spectral_engine = "base",
    laplacian = TRUE, regularize_tau = 0.1,
    regularize_subnetworks = FALSE
  )
  expect_true(subnet$parameters$regularize_subnetworks)
  expect_false(whole$parameters$regularize_subnetworks)
  expect_length(subnet$labels, nrow(A))
  expect_length(whole$labels, nrow(A))
})

test_that("ECV and NCV wrappers use new names and require no randnet calls", {
  expect_false(exists("ecv_stability_bm", asNamespace("netOP"),
                      inherits = FALSE))
  expect_false(exists("ncv_stability_bm", asNamespace("netOP"),
                      inherits = FALSE))
  ecv_text <- paste(deparse(body(ecv_stability_blockmodel)), collapse = " ")
  rdpg_text <- paste(deparse(body(ecv_stability_rdpg)), collapse = " ")
  expect_false(grepl("randnet", paste(ecv_text, rdpg_text), fixed = TRUE))

  A_dense <- small_block_network(24)
  A_sparse <- Matrix::Matrix(A_dense, sparse = TRUE)

  for (A in list(A_dense, A_sparse)) {
    ecv <- ecv_stability_blockmodel(
      A, 2, cv = 2, nrep = 1, ncores = 1, seed = 25,
      losses = "sse", verbose = FALSE
    )
    ncv <- ncv_stability_blockmodel(
      A, 2, cv = 2, nrep = 1, ncores = 1, seed = 26,
      losses = "sse", verbose = FALSE
    )
    expect_s3_class(ecv, "netcrop_blockmodel")
    expect_s3_class(ncv, "netcrop_blockmodel")
    expect_identical(ecv$algorithm, "ECV")
    expect_identical(ncv$algorithm, "NCV")
  }
})
