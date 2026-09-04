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
  expect_identical(names(formals(ecv_stability_blockmodel))[1:2],
                   c("A", "max_K"))
  expect_identical(names(formals(ecv_stability_rdpg))[1:2],
                   c("A", "max_d"))
  expect_identical(names(formals(ncv_stability_blockmodel))[1:2],
                   c("A", "max_K"))
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

test_that("ECV and NCV wrappers use new names and require no randnet calls", {
  expect_false(exists("ecv_stability_bm", asNamespace("netOP"),
                      inherits = FALSE))
  expect_false(exists("ncv_stability_bm", asNamespace("netOP"),
                      inherits = FALSE))
  ecv_text <- paste(deparse(body(ecv_stability_blockmodel)), collapse = " ")
  rdpg_text <- paste(deparse(body(ecv_stability_rdpg)), collapse = " ")
  expect_false(grepl("randnet", paste(ecv_text, rdpg_text), fixed = TRUE))

  A <- small_block_network(24)
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
})
