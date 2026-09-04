test_that("generators are deterministic and expose compact truth", {
  A_1 <- generate_sbm(n = 20, K = 2, alpha = 0.5, beta = 0.1,
                      seed = 10, ncores = 1)
  A_2 <- generate_sbm(n = 20, K = 2, alpha = 0.5, beta = 0.1,
                      seed = 10, ncores = 1)
  expect_identical(A_1, A_2)
  truth <- get_generator_parameters(A_1)
  expect_length(truth$g_true, 20)
  expect_false(any(names(attributes(A_1)) == "P"))
})

test_that("dense and sparse generation agree for fixed seeds", {
  skip_if_not_installed("Matrix")
  dense <- generate_er(n = 15, p = 0.2, representation = "dense",
                       seed = 11, ncores = 1)
  sparse <- generate_er(n = 15, p = 0.2, representation = "sparse",
                        seed = 11, ncores = 1)
  expect_equal(as.matrix(dense), as.matrix(sparse), ignore_attr = TRUE)
})

test_that("estimators and embedders return documented shapes", {
  A <- generate_sbm(n = 24, K = 2, alpha = 0.5, beta = 0.1,
                    seed = 12, ncores = 1)
  g <- get_generator_parameters(A)$g_true
  expect_equal(dim(estimate_sbm(A, g, K = 2)), c(2, 2))
  expect_equal(dim(estimate_sbm_P_hat(A, g, K = 2)), c(24, 24))
  embedding <- ase(A, d = 2)
  expect_equal(dim(embedding$Z_hat), c(24, 2))
  clustering <- spectral_cluster(A, K = 2, spectral_engine = "base",
                                 cluster_engine = "kmeans")
  expect_length(clustering$g_hat, 24)
})

test_that("network edge-list I/O round trips", {
  A <- generate_er(n = 10, p = 0.25, seed = 13, ncores = 1)
  path <- tempfile(fileext = ".csv")
  write_network(A, path)
  restored <- read_network(path, representation = "dense")
  expect_equal(as.matrix(restored), as.matrix(A), ignore_attr = TRUE)
  expect_error(write_network(A, path), "already exists")
  expect_invisible(write_network(A, path, overwrite = TRUE))
})
