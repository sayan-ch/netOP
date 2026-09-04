test_that("public namespace matches the explicit allowlist", {
  allowlist <- readLines(system.file("API", package = "netOP"), warn = FALSE)
  expect_setequal(getNamespaceExports("netOP"), allowlist)
})

test_that("common Matrix generics are re-exported for sparse networks", {
  A <- Matrix::Matrix(matrix(c(0, 1, 1, 0), 2, 2), sparse = TRUE)

  expect_equal(netOP::mean(A), 0.5)
  expect_equal(netOP::diag(A), c(0, 0))
  expect_equal(netOP::rowMeans(A), c(0.5, 0.5))
  expect_equal(netOP::rowSums(A), c(1, 1))
  expect_equal(netOP::colMeans(A), c(0.5, 0.5))
  expect_equal(netOP::colSums(A), c(1, 1))
  expect_equal(netOP::sum(A), 2)
})

test_that("implementation helpers remain internal", {
  exports <- getNamespaceExports("netOP")
  denylist <- c(
    "source_all", "ecv_stability_bm", "ncv_stability_bm",
    "holdout.evaluation.fast.all", "iter.SVD.core.fast.all", "ECV.BM",
    "ncv_bm", "plot_rdpg_comparison", "plot_blockmodel_comparison",
    "auc_cpp", "outer_add_cpp", "procrustes_translated_cpp",
    "connected_components_dense_cpp", "connected_components_sparse_cpp",
    "shortest_path_distances_dense_cpp",
    "shortest_path_distances_sparse_cpp", "lsm_pgd_cpp"
  )
  expect_false(any(denylist %in% exports))
  expect_false(any(grepl("^validate_|_cpp$", exports)))
})

test_that("S3 methods are registered and dispatch", {
  expect_true(is.function(getS3method("print", "sonnet")))
  expect_true(is.function(getS3method("summary", "netcrop_blockmodel")))
  expect_true(is.function(getS3method("plot", "netcrop_rdpg")))
})

test_that("randnet is not a declared dependency", {
  description <- packageDescription("netOP")
  fields <- c("Depends", "Imports", "Suggests", "Enhances")
  declared <- paste(unlist(description[fields]), collapse = " ")
  expect_false(grepl("randnet", declared, fixed = TRUE))
})
