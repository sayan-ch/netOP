test_that("losses and compiled fallbacks agree", {
  labels <- c(0, 1, 1, 0, 1)
  scores <- c(0.1, 0.8, 0.6, 0.2, 0.8)
  expect_equal(auc(labels, scores, use_cpp = TRUE),
               auc(labels, scores, use_cpp = FALSE))
  expect_equal(outer_add(1:3, 4:5, use_cpp = TRUE),
               outer_add(1:3, 4:5, use_cpp = FALSE))
  expect_equal(sse(1:3, c(1, 2, 4)), 1)
  expect_equal(mae(1:3, c(1, 2, 4)), 1 / 3)
})

test_that("graph compiled and R paths agree", {
  A <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
  expect_equal(connected_components(A, use_cpp = TRUE),
               connected_components(A, use_cpp = FALSE))
  expect_equal(shortest_path_distances(A, use_cpp = TRUE),
               shortest_path_distances(A, use_cpp = FALSE))
})

test_that("dense and sparse graph utilities agree", {
  skip_if_not_installed("Matrix")
  A <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), 3, 3)
  A_sparse <- methods::as(
    methods::as(Matrix::Matrix(A, sparse = TRUE), "generalMatrix"),
    "dgCMatrix"
  )
  expect_equal(connected_components(A), connected_components(A_sparse))
  expect_equal(shortest_path_distances(A),
               shortest_path_distances(A_sparse), ignore_attr = TRUE)
})

test_that("invalid inputs are rejected before fast paths", {
  expect_error(auc(c(0, 2), c(0.1, 0.9)), "only 0 and 1")
  expect_error(connected_components(matrix(1:6, 2, 3)), "square")
  expect_error(softmax(1:3, temperature = 0), "positive")
})
