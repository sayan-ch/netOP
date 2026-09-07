test_that("strict parallel backends execute workers and report task errors", {
  # These tiny backend probes intentionally use two workers even when a
  # constrained runner reports one CPU. Relax only the test-local soft limit;
  # keep the hard limit and restore the option after this test.
  withr::local_options(list(
    parallelly.maxWorkers.localhost = c(soft = 2, hard = 3)
  ))
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  for (force_windows in c(FALSE, TRUE)) {
    worker_pids <- uni_mclapply(
      1:2, function(x) Sys.getpid(), ncores = 2,
      force_windows = force_windows, stop_on_error = TRUE
    )
    expect_true(all(unlist(worker_pids) != Sys.getpid()))
    expect_identical(
      uni_mclapply(1:4, function(x) x * x, ncores = 2,
                   force_windows = force_windows, stop_on_error = TRUE),
      lapply(1:4, function(x) x * x)
    )
    expect_error(
      uni_mclapply(1:2, function(x) {
        if (x == 2L) stop("intentional worker failure")
        x
      }, ncores = 2, force_windows = force_windows, stop_on_error = TRUE),
      "Worker error in task\\(s\\) 2: intentional worker failure"
    )
  }
})

test_that("multisession restores the caller plan and environment on both paths", {
  # These tiny backend probes intentionally use two workers even when a
  # constrained runner reports one CPU. Relax only the test-local soft limit;
  # keep the hard limit and restore the option after this test.
  withr::local_options(list(
    parallelly.maxWorkers.localhost = c(soft = 2, hard = 3)
  ))
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  old_plan <- future::plan()
  old_environment <- Sys.getenv("RENV_CONFIG_SYNCHRONIZED_CHECK", unset = NA)
  uni_mclapply(1:2, identity, ncores = 2, force_windows = TRUE,
               stop_on_error = TRUE)
  expect_identical(future::plan(), old_plan)
  expect_identical(
    Sys.getenv("RENV_CONFIG_SYNCHRONIZED_CHECK", unset = NA), old_environment
  )
  expect_error(
    uni_mclapply(1:2, function(x) stop("restore after error"),
                 ncores = 2, force_windows = TRUE, stop_on_error = TRUE),
    "restore after error"
  )
  expect_identical(future::plan(), old_plan)
  expect_identical(
    Sys.getenv("RENV_CONFIG_SYNCHRONIZED_CHECK", unset = NA), old_environment
  )
})
