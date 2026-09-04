library(netOP)

set.seed(2026)
ncores <- max(1L, floor(parallel::detectCores() / 2))
A <- generate_dcbm(
  n = 1000,
  K = 5,
  average_degree = 60,
  ncores = ncores,
  seed = 10
)
tau_candidates <- seq(0, 2, by = 0.05)

print(system.time({
  netcrop_tau <- netcrop_tune_regularizer(
    A,
    K = 5,
    tau_candidates = tau_candidates,
    use_dcbm = TRUE,
    use_laplacian = TRUE,
    nrep = 5,
    ncores = ncores,
    seed = 11
  )
}))

print(system.time({
  dkest_tau <- dkest_tune_regularizer(
    A,
    K = 5,
    tau_candidates = tau_candidates,
    use_dcbm = TRUE,
    use_laplacian = TRUE,
    ncores = ncores,
    seed = 12
  )
}))

print(netcrop_tau)
print(dkest_tau)
print(sessionInfo())
