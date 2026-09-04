library(netOP)

set.seed(2026)
ncores <- max(1L, floor(parallel::detectCores() / 2))
A <- generate_dcbm(
  n = 3000,
  K = 5,
  average_degree = 80,
  ncores = ncores,
  seed = 1
)

print(system.time({
  sonnet_fit <- sonnet(A, K = 5, ncores = ncores, seed = 2)
}))

print(system.time({
  netcrop_fit <- netcrop_blockmodel(
    A,
    K_candidates = 1:10,
    nrep = 1,
    losses = "sse",
    ncores = ncores,
    seed = 3
  )
}))

print(sonnet_fit)
print(netcrop_fit)
print(sessionInfo())
