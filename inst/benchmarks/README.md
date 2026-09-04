# netOP benchmarks

These scripts preserve authentic large-workload patterns from the commented
programs in the authoritative `x*.R` files. They are not run by package checks.
Run them explicitly in a fresh R session after installing netOP, record
`sessionInfo()`, and adjust sizes and worker counts for the available machine.

- `benchmark_sonnet_netcrop.R`: large DCBM generation, spectral clustering,
  SONNET, and NETCROP selection.
- `benchmark_regularizers.R`: NETCROP and DKEST regularizer selection.
