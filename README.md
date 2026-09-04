# netOP

<!-- badges: start -->
[![R-CMD-check](https://github.com/sayan-ch/netOP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sayan-ch/netOP/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`netOP` provides network generation, estimation, embedding, clustering, and
model-selection tools for R. It brings together general graph utilities,
spectral and latent-space methods, SONNET, NETCROP, edge cross-validation
(ECV), node cross-validation (NCV), and regularization selection.

The development version is distributed through GitHub:

```r
install.packages("remotes")
remotes::install_github("sayan-ch/netOP")
```

## Quick start

```r
library(netOP)

A <- generate_sbm(
  n = 60,
  K = 2,
  alpha = 0.45,
  beta = 0.08,
  representation = "dense",
  seed = 100,
  ncores = 1
)

embedding <- ase(A, d = 2)
fit <- spectral_cluster(
  A,
  K = 2,
  spectral_engine = "base",
  cluster_engine = "kmeans"
)

selection <- netcrop_blockmodel(
  A,
  K_candidates = 1:3,
  num_subnetworks = 2,
  overlap_size = 20,
  nrep = 1,
  losses = "sse",
  ncores = 1,
  seed = 101,
  verbose = FALSE,
  sbm_est_options = list(spectral_cluster = list(spectral_engine = "base")),
  dcbm_est_options = list(spectral_cluster = list(spectral_engine = "base"))
)
summary(selection)
```

Public model-selection entry points consistently begin with the primary
network and model-size candidate set:

- `netcrop_blockmodel(A, K_candidates, ...)`
- `netcrop_rdpg(A, d_candidates, ...)`
- `netcrop_lsm(A, d_candidates, ...)`
- `netcrop_tune_regularizer(A, K, tau_candidates, ...)`
- `ecv_stability_blockmodel(A, max_K, ...)`
- `ecv_stability_rdpg(A, max_d, ...)`
- `ncv_stability_blockmodel(A, max_K, ...)`

The ECV maximum interfaces intentionally evaluate every candidate from one to
the maximum because that is part of the incorporated algorithm.

## Algorithm families

- Generators: ER, SBM, DCBM, RDPG, and latent-space models, with dense and
  sparse adjacency output where supported.
- Core methods: losses, graph utilities, decompositions, ASE, spectral
  clustering, latent-space projected-gradient fitting, and SBM/DCBM estimators.
- Scalable fitting: SONNET overlapping-subnetwork clustering.
- Model selection: NETCROP for block models, RDPGs, latent-space models, and
  spectral regularization; ECV and NCV stability wrappers; DKEST.

## Citations and implementation disclosures

For NETCROP, cite:

> Chakrabarty, S., Sengupta, S., and Chen, Y. (2026). Network
> Cross-Validation and Model Selection via Subsampling. arXiv:2504.06903.
> <https://doi.org/10.48550/arXiv.2504.06903>

NETCROP means “NETwork CRoss-Validation using Overlapping Partitions.”

For ECV, cite:

> Li, T., Levina, E., and Zhu, J. (2020). Network cross-validation by edge
> sampling. *Biometrika*, 107(2), 257–276.
> <https://doi.org/10.1093/biomet/asaa006>

netOP provides self-contained wrappers around an ECV implementation derived
from CRAN `randnet` 1.0. Installing or using netOP does not require `randnet`.
The ECV-specific implementation helpers are internal and are not part of the
netOP public API.

For NCV, cite:

> Chen, K. and Lei, J. (2018). Network Cross-Validation for Determining the
> Number of Communities in Network Data. *Journal of the American Statistical
> Association*, 113(521), 241–251.
> <https://doi.org/10.1080/01621459.2016.1246365>

The NCV wrapper is the netOP author's implementation of the exact algorithm in
Chen and Lei (2018). Numerical-stability measures and failsafes were added
without altering that algorithm. Run `citation("netOP")` for machine-readable
citations and see `inst/COPYRIGHTS` for file-level provenance.

## Support and development

Report problems at <https://github.com/sayan-ch/netOP/issues>. Contributions
are welcome under [CONTRIBUTING.md](CONTRIBUTING.md). The package is licensed
under GPL (>= 2).
