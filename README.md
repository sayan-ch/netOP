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
  n = 200,
  K = 3,
  alpha = 0.45,
  beta = 0.08,
  representation = "dense",
  seed = 100,
  ncores = 1
)

parameters <- get_generator_parameters(A)
table(parameters$g_true)

embedding <- ase(A, d = 3)
fit <- spectral_cluster(
  A,
  K = 3,
  spectral_engine = "base",
  cluster_engine = "kmeans"
)

selection <- netcrop_blockmodel(
  A,
  K_candidates = 1:5,
  num_subnetworks = 2,
  overlap_size = 50,
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

Generators use sparse output by default where supported. `netOP` re-exports
`mean()`, `sum()`, `diag()`, `rowMeans()`, `rowSums()`, `colMeans()`, and
`colSums()` from Matrix, so these familiar operations dispatch correctly for
sparse networks after `library(netOP)`. Use `representation = "dense"` only
when a dense matrix is required by a downstream workflow. Generator truth and
settings can be recovered with `get_generator_parameters()` as shown above.

For reproducible computations, supply `seed` explicitly. Documentation
examples use `ncores = 1` for portability and repeatability; supported routines
may use more workers in production while retaining their documented seeding
behavior.

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

See the package articles for a [method-selection
guide](vignettes/choosing-a-method.Rmd)
and a [getting-started
workflow](https://github.com/sayan-ch/netOP/blob/main/vignettes/getting-started.Rmd).

## Glossary

- **SBM**: stochastic block model.
- **DCBM**: degree-corrected stochastic block model.
- **RDPG**: random dot product graph.
- **LSM**: latent-space model.
- **ASE**: adjacency spectral embedding.
- **ECV**: edge cross-validation.
- **NCV**: node cross-validation.
- **NETCROP**: NETwork CRoss-Validation using Overlapping Partitions.

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
