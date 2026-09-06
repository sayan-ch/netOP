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
fit <- sonnet(
  A,
  K = 3,
  num_subnetworks = 2,
  overlap_size = 50,
  ncores = 1,
  seed = 101,
  verbose = FALSE,
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
plot(selection)
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

- **NETCROP**: NETwork CRoss-Validation using Overlapping Partitions.
- **SONNET**: Subsampling ON NETwork
- **SBM**: stochastic block model.
- **DCBM**: degree-corrected stochastic block model.
- **RDPG**: random dot product graph.
- **LSM**: latent-space model.
- **ASE**: adjacency spectral embedding.
- **ECV**: edge cross-validation.
- **NCV**: node cross-validation.


## Citations

For NETCROP, cite:

> Chakrabarty, S., Sengupta, S., and Chen, Y. (2026). Network
> Cross-Validation and Model Selection via Subsampling. arXiv:2504.06903.
> <https://doi.org/10.48550/arXiv.2504.06903>

NETCROP means “NETwork CRoss-Validation using Overlapping Partitions.”

For SONNET, cite:

> Chakrabarty, S., Sengupta, S., and Chen, Y. (2025). Subsampling
> Based Community Detection for Large Networks. Statistica Sinica (35)
> 1627 -- 1648. <https://doi.org/10.5705/ss.202022.0108>

Run `citation("netOP")` for machine-readable
citations and see `inst/COPYRIGHTS` for file-level provenance.

## Implementation Disclosures

* `netOP` provides self-contained wrappers around an ECV (<https://doi.org/10.1093/biomet/asaa006>)
implementation derived from CRAN `randnet` 1.0. Installing or using netOP does not require `randnet`.
The ECV-specific implementation helpers are internal and are not part of the netOP public API.

* The NCV wrapper is the `netOP` author's implementation of the exact algorithm in
<https://doi.org/10.1080/01621459.2016.1246365>. Numerical-stability measures and failsafes were added
without altering that algorithm.

## Support and development

Report problems at <https://github.com/sayan-ch/netOP/issues>. Contributions
are welcome under [CONTRIBUTING.md](CONTRIBUTING.md). The package is licensed
under GPL (>= 2).
