# netOP

[README](https://github.com/sayan-ch/netOP#readme) | [Dictionary](https://github.com/sayan-ch/netOP/blob/main/dictionary.md)

<!-- badges: start -->
[![R-CMD-check](https://github.com/sayan-ch/netOP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sayan-ch/netOP/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`netOP` provides network generation, estimation, embedding, clustering, and
model-selection tools for R. It brings together general graph utilities,
spectral and latent-space methods, SONNET and NETCROP for SBM, DCBM, RDPG, LSM and regularization selection.

## Installation

### Binary package (recommended for most macOS and Windows systems)

The 0.1.1 release provides binaries for R 4.4, R 4.5, and R 4.6 on
Apple Silicon and Intel macOS and on x86-64 Windows. The R 4.4 and R 4.5 Apple Silicon binaries and all Intel binaries target macOS
11 or newer. The R 4.6 Apple Silicon binary follows the official R 4.6 runtime
and requires macOS 14 or newer.

The following code selects the matching
asset and installs netOP without compiling it locally:

```r
local({
  version <- "0.1.1"
  r_series <- paste(
    R.version$major,
    sub("\\..*$", "", R.version$minor),
    sep = "."
  )
  if (!r_series %in% c("4.4", "4.5", "4.6")) {
    stop(
      "The current netOP binary release requires R 4.4.x, R 4.5.x, or R 4.6.x."
    )
  }

  r_architecture <- tolower(R.version$arch)
  asset <- if (.Platform$OS.type == "windows") {
    architecture <- if (identical(r_architecture, "x86_64")) {
      "x86_64"
    } else {
      stop(
        "No native netOP binary is available for this Windows R architecture. ",
        "On Windows ARM64, use the x86-64 build of R or install netOP from source."
      )
    }
    sprintf("netOP_%s_R-%s_%s.zip", version, r_series, architecture)
  } else if (identical(Sys.info()[["sysname"]], "Darwin")) {
    architecture <- if (grepl("arm64|aarch64", r_architecture)) {
      "arm64"
    } else if (identical(r_architecture, "x86_64")) {
      "x86_64"
    } else {
      stop("No netOP binary is available for this macOS architecture.")
    }
    sprintf("netOP_%s_R-%s_%s.tgz", version, r_series, architecture)
  } else {
    stop("Use the source installation instructions below on this platform.")
  }

  cran_repository <- getOption("repos")[["CRAN"]]
  if (is.null(cran_repository) || is.na(cran_repository) ||
      identical(cran_repository, "@CRAN@")) {
    cran_repository <- "https://cloud.r-project.org"
  }
  install.packages(
    c("cluster", "irlba", "Matrix", "Rcpp", "RcppEigen", "RSpectra", "tibble"),
    repos = cran_repository
  )

  binary_url <- sprintf(
    "https://github.com/sayan-ch/netOP/releases/download/v%s/%s",
    version,
    asset
  )
  binary_directory <- tempfile("netop-binary-")
  dir.create(binary_directory)
  on.exit(unlink(binary_directory, recursive = TRUE), add = TRUE)
  binary_file <- file.path(
    binary_directory,
    sprintf("netOP_%s.%s", version, tools::file_ext(asset))
  )
  download_status <- download.file(binary_url, binary_file, mode = "wb")
  if (!identical(download_status, 0L) || !file.exists(binary_file) ||
      file.info(binary_file)$size <= 0) {
    stop("The netOP binary could not be downloaded.")
  }
  install.packages(binary_file, repos = NULL, type = "binary")
})
```

Linux distributions do not share a portable R binary-package format. Each
GitHub release therefore includes a standard source tarball for Linux and
other Unix systems. Within the supported macOS and Windows versions, separate
assets are needed for each R major/minor series and processor architecture;
ordinary operating-system patch updates do not require another asset.

Install the released 0.1.1 source package on Linux or
another Unix-like system with:

```r
install.packages(
  c("cluster", "irlba", "Matrix", "Rcpp", "RcppEigen", "RSpectra", "tibble")
)
install.packages(
  "https://github.com/sayan-ch/netOP/releases/download/v0.1.1/netOP_0.1.1.tar.gz",
  repos = NULL,
  type = "source"
)
```

Compiling the source package requires the development tools described below.

### Development version

Install the latest development version from GitHub with:

```r
install.packages("remotes")
remotes::install_github("sayan-ch/netOP")
```

Because netOP contains C++ code, installing the development or source version
requires a working package-development toolchain:

- **macOS:** install Apple's Xcode Command Line Tools by running
  `xcode-select --install` in Terminal.
- **Windows:** install the version of
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R
  version.
- **Linux:** install GNU Make, a C++ compiler, and the R development headers.
  On Debian or Ubuntu, these are normally provided by `build-essential` and
  `r-base-dev`.

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

Generators use sparse output by default where supported. `netOP` exposes
`mean()`, `sum()`, `diag()`, `rowMeans()`, `rowSums()`, `colMeans()`, and
`colSums()` for sparse networks (six Matrix re-exports and a `sum()` wrapper
around `base::sum()`), so these familiar operations dispatch correctly for
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
are welcome under [CONTRIBUTING.md](https://github.com/sayan-ch/netOP/blob/main/CONTRIBUTING.md). The package is licensed
under GPL (>= 2).
