# netOP 0.1.1

- Revised package metadata and expanded the method and software citations for
  the initial CRAN submission.

# netOP 0.1.0

- Converted multiple network analysis methods and helper codes into an installable R
  package with registered Rcpp interfaces and curated exports.
- Added SONNET and NETCROP along with wrappers for ECV, NCV, and DKEST, documentation and citations.
- Standardized probability-estimator names as `estimate_sbm_P_hat()` and
  `estimate_dcbm_P_hat()`, and made network-file replacement opt-in through
  `write_network(..., overwrite = TRUE)`.
- Added tests, deterministic examples, benchmark scripts, a getting-started
  vignette, pkgdown configuration, and cross-platform R CMD check workflows.
- Licensed the package as GPL (>= 2) and documented randnet-derived ECV
  provenance.
- Made CLARA the configurable default clustering backend for NETCROP
  regularizer tuning. `oracle_plotter()` now derives a common candidate grid
  from tuner outcomes, diagnoses mismatched grids, and includes separate
  SONNET and spectral-clustering `tau = 0` baselines by default.
