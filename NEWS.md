# netOP 0.0.0.9000

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
