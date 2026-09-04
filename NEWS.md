# netOP 0.0.0.9000

- Converted the authoritative network-analysis sources into an installable R
  package with registered Rcpp interfaces and curated exports.
- Added SONNET, NETCROP, ECV, NCV, and DKEST documentation and citations.
- Renamed `ecv_stability_bm()` to `ecv_stability_blockmodel()` and exposed the
  stable NCV wrapper as `ncv_stability_blockmodel()`.
- Standardized probability-estimator names as `estimate_sbm_P_hat()` and
  `estimate_dcbm_P_hat()`, and made network-file replacement opt-in through
  `write_network(..., overwrite = TRUE)`.
- Added tests, deterministic examples, benchmark scripts, a getting-started
  vignette, pkgdown configuration, and cross-platform R CMD check workflows.
- Relicensed the package as GPL (>= 2) and documented randnet-derived ECV
  provenance.
