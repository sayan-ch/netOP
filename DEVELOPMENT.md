# Developing and releasing netOP

Run commands from the repository root. Release preparation keeps DESCRIPTION
at `0.1.1`; creating a release or release binaries is a separate manual action.

## Routine checks

Install package dependencies and development tools in R:

```r
install.packages(c("devtools", "rcmdcheck", "testthat", "urlchecker", "rhub"))
devtools::install_deps(dependencies = TRUE)
```

CI pins roxygen2 8.0.0 and Rcpp 1.1.2 for reproducible generation. Install those
versions with `remotes::install_version()` if your development library differs.
After changing native exports or roxygen comments:

```r
Rcpp::compileAttributes()
roxygen2::roxygenise()
```

Review and commit generated changes in DESCRIPTION, NAMESPACE, man/,
R/RcppExports.R, and src/RcppExports.cpp. CI regenerates them and rejects drift.

Fast functional checks:

```r
devtools::test(filter = "math-network", stop_on_failure = TRUE)
devtools::test(filter = "generators-estimators", stop_on_failure = TRUE)
devtools::test(filter = "model-selection", stop_on_failure = TRUE)
```

For the complete suite, test an installed package: documentation tests inspect
installed help indexes and the API allowlist. In Terminal:

```sh
R CMD build .
R CMD check --as-cran netOP_0.1.1.tar.gz
mkdir -p local-release-assets/test-library
R CMD INSTALL --preclean --install-tests \
  --library=local-release-assets/test-library netOP_0.1.1.tar.gz
```

Then in a fresh R session:

```r
.libPaths(c(normalizePath("local-release-assets/test-library"), .libPaths()))
library(netOP)
testthat::test_package("netOP", reporter = "summary", stop_on_failure = TRUE)
urlchecker::url_check()
```

Full source checks build both vignettes and the manual; install Pandoc and
LaTeX for these. All test workers are capped at two. Tiny parallel probes use a test-local
200% soft CPU limit when a runner reports one CPU, restoring the option after
each test. Parallel regression tests
exercise native and forced-multisession workers, task errors, restoration of
the caller's future plan/environment, and seeded generation/SONNET agreement.

## Other platforms

Every push to main and pull request runs the current Linux/macOS/Windows checks,
plus Linux R-devel. These jobs build/check source, run examples and vignettes,
and run the installed test suite.

Inspect `.github/workflows/R-CMD-check.yaml` for the exact reproducible setup.
R 4.4/4.5 and Intel macOS also receive the full installed suite when release
binaries are explicitly built. A deployment-target check alone does not replace
testing the complete dependency stack on the oldest supported macOS.

After committing/pushing the R-hub workflow, additional platform checks can be
started with `rhub::rhub_doctor()` and `rhub::rhub_check()`. Before CRAN submission,
also run `devtools::check_win_devel()` or upload the final source tarball to
https://win-builder.r-project.org/. Update cran-comments.md only from completed
checks of the final candidate. Submit the checked source tarball to CRAN;
CRAN produces its own binaries.

## Large-network performance

The programs in inst/benchmarks/ are deliberately outside package checks:

```sh
Rscript --vanilla inst/benchmarks/benchmark_sonnet_netcrop.R
Rscript --vanilla inst/benchmarks/benchmark_regularizers.R
```

Review problem sizes and worker counts before running. Record sessionInfo(),
elapsed time, memory where measured, seeds, dimensions, and worker counts when
comparing versions. Performance measurements do not replace correctness tests.

## Preparing release artifacts later

The release workflow is **manual-only**. Neither pushing main nor publishing a
release starts it. First commit all release changes, require passing checks,
and tag the tested commit with the matching stable version. Then:

```sh
git tag -a v0.1.1 -m "netOP 0.1.1"
git push origin v0.1.1
gh release create v0.1.1 --verify-tag --title "netOP 0.1.1" \
  --notes "See NEWS.md for release changes."
gh workflow run release-binaries.yaml --ref main -f tag=v0.1.1
```

The workflow requires an existing empty release and freezes the tag's commit.
It builds and checks one source tarball with rendered vignettes, transfers that
exact bundle to all nine binary jobs, and checks each installed binary's version,
compiled code, vignettes, and full test suite. Only after all jobs succeed does
one attachment job validate the ten-file inventory, generate SHA256SUMS, and
upload everything to the release. Existing assets are never overwritten.

Monitor with `gh run list --workflow release-binaries.yaml`, then
`gh run watch RUN_ID --exit-status`. The release must have nine binaries, the
package source tarball, and SHA256SUMS. Review with:

```sh
gh release view v0.1.1 --json assets --jq '.assets[].name'
```

If a build fails before attachment, use `gh run rerun RUN_ID --failed`. If an
upload is interrupted, the release may contain a partial set: inspect and remove
those assets explicitly before rerunning the failed attachment job. The
workflow refuses to overwrite them. Never move a published release tag.

After freezing/submitting 0.1.1, a separate development commit can advance main
to 0.1.1.9000 and update NEWS and the software citation. That does not require
another binary release.
