## Resubmission: netOP 0.1.2

In response to CRAN review of 0.1.1:

* Rephrased the opening of Description to start with "Implements methods".
* Removed single quotes around the method names SONNET and NETCROP.
* Explained both acronyms in Description.
* Updated the package version to 0.1.2.

## R CMD check results

Local macOS, R 4.6.1: 0 errors | 0 warnings | 1 note.

* This is a new submission. The single NOTE is "New submission".
* The full check included examples, tests, vignette rebuilds, and PDF/HTML manuals.
* The installed-package test suite also passed.

## Additional checks

GitHub Actions checks and release-binary checks are recorded in the release
handoff after completion.

The README URL checker reported HTTP 403 for two existing citation DOI links
and a timeout for the Statistica Sinica DOI destination; the same results
occurred on retry. The links were retained.
