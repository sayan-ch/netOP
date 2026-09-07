# Contributing to netOP

Please open an issue before proposing algorithmic changes. Mathematical
behavior must remain stable unless a defect has a minimal reproducer and the
maintainer approves the change.

Contributions must:

1. use underscore-separated identifiers except where R S3 dispatch or a
   preserved upstream boundary requires otherwise; keep helpers internal and
   list public exports explicitly in `inst/API`;
2. update `dictionary.md` for every function, signature, visibility,
   dependency, fallback, return-value, or source-tree change;
3. document public functions with roxygen2 and regenerate `NAMESPACE` and `man/`
   rather than editing generated files;
4. add deterministic, lightweight tests and examples, keeping large workloads
   in `inst/benchmarks/`;
5. run `Rcpp::compileAttributes()`, `roxygen2::roxygenise()`, the testthat suite,
   `R CMD build`, and `R CMD check --as-cran`.

Do not introduce `sourceCpp()`, runtime package installation, machine-specific
paths, credentials, or a dependency on large packages.

See [DEVELOPMENT.md](DEVELOPMENT.md) for reproducible checks, minimum-R testing,
and the draft-release procedure. Generated files must be committed and must
remain unchanged when regenerated with the CI-pinned generators.
