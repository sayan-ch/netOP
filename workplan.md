# netOP R Package Work Plan

## Purpose

Convert the current contents of `github/netOP` into a documented, installable,
and GitHub-ready R package named `netOP` while preserving the tested algorithms
in `library_functions/` as closely as possible.

This document is intended to be handed to a coding model as the implementation
specification. The implementing model must read `NAMING_CONVENTION.md` and
`dictionary.md` completely before editing any source file. Those two files are
binding project instructions. If this work plan conflicts with either file,
stop and ask the package owner which instruction takes precedence.

## Recommended Codex configuration

Use **GPT-5.6 Sol with high reasoning effort**. This task requires repository
restructuring, R and Rcpp package integration, API design, documentation,
dependency analysis, and repeated package checks. `high` is the recommended
quality/latency balance. Use `xhigh` if the package owner prefers a more
exhaustive API and security audit and accepts the additional latency. Reserve
`max` for difficult failures that remain after ordinary package diagnostics.

## Authoritative inputs and fixed decisions

- Work only inside `github/netOP` unless a temporary build directory is needed.
- `library_functions/` is the authoritative replacement for the deleted
  `SONNET_unified*` source trees.
- Preserve the mathematical algorithms and tested behavior. Change function
  bodies only when required for package installation, namespace operation,
  native-code registration, dependency correctness, documentation examples,
  consistent public argument ordering, or a demonstrated defect.
- The benchmark and test programs at the ends of the `x*.R` files are the
  primary source for documentation examples.
- Export useful general helpers, including loss functions, network utilities,
  decompositions, generators, estimators, and embedders. Do not limit exports
  to only SONNET, NETCROP, ECV, and NCV.
- Do not export validators, raw Rcpp functions, C++ helpers, package-loading
  helpers, source-order helpers, or similarly low-level implementation details.
- Export both ECV public wrappers:
  - rename `ecv_stability_bm()` to `ecv_stability_blockmodel()`;
  - retain `ecv_stability_rdpg()` under that name.
- Hide all other ECV-specific functions, including the legacy implementation
  functions used by the two wrappers.
- Treat consistent argument ordering as a desired API improvement. Apply it
  where it can be done systematically and safely; do not perform isolated
  reorderings that leave related functions inconsistent.
- The remote repository is `https://github.com/sayan-ch/netOP.git`.
- Use the initial GitHub repository description:
  `Network analysis and model selection tools for R.` The owner may replace it
  with a more detailed description later.
- Continue implementation directly on `main`; use small, coherent commits and
  push each validated checkpoint to `origin/main`. Never force-push.
- The initial distribution target is GitHub. Prepare for, but do not submit to,
  CRAN until the GitHub version is stable and the owner separately authorizes
  a CRAN submission.

## Non-negotiable safety rules

1. Preserve all pre-existing user changes. Do not use `git reset --hard`,
   `git checkout --`, or other destructive cleanup commands.
2. Before the package conversion, inspect the worktree and create one explicit
   baseline commit containing the current authoritative structure, then push
   it to `origin/main`. This baseline must make the starting state recoverable.
3. Do not include `.DS_Store`, `.Rhistory`, `.Rproj.user`, build products, shared
   libraries, or temporary check directories in the baseline or later commits.
4. Show the exact staged file list and staged diff summary before each commit.
5. Never rewrite Git history or force-push.
6. Do not silently fix a suspected algorithmic problem. Record it as a finding,
   add a reproducer when possible, and ask the owner before changing behavior.
7. Every function addition, removal, rename, relocation, signature change,
   dependency change, fallback change, or visibility change must be reflected
   in `dictionary.md` in the same commit.
8. Monitor the remaining model/context usage during implementation. When less
   than 5% remains, stop starting new work and create a clean continuation
   checkpoint as described in “Low-usage continuation checkpoint” below.

## Low-usage continuation checkpoint

When the agent has less than 5% of its available usage or context budget left,
it must prioritize recoverability over attempting one more implementation step:

1. Stop before beginning a new edit, migration, or long-running test.
2. Finish or safely roll back only the agent's own incomplete atomic edit; do
   not discard pre-existing user work.
3. Update or create `CHECKPOINT.md` with:
   - the current objective and completed phases;
   - the exact next action;
   - files changed and why;
   - unresolved decisions, failures, and risks;
   - validation commands already run and their outcomes;
   - commands the next agent should run first;
   - current branch, `HEAD`, upstream, and worktree status;
   - any process still running and how to inspect or stop it.
4. Run the fastest relevant non-destructive smoke check if enough usage remains.
5. Review `git diff --check`, `git status --short`, and the staged file list.
6. Create a clearly labeled checkpoint commit containing only safe, in-scope
   work when the tree is in a coherent state. Use a message such as
   `Checkpoint netOP package conversion after phase N`.
7. Push the checkpoint to the current authorized branch if pushing is already
   within scope and the commit contains no credentials, generated binaries, or
   known destructive breakage. Otherwise leave the local commit intact and
   state why it was not pushed.
8. End with a concise handoff pointing the next agent to `workplan.md` and
   `CHECKPOINT.md`. Do not claim completion unless all acceptance criteria pass.

## Phase 0: Establish the recoverable baseline

The repository currently contains many tracked deletions and untracked
replacement files. The implementing model must regard this as intentional user
work, not as damage to be reverted.

1. Run and save the results of:
   - `git status --short --branch`;
   - `git diff --stat`;
   - `git diff --cached --stat`;
   - `git ls-files --others --exclude-standard`;
   - `git remote -v`.
2. Confirm that `.gitignore` excludes `.DS_Store`, `.Rhistory`, and
   `.Rproj.user`.
3. Remove ignored metadata from the proposed commit only if it is already
   untracked or explicitly authorized; do not delete tracked material merely
   because it looks obsolete.
4. Stage the authoritative current structure, including
   `library_functions/`, `NAMING_CONVENTION.md`, and `dictionary.md`, together
   with the intended deletions of the replaced source trees.
5. Review the staged list for secrets, credentials, large archives, generated
   binaries, and personal machine state.
6. Commit with a message such as
   `Establish authoritative library_functions baseline`.
7. Push the commit to `origin/main` and verify that local `HEAD` matches
   `origin/main`.
8. Record the baseline commit hash in the implementation handoff.

Do not begin restructuring until this baseline push succeeds. If branch
protection prevents a direct push, stop and ask whether to create a branch and
pull request instead.

## Phase 0.5: Prepare the package-development environment

Install or update only tools and packages required for netOP; do not perform an
unscoped upgrade of every R package or system package on the machine.

1. Inventory and record versions of R, Git, `make`, C/C++ compilers, Fortran,
   Pandoc, qpdf, and the active R library paths.
2. Ensure the local toolchain includes Git, Apple/Xcode command-line build
   tools, `make`, Clang, a compatible Fortran compiler, Pandoc, and qpdf.
3. Ensure the R development stack includes `Rcpp`, `RcppEigen`, roxygen2,
   testthat, rcmdcheck, remotes, devtools, pkgdown, knitr, and rmarkdown.
4. Install the CRAN package `randnet` only in the local development environment
   so the ECV provenance and wrapper behavior can be compared against upstream.
   `randnet` must not become a netOP package dependency or be listed in
   `Depends`, `Imports`, `Suggests`, or `Enhances`.
5. After the source/dependency audit, install only dependencies that remain in
   the final `Imports`, `Suggests`, or `LinkingTo` fields. Do not install the
   current broad `DESCRIPTION` list blindly because several entries may be
   obsolete.
6. Prefer a project-local reproducibility record such as a version inventory or
   lockfile only after the final dependency set is known. Do not commit a large
   package cache or platform-specific binaries.
7. Run a minimal Rcpp compilation smoke test before restructuring all three C++
   sources, and record any compiler warnings separately from algorithm issues.

## Phase 1: Build a complete API and dependency inventory

### 1.1 Function inventory

Parse every top-level R function and every named/exported C++ function in
`library_functions/`. Reconcile the result against `dictionary.md`.

For every function, record:

- source file and line;
- current signature and defaults;
- intended public or internal status;
- callers and dependencies;
- whether it is an S3 method;
- whether it invokes optional compiled code;
- whether it is part of an ECV legacy implementation;
- whether it appears in an end-of-file benchmark block;
- proposed documentation topic and example source.

No function may disappear merely because it is omitted from an initial export
list. Resolve every discrepancy between code and dictionary explicitly.

### 1.2 Export policy

Export functions that are substantively useful to an R user, including these
categories:

- losses, clipping, alignment, and label matching;
- network input/output and graph utilities;
- eigendecomposition, singular decomposition, embeddings, and clustering;
- ER, SBM, DCBM, RDPG, and LSM generators;
- model estimators and probability reconstruction;
- simulation and resource-estimation utilities when they form a usable public
  feature rather than an implementation detail;
- SONNET entry points and parameter selection;
- NETCROP block-model, RDPG, LSM, and regularizer entry points;
- `ecv_stability_blockmodel()` and `ecv_stability_rdpg()`;
- the stable NCV public wrapper(s);
- useful result-comparison and plotting entry points.

Keep internal:

- every `validate_*()` function and equivalent validation-only helper;
- raw functions with `_cpp` suffix and named C++ implementation helpers;
- legacy dotted-name ECV functions such as `ECV.BM` and their internal cores;
- temporary compatibility environments and implementation closures;
- low-level split/evaluation workers used only by one public algorithm;
- `source_all()` and interactive source-order machinery, which package loading
  makes unnecessary;
- low-level seed-offset, formatting, or dispatch helpers with no independent
  user-facing purpose.

Do not infer visibility solely from a leading dot or current name. Use caller
analysis and the dictionary description. Register S3 methods in `NAMESPACE`
without exporting them as ordinary functions.

### 1.3 Dependency audit

Search all R and C++ sources for namespace calls, unqualified calls, S3/S4
requirements, and compiled headers. Classify each dependency as:

- `Imports`: required at installation/runtime;
- `Suggests`: optional plotting, examples, tests, benchmarks, or vignettes;
- `LinkingTo`: headers required for compilation;
- removable: currently listed but unused.

Prefer namespace-qualified calls over attaching packages. Optional packages
must be guarded with `requireNamespace()` and documented. Likely items requiring
special review include `Rcpp`, `RcppArmadillo`, `Matrix`, `RSpectra`, `irlba`,
`future`, `future.apply`, `ggplot2`, `randnet`, `peakRAM`, and the packages in
the current broad `Imports` field.

## Phase 2: Standardize the public argument API

### 2.1 Ordering rule

Use a consistent conceptual order across related public functions:

1. primary data (`A`, `P`, `Z`, or another object being operated on);
2. primary model size or candidate set (`K`, `K_candidates`, `d`,
   `d_candidates`, `tau_candidates`);
3. structural/model choices;
4. partition, fold, and repetition controls;
5. loss and estimator controls;
6. algorithm-specific option lists;
7. computational controls (`ncores`, backend/engine options);
8. reproducibility (`seed`);
9. diagnostics and behavior (`verbose`, `ram_check`, `failure_handling`);
10. memory/return controls such as `retain_intermediates`;
11. `...`, only where forwarding is deliberate and documented.

For the related model-selection APIs, use these canonical openings:

```r
netcrop_blockmodel(A, K_candidates, ...)
netcrop_rdpg(A, d_candidates, ...)
netcrop_lsm(A, d_candidates, ...)
netcrop_tune_regularizer(A, K, tau_candidates, ...)
ecv_stability_blockmodel(A, K_candidates, ...)
ecv_stability_rdpg(A, d_candidates, ...)
ncv_stability_blockmodel(A, K_candidates, ...)
```

Prefer `K_candidates` and `d_candidates` to inconsistent `max_K` and `max_d`
interfaces when the implementation can accept explicit vectors without
changing the published algorithm. If the exact ECV or NCV algorithm requires
the full sequence `1:max_K` or `1:max_d`, keep that mathematical restriction
but expose it consistently and document it; do not pretend arbitrary candidate
vectors are supported. In that case, use consistent names such as `max_K` and
`max_d` in equivalent positions across wrappers.

### 2.2 Compatibility and scope

- Inventory all benchmark calls and internal callers before changing a
  signature.
- Update all internal calls, examples, tests, documentation, and dictionary
  entries atomically.
- Prefer named arguments in examples and internal calls where signatures are
  long.
- Do not add deprecated aliases automatically. The package is at development
  version `0.0.0.9000`, so a clean initial API is preferable unless the owner
  reports external users depending on old names.
- Never change the order of S3 generics' required arguments: retain forms such
  as `print(x, ...)`, `summary(object, ...)`, and `plot(x, ...)`.
- After reordering, test both positional use for the first essential arguments
  and named use for all extended controls.

## Phase 3: Convert the source tree into an R package

1. Place R package source files in `R/`. Preserve meaningful grouping and
   deterministic load behavior. Do not rely on users calling `source()` in a
   particular order.
2. Place C++ source files in `src/`.
3. Add roxygen skeletons before generating namespace artifacts.
4. Generate and retain `RcppExports.R` and `RcppExports.cpp` using the package's
   chosen Rcpp workflow. Register native routines; do not use dynamic symbol
   lookup when avoidable.
5. Keep raw compiled functions internal. Public R wrappers must preserve the
   documented pure-R fallback behavior and `use_cpp` controls.
6. Add any necessary `src/Makevars` and `src/Makevars.win` conservatively and
   portably. Do not hard-code local compiler paths or architecture flags.
7. Replace source-order `exists()` checks that are inappropriate inside an
   installed namespace, but retain meaningful optional-backend checks.
8. Use roxygen2 as the single source of truth for `NAMESPACE` and `.Rd` files.
   Do not hand-edit generated files after generation.
9. Ensure package startup has no interactive prompts, benchmarks, examples,
   worker creation, compilation, or network access.

## Phase 4: Package metadata and repository conventions

### 4.1 `DESCRIPTION`

Use standard R package metadata:

- `Package: netOP`;
- a concise title in title case without a trailing period;
- an informative multi-sentence `Description` explaining network generation,
  estimation, embedding, SONNET, NETCROP, ECV, and NCV;
- `Version: 0.0.0.9000` until the owner chooses a first release;
- `Authors@R` with the confirmed creator and contributor roles;
- `License: MIT + file LICENSE`;
- `URL: https://github.com/sayan-ch/netOP`;
- `BugReports: https://github.com/sayan-ch/netOP/issues`;
- `Encoding: UTF-8`;
- current `Roxygen` and `RoxygenNote` fields;
- a defensible minimum R version if the code uses version-specific features;
- accurate `Imports`, `Suggests`, and `LinkingTo` fields;
- `Config/testthat/edition: 3` if testthat is adopted.

The usual practice is to use the repository URL and issue tracker above, keep
the package in a development version until the first release, and express MIT
as `MIT + file LICENSE`. Use Sayan Chakrabarty as the sole package author and
creator. The coauthors of the NETCROP and SONNET papers must be credited in the
relevant scholarly citations, but must not be listed as package contributors
because all package contributions are Sayan Chakrabarty's.

Use copyright year `2026`, as selected by the package owner.

### 4.2 Standard repository files

Add or update:

- `README.md` with purpose, development installation, a minimal quick start,
  main algorithm families, citations, and support links;
- `NEWS.md` beginning with `netOP 0.0.0.9000`;
- `CONTRIBUTING.md` with naming, dictionary, testing, and documentation rules;
- `CODE_OF_CONDUCT.md` if public outside a small known collaboration;
- `.Rbuildignore` and `.gitignore` for RStudio, check, compiled, and local
  artifacts;
- optional `cran-comments.md` only when preparing for CRAN submission.

Do not add a CRAN badge or claim CRAN availability before a release exists.
The initial README installation path must use GitHub, for example through
`remotes::install_github("sayan-ch/netOP")`. Prepare CRAN-compatible metadata
and checks during development, but treat actual CRAN submission as a later,
separately authorized release phase.

### 4.3 GitHub Actions

Add GitHub Actions R-CMD-check for Linux, macOS, and Windows, with an ordinary
release/devel matrix appropriate to dependency cost. At minimum, every push and
pull request should install dependencies, build documentation, and run
`R CMD check` with errors and warnings treated as failures.

This is usual practice for a public R package because platform-specific Rcpp,
parallel, path, and namespace problems are otherwise easy to miss. Keep the
workflow focused: do not add automatic publishing, releases, or credentials in
the initial package-conversion change.

### 4.4 Vignettes and pkgdown

Required for this package conversion:

- add one lightweight `vignettes/getting-started.Rmd` or `.qmd` covering data
  generation, fitting, model selection, and result inspection;
- keep large benchmarks out of vignette execution;
- make optional-package-dependent sections conditional;
- add pkgdown configuration only after the package installs and checks cleanly;
  do not enable automatic pkgdown deployment without separate authorization.

A README plus reference pages is the normal minimum. A vignette is worthwhile
here because several algorithms have multi-step workflows and nontrivial
outputs. Multiple vignettes and a hosted pkgdown site are useful later, but are
not required for the first installable GitHub version.

## Phase 5: Documentation and examples

### 5.1 Roxygen requirements

Every exported function must document:

- purpose and algorithm context;
- every argument and default;
- return value, including class and important list components;
- sparse/dense behavior and directed/undirected restrictions;
- reproducibility and parallel behavior where relevant;
- optional dependencies and fallback behavior;
- failure behavior and important validation constraints;
- references and related functions;
- at least one useful example unless execution is inherently unsuitable.

Use shared roxygen blocks, `@inheritParams`, and `@family` tags to keep related
interfaces consistent, but ensure the generated help page remains readable on
its own. Add package-level documentation in `R/netOP-package.R`.

### 5.2 Example conversion rules

The commented code at the ends of `x1_sonnet.R` through `x6_other_algos.R` is
the source material, not text to copy mechanically.

- Reduce network sizes, candidate ranges, repetitions, iterations, and worker
  counts so ordinary examples are deterministic and fast.
- Use `set.seed()` or explicit `seed` arguments.
- Prefer `ncores = 1L` in executable examples.
- Remove `system.time()` and `peakRAM::peakRAM()` from ordinary help examples.
- Put authentic large-scale benchmarks in `inst/benchmarks/` as runnable
  scripts, with environment/session information and expected usage notes.
- Use `\donttest{}` only for examples that are useful and reasonably bounded
  but too slow for routine checks.
- Use `\dontrun{}` only when execution truly requires substantial resources,
  unavailable services, or explicit user action. Do not hide broken examples
  inside `\dontrun{}`.
- Test every executable example through `R CMD check`.

### 5.3 Required method citations and wording

NETCROP documentation must cite:

> Chakrabarty, S., Sengupta, S., and Chen, Y. (2026). Network
> Cross-Validation and Model Selection via Subsampling. arXiv:2504.06903.
> https://doi.org/10.48550/arXiv.2504.06903

Include the citation on all primary NETCROP help topics and in the package
README citation section. The paper currently identifies NETCROP as “NETwork
CRoss-Validation using Overlapping Partitions.”

ECV documentation must cite:

> Li, T., Levina, E., and Zhu, J. (2020). Network cross-validation by edge
> sampling. *Biometrika*, 107(2), 257-276.
> https://doi.org/10.1093/biomet/asaa006

The help pages for `ecv_stability_blockmodel()` and
`ecv_stability_rdpg()` must state clearly, without implying authorship of the
underlying ECV method or code:

> netOP provides self-contained wrappers around an ECV implementation derived
> from the CRAN package `randnet`. Installing or using netOP does not require
> `randnet`. The ECV-specific implementation helpers are internal and are not
> part of the netOP public API.

The package owner has confirmed that the ECV code was taken from the CRAN
package `randnet`. CRAN `randnet` version 1.0 is licensed GPL (>= 2). During the
licensing audit, identify the exact upstream source files/functions
corresponding to the included implementation and record whether each portion
was copied or adapted. `randnet` is a development-only reference and must not
be a netOP dependency.

**Public-release blocker:** copied or adapted GPL code must not be distributed
as MIT-only code. Before making netOP public, obtain a defensible licensing
resolution. The expected choices are either (a) license netOP compatibly with
the upstream GPL terms and preserve required notices and attribution, or (b)
replace the copied/adapted ECV code with an independently written
implementation based on the published paper, with provenance documented. Do
not relabel copied GPL code as an independent implementation, and do not make
the repository public until this issue is resolved.

NCV documentation must cite:

> Chen, K. and Lei, J. (2018). Network Cross-Validation for Determining the
> Number of Communities in Network Data. *Journal of the American Statistical
> Association*, 113(521), 241-251.
> https://doi.org/10.1080/01621459.2016.1246365

The primary NCV wrapper documentation must say:

> This is the netOP authors' implementation of the exact algorithm described
> in Chen and Lei (2018). Numerical-stability measures and failsafes were added
> without altering the algorithm in the paper.

Do not broaden that claim to helpers or variants that deviate from the cited
algorithm. Document every stability measure and failsafe at a useful level,
and add regression tests showing that they preserve results on ordinary inputs.

## Phase 6: Testing strategy

### 6.1 Preserve and organize existing tests

Extract relevant checks from the end-of-file blocks into a standard
`tests/testthat/` suite. Do not execute benchmark-sized cases in routine tests.
Separate correctness tests from performance benchmarks.

Add test groups for:

- package installation and loading in a clean session;
- public export allowlist and internal-function denylist;
- absence of raw `_cpp` functions from `netOP::` completion;
- S3 method registration and dispatch;
- argument-order and named-argument behavior;
- ECV rename and absence of the old public name;
- ECV wrapper behavior without `randnet` installed, plus development-only
  equivalence comparisons against CRAN `randnet` when it is locally available;
- NCV ordinary cases and stability/failsafe cases;
- dense and sparse parity where supported;
- pure-R/compiled parity, including boundaries, invalid inputs, and ties;
- deterministic results across supported worker counts when promised;
- every generator, estimator, embedder, and high-level model-selection family;
- examples, README examples, and vignette code.

### 6.2 Package checks

Run, in order:

1. parse and source/load checks in a clean R session;
2. roxygen generation;
3. package installation into a temporary library;
4. testthat suite against the installed package;
5. `R CMD build`;
6. `R CMD check --as-cran` on the built tarball;
7. inspection of warnings, notes, compiled symbols, and undocumented exports;
8. optional `devtools::check()` or `rcmdcheck::rcmdcheck()` as a convenience,
   not as a substitute for the base commands;
9. GitHub Actions on all configured operating systems.

No warning should be waived merely because the code previously passed ad hoc
tests. Document any unavoidable NOTE and its rationale.

## Phase 7: Security, robustness, and licensing review

The implementing model must report findings even if no source change is made.
Inspect at least:

- file-reading and file-writing functions for path handling, unintended
  overwrites, unsafe object deserialization, and unbounded allocations;
- all system, shell, compiler, and subprocess calls for injection risks;
- parallel workers for cleanup, nested oversubscription, error propagation,
  package loading, and restoration of the caller's future plan;
- random seed validation and reproducibility across platforms/workers;
- RAM preflight calculations and denial-of-service-sized allocations;
- C++ index bounds, dimensions, finite-value handling, user interrupts,
  integer overflow, and exception translation;
- dynamic native-symbol lookup and visibility of internal compiled routines;
- dependency provenance and unnecessary high-risk dependencies;
- copied or adapted ECV/randnet code for license compatibility, attribution,
  and compliance with the upstream license;
- absence of `randnet` from all netOP dependency fields and clean-session
  operation of both ECV wrappers without `randnet` installed;
- benchmark code for machine-specific paths, personal data, or credentials;
- generated package artifacts for accidental disclosure of `.Rhistory`, IDE
  state, caches, archives, or local paths.

Potential vulnerabilities or correctness risks must be ranked by severity,
supported by file/line evidence, and separated from stylistic suggestions.

## Phase 8: Final API and documentation review

1. Generate the public function list from the installed namespace.
2. Confirm that useful general helpers requested by the owner remain exported.
3. Confirm that validators, source helpers, raw C++ functions, and ECV internals
   cannot be accessed with `netOP::name`.
4. Confirm both ECV wrappers are visible under their required names.
5. Confirm the old `ecv_stability_bm` name is absent unless the owner later
   requests a deprecated alias.
6. Compare signatures across all related function families and explain any
   remaining inconsistency that is required by generic conventions or the
   underlying algorithm.
7. Verify that every exported function has an indexed help topic and example.
8. Reconcile the final source tree and all signatures with `dictionary.md`.
9. Search the repository for stale dotted identifiers, old directory names,
   old function names, and source-order instructions.

## Phase 9: Commit and handoff

Use small, reviewable commits after the baseline, for example:

1. `Create standard R and Rcpp package structure`
2. `Define and document the public netOP API`
3. `Add package tests and benchmark scripts`
4. `Add README vignette and GitHub checks`
5. `Resolve R CMD check findings`

Do not push implementation commits until local checks pass, unless the owner
explicitly requests incremental remote checkpoints. Never create a release or
publish to CRAN without separate authorization.

The final handoff must include:

- baseline and final commit hashes;
- a concise public API list;
- renamed/reordered functions and migration notes;
- package dependency table and optional dependencies;
- commands run and their outcomes;
- R CMD check status on each tested platform;
- unresolved warnings/notes;
- security, robustness, and licensing findings;
- deliberate algorithmic changes, ideally none;
- links to the README, vignette, tests, and principal documentation files.

## Acceptance criteria

The work is complete only when all of the following hold:

- the original authoritative structure is recoverable from a pushed baseline
  commit;
- the package installs and loads in a clean R session;
- Rcpp code compiles on supported platforms and public wrappers retain working
  R fallbacks where required by `NAMING_CONVENTION.md`;
- `R CMD check --as-cran` completes without errors or warnings, and every NOTE
  is explained;
- the public API contains useful general helpers and high-level algorithms but
  excludes validators, raw C++ exports, and ECV internals;
- ECV exports are exactly the agreed wrappers, including
  `ecv_stability_blockmodel()` and `ecv_stability_rdpg()`;
- NETCROP, ECV, and NCV documentation contains the required citations and
  implementation disclosures;
- related public functions follow a documented, consistent argument-ordering
  scheme, with all callers and examples updated;
- benchmark-derived examples are deterministic, informative, and check-safe;
- `dictionary.md` accurately describes the final tree, function signatures,
  visibility, dependencies, fallbacks, and behavior;
- README, automated checks, and at least one lightweight workflow vignette are
  present unless the owner explicitly opts out;
- no user changes, credentials, local-history files, generated binaries, or
  unexplained algorithm changes have entered the repository.

## Adopted repository and metadata defaults

The package owner has approved the recommendations in sections 4.1 through
4.4. Therefore the implementation must:

- add the README, NEWS, contributing guidance, cross-platform GitHub Actions,
  and one lightweight getting-started vignette described above;
- use `https://github.com/sayan-ch/netOP` for `URL` and
  `https://sayan-ch.github.io` as the author website in `URL`, and use
  `https://github.com/sayan-ch/netOP/issues` for `BugReports`;
- use `License: MIT + file LICENSE`;
- list Sayan Chakrabarty as the sole package author and creator, retain the
  public maintainer email `sayanc@umich.edu`, and add ORCID
  `0000-0003-2372-2076` in the standard `Authors@R` comment field;
- the package owner has explicitly approved publishing `sayanc@umich.edu` in
  package metadata;
- do not list NETCROP or SONNET paper coauthors as package contributors unless
  their role in the package changes later;
- use copyright year `2026`, as approved by the package owner;
- set `Depends: R (>= 4.1.0)` and do not support R 3.x. The current source audit
  found no active syntax that would require R 4.6; base-pipe examples occur
  only in comments. Confirm dependency compatibility and run checks under R
  4.1.x or the closest available compatible environment before finalizing;
- continue development directly on `main`, using reviewable commits and
  pushing only coherent checkpoints that have received proportionate local
  validation;
- publish the development package through GitHub first. Prepare for CRAN, but
  do not submit or claim CRAN availability without separate owner approval;
- add pkgdown configuration only after clean local package checks, and obtain
  separate authorization before enabling deployment or publishing a site.

## Information still needed only at later release stages

No additional owner information is required to begin the package conversion.
Before publishing pkgdown or submitting to CRAN, ask for:

1. confirmation of the desired pkgdown URL and whether GitHub Pages should be
   enabled for this repository;
2. confirmation of the first public release version and release date;
3. explicit authorization to submit to CRAN;
4. any CRAN submission comments, downstream-use information, or additional
   maintainer contact details needed at that time.
