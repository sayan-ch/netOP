# Function Dictionary

This is the living function reference for the `netOP` package. It covers every
top-level named R function, every named C++ helper, and every Rcpp export.
Anonymous and lexically nested closures are implementation details of their
owning function and are not listed as separate APIs.

Update this file whenever a function is added, removed, renamed, moved, or has
its arguments, return value, dependencies, or behavior changed.

## File tree

```text
netOP/
├── DESCRIPTION, NAMESPACE   Package metadata and generated namespace
├── R/
│   ├── 01_basic_helpers.R through 06_estimators.R
│   ├── x0_helpers.R through x6_other_algos.R
│   ├── RcppExports.R         Generated internal compiled wrappers
│   └── netOP-package.R       Package documentation and native registration
├── src/
│   ├── 02_math_helpers.cpp
│   ├── 03_network_utils.cpp
│   ├── 05_embedders.cpp
│   └── RcppExports.cpp       Generated registered interfaces
├── tests/testthat/           Correctness, parity, API, and workflow tests
├── vignettes/                Lightweight getting-started workflow
├── inst/benchmarks/          Authentic non-check benchmark programs
├── inst/API                  Explicit public export allowlist
├── inst/CITATION             Package and method citations
├── inst/COPYRIGHTS           File-level licensing/provenance inventory
├── inst/LICENSE.note         Installed copy of the upstream-code notice
├── README.md, NEWS.md        User-facing package and release documentation
├── CONVENTIONS.md            Mandatory naming and implementation rules
└── dictionary.md             This living function reference
```

## Shared argument conventions

These meanings apply throughout the dictionary unless a function says
otherwise.

| Argument | Meaning |
|---|---|
| `A` | Adjacency matrix or another matrix being analyzed. |
| `P` | Edge-probability matrix, or predictions in loss functions. |
| `n` | Number of nodes. |
| `d` | Requested latent, eigen, or singular dimension. |
| `K` | Number of communities. |
| `g_true` | Ground-truth community labels. |
| `psi` | Node degree-correction parameters. |
| `Z` | Undirected latent-position matrix. |
| `Z_left`, `Z_right` | Source/left and target/right latent-position matrices for directed models. |
| `representation` | `"dense"` returns a base matrix; `"sparse"` returns a general `dgCMatrix`. |
| `directed` | Whether edge direction is meaningful. |
| `self_loops` | Generator Boolean controlling diagonal edges, or utility choice `"ignore"`/`"include"`. |
| `seed` | Optional nonnegative RNG seed. |
| `ncores` | Positive worker count; default is half the detected logical cores with safeguards. |
| `use_cpp` | Try an available compiled accelerator and fall back to R after absence or failure. |
| `lower_clip`, `upper_clip` | Inclusive numerical or probability bounds. |
| `average_degree` | Target expected row degree, not the realized random degree. |
| `average_degree_method` | Defaults to `"naive"`, the historical one-step or iterative approximation; `"calibrated"` targets the final expected degree. |
| `align_with` | Optional reference embedding used for rotation-only Procrustes alignment. |

## Visibility and installed-package behavior

Functions described as internal, validation helpers, splitters, raw Rcpp
wrappers, legacy ECV routines, and comparison dispatch helpers are namespace
internal. S3 methods are registered but are not ordinary exports. All other
user-facing rows are exported and are listed exactly in `inst/API`; generated
help topics contain the installed signatures and argument defaults. Package
loading resolves source relationships, so no user sourcing order or
`sourceCpp()` call is required. Compiled wrappers use registered symbols and
fall back through their documented R-facing functions.

## `01_basic_helpers.R`

Location: [`R/01_basic_helpers.R`](./R/01_basic_helpers.R)

### Functions

The RAM-query, formatting, estimation, and reporting helpers below are internal
implementation details. `uni_mclapply()`, `measure_peak_ram()`, and
`run_simulations()` remain public.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `uni_mclapply()` | Cross-platform `lapply` replacement. Caps workers at `length(X)`, uses `parallel::mclapply` on Unix and `future.apply::future_lapply` on Windows, and normally falls back to sequential execution after backend failure. Strict mode tags every worker result and stops with failed task indices and original messages instead of returning partial/error objects. High-level algorithms can retain one multisession pool across stages and preload packages needed for S4 method dispatch. | `X`; `FUN`; `...`; `ncores`; `force_windows`; `stop_on_error`; `manage_future_plan = TRUE`; `future_packages = NULL`; `suppress_worker_renv_sync_check = TRUE`: temporarily suppress identical renv synchronization messages during worker startup and restore the parent setting. | A list with one value per element of `X`; an empty list for empty `X`. |
| `available_ram()` | Queries currently available system memory through the optional lightweight `ps` package, accommodating platform/version-specific field names. | None. | Available bytes, or `NA_real_` when unavailable. |
| `format_bytes()` | Formats a non-negative byte count in binary units. | `x`: one byte count. | Human-readable string from bytes through PiB, or `"unknown"`. |
| `warn_if_insufficient_ram()` | Reports expected RAM out of currently available RAM. Warns when the estimate exceeds the conservative usable threshold; otherwise prints an informational message. | `estimated_bytes`; `max_fraction = 0.60`; `reserve_bytes = 2 GiB`; `operation`: description used in the report. | Invisibly `TRUE` when within the threshold, `FALSE` when above it, or `NA` when available RAM cannot be determined. |
| `estimate_rspectra_ram()` | Estimates Lanczos/Arnoldi RAM from the matrix dimension, requested rank, Krylov-basis size, and a safety factor. | `n`; `K`; `ncv`; `safety_factor`. | Estimated bytes. |
| `estimate_dense_eigen_ram()` | Estimates working-copy and eigenvector RAM for a dense symmetric eigendecomposition. | `n`; `input_already_counted`; `safety_factor`. | Estimated bytes. |
| `estimate_partial_svd_ram()` | Estimates working RAM for a partial rectangular singular decomposition. | `n`, `p`: dimensions; `K`: computed rank; `ncv`: basis size; `safety_factor`. | Estimated bytes. |
| `estimate_dense_svd_ram()` | Estimates input-copy, singular-vector, and workspace RAM for a dense SVD. | `n`, `p`; `nu`, `nv`; `input_already_counted`; `safety_factor`. | Estimated bytes. |
| `estimate_matrix_product_ram()` | Estimates result allocation and conservative BLAS packing workspace for one dense matrix product. It does not print or warn. | `nrow_left`; `shared_dimension`; `ncol_right`; `safety_factor`. | Estimated bytes. |
| `estimate_spectral_decomp_ram()` | Resolves the audited base/RSpectra/irlba fallback rules and estimates decomposition memory without executing a check. | `n`, `p`, `K`; `method`; `engine`; engine and vector controls; `dense_input`. | List with estimated bytes, resolved engine, computed dimension, and dense-input flag. |
| `report_ram_preflight()` | Reporting utility used only by high-level algorithms. Multiplies the largest per-operation estimate by the counted operations and concurrent cores, then warns when the deliberately overestimated total exceeds available RAM or availability is unknown. | `estimated_bytes`; `operation`; `operation_count`; `ncores`; optional `detail`. | Invisible diagnostic list. |
| `report_ram_formula()` | Reporting utility used only by high-level algorithms. Prints every additive RAM term as `(per-operation estimate x sequential operations) x parallel operations`, totals the terms conservatively, and warns when the total exceeds available RAM or availability is unknown. | `terms`: list of estimates, sequential counts, parallel counts, and labels; `operation`; optional `detail`. | Invisible diagnostic list containing normalized terms, total estimate, and available RAM. |
| `measure_peak_ram()` | Runs an expression or function exactly once through `peakRAM::peakRAM`, while capturing the original result. Use as `measure_peak_ram({ expression })` or `measure_peak_ram(FUN, ...)`. Requires `peakRAM`. | `process`: unevaluated expression or function; `...`: arguments used only when `process` resolves to a function. | List with `result` and one-row `metrics` containing `elapsed_seconds`, `total_ram_used_mib`, and `peak_ram_used_mib`. |
| `run_simulations()` | Repeats `one_simulation`, manages seeds, optional outer-loop parallelism, progress, errors, resource measurements, and resumable RDS results. The outer loop is sequential by default because simulation methods often parallelize internally. Resuming requires a caller-created or otherwise trusted RDS file. | `one_simulation`: function accepting named `simulation`; `nsim`: repetitions; `...`: forwarded inputs; `use_parallel_simulations`: parallelize outer loop; `ncores_outer`: outer workers; `seed` or `seeds`: base or explicit seeds; `results_file`: optional trusted RDS path; `action`: `"replace"`, `"resume"`, or `"archive"`; `retry_failed`: rerun failed records; `continue_on_error`: continue after failure; `measure_resources`: append elapsed time and peak RAM; `show_progress`: print completion messages; `force_windows`: select Windows backend. | List of simulation records; invisibly returned. Each record contains simulation number, seed, success, elapsed time, peak RAM, result, and error. |

## `02_math_helpers.R`

Location: [`02_math_helpers.R`](./R/02_math_helpers.R)

### Losses and scalar/vector transforms

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `modal()` | Finds the most frequent value; the first encountered value breaks frequency ties. | `x`: input vector; `na_rm`: remove missing values first. | One modal value, or a typed `NA` if removal leaves no values. |
| `validate_error_inputs()` | Internal validator shared by vector error metrics. | `x`, `y`: equal-length numeric vectors; `na_rm`: single logical. | `invisible(TRUE)` or an error. |
| `sse()` | Sum of squared errors. | `x`, `y`: equal-length numeric vectors; `na_rm`: remove missing differences; `validate_inputs = TRUE`: disable only after a high-level caller has validated both vectors. | Numeric scalar. |
| `sae()` | Sum of absolute errors. | Same as `sse()`. | Numeric scalar. |
| `mse()` | Mean squared error. | Same as `sse()`. | Numeric scalar. |
| `mae()` | Mean absolute error. | Same as `sse()`. | Numeric scalar. |
| `nmi()` | Computes the requested empirical normalized mutual information `MI(g_1, g_2) / H(g_1, g_2)` from the label contingency table. This joint-entropy normalization differs from marginal-entropy NMI definitions. Uses `entropy::mi.empirical()` and `entropy::entropy.empirical()` when the optional `entropy` package is installed and otherwise uses an equivalent dependency-free empirical plug-in calculation. | `g_1`, `g_2`: equal-length, non-empty atomic label vectors without missing values. | Numeric scalar; `NaN` when the joint entropy is zero. |
| `clip_values()` | Clips a numeric vector or matrix to inclusive bounds while preserving names and dimensions. | `x`: numeric vector/matrix; `lower_clip`, `upper_clip`: scalar bounds. | Clipped object with the shape of `x`. |
| `clip_probabilities()` | Protects logarithms and inverse links by clipping to `[eps, 1 - eps]`. | `P`: numeric probabilities; `eps`: scalar in `(0, 0.5)`. | Clipped object with the shape of `P`. |
| `bin_dev()` | Summed binary cross-entropy/negative log-likelihood. The historical name does not imply the extra factor of two. | `x`: binary response; `y`: probabilities; `epsilon`: clipping level; `na_rm`: missing-value removal; `validate_inputs = TRUE`: permit an audited high-level fast path. | Numeric loss scalar. |
| `bin_dev_mean()` | Mean binary cross-entropy over observed pairs. | Same as `bin_dev()`. | Numeric mean loss. |
| `auc()` | Binary ROC AUC with average ranks for ties. Tries `auc_cpp()` and otherwise uses a base-R rank implementation. | `A`: binary labels; `P`: finite scores; `use_cpp`: accelerator switch; `validate_inputs = TRUE`: permit an audited high-level fast path. | AUC in `[0,1]`, or `NaN` if only one class occurs. |
| `auc_as_loss()` | Converts AUC into a smaller-is-better loss. | Same as `auc()`. | `1 - auc(A, P)`. |
| `softplus()` | Stable `log(1 + exp(x))`. | `x`: numeric object. | Numeric object matching `x`. |
| `sigmoid()` | Stable logistic sigmoid using `stats::plogis`. | `x`: numeric object. | Probabilities matching `x`. |

### Extrema and selection helpers

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `validate_extrema_inputs()` | Internal validator for soft/hard extrema functions. | `x`: numeric vector; `na_rm`: logical; `temperature`: optional positive scalar. | `invisible(TRUE)` or an error. |
| `softmax()` | Normalized exponential weights favoring larger values, with stable shifting and explicit infinity handling. | `x`: scores; `temperature`: positive softness parameter; `na_rm`: missing handling. | Probability vector matching `x`. |
| `softmin()` | Softmax applied to `-x`, favoring smaller values. | Same as `softmax()`. | Probability vector. |
| `hardmax()` | Marks every maximal entry with one and all other observed entries with zero. | `x`: numeric vector; `na_rm`: missing handling. | Numeric 0/1 selector. |
| `hardmin()` | Marks every minimal entry. | Same as `hardmax()`. | Numeric 0/1 selector. |
| `which_softmax()` | Samples one index according to softmax probabilities. | `x`, `temperature`, `na_rm`: as in `softmax()`. | One integer index or `NA_integer_`. |
| `which_softmin()` | Samples one index according to softmin probabilities. | `x`, `temperature`, `na_rm`: as in `softmin()`. | One integer index or `NA_integer_`. |
| `which_hardmax()` | Returns every index tied for the maximum. | `x`, `na_rm`: as in `hardmax()`. | Integer index vector or `NA_integer_`. |
| `which_hardmin()` | Returns every index tied for the minimum. | `x`, `na_rm`: as in `hardmin()`. | Integer index vector or `NA_integer_`. |

### Matrix alignment and graph matrices

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `outer_add()` | Forms the additive outer combination with `y` as rows and `x` as columns. Uses `outer_add_cpp()` when possible. | `x`, `y`: numeric vectors; `operator`: currently only `"+"`; `use_cpp`: accelerator switch. | `length(y)` by `length(x)` numeric matrix. |
| `procrustes()` | Orthogonally aligns `X` to `X_star`, optionally translating, dilating, and reporting residual SSE. Pads `X` with zero columns when needed. | `X`, `X_star`: finite matrices with equal rows; `translate`: add translation; `dilate`: estimate scale; `sumsq`: return SSE; `use_cpp`: accelerate translation-only case; `validate_inputs = TRUE`: permit a validated high-level fast path. | List with `X_new`, `rotation`, and optional `translation`, `scale`, and `sse`. |
| `matrix_density()` | Fraction of matrix entries that are nonzero. Used by sparsity scaling in decompositions. | `A`: finite matrix-like object with positive dimensions. | Numeric scalar in `[0,1]`. |
| `is_symmetric_matrix()` | Numerical symmetry check supporting base and `Matrix` classes without attaching `Matrix`; base-matrix dimname differences are ignored. | `A`: matrix-like object. | Single logical. |
| `graph_laplacian()` | Builds `D - A` or the symmetric normalized Laplacian `I - D^{-1/2} A D^{-1/2}`; isolated-node diagonal entries remain zero when unregularized. Positive `tau` first adds `tau * mean(deg) / n` to every adjacency entry. | `A`: finite symmetric adjacency; `normalized`: choose normalized form; `tau`: dimensionless non-negative regularization strength. | Dense or sparse Laplacian; positive all-pairs regularization can make sparse input dense. |

### Decompositions and reconstruction

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `eig_decomp()` | Selected eigenpairs of a finite symmetric matrix, with safe fallbacks from RSpectra to base R. `irlba` is redirected because an SVD is not a signed eigendecomposition. This low-level helper performs no RAM preflight. Input validation remains enabled by default; an already-validated high-level algorithm may disable it to avoid repeated sparse finite-value and symmetry scans. | `A`: symmetric matrix; `d`: requested pairs; `only_values`: omit vectors; `scale_by`: `"none"`, `"dimension"`, or `"sparsity"`; `use_laplacian`: decompose graph Laplacian; `engine`: `"rspectra"`, `"base"`, or compatibility value `"irlba"`; `force_engine`: disable fallback; `order_by`: magnitude or signed value; `safe_d_multiplier`: extra partial components computed before truncation; `validate_inputs = TRUE`: perform matrix validation. | List with `values` and, unless omitted, `vectors`. |
| `singular_decomp()` | Selected singular values/vectors with RSpectra, irlba, and base fallbacks. This low-level helper performs no RAM preflight. Input validation remains enabled by default; an already-validated high-level algorithm may disable it to avoid repeated sparse finite-value scans. | `A`: finite matrix; `d`: requested values; `only_values`: omit vectors; `nu`, `nv`: left/right vector counts; `scale_by`, `use_laplacian`, `engine`, `force_engine`, `order_by`, `safe_d_multiplier`: decomposition controls; `validate_inputs = TRUE`: perform matrix validation. | List with `values` and optional `u` and `v`. |
| `truncated_svd_reconstruct()` | Rank-`d` SVD reconstruction using column scaling and `tcrossprod`, then clipping. Use infinite bounds for an unclipped result. | `A`: matrix; `d`: rank; clipping bounds; `engine`, `force_engine`, `safe_d_multiplier`: passed to `singular_decomp()`. | Reconstructed matrix with `dimnames(A)`. |
| `usvt()` | Historical USVT variant using rank `ceiling(n^(1/3))`, followed by truncated-SVD reconstruction and clipping. | `A`: square matrix; clipping bounds; decomposition engine controls. | Estimated probability matrix `P_hat`. |

## `02_math_helpers.cpp`

Location: [`02_math_helpers.cpp`](./src/02_math_helpers.cpp)

Compiled at package installation and registered through generated Rcpp
interfaces. All exports have R fallbacks in `R/02_math_helpers.R`.

| Function | Visibility | Description and arguments | Returns/use |
|---|---|---|---|
| `auc_cpp()` | Rcpp export | `group`: binary numeric labels; `predictions`: finite scores. Sorts compact label/score pairs and handles ties by average rank. | Numeric AUC; accelerator for `auc()`. |
| `outer_add_cpp()` | Rcpp export | `x`, `y`: numeric vectors. Allocates an uninitialized output and fills `y[row] + x[column]`. | Numeric matrix; accelerator for `outer_add()`. |
| `procrustes_translated_cpp()` | Rcpp export | `X`, `X_star`: equal finite dense matrices. Uses Eigen SVD without an `n` by `n` centering matrix. | List with `X_new`, `rotation`, and `translation`; accelerator for `procrustes()`. |

## `03_network_utils.R`

Location: [`03_network_utils.R`](./R/03_network_utils.R)

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `validate_adjacency()` | Internal square, numeric, finite adjacency validator supporting base and `Matrix` objects. | `A`: adjacency matrix. | `invisible(TRUE)` or an error. |
| `adjacency_neighbors()` | Builds an unweighted adjacency list without densifying sparse input. | `A`; `directed`: retain orientation or union both directions; `self_loops`: ignore/include diagonal neighbors. | List of integer neighbor vectors. |
| `connected_components()` | Breadth-first undirected/weak connected components. Edge weights are ignored. Tries dense/sparse C++ exports before the R implementation. | `A`; `self_loops`; `use_cpp`. | Named list of one-based node-index vectors. |
| `largest_connected_component()` | Extracts the first largest weak component and its induced adjacency matrix. | `A`; `self_loops`; `use_cpp`; `sort_nodes`: sort returned indices. | List with `nodes`, `size`, and `submatrix`. |
| `adjacency_weighted_neighbors()` | Builds node and weight adjacency lists for weighted paths without densifying sparse input. For undirected reciprocal weights, the smaller nonzero weight is used. | `A`; `directed`; `self_loops`. | List with parallel `nodes` and `weights` lists. |
| `shortest_path_distances()` | All-pairs shortest-path distance matrix. Uses BFS when unweighted and Dijkstra when weighted; has dense/sparse compiled paths and dependency-light R fallbacks. | `A`; `directed`; `weighted`: use nonnegative entries as distances; `self_loops`; `use_cpp`. | `n` by `n` matrix with zero diagonal and `Inf` for unreachable pairs. |

## `03_network_utils.cpp`

Location: [`03_network_utils.cpp`](./src/03_network_utils.cpp)

Compiled at package installation. The R-facing utilities validate inputs and
fall back to R when these registered internal routines are absent or fail.

| Function | Visibility | Description and arguments | Returns/use |
|---|---|---|---|
| `is_nonzero()` | Internal C++ | `value`: numeric adjacency entry; recognizes finite/nonmissing nonzero edge presence. | Boolean used by adjacency builders. |
| `undirected_weight()` | Internal C++ | `outgoing`, `incoming`: reciprocal edge entries. | The existing value, or smaller reciprocal nonzero weight. |
| `build_dense_adjacency()` | Internal C++ | Dense `A`; `directed`; `weighted`; `include_self_loops`. Validates weights and builds adjacency lists. | `adjacency_list`. |
| `build_sparse_adjacency()` | Internal C++ | Sparse S4 `A`; direction, weight, and loop controls. Deduplicates neighbors using minimum weight. | `adjacency_list`. |
| `connected_components_from_adjacency()` | Internal C++ | `neighbors`: undirected adjacency list. | R list of one-based BFS components. |
| `unweighted_distances_from_adjacency()` | Internal C++ | `neighbors`: adjacency list. | All-pairs BFS distance matrix. |
| `weighted_distances_from_adjacency()` | Internal C++ | `neighbors`: nonnegative weighted adjacency list. | All-pairs Dijkstra distance matrix. |
| `shortest_path_distances_from_adjacency()` | Internal C++ | `neighbors`; `weighted`: selects BFS or Dijkstra. | All-pairs distance matrix. |
| `connected_components_dense_cpp()` | Rcpp export | Dense `A`; `include_self_loops`. Always treats the network as undirected. | Accelerator for `connected_components()`. |
| `connected_components_sparse_cpp()` | Rcpp export | Sparse `A`; `include_self_loops`. | Sparse accelerator for `connected_components()`. |
| `shortest_path_distances_dense_cpp()` | Rcpp export | Dense `A`; `directed`; `weighted`; `include_self_loops`. | Accelerator for `shortest_path_distances()`. |
| `shortest_path_distances_sparse_cpp()` | Rcpp export | Sparse `A`; `directed`; `weighted`; `include_self_loops`. | Sparse accelerator for `shortest_path_distances()`. |

## `04_network_generator.R`

Location: [`04_network_generator.R`](./R/04_network_generator.R)

All public network generators return `A` directly. Generated truth is stored in
the lightweight `generator_parameters` attribute and retrieved with
`get_generator_parameters()`. Full `n` by `n` probability matrices are not
attached because that would defeat sparse-memory savings.

### Validation, metadata, and execution helpers

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `validate_generator_logical()` | Internal validator for single logical arguments. | `value`: candidate; `name`: argument name used in errors. | `invisible(TRUE)` or an error. |
| `validate_generator_count()` | Internal validator for integer-like counts. | `value`; `name`; `allow_zero`: permit zero. | Integer scalar. |
| `validate_generator_seed()` | Validates and normalizes optional seeds to R's integer range. | `seed`: `NULL` or nonnegative integer-like scalar. | `NULL` or integer seed. |
| `offset_generator_seed()` | Derives deterministic component-specific seeds without overflow. | `seed`; `offset`: numeric offset. | `NULL` or integer seed. |
| `validate_generator_ncores()` | Applies the standard positive-integer worker validation. | `ncores`. | Integer worker count. |
| `generator_lapply()` | Internal generator loop. Uses `uni_mclapply()` when available and requested; otherwise uses sequential `lapply`. | `X`, `FUN`, `...`, `ncores`. | List of loop results. |
| `set_generator_parameters()` | Attaches small generator truth/metadata without changing matrix class. | `A`; `parameters`: named list. | `A` with `generator_parameters` attribute. |
| `get_generator_parameters()` | Retrieves attached generator truth. | `A`: generated adjacency matrix. | Parameter list or `NULL`. |

### Network edge-list I/O

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `validate_network_for_io()` | Internal square, finite numeric matrix validator for I/O. | `A`. | `invisible(TRUE)` or an error. |
| `is_symmetric_network_matrix()` | Namespace-safe symmetry check for dense and sparse matrices. | `A`; `tolerance`: numerical comparison tolerance. | Single logical. |
| `resolve_network_format()` | Resolves `"auto"`, `"csv"`, or `"tsv"` from the explicit choice/file extension. | `file`; `format`. | `"csv"` or `"tsv"`. |
| `write_network()` | Writes an edge list with node pairs in the first two columns and optional weights in the third. Metadata preserves isolated nodes. Undirected files store one upper/lower triangle including loops. Existing paths are rejected unless replacement is explicit. | `A`; `file`; `directed`: explicit or inferred; `weighted`: explicit or inferred; `triangle`: `"upper"`/`"lower"`; `format`; `include_header`; `overwrite = FALSE`; `tolerance`. | Input file path, invisibly. |
| `read_network()` | Reads two- or three-column edge lists into dense/general sparse adjacency. Undirected matrices are formed with `A + t(A)` and loops added once. Rejects duplicate pairs. | `file`; `representation`; `directed`, `weighted`: explicit, metadata-derived, or default; `n`: optional dimension for isolates; `format`; `has_header`. | Dense matrix or general `dgCMatrix`. |

### Probability and average-degree helpers

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `validate_probability_matrix()` | Internal finite `[0,1]`, square, and optional symmetry validator. | `P`; `directed`; `tolerance`. | `invisible(TRUE)` or an error. |
| `clip_generator_probabilities()` | Internal probability clipping independent of math-helper source order. | `P`; `lower_clip`; `upper_clip`. | Clipped probability matrix. |
| `scale_to_average_degree()` | Calibrates a global probability multiplier by bisection against the final clipped, loop-adjusted expected degree. | `P`; `average_degree`; `self_loops`; clipping bounds; `tolerance`; `max_iterations`. | List with calibrated `P` and `multiplier`. |
| `scale_to_average_degree_naive()` | Historical one-step scaling: target divided by `mean(rowSums(P))`, followed by clipping and diagonal removal. | `P`; `average_degree`; `self_loops`; clipping bounds. | List with scaled `P` and `multiplier`; final degree may differ from target. |
| `apply_average_degree_scaling()` | Dispatches to naive or calibrated probability scaling. | `P`; `average_degree`; `average_degree_method = "naive"`; loop and clipping controls. | List with `P` and `multiplier`. |

### General, ER, RDPG, SBM, and DC-SBM generators

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `generate_adjacency()` | Samples independent Bernoulli edges from `P`. Undirected sampling draws one triangle and constructs a general symmetric matrix with `A + t(A)`; row seeds make results independent of worker count. Full probability validation remains enabled for direct matrix calls. A scalar `P` plus `n` avoids allocating a constant `n` by `n` probability matrix. | `P`: probability matrix, or scalar with `n`; `representation`; `directed`; `self_loops`; `seed`; `ncores`; `validate_inputs = TRUE`; `n = NULL`: dimension for scalar `P`. | Binary adjacency matrix. |
| `generate_er()` | Erdős–Rényi generator parameterized by exactly one of edge probability or expected average degree. It samples from a scalar probability and does not allocate a redundant dense probability matrix. | `n`; `p` or `average_degree`; representation, direction, loop, seed, and worker controls. | Adjacency `A`; metadata contains model, `p`, and expected degree. |
| `generate_latent_positions()` | Validates a supplied `Z` or generates an `n` by `d` matrix from a user-selected distribution. Default uniform bounds keep typical dot products controlled. | `n`, `d`: required when `Z` missing; `Z`: optional direct matrix; `latent_distribution`; named `latent_args`; `seed`. | Numeric latent-position matrix. |
| `generate_rdpg()` | Generates an undirected RDPG from `Z` or a directed generalized RDPG from `Z_left`, `Z_right`. The default legacy path generates missing positions from `Uniform(0,1)`, forms the raw Gram/cross-position matrix, divides it by its positive maximum, and only then multiplies by `sparsity_multiplier`, reproducing `P = rho * P0 / max(P0)`. Custom latent distributions/arguments and unscaled probabilities remain supported. | Dimensions and latent inputs; `latent_distribution = stats::runif`; `latent_args = NULL`, interpreted as `min = 0`, `max = 1` on the default RDPG path; `sparsity_multiplier = 1`; `scale_P = TRUE`; `average_degree`; `average_degree_method = "naive"`; clipping, representation, direction, loop, seed, and worker controls. | Adjacency `A`; metadata contains latent truth and degree-scaling details. |
| `generate_community_labels()` | Validates supplied labels and maps them to `1:K`, or samples labels from community probabilities. When neither labels nor probabilities are supplied, it uses the historical uniform default `rep(1 / K, K)`. Original label levels are retained after mapping. | `n`, `K`; `g_true`; `community_probabilities = NULL`; `seed`. | Integer label vector. |
| `generate_inverse_beta_degree_parameters()` | Generates the historical DCBM default `psi = 1 / Beta(4, 1)` with configurable positive beta-shape parameters. RNG seeding remains the responsibility of the owning generator/helper. | `n`: nonnegative count; `shape_1 = 4`; `shape_2 = 1`. | Length-`n` positive numeric degree-parameter vector. |
| `generate_degree_parameters()` | Validates direct `psi` or generates nonnegative degree parameters from a distribution; optionally normalizes within communities. Its default distribution is `generate_inverse_beta_degree_parameters`, reproducing `1 / Beta(4, 1)`; explicitly selecting `stats::runif` retains its existing `[0.5, 1]` default arguments. | `n`; `psi`; `degree_distribution = generate_inverse_beta_degree_parameters`; named `degree_args`; `degree_scale`: `"none"`, `"max_by_community"`, or `"mean_one_by_community"`; `g_true`; `seed`. | Numeric `psi` vector. |
| `generate_dcbm()` | Degree-corrected SBM using direct/generated labels and degree parameters. Missing labels default to uniform community sampling and missing `psi` defaults to `1 / Beta(4, 1)`. Uses direct `P_block` or constructs within/between probabilities from `alpha` and `beta`; all supplied overrides and degree-scaling choices remain supported. | `n`, `K`, labels/proportions; `P_block` or `alpha`, `beta`; `psi`; `degree_distribution = generate_inverse_beta_degree_parameters`; degree arguments/scaling; `sparsity_multiplier`; `average_degree`; `average_degree_method = "naive"`; clipping, representation, direction, loop, seed, and worker controls. | Adjacency `A`; metadata contains `g_true`, `psi`, `P_block`, and scaling details. |
| `generate_sbm()` | Ordinary SBM wrapper around `generate_dcbm()` with `psi = 1`. | Same block, label, sparsity, `average_degree`, `average_degree_method = "naive"`, clipping, representation, direction, loop, seed, and worker controls, excluding degree-correction inputs. | Adjacency `A` with model metadata changed to `"sbm"`. |

### Latent-space-model helpers and generator

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `generate_truncated_normal()` | Dependency-free inverse-CDF sampler for a standard normal truncated to finite bounds. | `n`: draws; `lower_bound`, `upper_bound`. | Numeric vector. |
| `normalize_lsm_positions()` | Centers columns and applies the Ma/Ma Frobenius-based latent scaling. Computes the exact identity `||ZZ' ||_F = ||Z'Z||_F` through the smaller `d` by `d` cross-product instead of materializing an `n` by `n` Gram matrix. | `Z`: finite nondegenerate latent matrix. | Centered, normalized `Z`. |
| `generate_lsm_positions()` | Generates mixture-component means and adds truncated-normal latent noise. | `n`, `d`, `K`, `g_true`; mean and noise bounds; `seed`. | `n` by `d` latent matrix. |
| `generate_lsm_alpha()` | Validates supplied scalar/vector intercepts or generates negative node effects. | `n`; `alpha`: `NULL`, scalar, or length `n`; `seed`. | Length-`n` numeric intercept vector. |
| `scale_lsm_to_average_degree()` | Adjusts the LSM logit intercept. Naive mode reproduces repeated `-log(current/target)` updates before diagonal removal; calibrated mode bisects against final probabilities. | `theta`; `average_degree`; `average_degree_method = "naive"`; `self_loops`; clipping bounds; `naive_iterations`; `tolerance`; `max_iterations`. | List with probability matrix `P` and `intercept_shift`. |
| `generate_lsm()` | Squared-distance latent-space generator. For an undirected graph it constructs `theta_ij = (alpha_i + alpha_j) / 2 - ||Z_i - Z_j||^2` by subtracting each squared row norm from its half-intercept and adding twice the latent inner product. A scalar supplied `alpha` therefore gives exactly `logit(P_ij) = alpha - ||Z_i - Z_j||^2`; the default `alpha_i = -runif(n, 1, 3) / 2` gives its node-intercept generalization. Missing positions use uniform `[-1,1]` community means plus standard-normal noise truncated to `[-2,2]`; generated or supplied positions are centered and Frobenius-normalized by default. Optional naive average-degree targeting performs ten global-intercept updates before diagonal removal. | `n`, `d`, `K = 1`, `g_true`, community probabilities; `alpha = NULL`: generate node-specific values, or supply a scalar/vector; latent inputs defaulting to `NULL`; `normalize_Z = TRUE`; `distance_adjustment = TRUE`; mean bounds `[-1,1]`; noise bounds `[-2,2]`; `average_degree`; `average_degree_method = "naive"`; `naive_iterations = 10L`; clipping, representation, direction, loop, seed, and worker controls. | Adjacency `A`; metadata includes normalized positions, original half-intercept inputs and distance-adjusted intercepts, labels, actual/target degree, method, and global intercept shift. |

## `05_embedders.R`

Location: [`05_embedders.R`](./R/05_embedders.R)

The R-facing PGD function works without the compiled routine and automatically
falls back after absence or runtime failure.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `ase()` | Adjacency spectral embedding. Its optional high-level preflight prints an explicit additive formula for decomposition memory and the largest explicit multiplication estimate times the number of reconstruction/alignment products. | `A`: finite square adjacency; `d`; `directed`; `reconstruct`; `align_with`; `ram_check = FALSE`. | Embedding(s), spectral values, optional alignment, and optional `A_hat`. |
| `spectral_cluster()` | Spectral clustering for undirected, non-negative networks, or direct clustering of a supplied representation. Its optional high-level preflight prints an explicit additive formula for decomposition, counted Laplacian operations, dense regularization, and PAM dissimilarities. Laplacian input is always estimated as dense even when stored and processed sparsely. | `A`; optional `U`; `K`; Laplacian, regularization, zero-degree, normalization, decomposition, and clustering controls; `ram_check = FALSE`; `validate_inputs = TRUE`: permit validated high-level callers to skip repeated scans. Legacy `spectral_options$ram_check` is promoted to this high-level argument. | List with full `g_hat`, clustered `U_hat`, spectral values, zero-degree and retained indices, clustering fit, and selected engines. |
| `lsm_pgd()` | Fits an undirected binary logistic latent-space model by projected gradient ascent. Its optional high-level preflight prints an explicit additive formula for spectral initialization, approximately six dense `n` by `n` working matrices, and the largest multiplication estimate times all `2 * niter + 2` explicit products. | `A`: symmetric binary zero-diagonal adjacency; `d`; `step_size = 0.3`; `niter = 100L`; `trace = FALSE`; optional initial values; `epsilon`; `use_cpp`; `ram_check = FALSE`. | List with `Z_hat`, `alpha_hat`, zero-diagonal `P_hat`, objective history, and resolved step sizes. |

## `05_embedders.cpp`

Location: [`05_embedders.cpp`](./src/05_embedders.cpp)

Compiled at package installation and registered through generated Rcpp
interfaces.

| Function | Visibility | Description and arguments | Returns/use |
|---|---|---|---|
| `stable_sigmoid()` | Internal C++ | `value`: scalar logit; uses sign-specific evaluation to avoid overflow. | Stable logistic probability. |
| `stable_softplus()` | Internal C++ | `value`: scalar logit. | Stable `log(1 + exp(value))`. |
| `lsm_log_likelihood_cpp()` | Internal C++ | `A`, `theta`: equal square matrices. Sums Bernoulli log likelihood only over `i < j`. | Undirected off-diagonal log-likelihood scalar. |
| `validate_lsm_pgd_cpp_inputs()` | Internal C++ | `A`, `Z`, `alpha`, the two step sizes, and `niter`. Checks dimensions, finiteness, symmetry, binary entries, zero diagonal, and positive controls. | Returns void or raises an R error. |
| `lsm_pgd_cpp()` | Rcpp export | Dense `A`; initial `Z`, `alpha`; `step_size_Z`, `step_size_alpha`; `niter`; `trace`. Uses stable likelihoods, excludes diagonal residuals, centers `Z` after every update, and checks interrupts. | Accelerator result with `Z_hat`, `alpha_hat`, `P_hat`, and `objective`; consumed by `lsm_pgd()`. |

## `06_estimators.R`

Location: [`06_estimators.R`](./R/06_estimators.R)

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `estimate_sbm()` | Estimates SBM block means from a full adjacency or a common rectangular NCV layout. The full undirected path computes only one block triangle. Rectangular denominators count the dyads actually present, genuine self-pairs are identified from original node indices, and undirected sufficient statistics pool both available orientations. | `A`; `g`; `K`; optional `fold_nodes`; `directed`; `self_loops`; `validate_inputs = TRUE`: skip repeated full matrix/label scans only after high-level validation. | Named `K` by `K` block-mean matrix `B_hat`. |
| `estimate_dcbm()` | Unified plug-in/spectral estimator of DCBM block and degree parameters for full networks and four common NCV matrix layouts. Supplied norms avoid recomputation. | `A`; `g`, `K`; `method`; `fold_nodes`; `row_norm`; `psi_omit`; `stabilizer`; `spectral_engine`, `spectral_options`; `validate_inputs = TRUE`: permit an audited high-level fast path. | List with named block matrix `B_hat` and node degree parameters `psi_hat`. |
| `estimate_sbm_P_hat()` | Calls `estimate_sbm()` and expands its block estimates into node-pair probabilities. Returns all-node probabilities normally and fold-by-fold probabilities in NCV mode. | Estimator arguments `A`, `g`, `K`, `fold_nodes`, `directed`, `self_loops`; `lower_clip`, `upper_clip`: finite output bounds. | Square probability matrix `P_hat` for the target nodes. |
| `estimate_dcbm_P_hat()` | Calls `estimate_dcbm()` and reconstructs `P_hat = tcrossprod(psi_hat) * B_hat[g, g]`. Returns all nodes normally, retained trailing nodes after `psi_omit`, or fold nodes in NCV mode. | All `estimate_dcbm()` arguments; `self_loops`: retain/zero the output diagonal; `lower_clip`, `upper_clip`: finite output bounds. | Square probability matrix `P_hat` for the estimated nodes. |
| `label_match_greedy()` | Relabels estimated communities to agree with a standard labeling. Identity and binary cases use exact shortcuts. The default deterministic greedy method repeatedly takes the largest remaining overlap; the optional dependency-free Hungarian method finds a globally optimal assignment. The audited mapping direction is source label to standard label. | `match_this`: labels to transform; `standard`: target labels; `K`: label range; `algorithm`: `"greedy"` (default) or `"hungarian"`; `return_mapping`: include assignment diagnostics. | List with `matched_labels` and proportional `mismatch_rate`; when requested, also source-to-standard `mapping`, achieved `agreement`, and overlap matrix. |
| `label_match_brute_force()` | Checks all `K!` source-to-standard label permutations, generated recursively without storing a factorial-size permutation matrix. For `K > 8`, interactive use prompts before proceeding; non-interactive use stops unless `confirm_large = TRUE`. Prefer `label_match_greedy()` for large `K`. | `match_this`, `standard`, `K`: matching inputs; `return_mapping`: include diagnostics; `confirm_large`: `NULL` prompts interactively, `TRUE` explicitly continues, and `FALSE` cancels when `K > 8`. | List with `matched_labels` and proportional `mismatch_rate`; when requested, also mapping, agreement, overlap table, and number of permutations evaluated. |

## `x1_sonnet.R`

Location: [`x1_sonnet.R`](./R/x1_sonnet.R)

Helper relationships are resolved inside the installed namespace.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `.sonnet_fit()` | Internal engine for both SONNET variants. When either partition parameter is `NULL`, it calls `sonnet_param_select()` with `n`, `K`, `theta = 3`, `q = 0.005` or the smallest feasible value, and the algorithm's `ncores`; a supplied non-`NULL` partner is preserved. When both parameters are supplied, they take priority and selection options are ignored. The remaining algorithm is unchanged. | Common SONNET arguments with `num_subnetworks = NULL`, `overlap_size = NULL`, `extra_nrep = 0L`, `verbose = TRUE`, `ram_check = FALSE`, internal `share_overlap`, and `parameter_select_options = list()` forwarded only during automatic selection; `n`, `K`, and `ncores` cannot be overridden; named `...` is forwarded to `spectral_cluster()`. | Full classed `sonnet` object. |
| `sonnet_shared_overlap()` | Explicit entry point for the historical shared-overlap SONNET variant. | `...`: arguments accepted by the SONNET engine. | Classed `sonnet` object with `share_overlap = TRUE`. |
| `sonnet_independent_overlap()` | Explicit entry point for independent overlaps. Every repetition yields all `n` labels and every extra labeling is globally matched to repetition zero before voting. | `...`: arguments accepted by the SONNET engine. | Classed `sonnet` object with `share_overlap = FALSE`. |
| `sonnet()` | Wrapper selecting the shared or independent SONNET flow. The default uses independent overlaps and automatically selects missing partition parameters through `sonnet_param_select()`. Explicitly supplied `num_subnetworks` and `overlap_size` always take priority, irrespective of `parameter_select_options`. Laplacian conversion and regularization are applied within each subnetwork by default; `regularize_subnetworks = FALSE` transforms the full network once before extracting subnetworks. Uses `spectral_cluster()` and strict `uni_mclapply()` worker-error propagation. | `A`; `K`; `num_subnetworks = NULL`; `overlap_size = NULL`; `extra_nrep = 0L`; `ncores`; `seed`; `matching_method`; `confirm_large`; `verbose = TRUE`; `force_windows = FALSE`; `ram_check = FALSE`; `share_overlap = FALSE`; `laplacian = FALSE`; `regularize_tau = 0`; `regularize_subnetworks = TRUE`; `parameter_select_options = list()`; named `...` forwarded to `spectral_cluster()`. | Classed `sonnet` object containing final labels, membership matrices, labels and alignment diagnostics, splits, resolved parameters, worker counts, timings, and call. |
| `print.sonnet()` | Compact print method for a fitted SONNET object. | `x`: `sonnet` object; `...`: method compatibility. | Input object invisibly. |
| `summary.sonnet()` | Builds SONNET diagnostics including community sizes, overlap augmentation, membership certainty, cores, and timings. | `object`: `sonnet` object; `...`: method compatibility. | Classed `summary.sonnet` object. |
| `print.summary.sonnet()` | Prints a SONNET summary. | `x`: `summary.sonnet`; `...`: method compatibility. | Input summary invisibly. |

## `x0_helpers.R`

Location: [`x0_helpers.R`](./R/x0_helpers.R)

These partition helpers are loaded into the installed namespace without a
user-visible source-order requirement.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `op_splitter()` | Repeatedly samples a common overlap and partitions the remaining nodes. When `s * m < n - o`, the default adds every remainder node to the common overlap; alternatively, extras are balanced randomly across the non-overlap pieces or ignored with a warning. Seed offsets are deterministic and overflow-safe. | `n`: nodes; `s`: splits per repetition; `o`: initial overlap size; `n_repetitions`: repeated partitions; `m`: base non-overlap nodes per split; `remainder_handling`: `"augment_overlap"` (default), `"distribute_randomly"`, or discouraged `"ignore"`; `seed`: optional non-negative integer seed. | List with repeated `splits`, effective overlap size `o`, `remainder_count`, and selected `remainder_handling`. |
| `sonnet_param_select()` | Uses a dependency-free integer grid search to maximize `1 - ((n - o) / n)^2 * ((s - 1) / s)` subject to `s * ((o + (n - o) / s) / n)^theta <= q * ncores`. It discards pairs violating `o <= n - s` or exact divisibility of `n - o` by `s`. If no pair meets the requested `q`, it uses the smallest feasible per-core runtime fraction and returns the corresponding pair. The default lower overlap is `min(K^3, n - s_lower)` and the per-candidate default upper overlap is `n - s`. Fractional bounds are normalized to integers and ties resolve in ascending grid order. | `n`: integer from 3 through R's maximum integer; `K = 5`; `theta = 3`; `q = 0.005`: requested per-core runtime fraction; `ncores`: standard positive-integer half-core default; `s_lower = 2`, `s_upper = 50`; `o_lower = NULL`: use the `K^3` default; `o_upper = NULL`: apply `n - s` separately for each `s`. | List with optimal `num_subnetworks`, `overlap_size`, exact integer `m`, objective, percentage `data_use`, aggregate `computational_fraction`, `parallel_runtime_fraction`, `K`, `theta`, effective `q`, `q_requested`, `ncores`, candidate counts, and normalized bounds. Errors only when the bounds contain no structurally admissible pair. |
| `netcrop_param_select()` | Computes NETCROP overlap proportions and subnetwork counts from the supplied formulas, with numerical tolerance before `ceiling()` to stabilize exact integer boundaries. Values of `o_range` interpolate between `1 - sqrt(2 * test_prop)` and `1 - sqrt(test_prop)`; endpoint values equal to one are replaced by `0.8` with a warning to avoid the singular upper endpoint. When `n` is provided, the initial overlap is rounded up and any division remainder augments it, ensuring `(n - overlap_size)` is divisible by `num_subnetworks` with a positive piece size. | `test_prop = 0.02`: finite fraction in `(0, 0.5]`; `n = NULL`: optional integer network size of at least 3; `o_range = c(0, 0.8)`: non-empty numeric vector in `[0, 1]`. | Data frame with `test_prop`, `overlap_proportion`, and `num_subnetworks`; when `n` is supplied, also `n` and adjusted `overlap_size`. |
| `sonnet_splitter()` | SONNET-specific splitter with one overlap shared across all repetitions. Remainders always augment that overlap. First-repetition subnetworks include the overlap; extra repetitions freshly permute and split only the non-overlap nodes, preserving the original SONNET flow. | `n`; `num_subnetworks`; `overlap_size`; `extra_nrep`; `m`; `seed`. | List with per-repetition `subnetworks`, overlap/non-overlap indices, effective/requested overlap, remainder, and piece size. |
| `sonnet_splitter_independent()` | Samples a separate augmented overlap and complete overlap-plus-partition collection for every repetition. Unlike `sonnet_splitter()`, each repetition covers all `n` nodes. | `n`; `num_subnetworks`; `overlap_size`; `extra_nrep`; `m`; `seed`. | List with per-repetition subnetworks, lists of overlap/non-overlap indices, effective/requested overlap, remainder, and piece size. |
| `netcrop_splitter()` | Shared NETCROP splitter that generates a fresh overlap for every repetition and partitions all remaining nodes into equal non-overlap pieces. Remainders always augment the overlap, so no nodes are dropped and the effective overlap is explicitly tracked. | `n`; `num_subnetworks >= 2`; `overlap_size`; `nrep = 1`; `seed`. | Subnetworks, repetition-specific overlap nodes and non-overlap pieces, requested/effective overlap, remainder count, piece size, and repetition metadata. |

## `x2_netcrop_bm.R`

Location: [`x2_netcrop_bm.R`](./R/x2_netcrop_bm.R)

The core algorithm uses base R plus the dependencies required by the selected
decomposition, clustering, and matrix representation. Plotting optionally
requires `ggplot2`.
Periods in the S3 method names below are required by R method dispatch and are
external-interface exceptions to the underscore naming rule.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `netcrop_blockmodel()` | Selects SBM/DCBM and `K` by the original three-stage overlapping-subnetwork flow. When either partition parameter is `NULL`, it calls `netcrop_param_select(test_prop = 0.02, n = n, o_range = 0)` and preserves any supplied non-`NULL` partner. Explicitly supplied `num_subnetworks` and `overlap_size` always take priority and cause selection options to be ignored. Direct Laplacian and regularization controls are shared by SBM and DCBM; matching legacy nested values remain accepted when the direct arguments are omitted. The remaining model/loss algorithm is unchanged; the result carries `algorithm = "NETCROP"` for shared method dispatch with ECV. | Core NETCROP arguments with `num_subnetworks = NULL`, `overlap_size = NULL`, `nrep = 1L`, `losses = "sse"`, `laplacian = FALSE`, `regularize_tau = 0`, `parameter_select_options = list()` forwarded only during automatic selection, nested estimator/clustering options, matching controls, `ncores`, `seed`, `verbose = TRUE`, `ram_check = FALSE`, and retention choice. Automatic options must select exactly one pair and cannot override `n`. | Classed `netcrop_blockmodel` object with unchanged retention behavior and an explicit algorithm label. |
| `print.netcrop_blockmodel()` | Prints the overall selected model and `K` for every loss. Uses the result's `algorithm` field to label NETCROP and ECV results accurately. | `x`; `...`: S3 compatibility. | Input object invisibly. |
| `summary.netcrop_blockmodel()` | Builds an algorithm-aware durable summary. NETCROP summaries contain overlap and matching diagnostics; ECV summaries contain edge-sampling diagnostics; NCV summaries contain fold sizes, DCBM method, Laplacian choice, failures, workers, and timings. | `object`; `...`: S3 compatibility. | Classed `summary.netcrop_blockmodel` object. |
| `print.summary.netcrop_blockmodel()` | Prints an algorithm-aware NETCROP, ECV, or NCV block-model summary and the per-repetition and overall selections. | `x`; `...`: S3 compatibility. | Input summary invisibly. |
| `plot.netcrop_blockmodel()` | Plots finite candidate CV losses against `K` and labels the plot from the result's `algorithm` field. Passing two or three distinct NETCROP, ECV, and NCV results in any order produces a comparison, including `plot(netcrop_output, ecv_output, ncv_output)`. Irregular and singleton `K` sets are supported. | `x`; `aggregate = TRUE`: logical aggregation choice or a second block-model result; `...`: additional block-model results and comparison options such as `loss_scale`. | `ggplot` object; requires `ggplot2`. |

## `x3_netcrop_rdpg.R`

Location: [`x3_netcrop_rdpg.R`](./R/x3_netcrop_rdpg.R)

This variant supports symmetric, non-negative, loop-free adjacency matrices.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `netcrop_rdpg()` | Selects a symmetric-RDPG latent dimension using the original four-stage flow. When either partition parameter is `NULL`, it calls `netcrop_param_select(test_prop = 0.02, n = n, o_range = 0)` and preserves any supplied non-`NULL` partner. Explicitly supplied `num_subnetworks` and `overlap_size` always take priority and cause selection options to be ignored. The remaining decomposition, alignment, and loss algorithm is unchanged; the result carries `algorithm = "NETCROP"` for shared dispatch with ECV. | `A`; `d_candidates`; `num_subnetworks = NULL`; `overlap_size = NULL`; `nrep = 1L`; `parameter_select_options = list()` forwarded only during automatic selection; losses; eigendecomposition options; `ncores`; `seed`; `verbose = TRUE`; `ram_check = FALSE`; retention choice. Automatic options must select exactly one pair and cannot override `n`. | Classed `netcrop_rdpg` object with unchanged retention behavior and an explicit algorithm label. |
| `print.netcrop_rdpg()` | Prints the overall selected dimension for every loss, labeled from the NETCROP or ECV algorithm field. | `x`; `...`: S3 compatibility. | Input object invisibly. |
| `summary.netcrop_rdpg()` | Builds an algorithm-aware durable summary. NETCROP summaries contain overlap, sparsity, and eigendecomposition diagnostics; ECV summaries contain edge-sampling, repetition-failure, worker, and timing diagnostics. | `object`; `...`: S3 compatibility. | Classed `summary.netcrop_rdpg` object. |
| `print.summary.netcrop_rdpg()` | Prints the algorithm-aware RDPG summary and per-repetition and overall dimension selections. | `x`; `...`: S3 compatibility. | Input summary invisibly. |
| `plot.netcrop_rdpg()` | Plots finite CV loss against `d`, labeled by algorithm. The default aggregates repetitions with mean curves and one-standard-deviation ribbons; the unaggregated view distinguishes repetitions. Passing one NETCROP and one ECV RDPG result in either order produces their comparison through ordinary `plot(result_1, result_2)`. | `x`; `aggregate = TRUE`: logical aggregation choice or a second RDPG result; `...`: another RDPG result or comparison options such as `loss_scale`. | `ggplot` object; requires `ggplot2`. |

## `x4_netcrop_lsm.R`

Location: [`x4_netcrop_lsm.R`](./R/x4_netcrop_lsm.R)

This variant requires a symmetric, binary, loop-free adjacency matrix. Its LSM fits
are dense even when the supplied adjacency is sparse. Periods in S3 method
names are required by R method dispatch and are external-interface exceptions
to the underscore naming rule.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `netcrop_lsm()` | Selects an LSM latent dimension through the established three-stage NETCROP flow. It leaves `lsm_pgd()` unchanged and exactly rewrites each fitted `beta_i + beta_j + X_i'X_j` as `gamma_i + gamma_j - ||Z_i-Z_j||^2`, with `Z_i = X_i / sqrt(2)` and `gamma_i = beta_i + ||X_i||^2 / 2`. It rigidly aligns distance coordinates through overlap nodes without changing invariant `gamma`, then evaluates held-out probabilities from squared cross-distances. The penalty `2 * (effective_overlap_size + piece_size + d)` uses the actual candidate dimension. Candidates are an exact set, including singletons. When either partition parameter is `NULL`, automatic selection uses `netcrop_param_select(test_prop = 0.02, n = n, o_range = 0)` while preserving a supplied partner; supplying both makes explicit values take priority. | `A`; `d_candidates`; `num_subnetworks = NULL`; `overlap_size = NULL`; `nrep = 1L`; `losses = "sse"`; `lsm_options = list()` may override `step_size = 0.3`, `niter = 100L`, `trace = FALSE`, `epsilon = 1e-6`, and `use_cpp = TRUE`; standard `ncores`; `seed`; `verbose = TRUE`; `force_windows = FALSE`; `ram_check = FALSE`; `parameter_select_options = list()`; `retain_intermediates = c("all", "minimal")`. | Classed `netcrop_lsm` object containing detailed and aggregated penalized CV losses, per-repetition and overall dimensions, partitions, task grids, optional raw/reparameterized fits and rigid alignments, resolved options, worker counts, timings, and RAM diagnostics. |
| `print.netcrop_lsm()` | Prints the overall selected LSM dimension for every loss. | `x`; `...`: S3 compatibility. | Input object invisibly. |
| `summary.netcrop_lsm()` | Builds a durable summary with candidate dimensions, effective overlap, held-out coverage, selections, resolved LSM options, workers, and timings. | `object`; `...`: S3 compatibility. | Classed `summary.netcrop_lsm` object. |
| `print.summary.netcrop_lsm()` | Prints the LSM NETCROP summary and its per-repetition and overall selections. | `x`; `...`: S3 compatibility. | Input summary invisibly. |
| `plot.netcrop_lsm()` | Plots penalized CV loss against `d`; the default aggregates repetitions with mean curves and one-standard-deviation ribbons, while the unaggregated view distinguishes repetitions. Arbitrary and singleton candidate sets are supported. | `x`; `aggregate = TRUE`; `...`: S3 compatibility. | `ggplot` object; requires `ggplot2`. |

## `x5_netcrop_regularizer.R`

Location: [`x5_netcrop_regularizer.R`](./R/x5_netcrop_regularizer.R)

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `pair_nmi_loss()` | Forms the same Cartesian-product community labels as the legacy routine and returns their negative joint-entropy-normalized mutual information. Constant-versus-constant partitions have loss `-1`. | Positive-integer label vectors `g_1`, `g_2`: first pair; `g_3`, `g_4`: reference pair; the Cartesian products must have equal sizes. | Finite numeric scalar when NMI is defined, otherwise `NA_real_`. |
| `pair_hamming_loss()` | Forms legacy Cartesian-product community labels, matches them globally with the Hungarian algorithm, and returns their mismatch proportion. | Positive-integer label vectors `g_1`, `g_2`: first pair; `g_3`, `g_4`: reference pair; the Cartesian products must have equal sizes. | Numeric scalar in `[0, 1]`. |
| `netcrop_tune_regularizer()` | Cross-validates candidate spectral regularization values for a fixed SBM or DCBM `K` using overlapping subnetworks. The default SSE path preserves the legacy probability-estimation and held-out-edge logic. Label losses preserve the legacy full-network reference by default; `label_reference = "leave_pair_out"` optionally masks the evaluated cross-piece edges. Inputs, option conflicts, folds, worker errors, non-finite losses, and candidate selection are validated or diagnosed. Remainders augment the overlap so no node is omitted. Overall tau is an evaluated candidate minimizing loss averaged across repetitions, with the smallest tau breaking ties. `K = 1` emits a message and returns invisibly because tuning is unnecessary. The output carries `algorithm = "NETCROP"` for shared dispatch with DKEST. | `A`; `K`; unique non-negative `tau_candidates`; `use_dcbm = FALSE`; `num_subnetworks = NULL`; `overlap_size = NULL`; `nrep = 1L`; `use_laplacian = FALSE`; `dcbm_est_method`; `losses = "sse"`; optional named `loss_types`; `label_reference = c("full_network", "leave_pair_out")`; spectral, clustering, estimator, matching, worker, seed, verbosity, RAM, automatic-parameter-selection, and retention controls. | Classed `netcrop_regularizer` result containing compatibility `loss` and `selection` objects, tidy CV losses, per-repetition and overall selections, diagnostics, split metadata, options, worker counts, timings, RAM preflight, optional intermediates, and an explicit algorithm label. |
| `print.netcrop_regularizer()` | Prints the model, fixed `K`, and overall selected tau, labeled for NETCROP or DKEST. | `x`; `...`. | Input invisibly. |
| `summary.netcrop_regularizer()` | Builds an algorithm-aware summary: split and CV diagnostics for NETCROP, or candidate validity and DK diagnostics for DKEST. | `object`; `...`. | Classed `summary.netcrop_regularizer` object. |
| `print.summary.netcrop_regularizer()` | Prints the algorithm-aware regularizer-selection summary. | `x`; `...`. | Input invisibly. |
| `plot.netcrop_regularizer()` | Plots finite NETCROP CV losses against tau, aggregated across repetitions or separated by repetition. DKEST results explicitly decline this method because DKEST estimates `tau_hat` directly and has no CV loss. | `x`; `aggregate = TRUE`; `...`: S3 compatibility. | `ggplot` object for NETCROP; informative error for DKEST. |

## `x6_other_algos.R`

Location: [`x6_other_algos.R`](./R/x6_other_algos.R)

The three low-level ECV functions preserve the supplied randnet-derived code,
including its legacy dotted external identifiers. The only changes within
those functions are the canonical selectable loss names `sse`, `bin_dev`, and
`auc_as_loss`, plus using the plug-in DCBM probability estimate for the
`K = 1` DCBM SSE. The high-level stability wrapper performs validation, dependency
injection, cross-platform outer parallelism, result auditing, and conversion to
the shared `netcrop_blockmodel` result structure. The file also contains the
validated RDPG ECV wrapper and the oracle-free DKEST regularizer selector, both
using the corresponding NETCROP result classes and S3 methods. Multi-tau SONNET
and spectral-clustering runners support single networks or matched lists, and an
oracle accuracy plotter can evaluate NETCROP and DKEST selections under either
clustering engine independently of which tuner produced each selected tau.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `holdout.evaluation.fast.all()` | Preserved low-level ECV fold evaluator. It removes selected upper-triangle dyads, performs successive truncated-SVD reconstructions, fits SBM and degree-corrected candidates, and evaluates the requested canonical losses. At DCBM `K = 1`, SSE now uses the same plug-in `P.hat` used by binary deviance and AUC-as-loss. Legacy identifiers and raw component names are retained. | Legacy arguments `holdout.index`, `A`, `max.K`, `tau = 0`, `dc.est = 2`, `p.sample = 1`; `loss = c("sse", "bin_dev", "auc_as_loss")`. | Legacy list of per-`K` SBM/DCBM loss vectors and diagnostics. |
| `iter.SVD.core.fast.all()` | Preserved low-level successive truncated-SVD reconstruction used by ECV. It optionally converts to sparse storage, treats missing entries as zero, rescales by the sampling probability, and retains raw and clipped reconstructions through `Kmax`. | Legacy arguments `A`, `Kmax`, `tol = 1e-5`, `max.iter = 100`, `sparse = TRUE`, `tau = 0`, `p.sample = 1`. | List of SVD components plus raw and thresholded reconstructions for every rank. |
| `ECV.BM()` | Preserved single-repetition ECV driver. It samples `cv` holdouts, calls the low-level evaluator, aggregates legacy loss matrices, and returns its original raw selection structure. Only canonical loss names replace the pasted names. | Legacy arguments `A`, `max.K`, `cv = 3`, `holdout.p = 0.1`, `tau = 0`, `dc.est = 2`, `loss = c("sse", "bin_dev", "auc_as_loss")`, `ncore = 1`, `seed = 100`. | Original ECV list of raw fold matrices, aggregated loss vectors, and legacy best-model strings. |
| `ecv_stability_blockmodel()` | Validated high-level repeated ECV block-model selector using only the plug-in DCBM estimator (`dc.est = 2`). It checks binary undirected loop-free input, candidate and sampling feasibility, dependencies, seeds, workers, and requested losses before calling the preserved routines. Repetitions use `uni_mclapply()` while each nested `ECV.BM()` call is sequential, preventing worker oversubscription. A private runtime environment supplies namespace-safe `mclapply`, `Matrix`, transpose, sparse-`which`, and clustering behavior required by the legacy code. It performs conservative RAM preflight, validates every raw result, optionally stops on or omits failed repetitions, converts raw fold matrices to the shared tidy loss schema, and minimizes canonical losses. When `seed = NULL`, independent valid repetition seeds are generated before parallel execution and retained in the result. | `A`; `max_K`; `train_proportion = 0.9`; `cv = 3L`; `nrep = 20L`; `tau = 0`; `losses = "sse"`; standard `ncores`; `seed = NULL`; `verbose = TRUE`; `force_windows = FALSE`; `ram_check = FALSE`; `failure_handling = c("stop", "omit")`; `retain_intermediates = c("all", "minimal")`. | A classed `netcrop_blockmodel` object with `algorithm = "ECV"`, common tidy losses/selections/candidates, ECV sampling metadata and realized repetition seeds, diagnostics, worker allocation, timing, RAM report, and optional legacy raw output. |
| `ncv_bm()` | Runs one corrected Chen-Lei node-CV repetition without changing its stage order: permute nodes, form rectangular training matrices and held-out fold blocks, obtain right singular vectors, cluster SBM and row-normalized DCBM representations, estimate rectangular SBM/DCBM parameters, evaluate held-out dyads, and average folds. It uses `estimate_sbm()` and `estimate_dcbm()` directly; spectral DCBM is the default. Undirected losses use one off-diagonal triangle, and probabilities are clipped consistently. | `A`; `max_K`; `cv = 3L`; `dc_est = c("spectral", "plugin")`; `tau = 0`; `use_laplacian = FALSE`; `losses = "sse"`; standard `ncores`; `seed = NULL`; `force_windows = FALSE`; `validate_inputs = TRUE`. | List with tidy fold losses, fold-averaged candidate losses, best candidates, fold membership/sizes, permutation, and realized fold seeds. |
| `ncv_stability_blockmodel()` | Repeats `ncv_bm()` using `uni_mclapply()` at the outer level with sequential fold stages, spectral DCBM by default, preflight validation, deterministic or generated repetition seeds, conservative RAM reporting, strict or omission-based failure handling, and optional raw-output retention. | `A`; `max_K`; `cv = 3L`; `nrep = 20L`; `dc_est = c("spectral", "plugin")`; `tau = 0`; `use_laplacian = FALSE`; canonical `losses`; standard `ncores`; `seed = NULL`; `verbose = TRUE`; `force_windows = FALSE`; `ram_check = FALSE`; `failure_handling = c("stop", "omit")`; `retain_intermediates = c("all", "minimal")`. | A classed `netcrop_blockmodel` object with `algorithm = "NCV"`, common tidy losses/selections/candidates, fold and estimator metadata, failures, workers, timing, RAM report, seeds, and optional raw repetition output. |
| `ecv_stability_rdpg()` | Validated repeated edge-CV selector for symmetric RDPG dimension. The three pasted lower-level functions are preserved unchanged inside an isolated environment, preventing their shared legacy SVD name from overwriting block-model ECV. The wrapper checks graph shape, symmetry, loops, values, binary-loss compatibility, rank and holdout feasibility, dependencies, seeds, workers, and requested losses; supplies namespace-safe legacy dependencies and canonical loss adapters; runs repetitions through `uni_mclapply()` with sequential folds; performs conservative dense-RAM preflight; audits each result; and supports strict or omission-based failures. Independent realized seeds are retained when `seed = NULL`. The canonical `mse` label reflects the legacy helper's mean squared error; relative comparison treats it as the same family as NETCROP `sse`. | `A`; `max_d`; `cv = 3L`; `nrep = 1L`; `train_proportion = 0.9`; `losses = c("mse", "bin_dev", "auc_as_loss")`; standard `ncores`; `seed = NULL`; `verbose = TRUE`; `force_windows = FALSE`; `ram_check = FALSE`; `failure_handling = c("stop", "omit")`; `retain_intermediates = c("all", "minimal")`. | A classed `netcrop_rdpg` object with `algorithm = "ECV"`, tidy repetition-averaged loss curves, selections, edge-holdout metadata, realized seeds, failures, worker allocation, timing, RAM diagnostics, and optional legacy raw output. |
| `plot_rdpg_comparison()` | Internal implementation used by standard `plot()` dispatch for one NETCROP and one ECV RDPG output in either order. Relative mode matches compatible loss families and min-max normalizes within algorithm; raw mode requires identical loss names. | `...`: two classed `netcrop_rdpg` results; `loss_scale = c("relative", "raw")`. | `ggplot` object; requires `ggplot2`; call through `plot(netcrop_result, ecv_result)`. |
| `dkest_tune_regularizer()` | Directly estimates tau by minimizing the DK eigenspectral ratio while preserving the supplied calculation order: regularize the observed matrix, optionally form its normalized Laplacian, cluster the fixed-`K` embedding, estimate SBM/DCBM probabilities, construct the fitted regularized matrix, and compute the residual-to-fitted eigenvalue ratio. It does not estimate or expose CV loss. The oracle `true.g` input and accuracy calculation are removed. Preflight checks cover graph validity, ranks, tau candidates, zero-degree feasibility, dependencies, seeds, workers, sparse densification, and RAM. Candidate tasks use `uni_mclapply()` and independently realized seeds; malformed clustering, eigensolver failures, zero denominators, and non-finite statistics are stopped or omitted according to policy. | `A`; fixed `K`; unique non-negative `tau_candidates`; `use_laplacian = TRUE`; `use_dcbm = TRUE`; `dcbm_est_method = c("plugin", "spectral")`; standard `ncores`; `seed = NULL`; `verbose = TRUE`; `force_windows = FALSE`; `ram_check = FALSE`; `failure_handling = c("stop", "omit")`; `retain_intermediates = c("all", "minimal")`. | Classed `netcrop_regularizer` object with `algorithm = "DKEST"`, direct `tau_hat`, selected DK statistic, per-candidate numerator/denominator/statistic table, diagnostics, realized candidate seeds, workers, timing, RAM report, and optional raw candidate output; no `cv_loss` field. |
| `mult_reg_spectral_cluster()` | Runs `spectral_cluster()` for every network-by-tau combination after validating common network dimensions, candidate uniqueness, graph structure, options, seeds, worker settings, fitted labels, failures, and optional RAM demand. Candidate tasks are parallelized with `uni_mclapply()`; the same realized seed is used across tau candidates within a network to make stochastic clustering comparisons fair. | `A`: one adjacency matrix or a non-empty list with shared dimensions; `K`; unique non-negative `tau_candidates`; all applicable `spectral_cluster()` controls with their current defaults; standard `ncores`; `seed = NULL`; `verbose = TRUE`; `force_windows = FALSE`; `ram_check = FALSE`; `failure_handling = c("stop", "omit")`; `retain_fits = TRUE`. | Classed `mult_reg_clustering` object containing nested labels and optional fits by network and tau, resolved parameters, task grid, diagnostics, realized network seeds, worker allocation, timing, and RAM preflight metadata. |
| `mult_reg_sonnet()` | Runs `sonnet()` at every requested tau for one network or a matched list. It retains SONNET's current defaults, validates common dimensions, candidates, named forwarded options, protected arguments, seeds, output labels, and failures, and deliberately keeps the network-by-tau loop sequential so `ncores` parallelism remains inside SONNET. The same realized seed is reused across tau candidates within each network. | `A`: one adjacency matrix or a non-empty list with shared dimensions; `K`; unique non-negative `tau_candidates`; all SONNET controls with their current defaults; named spectral-clustering arguments through `...`; `failure_handling = c("stop", "omit")`; `retain_fits = TRUE`. | Classed `mult_reg_clustering` object containing nested labels and optional SONNET fits by network and tau, resolved parameters, task grid, diagnostics, realized network seeds, timing, and the call. |
| `oracle_plotter()` | Computes accuracy as one minus the validated label-matching mismatch rate and plots mean accuracy. A single network produces points without error bars; multiple networks add horizontal plus/minus-one-SD bars. Without tuner outcomes it compares every candidate tau. With optional NETCROP and/or DKEST outcomes it compares engine-specific oracle accuracy with every available tuner selection: NETCROP's first repetition and optionally its repetition mean and mode for each chosen loss, plus DKEST's direct `tau_hat`. Off-grid selected values trigger an extra fit. `engines` independently selects SONNET, spectral clustering, or both, so every supplied tuner's tau is evaluated under every requested engine regardless of provenance. For outcome-derived SONNET settings, NETCROP supplies its subnetwork count and overlap but `extra_nrep` defaults to zero; explicit `sonnet_options` override all recovered defaults. | `A` and required `g_true`: one matched matrix/vector pair or same-length lists with shared dimensions; `tau_candidates`; optional `K`; `netcrop_outcomes` and `dkest_outcomes` may be single classed outputs for scalar inputs, while list inputs require matching outcome lists; `include_netcrop_mean = TRUE`; `include_netcrop_mode = TRUE`; optional `losses`; `engines = c("sonnet", "spectral_cluster")`; validated `matching_method`; brute-force confirmation; named fitter option lists; standard worker, seed, verbosity, cross-platform, and RAM controls. | Invisibly returns the rendered `ggplot`, faceted by requested engine. Raw accuracies, aggregated plotting data, effective fit metadata, and diagnostics are attached as `accuracy_data`, `plot_data`, `metadata`, and `diagnostics` attributes. |
| `plot_blockmodel_comparison()` | Internal implementation used by standard `plot()` dispatch for two or three distinct NETCROP, ECV, and NCV outputs in any argument order. Relative mode matches compatible loss families and min-max normalizes within algorithm; raw mode requires identical loss names. | `...`: two or three classed block-model results; `loss_scale = c("relative", "raw")`. | `ggplot` object; requires `ggplot2`; call through `plot(result_1, result_2, result_3)`. |

## Generated Rcpp interfaces

Locations: [`R/RcppExports.R`](./R/RcppExports.R) and
[`src/RcppExports.cpp`](./src/RcppExports.cpp).

`Rcpp::compileAttributes()` generates internal R wrappers and registered C++
entry points for `auc_cpp()`, `outer_add_cpp()`,
`procrustes_translated_cpp()`, `connected_components_dense_cpp()`,
`connected_components_sparse_cpp()`, `shortest_path_distances_dense_cpp()`,
`shortest_path_distances_sparse_cpp()`, and `lsm_pgd_cpp()`. Their signatures
and return values are recorded in the corresponding C++ sections above. None
is exported through `netOP::`; each is reached only by a validated R wrapper.

## Dependency classification

| Classification | Packages and use |
|---|---|
| `Imports` | `Rcpp` for registered interfaces; `Matrix` for sparse matrices; `RSpectra` and `irlba` for partial decompositions; `cluster` for PAM/CLARA; `tibble` for regularizer output. |
| `LinkingTo` | `Rcpp`, `RcppEigen`. |
| `Suggests` | `future` and `future.apply` for Windows/multisession work; `ggplot2` for plots; `ps` for available RAM; `entropy` for an optional NMI implementation; `peakRAM` for explicit resource measurement; `testthat`, `knitr`, `rmarkdown`, and `pkgdown` for development and documentation. |
| Development-only reference | CRAN `randnet` 1.0, used for provenance/equivalence auditing and deliberately absent from every dependency field. |
| Removed legacy declarations | `dplyr`, `Rfast`, `data.table`, `readr`, `softImpute`, `pROC`, `IMIFA`, `rlist`, and `latentnet` were not required by active installed-package code. |

The file-level license, copyright, origin, modification, and dependency-license
inventory is maintained in `inst/COPYRIGHTS`.

## Maintenance checklist

Whenever code changes:

1. Update the file tree if a file was added, removed, or renamed.
2. Add, remove, or revise the corresponding function row.
3. Keep argument names/default behavior and return descriptions synchronized.
4. Mark internal helpers and compiled exports clearly.
5. Record new optional package dependencies and fallback behavior.
6. Check that every top-level R function and every named C++ function appears here.
