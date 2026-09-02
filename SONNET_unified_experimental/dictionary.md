# Function Dictionary

This is the living function reference for `SONNET_unified_experimental`. It covers every
top-level named R function, every named C++ helper, and every Rcpp export.
Anonymous and lexically nested closures are implementation details of their
owning function and are not listed as separate APIs.

Update this file whenever a function is added, removed, renamed, moved, or has
its arguments, return value, dependencies, or behavior changed.

## File tree

```text
SONNET_unified_experimental/
├── 01_basic_helpers.R       Cross-platform execution, RAM, and simulation helpers
├── 02_math_helpers.R        Losses, transforms, alignment, and decompositions
├── 02_math_helpers.cpp      Optional math accelerators exposed through Rcpp
├── 03_network_utils.R       Components and shortest-path utilities
├── 03_network_utils.cpp     Optional network-algorithm accelerators
├── 04_network_generator.R   Network I/O and random graph generators
├── 05_embedders.R           Spectral embedding and latent-space fitting
├── 05_embedders.cpp         Optional latent-space PGD accelerator
├── 06_estimators.R          Network-model parameter estimators
├── x1_sonnet.R              SONNET splitting, estimation, and S3 methods
├── xx_helpers.R             Helpers shared by NETCROP variants
├── x2_netcrop_bm.R          NETCROP block-model selection and S3 methods
├── x3_netcrop_rdpg.R        NETCROP symmetric-RDPG selection and S3 methods
├── source_all.R              Interactive ordered selector for sourcing R files
├── Archive.zip               Preserved source archive from the copied directory
├── NAMING_CONVENTION.md     Mandatory naming and implementation rules
└── dictionary.md            This living function reference
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
| `average_degree_method` | `"calibrated"` targets the final expected degree; `"naive"` uses the historical one-step or iterative approximation. |
| `align_with` | Optional reference embedding used for rotation-only Procrustes alignment. |

## `source_all.R`

Location: [`source_all.R`](./source_all.R)

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `source_all()` | Detects the directory containing `source_all.R`, lists every other `.R` file alphabetically, and interactively accepts serial numbers in the desired source order. Enter `0` alone for every file alphabetically or `-1` alone for none. The selector excludes itself to prevent recursion. It runs automatically when `source_all.R` is sourced in an interactive R session and can be called programmatically in non-interactive sessions. | `selection = NULL`: prompt, ordered integer vector, or comma/space-separated string; `envir = parent.frame()`: sourcing environment; `echo = FALSE`; `verbose = TRUE`. | Invisibly, the ordered absolute paths that were sourced, or an empty character vector. |

## `01_basic_helpers.R`

Location: [`01_basic_helpers.R`](./01_basic_helpers.R)

### Public functions

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
| `run_simulations()` | Repeats `one_simulation`, manages seeds, optional outer-loop parallelism, progress, errors, resource measurements, and resumable RDS results. The outer loop is sequential by default because simulation methods often parallelize internally. | `one_simulation`: function accepting named `simulation`; `nsim`: repetitions; `...`: forwarded inputs; `use_parallel_simulations`: parallelize outer loop; `ncores_outer`: outer workers; `seed` or `seeds`: base or explicit seeds; `results_file`: optional RDS path; `action`: `"replace"`, `"resume"`, or `"archive"`; `retry_failed`: rerun failed records; `continue_on_error`: continue after failure; `measure_resources`: append elapsed time and peak RAM; `show_progress`: print completion messages; `force_windows`: select Windows backend. | List of simulation records; invisibly returned. Each record contains simulation number, seed, success, elapsed time, peak RAM, result, and error. |

## `02_math_helpers.R`

Location: [`02_math_helpers.R`](./02_math_helpers.R)

### Losses and scalar/vector transforms

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `modal()` | Finds the most frequent value; the first encountered value breaks frequency ties. | `x`: input vector; `na_rm`: remove missing values first. | One modal value, or a typed `NA` if removal leaves no values. |
| `validate_error_inputs()` | Internal validator shared by vector error metrics. | `x`, `y`: equal-length numeric vectors; `na_rm`: single logical. | `invisible(TRUE)` or an error. |
| `sse()` | Sum of squared errors. | `x`, `y`: equal-length numeric vectors; `na_rm`: remove missing differences; `validate_inputs = TRUE`: disable only after a high-level caller has validated both vectors. | Numeric scalar. |
| `sae()` | Sum of absolute errors. | Same as `sse()`. | Numeric scalar. |
| `mse()` | Mean squared error. | Same as `sse()`. | Numeric scalar. |
| `mae()` | Mean absolute error. | Same as `sse()`. | Numeric scalar. |
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

Location: [`02_math_helpers.cpp`](./02_math_helpers.cpp)

Compile with `Rcpp::sourceCpp("02_math_helpers.cpp")`. All exports have R
fallbacks in `02_math_helpers.R`.

| Function | Visibility | Description and arguments | Returns/use |
|---|---|---|---|
| `auc_cpp()` | Rcpp export | `group`: binary numeric labels; `predictions`: finite scores. Sorts compact label/score pairs and handles ties by average rank. | Numeric AUC; accelerator for `auc()`. |
| `outer_add_cpp()` | Rcpp export | `x`, `y`: numeric vectors. Allocates an uninitialized output and fills `y[row] + x[column]`. | Numeric matrix; accelerator for `outer_add()`. |
| `procrustes_translated_cpp()` | Rcpp export | `X`, `X_star`: equal finite dense matrices. Uses Eigen SVD without an `n` by `n` centering matrix. | List with `X_new`, `rotation`, and `translation`; accelerator for `procrustes()`. |

## `03_network_utils.R`

Location: [`03_network_utils.R`](./03_network_utils.R)

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `validate_adjacency()` | Internal square, numeric, finite adjacency validator supporting base and `Matrix` objects. | `A`: adjacency matrix. | `invisible(TRUE)` or an error. |
| `adjacency_neighbors()` | Builds an unweighted adjacency list without densifying sparse input. | `A`; `directed`: retain orientation or union both directions; `self_loops`: ignore/include diagonal neighbors. | List of integer neighbor vectors. |
| `connected_components()` | Breadth-first undirected/weak connected components. Edge weights are ignored. Tries dense/sparse C++ exports before the R implementation. | `A`; `self_loops`; `use_cpp`. | Named list of one-based node-index vectors. |
| `largest_connected_component()` | Extracts the first largest weak component and its induced adjacency matrix. | `A`; `self_loops`; `use_cpp`; `sort_nodes`: sort returned indices. | List with `nodes`, `size`, and `submatrix`. |
| `adjacency_weighted_neighbors()` | Builds node and weight adjacency lists for weighted paths without densifying sparse input. For undirected reciprocal weights, the smaller nonzero weight is used. | `A`; `directed`; `self_loops`. | List with parallel `nodes` and `weights` lists. |
| `shortest_path_distances()` | All-pairs shortest-path distance matrix. Uses BFS when unweighted and Dijkstra when weighted; has dense/sparse compiled paths and dependency-light R fallbacks. | `A`; `directed`; `weighted`: use nonnegative entries as distances; `self_loops`; `use_cpp`. | `n` by `n` matrix with zero diagonal and `Inf` for unreachable pairs. |

## `03_network_utils.cpp`

Location: [`03_network_utils.cpp`](./03_network_utils.cpp)

Compile with `Rcpp::sourceCpp("03_network_utils.cpp")`. The R-facing utilities
validate inputs and fall back to R when these exports are absent or fail.

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

Location: [`04_network_generator.R`](./04_network_generator.R)

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
| `write_network()` | Writes an edge list with node pairs in the first two columns and optional weights in the third. Metadata preserves isolated nodes. Undirected files store one upper/lower triangle including loops. | `A`; `file`; `directed`: explicit or inferred; `weighted`: explicit or inferred; `triangle`: `"upper"`/`"lower"`; `format`; `include_header`; `tolerance`. | Input file path, invisibly. |
| `read_network()` | Reads two- or three-column edge lists into dense/general sparse adjacency. Undirected matrices are formed with `A + t(A)` and loops added once. Rejects duplicate pairs. | `file`; `representation`; `directed`, `weighted`: explicit, metadata-derived, or default; `n`: optional dimension for isolates; `format`; `has_header`. | Dense matrix or general `dgCMatrix`. |

### Probability and average-degree helpers

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `validate_probability_matrix()` | Internal finite `[0,1]`, square, and optional symmetry validator. | `P`; `directed`; `tolerance`. | `invisible(TRUE)` or an error. |
| `clip_generator_probabilities()` | Internal probability clipping independent of math-helper source order. | `P`; `lower_clip`; `upper_clip`. | Clipped probability matrix. |
| `scale_to_average_degree()` | Calibrates a global probability multiplier by bisection against the final clipped, loop-adjusted expected degree. | `P`; `average_degree`; `self_loops`; clipping bounds; `tolerance`; `max_iterations`. | List with calibrated `P` and `multiplier`. |
| `scale_to_average_degree_naive()` | Historical one-step scaling: target divided by `mean(rowSums(P))`, followed by clipping and diagonal removal. | `P`; `average_degree`; `self_loops`; clipping bounds. | List with scaled `P` and `multiplier`; final degree may differ from target. |
| `apply_average_degree_scaling()` | Dispatches to calibrated or naive probability scaling. | `P`; `average_degree`; `average_degree_method`; loop and clipping controls. | List with `P` and `multiplier`. |

### General, ER, RDPG, SBM, and DC-SBM generators

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `generate_adjacency()` | Samples independent Bernoulli edges from `P`. Undirected sampling draws one triangle and constructs a general symmetric matrix with `A + t(A)`; row seeds make results independent of worker count. Full probability validation remains enabled for direct matrix calls. A scalar `P` plus `n` avoids allocating a constant `n` by `n` probability matrix. | `P`: probability matrix, or scalar with `n`; `representation`; `directed`; `self_loops`; `seed`; `ncores`; `validate_inputs = TRUE`; `n = NULL`: dimension for scalar `P`. | Binary adjacency matrix. |
| `generate_er()` | Erdős–Rényi generator parameterized by exactly one of edge probability or expected average degree. It samples from a scalar probability and does not allocate a redundant dense probability matrix. | `n`; `p` or `average_degree`; representation, direction, loop, seed, and worker controls. | Adjacency `A`; metadata contains model, `p`, and expected degree. |
| `generate_latent_positions()` | Validates a supplied `Z` or generates an `n` by `d` matrix from a user-selected distribution. Default uniform bounds keep typical dot products controlled. | `n`, `d`: required when `Z` missing; `Z`: optional direct matrix; `latent_distribution`; named `latent_args`; `seed`. | Numeric latent-position matrix. |
| `generate_rdpg()` | Generates undirected RDPG from `Z` or directed generalized RDPG from `Z_left`, `Z_right`; any missing positions are generated. Supports probability scaling, clipping, sparsity, and average-degree targeting. | Dimensions and latent inputs; `latent_distribution`, `latent_args`; `sparsity_multiplier`; `scale_P`; `average_degree`, method; clipping; representation, direction, loops, seed, workers. | Adjacency `A`; metadata contains latent truth and degree-scaling details. |
| `generate_community_labels()` | Validates supplied labels and maps them to `1:K`, or samples labels from community probabilities. Original label levels are retained as an attribute after mapping. | `n`, `K`; `g_true`; `community_probabilities`; `seed`. | Integer label vector. |
| `generate_degree_parameters()` | Validates direct `psi` or generates nonnegative degree parameters from a distribution; optionally normalizes within communities. | `n`; `psi`; `degree_distribution`; named `degree_args`; `degree_scale`: `"none"`, `"max_by_community"`, or `"mean_one_by_community"`; `g_true`; `seed`. | Numeric `psi` vector. |
| `generate_dcbm()` | Degree-corrected SBM using direct/generated labels and degree parameters. Uses direct `P_block` or constructs within/between probabilities from `alpha` and `beta`. Supports directed blocks and both degree-targeting methods. | `n`, `K`, labels/proportions; `P_block` or `alpha`, `beta`; `psi` or degree distribution controls; `sparsity_multiplier`; average-degree and clipping controls; representation, direction, loops, seed, workers. | Adjacency `A`; metadata contains `g_true`, `psi`, `P_block`, and scaling details. |
| `generate_sbm()` | Ordinary SBM wrapper around `generate_dcbm()` with `psi = 1`. | Same block, label, sparsity, average-degree, clipping, representation, direction, loop, seed, workers, excluding degree-correction inputs. | Adjacency `A` with model metadata changed to `"sbm"`. |

### Latent-space-model helpers and generator

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `generate_truncated_normal()` | Dependency-free inverse-CDF sampler for a standard normal truncated to finite bounds. | `n`: draws; `lower_bound`, `upper_bound`. | Numeric vector. |
| `normalize_lsm_positions()` | Centers columns and applies the Ma/Ma Frobenius-based latent scaling. Computes the exact identity `||ZZ' ||_F = ||Z'Z||_F` through the smaller `d` by `d` cross-product instead of materializing an `n` by `n` Gram matrix. | `Z`: finite nondegenerate latent matrix. | Centered, normalized `Z`. |
| `generate_lsm_positions()` | Generates mixture-component means and adds truncated-normal latent noise. | `n`, `d`, `K`, `g_true`; mean and noise bounds; `seed`. | `n` by `d` latent matrix. |
| `generate_lsm_alpha()` | Validates supplied scalar/vector intercepts or generates negative node effects. | `n`; `alpha`: `NULL`, scalar, or length `n`; `seed`. | Length-`n` numeric intercept vector. |
| `scale_lsm_to_average_degree()` | Adjusts the LSM logit intercept. Calibrated mode bisects against final probabilities; naive mode reproduces repeated `-log(current/target)` updates before diagonal removal. | `theta`; `average_degree`; method; `self_loops`; clipping bounds; `naive_iterations`; `tolerance`; `max_iterations`. | List with probability matrix `P` and `intercept_shift`. |
| `generate_lsm()` | Hoff/Ma-style latent-space generator. Accepts/generates `alpha`, undirected `Z`, or directed `Z_left`/`Z_right`; optionally generates mixture labels, normalizes positions, applies the corrected squared-norm distance identity, and targets average degree. | `n`, `d`, `K`, `g_true`, community probabilities; `alpha`; latent inputs; `normalize_Z`; `distance_adjustment`; mixture/noise bounds; average-degree controls; clipping; representation, direction, loops, seed, workers. | Adjacency `A`; metadata includes positions, intercepts, labels, actual/target degree, method, and intercept shift. |

## `05_embedders.R`

Location: [`05_embedders.R`](./05_embedders.R)

Source `02_math_helpers.R` first. The R-facing PGD function works without the
compiled file and automatically falls back after absence or runtime failure.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `ase()` | Adjacency spectral embedding. Its high-level preflight prints an explicit additive formula for decomposition memory and the largest explicit multiplication estimate times the number of reconstruction/alignment products. | `A`: finite square adjacency; `d`; `directed`; `reconstruct`; `align_with`; `ram_check = TRUE`. | Embedding(s), spectral values, optional alignment, and optional `A_hat`. |
| `spectral_cluster()` | Spectral clustering for undirected, non-negative networks, or direct clustering of a supplied representation. Its high-level preflight prints an explicit additive formula for decomposition, counted Laplacian operations, dense regularization, and PAM dissimilarities. Laplacian input is always estimated as dense even when stored and processed sparsely. | `A`; optional `U`; `K`; Laplacian, regularization, zero-degree, normalization, decomposition, and clustering controls; `ram_check = TRUE`; `validate_inputs = TRUE`: permit validated high-level callers to skip repeated scans. Legacy `spectral_options$ram_check` is promoted to this high-level argument. | List with full `g_hat`, clustered `U_hat`, spectral values, zero-degree and retained indices, clustering fit, and selected engines. |
| `lsm_pgd()` | Fits an undirected binary logistic latent-space model by projected gradient ascent. Its high-level preflight prints an explicit additive formula for spectral initialization, approximately six dense `n` by `n` working matrices, and the largest multiplication estimate times all `2 * niter + 2` explicit products. | `A`: symmetric binary zero-diagonal adjacency; `d`; step and iteration controls; optional initial values; `epsilon`; `use_cpp`; `ram_check = TRUE`. | List with `Z_hat`, `alpha_hat`, zero-diagonal `P_hat`, objective history, and resolved step sizes. |

## `05_embedders.cpp`

Location: [`05_embedders.cpp`](./05_embedders.cpp)

Compile with `Rcpp::sourceCpp("05_embedders.cpp")`.

| Function | Visibility | Description and arguments | Returns/use |
|---|---|---|---|
| `stable_sigmoid()` | Internal C++ | `value`: scalar logit; uses sign-specific evaluation to avoid overflow. | Stable logistic probability. |
| `stable_softplus()` | Internal C++ | `value`: scalar logit. | Stable `log(1 + exp(value))`. |
| `lsm_log_likelihood_cpp()` | Internal C++ | `A`, `theta`: equal square matrices. Sums Bernoulli log likelihood only over `i < j`. | Undirected off-diagonal log-likelihood scalar. |
| `validate_lsm_pgd_cpp_inputs()` | Internal C++ | `A`, `Z`, `alpha`, the two step sizes, and `niter`. Checks dimensions, finiteness, symmetry, binary entries, zero diagonal, and positive controls. | Returns void or raises an R error. |
| `lsm_pgd_cpp()` | Rcpp export | Dense `A`; initial `Z`, `alpha`; `step_size_Z`, `step_size_alpha`; `niter`; `trace`. Uses stable likelihoods, excludes diagonal residuals, centers `Z` after every update, and checks interrupts. | Accelerator result with `Z_hat`, `alpha_hat`, `P_hat`, and `objective`; consumed by `lsm_pgd()`. |

## `06_estimators.R`

Location: [`06_estimators.R`](./06_estimators.R)

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `estimate_sbm()` | Estimates SBM block means from a full adjacency or a common rectangular NCV layout. The full undirected path computes only one block triangle. Rectangular denominators count the dyads actually present, genuine self-pairs are identified from original node indices, and undirected sufficient statistics pool both available orientations. | `A`; `g`; `K`; optional `fold_nodes`; `directed`; `self_loops`; `validate_inputs = TRUE`: skip repeated full matrix/label scans only after high-level validation. | Named `K` by `K` block-mean matrix `B_hat`. |
| `estimate_dcbm()` | Unified plug-in/spectral estimator of DCBM block and degree parameters for full networks and four common NCV matrix layouts. Supplied norms avoid recomputation. | `A`; `g`, `K`; `method`; `fold_nodes`; `row_norm`; `psi_omit`; `stabilizer`; `spectral_engine`, `spectral_options`; `validate_inputs = TRUE`: permit an audited high-level fast path. | List with named block matrix `B_hat` and node degree parameters `psi_hat`. |
| `estimate_sbm_Phat()` | Calls `estimate_sbm()` and expands its block estimates into node-pair probabilities. Returns all-node probabilities normally and fold-by-fold probabilities in NCV mode. | Estimator arguments `A`, `g`, `K`, `fold_nodes`, `directed`, `self_loops`; `lower_clip`, `upper_clip`: finite output bounds. | Square probability matrix `P_hat` for the target nodes. |
| `estimate_dcbm_Phat()` | Calls `estimate_dcbm()` and reconstructs `P_hat = tcrossprod(psi_hat) * B_hat[g, g]`. Returns all nodes normally, retained trailing nodes after `psi_omit`, or fold nodes in NCV mode. | All `estimate_dcbm()` arguments; `self_loops`: retain/zero the output diagonal; `lower_clip`, `upper_clip`: finite output bounds. | Square probability matrix `P_hat` for the estimated nodes. |
| `label_match_greedy()` | Relabels estimated communities to agree with a standard labeling. Identity and binary cases use exact shortcuts. The default deterministic greedy method repeatedly takes the largest remaining overlap; the optional dependency-free Hungarian method finds a globally optimal assignment. The audited mapping direction is source label to standard label. | `match_this`: labels to transform; `standard`: target labels; `K`: label range; `algorithm`: `"greedy"` (default) or `"hungarian"`; `return_mapping`: include assignment diagnostics. | List with `matched_labels` and proportional `mismatch_rate`; when requested, also source-to-standard `mapping`, achieved `agreement`, and overlap matrix. |
| `label_match_brute_force()` | Checks all `K!` source-to-standard label permutations, generated recursively without storing a factorial-size permutation matrix. For `K > 8`, interactive use prompts before proceeding; non-interactive use stops unless `confirm_large = TRUE`. Prefer `label_match_greedy()` for large `K`. | `match_this`, `standard`, `K`: matching inputs; `return_mapping`: include diagnostics; `confirm_large`: `NULL` prompts interactively, `TRUE` explicitly continues, and `FALSE` cancels when `K > 8`. | List with `matched_labels` and proportional `mismatch_rate`; when requested, also mapping, agreement, overlap table, and number of permutations evaluated. |

## `x1_sonnet.R`

Location: [`x1_sonnet.R`](./x1_sonnet.R)

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `op_splitter()` | Repeatedly samples a common overlap and partitions the remaining nodes. When `s * m < n - o`, the default adds every remainder node to the common overlap; alternatively, extras are balanced across randomly selected non-overlap pieces or ignored with a warning. Seed offsets are deterministic and overflow-safe. | `n`: nodes; `s`: splits per repetition; `o`: initial overlap size; `n_repetitions`: repeated partitions; `m`: base non-overlap nodes per split; `remainder_handling`: `"augment_overlap"` (default), `"distribute_randomly"`, or discouraged `"ignore"`; `seed`: optional non-negative integer seed. | List with repeated `splits`, effective overlap size `o`, `remainder_count`, and selected `remainder_handling`. |
| `sonnet_splitter()` | SONNET-specific splitter with one overlap shared across all repetitions. Remainders always augment that overlap. First-repetition subnetworks include the overlap; extra repetitions freshly permute and split only the non-overlap nodes, preserving the original SONNET flow. | `n`; `num_subnetworks`; `overlap_size`; `extra_nrep`; `m`; `seed`. | List with per-repetition `subnetworks`, overlap/non-overlap indices, effective/requested overlap, remainder, and piece size. |
| `sonnet_splitter_independent()` | Samples a separate augmented overlap and complete overlap-plus-partition collection for every repetition. Unlike `sonnet_splitter()`, each repetition covers all `n` nodes. | `n`; `num_subnetworks`; `overlap_size`; `extra_nrep`; `m`; `seed`. | List with per-repetition subnetworks, lists of overlap/non-overlap indices, effective/requested overlap, remainder, and piece size. |
| `.sonnet_fit()` | Internal engine for both SONNET variants. Before parallel clustering, it prints an explicit additive formula for the largest scheduled decomposition, counted Laplacian operations, a full adjacency-copy allowance per worker, and retained label/membership storage. Each concurrent term uses the actual smart worker count. Laplacian input is estimated as dense. | Common SONNET arguments plus `force_windows`, `ram_check = TRUE`, and internal `share_overlap`; named `...` is forwarded to `spectral_cluster()`. | Full classed `sonnet` object. |
| `sonnet_shared_overlap()` | Explicit entry point for the historical shared-overlap SONNET variant. | `...`: arguments accepted by the SONNET engine. | Classed `sonnet` object with `share_overlap = TRUE`. |
| `sonnet_independent_overlap()` | Explicit entry point for independent overlaps. Every repetition yields all `n` labels and every extra labeling is globally matched to repetition zero before voting. | `...`: arguments accepted by the SONNET engine. | Classed `sonnet` object with `share_overlap = FALSE`. |
| `sonnet()` | Wrapper selecting the shared or independent SONNET flow. The default uses independent overlaps. Uses `spectral_cluster()` and strict `uni_mclapply()` worker-error propagation; spectral defaults, including CLARA, are retained. | `A`; `K`; `num_subnetworks >= 2`; positive `overlap_size`; `extra_nrep`; `ncores`; `seed`; `matching_method`; `confirm_large`; `verbose`; `force_windows = FALSE`; `ram_check = TRUE`: print the explicitly factored decomposition, matrix-operation, worker-copy, and retained-result estimate and warn when its conservative sum exceeds available RAM; `share_overlap = FALSE`; named `...` forwarded to `spectral_cluster()`. | Classed `sonnet` object containing final labels, `extra_nrep + 1` proportional `n` by `K` membership matrices, complete repetition labels, raw/aligned subnetwork labels, alignment diagnostics, splits, parameters, worker counts, timings, and call. |
| `print.sonnet()` | Compact print method for a fitted SONNET object. | `x`: `sonnet` object; `...`: method compatibility. | Input object invisibly. |
| `summary.sonnet()` | Builds SONNET diagnostics including community sizes, overlap augmentation, membership certainty, cores, and timings. | `object`: `sonnet` object; `...`: method compatibility. | Classed `summary.sonnet` object. |
| `print.summary.sonnet()` | Prints a SONNET summary. | `x`: `summary.sonnet`; `...`: method compatibility. | Input summary invisibly. |

## `xx_helpers.R`

Location: [`xx_helpers.R`](./xx_helpers.R)

Source `x1_sonnet.R` first.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `netcrop_splitter()` | Shared NETCROP splitter that generates a fresh overlap for every repetition and partitions all remaining nodes into equal non-overlap pieces. Remainders always augment the overlap, so no nodes are dropped and the effective overlap is explicitly tracked. | `n`; `num_subnetworks >= 2`; `overlap_size`; `nrep = 1`; `seed`. | Subnetworks, repetition-specific overlap nodes and non-overlap pieces, requested/effective overlap, remainder count, piece size, and repetition metadata. |

## `x2_netcrop_bm.R`

Location: [`x2_netcrop_bm.R`](./x2_netcrop_bm.R)

Source the numbered helpers, `x1_sonnet.R`, and `xx_helpers.R` first. The core algorithm uses
base R plus the dependencies already required by the selected decomposition,
clustering, and matrix representation. Plotting optionally requires `ggplot2`.
Periods in the S3 method names below are required by R method dispatch and are
external-interface exceptions to the underscore naming rule.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `netcrop_blockmodel()` | Selects SBM/DCBM and `K` by the original three-stage overlapping-subnetwork flow. The experimental engine preserves the model/loss calculations while computing all candidates for one subnetwork in one worker task, extracting that subgraph once, bypassing clustering/estimation wrappers for `K = 1`, reusing DCBM row norms, using direct task-index tables, caching numeric loss inputs, returning compact worker records, and skipping validations already completed by this high-level function. On Windows all stages reuse one multisession pool and preload `Matrix` when required. | Core NETCROP arguments; nested estimator/clustering options; matching controls; `ncores`; `seed`; `verbose`; `force_windows`; `ram_check`; `retain_intermediates = c("all", "minimal")`. CLARA defaults to `samples = 5`; engine-specific options are filtered before dispatch. | Classed `netcrop_blockmodel` object. With `retain_intermediates = "minimal"`, large raw, estimation, averaged-estimate, and detailed matching objects are omitted while selections, CV losses, matching summary, timing, options, and diagnostics remain. |
| `print.netcrop_blockmodel()` | Prints the overall selected model and `K` for every loss. | `x`; `...`: S3 compatibility. | Input object invisibly. |
| `summary.netcrop_blockmodel()` | Builds a durable summary with candidates, effective overlap, held-out-pair coverage, selections, underidentified-match count, workers, and timings. | `object`; `...`: S3 compatibility. | Classed `summary.netcrop_blockmodel` object. |
| `print.summary.netcrop_blockmodel()` | Prints a NETCROP block-model summary and the per-repetition and overall selections. | `x`; `...`: S3 compatibility. | Input summary invisibly. |
| `plot.netcrop_blockmodel()` | Plots candidate CV loss against `K`. The default aggregates repetitions using mean curves and one-standard-deviation ribbons; the unaggregated view distinguishes repetitions. Irregular and singleton `K` sets are supported. | `x`; `aggregate = TRUE`; `...`: S3 compatibility. | `ggplot` object; requires `ggplot2`. |

## `x3_netcrop_rdpg.R`

Location: [`x3_netcrop_rdpg.R`](./x3_netcrop_rdpg.R)

Source the numbered helpers, `x1_sonnet.R`, and `xx_helpers.R` first. This
variant supports symmetric, non-negative, loop-free adjacency matrices.

| Function | Description and use | Arguments | Returns |
|---|---|---|---|
| `netcrop_rdpg()` | Selects a symmetric-RDPG latent dimension using the original four-stage flow: magnitude-ordered subnetwork eigendecompositions of unscaled `A`, candidate embeddings, overlap Procrustes alignment, and pairwise held-out losses. Predictions stay on the same unscaled adjacency scale, so no extra `rho_hat` multiplier is applied; `rho_hat` is diagnostic only. The experimental engine replaces repeated grid scans with direct index arrays, prevalidates losses once, caches numeric blocks, emits compact worker records, skips repeated Procrustes/loss validation, supports sparse diagonal/degree access, and reuses one Windows worker pool with `Matrix` preloaded. | `A`; `d_candidates`; `num_subnetworks`; `overlap_size`; `nrep = 1`; `losses`; `eig_options`; `ncores`; `seed`; `verbose`; `force_windows`; `ram_check`; `retain_intermediates = c("all", "minimal")`. | Classed `netcrop_rdpg` object. Minimal retention omits raw decompositions and candidate/aligned embeddings but keeps CV results, selections, splitter, diagnostics, timing, and options. |
| `print.netcrop_rdpg()` | Prints the overall selected dimension for every loss. | `x`; `...`: S3 compatibility. | Input object invisibly. |
| `summary.netcrop_rdpg()` | Builds a durable summary with candidates, effective overlap, held-out coverage, sparsity estimate, selections, negative retained eigenvalue count, workers, and timings. | `object`; `...`: S3 compatibility. | Classed `summary.netcrop_rdpg` object. |
| `print.summary.netcrop_rdpg()` | Prints the RDPG NETCROP summary and per-repetition and overall dimension selections. | `x`; `...`: S3 compatibility. | Input summary invisibly. |
| `plot.netcrop_rdpg()` | Plots CV loss against `d`; the default aggregates repetitions with mean curves and one-standard-deviation ribbons, while the unaggregated view distinguishes repetitions. Arbitrary and singleton candidate sets are supported. | `x`; `aggregate = TRUE`; `...`: S3 compatibility. | `ggplot` object; requires `ggplot2`. |

## Maintenance checklist

Whenever code changes:

1. Update the file tree if a file was added, removed, or renamed.
2. Add, remove, or revise the corresponding function row.
3. Keep argument names/default behavior and return descriptions synchronized.
4. Mark internal helpers and compiled exports clearly.
5. Record new optional package dependencies and fallback behavior.
6. Check that every top-level R function and every named C++ function appears here.
