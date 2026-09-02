# Naming Convention

## Mandatory rule

All variable, function, argument, method, class, constant, data-field, and other code identifiers must use an underscore (`_`) as the separator between words.

Do not use a period (`.`) as a word separator in any identifier. This rule applies to all existing and future code so that names can be transferred to and used directly in Python code.

## Examples

| Do not use | Use instead |
|---|---|
| `old.spectral.cluster` | `old_spectral_cluster` |
| `label.match` | `label_match` |
| `make.basis` | `make_basis` |
| `tmp.std` | `tmp_std` |
| `cut.over` | `cut_over` |
| `n.core` | `n_core` |

## Scope

This convention applies to:

- Variable and object names
- Function and method names
- Function arguments and parameters
- Class and module names
- Constants
- List keys, named-vector elements, and internally created data fields
- File-derived identifiers and any other names created in code

Names required by an external package, API, file format, or dataset may retain their original spelling at the external boundary. Convert them to underscore-separated names before using them internally whenever practical.

## Living function dictionary

`dictionary.md` is the mandatory living reference for this directory. Every top-level named R function, named C++ helper, and exported compiled function must be listed there with its description, argument details, intended use, return value, and source location.

Any change that adds, removes, renames, moves, or modifies a function must update `dictionary.md` in the same task. Changes to arguments, defaults, return structure, dependencies, fallback behavior, or intended use also require a corresponding dictionary update. The file tree in `dictionary.md` must remain synchronized with the directory contents.

Anonymous and lexically nested implementation closures may be documented through their owning top-level function rather than as separate callable APIs.

## Universal names

Use the following names consistently throughout the codebase. Preserve the capitalization shown here.

| Name | Meaning |
|---|---|
| `A` | Adjacency matrix |
| `P` | Edge-probability or block-probability matrix |
| `L` | Graph Laplacian or normalized adjacency/Laplacian matrix |
| `K` | Number of communities |
| `n` | Number of nodes |
| `d` | Latent dimension of an RDPG |
| `Z` | Latent-position matrix for an undirected latent-position model |
| `Z_left` | Source/left latent-position matrix for a directed model |
| `Z_right` | Target/right latent-position matrix for a directed model |
| `alpha` | Latent-space-model intercept or node-intercept vector |
| `deg` | Node-degree vector |
| `average_degree` | Target or observed average node degree |
| `ncores` | Number of processor cores |
| `g` | Vector of estimated or observed community labels |
| `g_true` | Ground-truth community-label vector |
| `g_one_hot` | n x K matrix of one-hot community matrix |
| `nreps` | Number of repetitions of an algorithm |
| `seed` | Random-number-generator seed |
| `tau` | Regularization parameter |
| `nsim` | Number of simulations |
| `nboot` | Number of bootstraps or resamples or subsamples|
| `psi` | Degree correction parameters |
| `psi_one_hot` | n x K degree parameter matrix |

Do not introduce synonyms for these concepts. For example, use `K` rather than alternating among `k`, `n_communities`, and `num_clusters`, and use `ncores` rather than `ncore`, `n_cores`, or `cores`.

### Default number of cores

Unless a function has a documented reason to behave differently, its default value for `ncores` must be half of the detected logical cores, rounded down:

```r
ncores = max(
  floor(parallel::detectCores() / 2),
  1L,
  na.rm = TRUE
)
```

The safeguards are mandatory:

- `floor()` ensures that the worker count is an integer.
- `max(..., 1L)` ensures that at least one core is used.
- `na.rm = TRUE` ensures that a failed `parallel::detectCores()` call, which may return `NA`, produces the safe default of one core.
- Functions accepting `ncores` must also validate user-supplied values and reject missing, non-finite, non-integer, or non-positive values.

## Estimated quantities

Append `_hat` directly to the base variable name for an estimated version of a mathematical quantity:

| Quantity | Estimate |
|---|---|
| `P` | `P_hat` |
| `A` | `A_hat` |
| `L` | `L_hat` |
| `deg` | `deg_hat` |
| `labels` | `labels_hat` |

Additional descriptors and cases must follow `_hat`, such as `P_hat_regularized` or `P_hat_from_sub1`.

## Composite names

Construct composite identifiers in this order:

```text
[varname]_[descriptor]_[case]
```

- `varname` identifies the underlying quantity.
- `descriptor` states its transformation, estimator, role, or property.
- `case` identifies the source, dataset, method, scenario, or indexed instance.

Use only the components needed to make a name unambiguous. Do not add empty, redundant, or repeated components.

Examples:

| Meaning | Name |
|---|---|
| Estimated probability matrix from subgraph 1 | `P_hat_from_sub1` |
| Regularized Laplacian for the full graph | `L_regularized_full` |
| Degree vector for subgraph 2 | `deg_sub2` |
| Estimated labels from SONNET | `labels_hat_sonnet` |
| Adjacency matrix for the training sample | `A_train` |
| Probability estimate for the DCBM case | `P_hat_dcbm` |

Keep common mathematical suffixes in a stable position:

```text
[varname]_hat_[transformation]_[source_or_case]
```

For example, use `P_hat_regularized_sub1`, not `regularized_sub1_P_hat`.

## Additional conventions

- Use `snake_case` for descriptive identifiers and retain uppercase only for the canonical mathematical objects listed above.
- Use nouns for data objects (`degree_weights`) and verbs for functions (`estimate_probability`, `match_labels`).
- Prefix Boolean variables with `is_`, `has_`, `use_`, or `should_`, such as `is_sparse` or `use_regularization`.
- Use plural names for collections and singular names for individual elements, such as `subgraphs` and `subgraph`.
- Prefix counts with `n`, using the established compact forms where specified: `ncores`, `nreps`, `n_nodes`, and `n_edges`.
- Use descriptive loop indices when their role is known (`node_idx`, `subgraph_idx`, `rep_idx`). Reserve `i`, `j`, and `r` for short, local mathematical loops.
- Use `_idx` for a single index, `_indices` for multiple indices, and `_mask` for a logical selection vector.
- Use `_mat`, `_vec`, or `_list` only when the type is not otherwise evident and the suffix prevents ambiguity.
- Use established abbreviations consistently. Prefer `subgraph` in public names; use compact forms such as `sub1` only for clearly indexed cases.

## RAM preflight convention

RAM checks belong in high-level algorithms, not in low-level decomposition or
matrix-operation helpers. A preflight message must label each conservative
term explicitly in this form:

```text
(per-operation RAM x sequential operations) x parallel operations + ...
```

Use the actual worker count for parallel terms. Count repeated matrix products
using the largest relevant product estimate, and treat Laplacian construction
as dense for estimation even when the implementation stores or processes a
sparse matrix. Prefer a documented overestimate to a possible underestimate.
The message becomes a warning when the summed estimate exceeds available RAM
or available RAM cannot be determined.

### Sparse-matrix preservation

When an algorithm receives a strict sparse adjacency matrix, its default path
must preserve sparse adjacency-derived matrices. Do not evaluate `A + 0`,
`A + tau * value`, or an equivalent scalar addition unless `tau > 0`; positive
all-pairs scalar regularization is the explicit case where densification is
expected. Use sparse diagonal scaling instead of materializing a dense outer
product for normalized adjacency or Laplacian construction. A high-level
algorithm must not silently fall back from a sparse partial decomposition to a
dense base decomposition: it should use a sparse-capable engine or stop with an
informative error. Dense embeddings, probability matrices, and loss outputs
that are mathematically dense are not adjacency-storage representations and
are outside this preservation rule.
- Avoid names that shadow common functions or language keywords, such as `mean`, `sum`, `list`, `function`, `class`, or `return`.
- Keep identical concepts identically named across R and Python. Do not translate a name merely because it appears in a different language.

## Compiled-code fallbacks

Every function accelerated by C++, Rcpp, or another compiled extension must have a functionally equivalent implementation in base R or in non-compiled R package code.

The R-facing wrapper must:

- Work when the compiled package, compiler, shared library, or exported function is unavailable.
- Check whether the compiled function exists before calling it.
- Catch failures raised while calling the compiled function and continue with the R fallback.
- Validate inputs before selecting either implementation so both paths accept and reject the same inputs.
- Return the same type, shape, missing-value behavior, and numerical interpretation from both paths.
- Preserve a `use_cpp`-style option when users may need to select or test the R implementation explicitly.
- Avoid automatically installing compilers or packages. Missing compiled dependencies must not prevent the R file from being sourced.
- Emit a concise diagnostic when a loaded compiled implementation fails at runtime; absence of an optional compiled implementation may fall back silently.
- Include parity tests covering ordinary inputs, boundary cases, missing or invalid inputs, and ties where relevant.

Compiled functions should use the same base name as their R wrapper with the `_cpp` suffix, such as `auc` and `auc_cpp`.

## Matrix operations

Use the following order of preference for matrix computations:

1. Use vector operations, broadcasting, or row/column scaling when they express the same computation clearly and do not require significantly more peak memory.
2. Use `crossprod()` and `tcrossprod()` instead of `%*%` whenever the required transposition and multiplication match their semantics.
3. Use `%*%` only when neither vector operations nor `crossprod()`/`tcrossprod()` express the operation correctly and clearly.

Additional requirements:

- Do not construct a dense diagonal matrix merely to scale matrix rows or columns. Use vectorized row/column scaling for dense matrices.
- Preserve sparse matrix structure. Sparse diagonal multiplication may be preferable to vector operations when vectorization would densify the matrix or require a large dense intermediate.
- Do not construct undirected adjacency matrices with `Matrix::sparseMatrix(..., symmetric = TRUE)` or coerce them to symmetric sparse-matrix classes. These representations can interact poorly with downstream eigensolvers such as RSpectra.
- Construct an undirected adjacency matrix from one off-diagonal triangle with `A <- A + t(A)`. Add diagonal self-loop entries separately afterward so their weights are not doubled. Sparse implementations may use the namespace-explicit equivalent `A <- A + Matrix::t(A)` when the Matrix package is loaded but not attached.
- Do not construct a dense outer-product matrix solely to perform elementwise row-and-column scaling when sequential scaling uses less peak memory.
- Use `crossprod(A, B)` for `t(A) %*% B` and `tcrossprod(A, B)` for `A %*% t(B)`.
- Confirm dimensions and orientation before replacing a product; performance preferences must never change the mathematical result.
- Prefer operations that avoid unnecessary transposes, copies, and intermediate matrices, especially for large network and probability matrices.
- Include dense and sparse parity tests when optimizing a matrix operation that supports both representations.
