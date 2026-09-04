#' Resource estimation and reproducible execution utilities
#'
#' Utilities for cross-platform mapping, conservative RAM estimation and
#' reporting, and reproducible simulation studies. Parallel helpers restore
#' caller state and high-level RAM checks report conservative estimates rather
#' than reserving memory.
#'
#' @param X A collection of tasks, or a numeric value to format where noted.
#' @param FUN A function applied to each element of `X`.
#' @param ... Additional arguments passed to a supplied function.
#' @param ncores,ncores_outer Positive worker counts.
#' @param force_windows Use the Windows-compatible future backend.
#' @param stop_on_error Stop instead of falling back or returning failures.
#' @param manage_future_plan Whether the helper creates and restores a future plan.
#' @param future_packages Character vector of packages required by workers.
#' @param suppress_worker_renv_sync_check Suppress repeated worker startup messages.
#' @param x A value to format.
#' @param estimated_bytes A nonnegative byte estimate.
#' @param max_fraction Maximum fraction of available RAM considered usable.
#' @param reserve_bytes Bytes reserved for the operating system and other work.
#' @param operation,detail Descriptive text included in reports.
#' @param n,p Matrix dimensions.
#' @param K,ncv Requested rank and Krylov subspace size.
#' @param safety_factor Multiplicative conservatism factor.
#' @param input_already_counted Whether input allocation was counted elsewhere.
#' @param nu,nv Numbers of left and right singular vectors.
#' @param nrow_left,shared_dimension,ncol_right Matrix-product dimensions.
#' @param method,engine,force_engine,safe_d_multiplier,dense_input Decomposition controls.
#' @param operation_count Number of sequential operations.
#' @param terms Named RAM formula terms.
#' @param process An expression or function to measure.
#' @param one_simulation Function implementing one simulation.
#' @param nsim Number of simulations.
#' @param use_parallel_simulations Whether to parallelize simulations.
#' @param seed,seeds A base seed or explicit per-simulation seeds.
#' @param results_file Optional resumable RDS path. Resume only files created
#'   by or otherwise trusted by the caller; RDS is not an untrusted interchange
#'   format.
#' @param action Replace, resume, or archive existing results.
#' @param retry_failed,continue_on_error,measure_resources,show_progress Simulation controls.
#' @return The utility-specific scalar, diagnostic list, mapped list, or
#'   simulation records described in each function's usage.
#' @examples
#' format_bytes(1024)
#' estimate_dense_eigen_ram(20)
#' @name resource_utilities
NULL

#' Mathematical losses, transforms, alignment, and decompositions
#'
#' General-purpose numerical helpers used by netOP algorithms. Compiled fast
#' paths validate inputs and fall back to the corresponding R implementation
#' when unavailable or unsuccessful.
#'
#' @param x,y Numeric vectors or objects operated on by a scalar loss.
#' @param na_rm Remove missing comparisons.
#' @param validate_inputs Validate inputs; only audited callers should disable this.
#' @param g_1,g_2 Label vectors.
#' @param lower_clip,upper_clip Inclusive clipping bounds.
#' @param P Probabilities or prediction scores.
#' @param eps,epsilon Probability clipping constants.
#' @param A A matrix or binary response vector.
#' @param use_cpp Try the registered compiled accelerator.
#' @param temperature Positive soft-selection temperature.
#' @param operator Outer operation; currently `"+"`.
#' @param X,X_star Matrices to align.
#' @param translate,dilate,sumsq Procrustes output controls.
#' @param normalized Whether to construct a normalized graph Laplacian.
#' @param tau Nonnegative graph regularization strength.
#' @param d Requested decomposition or reconstruction rank.
#' @param only_values Omit decomposition vectors.
#' @param scale_by Matrix scaling convention.
#' @param use_laplacian Decompose a graph Laplacian.
#' @param engine,force_engine Decomposition backend and fallback control.
#' @param order_by Eigenvalue ordering convention.
#' @param safe_d_multiplier Extra partial components computed before truncation.
#' @param nu,nv Requested left and right singular vectors.
#' @return A scalar loss, transformed object, alignment result, graph matrix,
#'   decomposition, or reconstructed matrix, as appropriate.
#' @examples
#' sse(c(1, 2), c(1, 3))
#' auc(c(0, 1, 1), c(0.1, 0.7, 0.8), use_cpp = FALSE)
#' eig_decomp(diag(c(3, 2, 1)), d = 2, engine = "base")
#' @name math_utilities
NULL

#' Network graph utilities
#'
#' Dependency-light adjacency-list, connected-component, and shortest-path
#' operations for dense and sparse network matrices.
#'
#' @param A A finite square adjacency matrix.
#' @param directed Preserve edge direction.
#' @param self_loops Include or ignore diagonal edges.
#' @param use_cpp Try a registered compiled accelerator before the R fallback.
#' @param sort_nodes Sort returned component indices.
#' @param weighted Treat nonzero entries as nonnegative path lengths.
#' @return An adjacency list, component result, induced subnetwork, or distance matrix.
#' @examples
#' A <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), 3, 3)
#' connected_components(A, use_cpp = FALSE)
#' shortest_path_distances(A, use_cpp = FALSE)
#' @name network_utilities
NULL

#' Network input, generation, and probability calibration
#'
#' Read and write edge lists and generate ER, SBM, DCBM, RDPG, and latent-space
#' networks. Generators return an adjacency matrix with compact truth metadata
#' available through [get_generator_parameters()].
#'
#' @param A A network matrix.
#' @param tolerance Numerical comparison or calibration tolerance.
#' @param file Edge-list path.
#' @param directed,weighted Direction and edge-weight controls.
#' @param triangle Triangle written for an undirected network.
#' @param format File format (`"auto"`, `"csv"`, or `"tsv"`).
#' @param include_header,has_header Header writing or reading controls.
#' @param overwrite Permit [write_network()] to replace an existing file.
#' @param representation Dense or sparse output representation.
#' @param n,d,K Numbers of nodes, latent dimensions, and communities.
#' @param P,P_block Edge-probability or block-probability matrix.
#' @param p ER edge probability.
#' @param average_degree,average_degree_method Expected-degree target and calibration method.
#' @param self_loops Include diagonal edges.
#' @param lower_clip,upper_clip Probability bounds.
#' @param max_iterations,naive_iterations Calibration iteration limits.
#' @param seed,ncores Reproducibility seed and worker count.
#' @param validate_inputs Validate direct probability input.
#' @param Z,Z_left,Z_right Latent-position matrices.
#' @param latent_distribution,latent_args Position distribution and named arguments.
#' @param sparsity_multiplier,scale_P Probability sparsity and scaling controls.
#' @param g_true,community_probabilities Labels or community proportions.
#' @param shape_1,shape_2 Beta distribution shapes.
#' @param psi,degree_distribution,degree_args,degree_scale Degree-correction controls.
#' @param alpha,beta Block or latent-space intercept parameters.
#' @param lower_bound,upper_bound Truncation bounds.
#' @param mean_lower,mean_upper,noise_lower,noise_upper Latent-mixture bounds.
#' @param theta A latent-space logit matrix.
#' @param normalize_Z,distance_adjustment Latent-space normalization controls.
#' @return A generated adjacency matrix, edge-list path, parsed network,
#'   calibrated probabilities, latent values, or generator metadata.
#' @examples
#' A <- generate_er(n = 200, average_degree = 5, seed = 1, ncores = 1)
#' get_generator_parameters(A)
#' generate_sbm(n = 200, K = 3, alpha = 0.4, beta = 0.1,
#'              seed = 2, ncores = 1)
#' @name network_generators
NULL

#' Network embedding, clustering, estimation, and label alignment
#'
#' Spectral and latent-space embedders, SBM/DCBM probability estimators, and
#' deterministic or exact label-alignment helpers.
#'
#' @param A A network or rectangular training matrix.
#' @param d,K Embedding dimension or number of communities.
#' @param directed,self_loops Direction and diagonal controls.
#' @param reconstruct Return a rank-reconstructed adjacency matrix.
#' @param align_with Optional reference embedding.
#' @param ram_check Report a conservative high-level RAM estimate.
#' @param U Optional supplied spectral representation.
#' @param laplacian,normalize_laplacian,regularize_tau Spectral graph controls.
#' @param handle_zero_degree_nodes Zero-degree node policy.
#' @param row_normalize Normalize representation rows before clustering.
#' @param spectral_method,spectral_engine,spectral_options Spectral controls.
#' @param cluster_engine,cluster_options Clustering backend and named options.
#' @param validate_inputs Validate inputs; only audited callers should disable this.
#' @param step_size,niter,trace Latent-space optimizer controls.
#' @param Z_init,alpha_init Optional optimizer initial values.
#' @param epsilon Probability clipping constant.
#' @param use_cpp Try the compiled latent-space optimizer.
#' @param g Community labels.
#' @param fold_nodes Original NCV fold-node indices.
#' @param method DCBM estimation method.
#' @param row_norm,psi_omit,stabilizer DCBM estimation controls.
#' @param lower_clip,upper_clip Probability bounds.
#' @param match_this,standard Label vectors to align and target.
#' @param algorithm Greedy or Hungarian assignment algorithm.
#' @param return_mapping Include mapping diagnostics.
#' @param confirm_large Confirm factorial exact matching for large `K`.
#' @return An embedding, fitted clustering/model, reconstructed probability
#'   matrix, or label-matching result.
#' @examples
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.5, beta = 0.1,
#'                   representation = "dense", seed = 3, ncores = 1)
#' spectral_cluster(A, K = 3, spectral_engine = "base",
#'                  cluster_engine = "kmeans")
#' label_match_greedy(c(2, 2, 1, 1), c(1, 1, 2, 2), K = 2)
#' @name embedding_estimating
NULL

#' SONNET, NETCROP, and regularization model selection
#'
#' High-level scalable fitting and model-selection workflows. Related public
#' functions place the primary network and candidate set first, followed by
#' structural, resampling, loss, computational, seed, diagnostic, and memory
#' controls. Plotting methods require the suggested `ggplot2` package.
#'
#' NETCROP refers to NETwork CRoss-Validation using Overlapping Partitions.
#'
#' @param n,K Network size or fixed community count.
#' @param theta,q,s_lower,s_upper,o_lower,o_upper SONNET parameter-grid controls.
#' @param ncores Positive worker count.
#' @param test_prop,o_range NETCROP automatic partition controls.
#' @param ... Named arguments forwarded deliberately, or S3 method arguments.
#' @param A Primary adjacency matrix, or a matched list where documented.
#' @param num_subnetworks,overlap_size,extra_nrep Partition and repetition controls.
#' @param seed Reproducibility seed.
#' @param matching_method,confirm_large Label-matching controls.
#' @param verbose,force_windows,ram_check Diagnostic and backend controls.
#' @param share_overlap Select the historical shared-overlap SONNET variant.
#' @param parameter_select_options Named automatic-partition options.
#' @param K_candidates,d_candidates,tau_candidates Candidate sets.
#' @param nrep Number of model-selection repetitions.
#' @param losses Requested canonical losses.
#' @param model_candidates Candidate block-model families.
#' @param sbm_est_options,dcbm_est_options,eig_options,lsm_options Nested fit options.
#' @param retain_intermediates Retain all or minimal intermediate results.
#' @param g_1,g_2,g_3,g_4 Community labels used by pairwise losses.
#' @param use_dcbm,use_laplacian,dcbm_est_method Model and estimator choices.
#' @param loss_types Optional named legacy loss mapping.
#' @param label_reference Reference labeling used for label losses.
#' @param spectral_options,cluster_options,estimator_options Nested algorithm options.
#' @param failure_handling Stop or omit failed tasks.
#' @param laplacian,normalize_laplacian,handle_zero_degree_nodes Spectral controls.
#' @param row_normalize,spectral_method,spectral_engine,cluster_engine Backend controls.
#' @param retain_fits Retain fitted candidate models.
#' @param g_true Ground-truth labels for oracle comparison.
#' @param netcrop_outcomes,dkest_outcomes Optional tuner results.
#' @param include_netcrop_mean,include_netcrop_mode Include NETCROP summaries.
#' @param engines Clustering engines evaluated by the oracle plotter.
#' @param sonnet_options,spectral_cluster_options Named fitter options.
#' @param x,object An object handled by an S3 method.
#' @param aggregate Aggregate repetitions in a plot.
#' @return A parameter recommendation, fitted SONNET object, model-selection
#'   result, summary, plot, pairwise loss, or multi-regularization result.
#' @references
#' Chakrabarty, S., Sengupta, S., and Chen, Y. (2026). Network
#' Cross-Validation and Model Selection via Subsampling. arXiv:2504.06903.
#' \doi{10.48550/arXiv.2504.06903}
#' @examples
#' netcrop_param_select(test_prop = 0.05, n = 200, o_range = 0)
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.5, beta = 0.1,
#'                   seed = 4, ncores = 1)
#' fit <- sonnet(A, K = 3, num_subnetworks = 2, overlap_size = 20,
#'               ncores = 1, seed = 5, spectral_engine = "base")
#' fit
#' @name model_selection
NULL

#' Edge cross-validation stability selection
#'
#' Select a block-model community count or RDPG dimension by repeated edge
#' sampling. netOP provides self-contained wrappers around an ECV
#' implementation derived from the CRAN package `randnet`. Installing or using
#' netOP does not require `randnet`. The ECV-specific implementation helpers are
#' internal and are not part of the netOP public API.
#'
#' @param A A binary, symmetric, loop-free adjacency matrix.
#' @param max_K,max_d Maximum community count or latent dimension. The exact
#'   ECV implementation evaluates the full sequence from one to this maximum.
#' @param train_proportion Proportion of dyads retained for training.
#' @param cv Number of edge-sampling folds.
#' @param nrep Number of repeated ECV runs.
#' @param tau Nonnegative block-model regularization.
#' @param losses Canonical loss names.
#' @param ncores Positive outer worker count; inner folds run sequentially.
#' @param seed Optional reproducibility seed.
#' @param verbose Print progress information.
#' @param force_windows Use the Windows-compatible parallel backend.
#' @param ram_check Report conservative RAM demand.
#' @param failure_handling Stop or omit failed repetitions.
#' @param retain_intermediates Retain all or minimal legacy results.
#' @return A `netcrop_blockmodel` or `netcrop_rdpg` result with
#'   `algorithm = "ECV"`, tidy losses, selections, seeds, and diagnostics.
#' @references
#' Li, T., Levina, E., and Zhu, J. (2020). Network cross-validation by edge
#' sampling. *Biometrika*, 107(2), 257-276. \doi{10.1093/biomet/asaa006}
#' @examples
#' \donttest{
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.5, beta = 0.1,
#'                   seed = 6, ncores = 1)
#' ecv_stability_blockmodel(A, max_K = 5, cv = 2, nrep = 1,
#'                          ncores = 1, seed = 7, verbose = FALSE)
#' }
#' @name ecv_stability
NULL

#' Node cross-validation stability selection
#'
#' Repeated node cross-validation for selecting the number of block-model
#' communities. This is the netOP authors' implementation of the exact
#' algorithm described in Chen and Lei (2018). Numerical-stability measures
#' and failsafes were added without altering the algorithm in the paper.
#' Probabilities are clipped before logarithmic losses, zero-degree cases are
#' checked before normalization, non-finite candidate results are audited, and
#' failed repetitions can either stop or be omitted explicitly.
#'
#' @param A A binary, symmetric, loop-free adjacency matrix.
#' @param max_K Maximum community count; candidates are `1:max_K`.
#' @param cv,nrep Numbers of folds and repetitions.
#' @param dc_est DCBM estimator (`"spectral"` or `"plugin"`).
#' @param tau,use_laplacian Spectral regularization controls.
#' @param losses Canonical loss names.
#' @param ncores Positive outer worker count.
#' @param seed Optional reproducibility seed.
#' @param verbose,force_windows,ram_check Diagnostic and backend controls.
#' @param failure_handling Stop or omit failed repetitions.
#' @param retain_intermediates Retain all or minimal repetition output.
#' @return A `netcrop_blockmodel` result with `algorithm = "NCV"`, tidy
#'   losses, selections, fold metadata, seeds, and diagnostics.
#' @references
#' Chen, K. and Lei, J. (2018). Network Cross-Validation for Determining the
#' Number of Communities in Network Data. *Journal of the American Statistical
#' Association*, 113(521), 241-251. \doi{10.1080/01621459.2016.1246365}
#' @examples
#' \donttest{
#' A <- generate_sbm(n = 200, K = 3, alpha = 0.55, beta = 0.08,
#'                   seed = 24, ncores = 1)
#' ncv_stability_blockmodel(A, max_K = 5, cv = 2, nrep = 1,
#'                          ncores = 1, seed = 26, verbose = FALSE,
#'                          losses = "sse")
#' }
#' @name ncv_stability
NULL
