# This file supplies concise topic-level metadata for public functions whose
# implementation comments remain deliberately close to the algorithm source.

#' Apply a function consistently across platforms
#'
#' Applies a function sequentially or with supported parallel workers while
#' preserving a consistent return structure across operating systems.
#' @inheritParams resource_utilities
#' @param ncores Positive worker count.
#' @return A list containing one result per input element.
#' @seealso [run_simulations()], [sonnet()]
#' @name uni_mclapply
NULL

#' Measure peak memory use
#'
#' Measures peak memory consumption for a supplied expression or function.
#' @inheritParams resource_utilities
#' @return A memory-measurement result in bytes with execution metadata.
#' @seealso [available_ram()], [report_ram_preflight()]
#' @name measure_peak_ram
NULL

#' Run reproducible simulation jobs
#'
#' Executes, saves, and optionally resumes a collection of reproducible
#' simulation jobs.
#' @inheritParams resource_utilities
#' @param ncores_outer Positive outer worker count.
#' @return A list of simulation results and execution metadata.
#' @seealso [uni_mclapply()], [get_generator_parameters()]
#' @name run_simulations
NULL

#' Compute the statistical mode
#'
#' Returns the most frequent non-missing value in a vector.
#' @inheritParams math_utilities
#' @param x Input vector.
#' @return A scalar value with the same basic type as the input.
#' @seealso [nmi()]
#' @name modal
NULL

#' Compute normalized mutual information
#'
#' Measures agreement between two community-label vectors using normalized
#' mutual information.
#' @inheritParams math_utilities
#' @return A numeric agreement score.
#' @seealso [label_match_greedy()], [pair_nmi_loss()]
#' @name nmi
NULL

#' Clip numeric values to an interval
#'
#' Restricts numeric values to supplied lower and upper bounds.
#' @inheritParams math_utilities
#' @param x Numeric vector or matrix to clip.
#' @return An object shaped like the input with clipped values.
#' @seealso [clip_probabilities()]
#' @name clip_values
NULL

#' Clip probability values safely
#'
#' Restricts probabilities away from zero and one for stable likelihood
#' calculations.
#' @inheritParams math_utilities
#' @param eps Positive probability-clipping constant.
#' @return An object shaped like the input with clipped probabilities.
#' @seealso [clip_values()], [bin_dev()]
#' @name clip_probabilities
NULL

#' Compute the softplus transform
#'
#' Evaluates the numerically stable softplus transform elementwise.
#' @param x Numeric vector or matrix.
#' @return A numeric vector or matrix matching the input shape.
#' @seealso [sigmoid()]
#' @name softplus
NULL

#' Compute the sigmoid transform
#'
#' Evaluates the numerically stable logistic sigmoid elementwise.
#' @param x Numeric vector or matrix.
#' @return A numeric vector or matrix matching the input shape.
#' @seealso [softplus()]
#' @name sigmoid
NULL

#' Add vectors by outer expansion
#'
#' Forms the matrix of all pairwise sums of two numeric vectors.
#' @inheritParams math_utilities
#' @return A numeric matrix.
#' @seealso [procrustes()]
#' @name outer_add
NULL

#' Align point configurations by Procrustes transformation
#'
#' Aligns two point configurations using an orthogonal transformation and
#' optional translation.
#' @inheritParams math_utilities
#' @return An aligned configuration with transformation details.
#' @seealso [ase()], [generate_latent_positions()]
#' @name procrustes
NULL

#' Compute matrix density
#'
#' Computes the proportion of nonzero entries in a dense or sparse matrix.
#' @inheritParams math_utilities
#' @return A numeric scalar between zero and one.
#' @seealso [generate_adjacency()]
#' @name matrix_density
NULL

#' Test matrix symmetry
#'
#' Tests whether a dense or sparse matrix is symmetric within tolerance.
#' @inheritParams math_utilities
#' @return A logical scalar.
#' @seealso [is_symmetric_network_matrix()]
#' @name is_symmetric_matrix
NULL

#' Construct a graph Laplacian
#'
#' Constructs an unnormalized or normalized graph Laplacian from a dense or
#' sparse adjacency matrix.
#' @inheritParams math_utilities
#' @return A dense or sparse Laplacian matrix.
#' @seealso [spectral_cluster()], [ase()]
#' @name graph_laplacian
NULL

#' Compute an eigendecomposition
#'
#' Computes a full or truncated eigendecomposition with a selected backend.
#' @inheritParams math_utilities
#' @return Eigenvalues, eigenvectors, and decomposition metadata.
#' @seealso [singular_decomp()], [ase()]
#' @name eig_decomp
NULL

#' Compute a singular-value decomposition
#'
#' Computes a full or truncated singular-value decomposition with a selected
#' backend.
#' @inheritParams math_utilities
#' @return Singular values, vectors, and decomposition metadata.
#' @seealso [eig_decomp()], [truncated_svd_reconstruct()]
#' @name singular_decomp
NULL

#' Reconstruct a truncated singular-value approximation
#'
#' Reconstructs a matrix from a selected number of singular components.
#' @inheritParams math_utilities
#' @return A numeric low-rank matrix approximation.
#' @seealso [singular_decomp()], [usvt()]
#' @name truncated_svd_reconstruct
NULL

#' Estimate a matrix with universal singular-value thresholding
#'
#' Applies universal singular-value thresholding and optional probability
#' clipping to a dense or sparse matrix.
#' @inheritParams math_utilities
#' @return A dense thresholded matrix estimate.
#' @seealso [singular_decomp()], [ase()]
#' @name usvt
NULL

#' Find adjacent nodes
#'
#' Builds unweighted or weighted adjacency lists from dense or sparse networks.
#' @inheritParams network_utilities
#' @return A list with one neighbor vector per node.
#' @seealso [connected_components()], [shortest_path_distances()]
#' @name adjacency_neighbors
NULL

#' Find connected components
#'
#' Finds every component or extracts the largest induced component of a dense
#' or sparse network.
#' @inheritParams network_utilities
#' @return Component memberships or an induced adjacency matrix, as documented.
#' @seealso [adjacency_neighbors()], [shortest_path_distances()]
#' @name connected_components
NULL

#' Compute shortest-path distances
#'
#' Computes all-pairs shortest-path distances in a dense or sparse network.
#' @inheritParams network_utilities
#' @return A numeric distance matrix.
#' @seealso [adjacency_neighbors()], [connected_components()]
#' @name shortest_path_distances
NULL

#' Compute an adjacency spectral embedding
#'
#' Embeds a dense or sparse adjacency matrix in `d` dimensions.
#' @inheritParams embedding_estimating
#' @param d Positive embedding dimension.
#' @param directed Whether to treat `A` as directed.
#' @return An embedding object containing coordinates and decomposition details.
#' @seealso [spectral_cluster()], [generate_rdpg()], [netcrop_rdpg()]
#' @name ase
NULL

#' Cluster a network spectrally
#'
#' Clusters a dense or sparse network into `K` communities using a configurable
#' spectral representation and clustering backend.
#' @inheritParams embedding_estimating
#' @param K Positive number of communities.
#' @return A fitted spectral-clustering object with labels and diagnostics.
#' @seealso [ase()], [estimate_sbm()], [sonnet()]
#' @name spectral_cluster
NULL

#' Fit a latent-space model by projected gradient descent
#'
#' Estimates latent positions and intercepts for a dense or sparse network.
#' @inheritParams embedding_estimating
#' @param d Positive latent dimension.
#' @return A fitted latent-space-model object.
#' @seealso [generate_lsm()], [netcrop_lsm()]
#' @name lsm_pgd
NULL

#' Estimate stochastic block models
#'
#' Estimates SBM or DCBM parameters from dense, sparse, full, or supported NCV
#' adjacency layouts.
#' @inheritParams embedding_estimating
#' @param K Positive number of communities.
#' @param spectral_engine Spectral-decomposition backend.
#' @param spectral_options Named options passed to the spectral decomposition.
#' @return A fitted block-model object.
#' @seealso [estimate_sbm_P_hat()], [estimate_dcbm_P_hat()], [spectral_cluster()]
#' @name estimate_blockmodel
NULL

#' Reconstruct block-model probability matrices
#'
#' Reconstructs fitted SBM or DCBM edge probabilities.
#' @inheritParams embedding_estimating
#' @param K Positive number of communities.
#' @param spectral_engine Spectral-decomposition backend.
#' @param spectral_options Named options passed to the spectral decomposition.
#' @return A fitted edge-probability matrix.
#' @seealso [estimate_sbm()], [estimate_dcbm()], [clip_probabilities()]
#' @name estimate_blockmodel_P_hat
NULL

#' Match community labels
#'
#' Aligns labels using greedy assignment or exact brute-force matching.
#' @inheritParams embedding_estimating
#' @param K Positive number of label classes.
#' @return A relabeled vector, optionally with mapping diagnostics.
#' @seealso [spectral_cluster()], [sonnet()]
#' @name label_match
NULL

#' Compute squared-error losses
#'
#' Computes summed or mean squared error for equal-length numeric vectors.
#' @inheritParams math_utilities
#' @return A nonnegative numeric scalar.
#' @seealso [sae()], [mae()]
#' @name squared_error_losses
NULL

#' Compute absolute-error losses
#'
#' Computes summed or mean absolute error for equal-length numeric vectors.
#' @inheritParams math_utilities
#' @return A nonnegative numeric scalar.
#' @seealso [sse()], [mse()]
#' @name absolute_error_losses
NULL

#' Compute binary deviance losses
#'
#' Computes summed or mean binary deviance after stable probability clipping.
#' @inheritParams math_utilities
#' @param epsilon Positive probability-clipping constant.
#' @return A nonnegative numeric scalar.
#' @seealso [auc()], [clip_probabilities()]
#' @name binary_deviance_losses
NULL

#' Compute AUC scores and losses
#'
#' Computes area under the ROC curve or its loss transformation.
#' @inheritParams math_utilities
#' @return A numeric scalar.
#' @seealso [bin_dev()], [bin_dev_mean()]
#' @name auc_losses
NULL

#' Compute hard and soft extrema
#'
#' Computes hard or temperature-smoothed extrema and their selected indices.
#' @inheritParams math_utilities
#' @param x Numeric vector whose extrema are selected.
#' @return A numeric weight vector or selected index, according to the function.
#' @seealso [softplus()], [sigmoid()]
#' @name extrema_selectors
NULL

#' Estimate RAM requirements
#'
#' Estimates memory for dense and partial eigenvalue decompositions, SVDs,
#' spectral decompositions, and matrix products.
#' @inheritParams resource_utilities
#' @return A named RAM-estimate object expressed in bytes.
#' @seealso [available_ram()], [report_ram_preflight()], [measure_peak_ram()]
#' @name ram_estimators
NULL

#' Inspect and report RAM availability
#'
#' Queries available memory, formats byte counts, and reports preflight formulas.
#' @inheritParams resource_utilities
#' @param ncores Positive worker count used by the availability report.
#' @return A byte count, formatted string, or diagnostic report.
#' @seealso [estimate_spectral_decomp_ram()], [measure_peak_ram()]
#' @name ram_reporting
NULL

#' Select SONNET partition parameters
#'
#' Recommends SONNET partition parameters for a requested network size.
#' @inheritParams model_selection
#' @param n Positive network size.
#' @param K Positive number of communities.
#' @return A named parameter recommendation.
#' @seealso [sonnet()], [netcrop_param_select()]
#' @name sonnet_param_select
NULL

#' Select NETCROP partition parameters
#'
#' Recommends NETCROP partition parameters for a requested network size.
#' @inheritParams model_selection
#' @param n Positive network size.
#' @return A named parameter recommendation.
#' @seealso [netcrop_blockmodel()], [netcrop_rdpg()], [netcrop_lsm()]
#' @name netcrop_param_select
NULL

#' Fit SONNET overlap variants
#'
#' Fits the shared-overlap or independent-overlap SONNET variant.
#' @inheritParams model_selection
#' @return A fitted `sonnet` object.
#' @seealso [sonnet()], [sonnet_param_select()]
#' @name sonnet_overlap_variants
NULL

#' Fit a SONNET network clustering
#'
#' Fits scalable overlapping-partition clustering to a dense or sparse network.
#' @inheritParams model_selection
#' @param K Positive number of communities.
#' @return A fitted `sonnet` object with print and summary methods.
#' @seealso [sonnet_shared_overlap()], [sonnet_independent_overlap()],
#'   [spectral_cluster()]
#' @name sonnet
NULL

#' Select a block model with NETCROP
#'
#' Selects community count and block-model family by overlapping-subnetwork
#' cross-validation.
#' @inheritParams model_selection
#' @param K_candidates Candidate community counts, conventionally `1:5`.
#' @param num_subnetworks Positive number of subnetworks.
#' @param overlap_size Nonnegative overlap size.
#' @param sbm_est_options,dcbm_est_options Named estimator options.
#' @return A `netcrop_blockmodel` result with print, summary, and plot methods.
#' @seealso [ecv_stability_blockmodel()], [ncv_stability_blockmodel()]
#' @name netcrop_blockmodel
NULL

#' Select RDPG dimension with NETCROP
#'
#' Selects latent dimension by overlapping-subnetwork cross-validation.
#' @inheritParams model_selection
#' @param d_candidates Candidate latent dimensions, conventionally `1:5`.
#' @param num_subnetworks Positive number of subnetworks.
#' @param overlap_size Nonnegative overlap size.
#' @param eig_options Named eigendecomposition options.
#' @return A `netcrop_rdpg` result with print, summary, and plot methods.
#' @seealso [ecv_stability_rdpg()], [ase()]
#' @name netcrop_rdpg
NULL

#' Select latent-space dimension with NETCROP
#'
#' Selects latent-space-model dimension by overlapping-subnetwork validation.
#' @inheritParams model_selection
#' @param d_candidates Candidate latent dimensions, conventionally `1:5`.
#' @param num_subnetworks Positive number of subnetworks.
#' @param overlap_size Nonnegative overlap size.
#' @param lsm_options Named latent-space fitting options.
#' @return A `netcrop_lsm` result with print, summary, and plot methods.
#' @seealso [generate_lsm()], [lsm_pgd()]
#' @name netcrop_lsm
NULL

#' Compare pairs of clusterings
#'
#' Computes NMI- or Hamming-based disagreement between two pairs of labelings.
#' @inheritParams model_selection
#' @return A numeric loss.
#' @seealso [nmi()], [label_match_greedy()]
#' @name pairwise_clustering_losses
NULL

#' Tune a network regularizer with NETCROP
#'
#' Selects a regularization parameter by overlapping-subnetwork validation.
#' @inheritParams model_selection
#' @param K Positive number of communities.
#' @param tau_candidates Candidate regularization values.
#' @param num_subnetworks Positive number of subnetworks.
#' @param overlap_size Nonnegative overlap size.
#' @return A `netcrop_regularizer` result with print, summary, and plot methods.
#' @seealso [dkest_tune_regularizer()], [mult_reg_spectral_cluster()]
#' @name netcrop_tune_regularizer
NULL

#' Select block-model size with edge cross-validation
#'
#' Runs repeated ECV over candidate community counts from 1 through `max_K`.
#' @inheritParams ecv_stability
#' @param max_K Maximum candidate community count; candidates are `1:max_K`.
#' @return A `netcrop_blockmodel` result with `algorithm = "ECV"`.
#' @seealso [ecv_stability_rdpg()], [ncv_stability_blockmodel()], [netcrop_blockmodel()]
#' @name ecv_stability_blockmodel
NULL

#' Select RDPG dimension with edge cross-validation
#'
#' Runs repeated ECV over candidate dimensions from 1 through `max_d`.
#' @inheritParams ecv_stability
#' @param max_d Maximum candidate dimension; candidates are `1:max_d`.
#' @return A `netcrop_rdpg` result with `algorithm = "ECV"`.
#' @seealso [ecv_stability_blockmodel()], [netcrop_rdpg()]
#' @name ecv_stability_rdpg
NULL

#' Select block-model size with node cross-validation
#'
#' Runs repeated NCV over candidate community counts from 1 through `max_K`.
#' @inheritParams ncv_stability
#' @return A `netcrop_blockmodel` result with `algorithm = "NCV"`.
#' @seealso [ecv_stability_blockmodel()], [netcrop_blockmodel()]
#' @name ncv_stability_blockmodel
NULL

#' Tune a degree-regularized spectral model
#'
#' Selects the spectral regularizer using the DKEST workflow.
#' @inheritParams model_selection
#' @param K Positive number of communities.
#' @param tau_candidates Candidate regularization values.
#' @return A regularizer-selection result.
#' @seealso [netcrop_tune_regularizer()], [mult_reg_spectral_cluster()]
#' @name dkest_tune_regularizer
NULL

#' Fit spectral clusterings over multiple regularizers
#'
#' Fits spectral clustering across a requested regularization grid.
#' @inheritParams model_selection
#' @param K Positive number of communities.
#' @param tau_candidates Candidate regularization values.
#' @param spectral_options,cluster_options Named spectral and clustering options.
#' @return A multi-regularization clustering result.
#' @seealso [dkest_tune_regularizer()], [netcrop_tune_regularizer()]
#' @name mult_reg_spectral_cluster
NULL

#' Fit SONNET over multiple regularizers
#'
#' Fits SONNET across a requested regularization grid.
#' @inheritParams model_selection
#' @param K Positive number of communities.
#' @param tau_candidates Candidate regularization values.
#' @return A multi-regularization SONNET result.
#' @seealso [sonnet()], [mult_reg_spectral_cluster()]
#' @name mult_reg_sonnet
NULL

#' Compare clustering methods with known labels
#'
#' Evaluates candidate clustering engines against supplied ground-truth labels.
#' @inheritParams model_selection
#' @param K Positive number of communities.
#' @param tau_candidates Candidate regularization values.
#' @return An oracle comparison plot or its underlying evaluation result.
#' @seealso [spectral_cluster()], [sonnet()]
#' @name oracle_plotter
NULL
