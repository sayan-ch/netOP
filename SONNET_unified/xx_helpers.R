# Shared NETCROP helpers
#
# Source x1_sonnet.R before this file. Both block-model and RDPG NETCROP use
# the same independent-overlap partitioning rule.

# Generate independent overlap partitions for NETCROP repetitions.
#
# Remainder nodes augment the overlap in every repetition. Consequently every
# node is retained, every non-overlap piece has the same size, and the returned
# effective overlap may be larger than overlap_size.
netcrop_splitter <- function(
    n,
    num_subnetworks,
    overlap_size,
    nrep = 1L,
    seed = NULL) {
  if (!exists("op_splitter", mode = "function", inherits = TRUE)) {
    stop("Source x1_sonnet.R before calling netcrop_splitter().",
         call. = FALSE)
  }
  split <- op_splitter(
    n = n,
    s = num_subnetworks,
    o = overlap_size,
    n_repetitions = nrep,
    m = floor((n - overlap_size) / num_subnetworks),
    remainder_handling = "augment_overlap",
    seed = seed
  )
  overlap_nodes <- lapply(split$splits, function(subnetworks) {
    Reduce(intersect, subnetworks)
  })
  nonoverlap_pieces <- lapply(seq_along(split$splits), function(repetition) {
    lapply(split$splits[[repetition]], function(nodes) {
      setdiff(nodes, overlap_nodes[[repetition]])
    })
  })
  names(nonoverlap_pieces) <- names(split$splits)
  list(
    subnetworks = split$splits,
    overlap_nodes = overlap_nodes,
    nonoverlap_pieces = nonoverlap_pieces,
    overlap_size = split$o,
    requested_overlap_size = as.integer(overlap_size),
    remainder_count = split$remainder_count,
    piece_size = as.integer(
      floor((n - overlap_size) / num_subnetworks)
    ),
    nrep = as.integer(nrep),
    remainder_handling = "augment_overlap"
  )
}
