rd_text <- function(node) {
  paste(unlist(node, recursive = TRUE, use.names = FALSE), collapse = "")
}

rd_sections <- function(rd, tag) {
  Filter(function(x) identical(attr(x, "Rd_tag"), tag), rd)
}

rd_aliases <- function(rd) {
  trimws(vapply(rd_sections(rd, "\\alias"), rd_text, character(1)))
}

test_that("every public symbol is indexed with useful help metadata", {
  exports <- getNamespaceExports("netOP")
  db <- tools::Rd_db("netOP")
  aliases <- unlist(lapply(db, rd_aliases), use.names = FALSE)

  index_file <- system.file("help", "AnIndex", package = "netOP")
  index <- utils::read.table(
    index_file, sep = "\t", quote = "", comment.char = "",
    stringsAsFactors = FALSE, fill = TRUE
  )

  expect_setequal(intersect(aliases, exports), exports)
  expect_setequal(intersect(index[[1]], exports), exports)

  for (symbol in exports) {
    topic <- db[[which(vapply(db, function(rd) symbol %in% rd_aliases(rd), logical(1)))]]
    title <- trimws(rd_text(rd_sections(topic, "\\title")))
    description <- trimws(rd_text(rd_sections(topic, "\\description")))

    expect_true(nzchar(title), info = paste(symbol, "has no help title"))
    expect_true(nchar(title) >= 8L,
                info = paste(symbol, "has an uninformative help title"))
    expect_true(nzchar(description),
                info = paste(symbol, "has no help description"))
    expect_true(nchar(description) >= 20L,
                info = paste(symbol, "has an uninformative help description"))
  }
})

test_that("library(help) lists every public function", {
  exports <- getNamespaceExports("netOP")
  listing <- library(help = "netOP", character.only = TRUE)$info[[2]]
  listed_symbols <- sub("\\s.*$", "", trimws(listing[nzchar(trimws(listing))]))

  expect_setequal(intersect(listed_symbols, exports), exports)
})

test_that("only approved function families share public help topics", {
  approved_groups <- list(
    c("hardmin", "hardmax", "softmin", "softmax", "which_hardmin",
      "which_hardmax", "which_softmin", "which_softmax"),
    c("label_match_brute_force", "label_match_greedy"),
    c("generate_sbm", "generate_dcbm"),
    c("estimate_sbm", "estimate_dcbm"),
    c("estimate_sbm_P_hat", "estimate_dcbm_P_hat"),
    c("mse", "sse"),
    c("mae", "sae"),
    c("bin_dev", "bin_dev_mean"),
    c("auc", "auc_as_loss"),
    c("connected_components", "largest_connected_component"),
    c("adjacency_neighbors", "adjacency_weighted_neighbors"),
    c("read_network", "write_network"),
    c("scale_to_average_degree", "scale_to_average_degree_naive"),
    c("sonnet_shared_overlap", "sonnet_independent_overlap"),
    c("pair_hamming_loss", "pair_nmi_loss"),
    c("estimate_dense_eigen_ram", "estimate_dense_svd_ram",
      "estimate_partial_svd_ram", "estimate_rspectra_ram",
      "estimate_spectral_decomp_ram", "estimate_matrix_product_ram"),
    c("available_ram", "format_bytes", "warn_if_insufficient_ram",
      "report_ram_preflight", "report_ram_formula"),
    c("mean", "sum", "diag", "rowMeans", "rowSums", "colMeans", "colSums")
  )
  exports <- getNamespaceExports("netOP")
  db <- tools::Rd_db("netOP")
  shared <- Filter(function(x) length(x) > 1L, lapply(db, function(rd) {
    intersect(rd_aliases(rd), exports)
  }))

  is_approved <- function(group) {
    any(vapply(approved_groups, setequal, logical(1), y = group))
  }
  unexpected <- Filter(function(group) !is_approved(group), shared)

  expect_equal(
    unname(unexpected),
    list(),
    info = paste(
      "Unexpected grouped public help topics:",
      paste(vapply(unexpected, paste, character(1), collapse = ", "),
            collapse = "; ")
    )
  )

  for (symbol in c("generate_lsm", "generate_lsm_alpha",
                   "generate_lsm_positions")) {
    containing <- Filter(function(rd) symbol %in% rd_aliases(rd), db)
    expect_length(containing, 1L)
    expect_setequal(intersect(rd_aliases(containing[[1]]), exports), symbol)
  }
})

test_that("local help links resolve to installed aliases", {
  db <- tools::Rd_db("netOP")
  installed_aliases <- unique(unlist(lapply(db, rd_aliases), use.names = FALSE))

  collect_links <- function(x) {
    found <- character()
    if (identical(attr(x, "Rd_tag"), "\\link")) {
      option <- attr(x, "Rd_option")
      if (is.null(option) || !grepl(":", option, fixed = TRUE)) {
        target <- if (is.null(option) || !nzchar(option)) rd_text(x) else option
        target <- sub("^=", "", target)
        found <- trimws(target)
      }
    }
    if (is.list(x)) {
      found <- c(found, unlist(lapply(x, collect_links), use.names = FALSE))
    }
    found[nzchar(found)]
  }
  local_links <- unique(unlist(lapply(db, collect_links), use.names = FALSE))

  expect_setequal(intersect(local_links, installed_aliases), local_links)
})

test_that("network examples use the documented standard problem sizes", {
  db <- tools::Rd_db("netOP")
  examples <- paste(unlist(lapply(db, function(rd) {
    vapply(rd_sections(rd, "\\examples"), rd_text, character(1))
  }), use.names = FALSE), collapse = "\n")

  expect_false(grepl("generate_(sbm|dcbm|lsm)\\([^)]*\\bn\\s*=\\s*+(?!200\\b)",
                     examples, perl = TRUE))
  expect_false(grepl("generate_(sbm|dcbm|lsm)\\([^)]*\\bK\\s*=\\s*+(?!3\\b)",
                     examples, perl = TRUE))
  expect_false(grepl("generate_rdpg\\([^)]*\\bn\\s*=\\s*+(?!200\\b)",
                     examples, perl = TRUE))
  expect_false(grepl("generate_rdpg\\([^)]*\\bd\\s*=\\s*+(?!3\\b)",
                     examples, perl = TRUE))
  expect_false(grepl("max_K\\s*=\\s*+(?!5(?:L)?\\b)", examples, perl = TRUE))
  expect_false(grepl("max_d\\s*=\\s*+(?!5(?:L)?\\b)", examples, perl = TRUE))
})
