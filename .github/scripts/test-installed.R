# Called in a fresh session after installing the candidate binary. The caller
# supplies library_path and version; never load the development namespace here.
library(netOP, lib.loc = library_path, character.only = FALSE)
stopifnot(
  normalizePath(find.package("netOP")) ==
    normalizePath(file.path(library_path, "netOP")),
  as.character(packageVersion("netOP")) == version
)

compiled_result <- getFromNamespace("outer_add_cpp", "netOP")(
  c(1, 2), c(3, 4)
)
stopifnot(isTRUE(all.equal(compiled_result, matrix(c(4, 5, 5, 6), 2))))

articles <- vignette(package = "netOP")$results
stopifnot(all(
  c("getting-started", "choosing-a-method") %in% articles[, "Item"]
))

# Plotting tests need a device but must not leave Rplots.pdf in the checkout.
grDevices::pdf(file = NULL)
tryCatch(
  testthat::test_dir(
    "tests/testthat", reporter = "summary",
    stop_on_failure = TRUE, stop_on_warning = TRUE
  ),
  finally = grDevices::dev.off()
)
print(sessionInfo())
