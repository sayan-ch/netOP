# Interactive source-file selector
#
# When this file is sourced interactively, it locates its own directory,
# lists the other .R files alphabetically, and asks which files to source and
# in which order. Enter 0 to source every listed file alphabetically or -1 to
# source none. The current file is excluded to prevent recursive sourcing.

source_all <- local({
  source_file <- NULL
  source_frames <- sys.frames()
  if (length(source_frames) > 0L) {
    source_candidates <- vapply(
      source_frames,
      function(frame) {
        file <- frame$ofile
        if (is.null(file) || length(file) != 1L) "" else as.character(file)
      },
      character(1)
    )
    source_candidates <- source_candidates[nzchar(source_candidates)]
    if (length(source_candidates) > 0L) {
      source_file <- source_candidates[length(source_candidates)]
    }
  }

  if (is.null(source_file)) {
    file_argument <- grep(
      "^--file=",
      commandArgs(trailingOnly = FALSE),
      value = TRUE
    )
    if (length(file_argument) > 0L) {
      source_file <- sub("^--file=", "", file_argument[length(file_argument)])
    }
  }

  source_directory <- if (is.null(source_file)) {
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  } else {
    dirname(normalizePath(source_file, winslash = "/", mustWork = TRUE))
  }
  source_file_normalized <- if (is.null(source_file)) {
    file.path(source_directory, "source_all.R")
  } else {
    normalizePath(source_file, winslash = "/", mustWork = TRUE)
  }

  function(
      selection = NULL,
      envir = parent.frame(),
      echo = FALSE,
      verbose = TRUE) {
    if (!is.environment(envir)) {
      stop("envir must be an environment.", call. = FALSE)
    }
    logical_inputs <- list(echo = echo, verbose = verbose)
    invalid_logical <- vapply(
      logical_inputs,
      function(value) {
        length(value) != 1L || !is.logical(value) || is.na(value)
      },
      logical(1)
    )
    if (any(invalid_logical)) {
      stop(
        paste(names(logical_inputs)[invalid_logical], collapse = ", "),
        " must be TRUE or FALSE.",
        call. = FALSE
      )
    }

    r_files <- list.files(
      source_directory,
      pattern = "[.]R$",
      full.names = TRUE
    )
    if (length(r_files) > 0L) {
      r_files <- normalizePath(r_files, winslash = "/", mustWork = TRUE)
      r_files <- r_files[r_files != source_file_normalized]
      r_files <- r_files[order(basename(r_files), method = "radix")]
    }

    if (length(r_files) == 0L) {
      if (verbose) message("No other .R files were found in ", source_directory, ".")
      return(invisible(character()))
    }

    if (is.null(selection)) {
      cat("R files in ", source_directory, ":\n", sep = "")
      cat(sprintf("%d: %s\n", seq_along(r_files), basename(r_files)))
      cat("0: source all files in alphabetical order\n")
      cat("-1: source none\n")
      selection <- readline(
        paste0(
          "Enter serial numbers in the desired order ",
          "(separated by spaces or commas): "
        )
      )
    }

    if (is.character(selection)) {
      selection <- trimws(paste(selection, collapse = " "))
      if (!nzchar(selection)) {
        stop("No selection was entered.", call. = FALSE)
      }
      selection_tokens <- strsplit(
        selection,
        "[,[:space:]]+",
        perl = TRUE
      )[[1L]]
      if (any(!grepl("^-?[0-9]+$", selection_tokens))) {
        stop("Selections must be integer serial numbers.", call. = FALSE)
      }
      selection <- as.integer(selection_tokens)
    } else if (is.numeric(selection)) {
      if (length(selection) < 1L || anyNA(selection) ||
          any(!is.finite(selection)) || any(selection != floor(selection))) {
        stop("selection must contain integer serial numbers.", call. = FALSE)
      }
      selection <- as.integer(selection)
    } else {
      stop(
        "selection must be NULL, an integer vector, or a character entry.",
        call. = FALSE
      )
    }

    if (-1L %in% selection) {
      if (length(selection) != 1L) {
        stop("-1 must be entered by itself.", call. = FALSE)
      }
      if (verbose) message("No files were sourced.")
      return(invisible(character()))
    }
    if (0L %in% selection) {
      if (length(selection) != 1L) {
        stop("0 must be entered by itself.", call. = FALSE)
      }
      selection <- seq_along(r_files)
    }
    if (any(selection < 1L | selection > length(r_files))) {
      stop(
        "Every serial number must be between 1 and ", length(r_files),
        ", or use 0/-1 as a standalone choice.",
        call. = FALSE
      )
    }
    if (anyDuplicated(selection)) {
      stop("Each file may be selected at most once.", call. = FALSE)
    }

    selected_files <- r_files[selection]
    for (file in selected_files) {
      if (verbose) message("Sourcing ", basename(file), " ...")
      source(
        file,
        local = envir,
        echo = echo,
        chdir = TRUE,
        encoding = "UTF-8"
      )
    }
    invisible(selected_files)
  }
})

if (interactive()) {
  invisible(source_all())
}
