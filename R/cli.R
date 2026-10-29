# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Bogdan Oancea

#' RADAR command-line interface
#'
#' Run with the installed exec/radar.R script, or call this function directly.
#' Commands are decompose, convert, generate, train, evaluate, formats and info.
#' Use --help to display usage.
#' @param args Character vector of command-line arguments.
#' @return Invisibly, the command result.
#' @export
#' @examples
#' # The same interface as the exec/radar.R Rscript entry point.
#' radar_cli(c("formats"))
#'
#' out <- file.path(tempdir(), "generated.rds")
#' radar_cli(c("generate", "-n", "2", "--length", "64", "-o", out))
#' unlink(out)
radar_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  help <- paste(
    "RADAR: pure R Raman spectral decomposition",
    "  decompose INPUT [INPUT ...] -o OUTPUT [-m ensemble|standard|MODEL.rds]",
    "  convert INPUT OUTPUT",
    "  generate -n N -o DATASET [--length 1024] [--seed 0]",
    "  train -o MODEL.rds [-a ensemble|standard] [--init pretrained|MODEL.rds]",
    "  evaluate DATASET [-m ensemble|standard|MODEL.rds] [--statistics]",
    "  formats",
    "  info [MODEL.rds|standard|ensemble]",
    "Read options: --format, --layout, --key, --sheet, --table, --query,",
    "              --delimiter, --decimal",
    "Training options: --data, --epochs, --steps, --batch-size, --length,",
    "                  --learning-rate, --validation-size, --checkpoint, --seed",
    sep = "\n")
  if (!length(args) || any(args %in% c("-h", "--help"))) {
    cat(help, "\n")
    return(invisible(NULL))
  }
  cmd <- args[1L]
  aliases <- c(o = "output", m = "model", n = "n", a = "architecture")
  flags <- c("no-normalize", "statistics", "train")
  values <- c("output", "model", "n", "architecture", "units", "batch-size",
              "format", "layout", "key", "sheet", "table", "query", "delimiter",
              "decimal", "out-layout", "length", "min-width", "max-width",
              "seed", "init", "data", "epochs", "steps", "learning-rate",
              "validation-size", "checkpoint", "verbose")
  opt <- list()
  pos <- character()
  i <- 2L
  while (i <= length(args)) {
    token <- args[i]
    if (startsWith(token, "-")) {
      key <- sub("^-+", "", token)
      if (key %in% names(aliases)) key <- unname(aliases[key])
      if (key %in% flags) opt[[key]] <- TRUE
      else {
        if (!key %in% values) stop("Unknown option: ", token, call. = FALSE)
        if (i == length(args)) stop("Missing value for ", token, call. = FALSE)
        i <- i + 1L
        opt[[key]] <- args[i]
      }
    } else pos <- c(pos, token)
    i <- i + 1L
  }
  get_opt <- function(key, default = NULL) {
    if (is.null(opt[[key]])) default else opt[[key]]
  }
  num <- function(key, default) {
    val <- suppressWarnings(as.numeric(get_opt(key, default)))
    if (length(val) != 1L || !is.finite(val))
      stop("--", key, " must be a number", call. = FALSE)
    val
  }
  required_output <- function() {
    p <- get_opt("output")
    if (is.null(p)) stop("Specify -o OUTPUT", call. = FALSE)
    p
  }
  read_opt <- opt[intersect(names(opt),
    c("format", "layout", "key", "sheet", "table", "query", "delimiter", "decimal"))]
  if (!is.null(read_opt$sheet) && grepl("^[0-9]+$", read_opt$sheet))
    read_opt$sheet <- as.integer(read_opt$sheet)
  out_opt <- if (is.null(opt[["out-layout"]])) list() else list(layout = opt[["out-layout"]])
  model <- get_opt("model", "ensemble")
  result <- switch(cmd,
    decompose = {
      if (!length(pos)) stop("Specify an input file", call. = FALSE)
      output <- required_output()
      net <- load_model(model, architecture = get_opt("architecture"))
      outputs <- lapply(pos, function(path) {
        stem <- tools::file_path_sans_ext(basename(path))
        out <- if (grepl("{stem}", output, fixed = TRUE))
          gsub("{stem}", stem, output, fixed = TRUE)
        else if (length(pos) > 1L || dir.exists(output))
          file.path(output, paste0(stem, "_components.rds")) else output
        data <- do.call(read, c(list(path = path), read_opt))
        parts <- decompose(data, model = net, normalize = !isTRUE(opt[["no-normalize"]]),
          units = get_opt("units", "original"), batch_size = num("batch-size", 32L))
        do.call(write, c(list(parts, out), out_opt))
        cat(path, "->", out, "\n")
        out
      })
      unlist(outputs, use.names = FALSE)
    },
    convert = {
      if (length(pos) != 2L) stop("Usage: convert INPUT OUTPUT", call. = FALSE)
      obj <- tryCatch(load_components(pos[1L]), error = function(e) NULL)
      if (is.null(obj)) obj <- do.call(read, c(list(path = pos[1L]), read_opt))
      do.call(write, c(list(obj, pos[2L]), out_opt))
    },
    generate = {
      n <- num("n", NA_real_)
      kw <- list(n = n, length = num("length", 1024), seed = num("seed", 0))
      if (isTRUE(opt$train)) kw <- c(kw, list(min_width = 5, max_width = 256, max_num_peaks = 5,
                                             increase_noise = FALSE))
      if (!is.null(opt[["min-width"]])) kw$min_width <- num("min-width", 10)
      if (!is.null(opt[["max-width"]])) kw$max_width <- num("max-width", 300)
      ds <- do.call(generate_dataset, kw)
      write(ds, required_output())
    },
    train = {
      val_n <- num("validation-size", 1000)
      lr <- if (is.null(opt[["learning-rate"]])) "original" else num("learning-rate", 1e-5)
      train(data = get_opt("data"), model = get_opt("architecture", "ensemble"),
        init = get_opt("init"), save_to = required_output(),
        epochs = num("epochs", 10), steps_per_epoch = num("steps", 1000),
        batch_size = num("batch-size", 1), length = num("length", 750),
        learning_rate = lr, validation_data = if (val_n > 0) "synthetic" else NULL,
        validation_size = val_n, seed = num("seed", 0),
        checkpoint_path = get_opt("checkpoint"), verbose = num("verbose", 1))
    },
    evaluate = {
      if (length(pos) != 1L) stop("Specify a dataset file", call. = FALSE)
      ans <- evaluate(model, pos[1L], statistics = isTRUE(opt$statistics))
      print(ans)
      ans
    },
    formats = formats(show = TRUE),
    info = {
      if (length(pos)) model_info(pos[1L], show = TRUE)
      else {
        cat("RADAR", as.character(utils::packageVersion("RADAR")),
            "\n", R.version.string, "\nNative R inference and training\n")
        invisible(lapply(c("standard", "ensemble"), model_info, show = TRUE))
      }
    },
    stop("Unknown command: ", cmd, "\n", help, call. = FALSE))
  invisible(result)
}
