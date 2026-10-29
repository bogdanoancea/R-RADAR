# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjoberg; Copyright (c) 2026 Bogdan Oancea

.radar_matrix <- function(x, name = "spectra", min_length = 2L) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
  if (!is.numeric(x) || length(dim(x)) != 2L ||
      nrow(x) < 1L || ncol(x) < min_length)
    stop(name, " must be a non-empty numeric matrix with at least ",
         min_length, " spectral points", call. = FALSE)
  if (any(!is.finite(x))) stop(name, " contains missing or infinite values", call. = FALSE)
  storage.mode(x) <- "double"
  x
}

.radar_model_cache <- new.env(parent = emptyenv())

.cached_model <- function(model) {
  if (inherits(model, "radar_model")) return(model)
  if (is.character(model) && length(model) == 1L &&
      model %in% c("standard", "ensemble") && !file.exists(model)) {
    if (!exists(model, .radar_model_cache, inherits = FALSE))
      assign(model, load_model(model), .radar_model_cache)
    return(get(model, .radar_model_cache, inherits = FALSE))
  }
  load_model(model)
}

#' Decompose Raman spectra with a native R neural network
#'
#' Predict baseline, cosmic rays, noise and Raman peaks using the original
#' standard or ensemble network. Every numerical operation is implemented in R.
#' Rows of a matrix are spectra; the last dimension of a map array is the
#' spectral axis. Lists may contain spectra of different lengths.
#'
#' @param spectra Numeric vector, matrix, map array, list of vectors,
#'   spectral_data object, or file name.
#' @param model Bundled model name, saved model file, or radar_model object.
#' @param output Optional file name to which components are written.
#' @param wavenumbers Optional spectral axis of matching length.
#' @param normalize Min-max normalize each input spectrum before inference.
#' @param units Return components in original or normalized intensity units.
#' @param batch_size Number of spectra processed in a prediction batch.
#' @param verbose Print progress when positive.
#' @param ... Options passed to read() when spectra is a file name.
#' @return A radar_components object, or a list of objects for unequal lengths.
#' @export
decompose <- function(spectra, model = "ensemble", output = NULL,
                      wavenumbers = NULL, normalize = TRUE,
                      units = c("original", "normalized"), batch_size = 32L,
                      verbose = 0L, ...) {
  units <- match.arg(units)
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize))
    stop("normalize must be TRUE or FALSE", call. = FALSE)
  net <- .cached_model(model)
  info <- list()
  read_options <- list(...)
  if (is.character(spectra) && length(spectra) == 1L) {
    spectra <- do.call(read, c(list(path = spectra), read_options))
  } else if (length(read_options)) {
    stop("File-reading options require a file name", call. = FALSE)
  }
  if (inherits(spectra, "spectral_data") || inherits(spectra, "radar_spectral_data")) {
    if (is.null(wavenumbers)) wavenumbers <- spectra$wavenumbers
    for (key in c("names", "positions", "position_names", "map_shape", "x_units"))
      info[key] <- list(spectra[[key]])
    info$metadata <- spectra$metadata
    spectra <- if (!is.null(spectra$spectra)) spectra$spectra else spectra$intensities
  }
  if (inherits(spectra, "synthetic_dataset")) spectra <- spectra$X
  if (is.list(spectra) && !is.data.frame(spectra)) {
    if (!length(spectra)) stop("spectra must not be empty", call. = FALSE)
    if (!all(vapply(spectra, is.numeric, logical(1L))))
      stop("Each list element must be a numeric spectrum", call. = FALSE)
    lens <- lengths(spectra)
    if (length(unique(lens)) > 1L) {
      if (!is.null(output))
        stop("Save spectra of different lengths separately", call. = FALSE)
      if (is.list(wavenumbers) && length(wavenumbers) != length(spectra))
        stop("Supply one spectral axis for each spectrum", call. = FALSE)
      return(lapply(seq_along(spectra), function(i) {
        axis <- if (is.list(wavenumbers)) wavenumbers[[i]] else wavenumbers
        decompose(spectra[[i]], model = net, wavenumbers = axis,
                  normalize = normalize, units = units, batch_size = batch_size,
                  verbose = verbose)
      }))
    }
    spectra <- do.call(rbind, lapply(spectra, as.numeric))
    if (is.list(wavenumbers)) wavenumbers <- wavenumbers[[1L]]
  }
  dims <- dim(spectra)
  single <- is.null(dims)
  if (length(dims) > 2L) {
    if (is.null(info$map_shape)) info$map_shape <- dims[-length(dims)]
    spectra <- matrix(spectra, nrow = prod(dims[-length(dims)]),
                      ncol = dims[length(dims)])
  }
  arr <- .radar_matrix(spectra)
  if (!is.null(wavenumbers) &&
      (!is.numeric(wavenumbers) || length(wavenumbers) != ncol(arr) ||
       any(!is.finite(wavenumbers))))
    stop("wavenumbers must be a finite numeric vector matching the spectrum length",
         call. = FALSE)
  if (normalize) {
    norm <- minmax_normalize(arr)
    x <- norm$normalized
    scaling <- norm$scaling
  } else {
    x <- arr
    scaling <- list(offset = rep(0, nrow(arr)), scale = rep(1, nrow(arr)))
  }
  parts <- stats::predict(net, x, batch_size = batch_size, verbose = verbose)
  if (units == "original") {
    parts <- denormalize_components(parts, scaling)
    inp <- arr
  } else inp <- x
  if (single) {
    inp <- as.numeric(inp)
    parts <- lapply(parts, as.numeric)
  }
  meta <- list(offset = scaling$offset, scale = scaling$scale,
               normalized_input = normalize)
  if (!is.null(info$metadata)) meta$source_metadata <- info$metadata
  if (!is.null(info$metadata$map_order)) meta$map_order <- info$metadata$map_order
  info$metadata <- NULL
  result <- do.call(radar_components, c(list(spectra = inp), parts,
    list(wavenumbers = wavenumbers, units = if (normalize) units else "original",
         model = net$architecture, meta = meta), info))
  if (!is.null(output)) write(result, output)
  result
}

#' Evaluate a model on spectra with known components
#' @param model Model name, file or radar_model object.
#' @param dataset Synthetic dataset or a file containing one.
#' @param statistics Return summaries of per-spectrum errors.
#' @param normalize Normalize inputs before predicting. Synthetic data are
#'   already scaled, so the default is FALSE.
#' @param batch_size Prediction batch size.
#' @return Named lists of error metrics for each component.
#' @export
evaluate <- function(model, dataset, statistics = FALSE, normalize = FALSE,
                     batch_size = 32L) {
  if (is.character(dataset) && length(dataset) == 1L) dataset <- load_dataset(dataset)
  if (!inherits(dataset, "synthetic_dataset"))
    stop("dataset must be a synthetic_dataset or a saved dataset file", call. = FALSE)
  parts <- decompose(dataset$X, model = model, normalize = normalize,
                     units = "original", batch_size = batch_size)
  if (isTRUE(statistics)) error_statistics(dataset, parts)
  else evaluate_components(dataset, parts)
}

#' Reference for the RADAR method
#' @return The article citation as a character string.
#' @export
radar_citation <- function() {
  paste("J. Sjoberg, N. Siminea, A. Paun, A. Lita, M. Larion, I. Petre,",
        "RADAR: Raman Spectral Analysis Using Deep Learning for Artifact Removal,",
        "Advanced Optical Materials (2025) 2500736.",
        "https://doi.org/10.1002/adom.202500736")
}
