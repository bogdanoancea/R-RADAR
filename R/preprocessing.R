# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjoberg; Copyright (c) 2026 Bogdan Oancea

#' Per-spectrum min-max normalization
#' @param spectra Numeric matrix, one spectrum per row.
#' @return A list with normalized spectra and affine scaling parameters.
#' @export
minmax_normalize <- function(spectra) {
  if (!is.matrix(spectra) || !is.numeric(spectra) || !nrow(spectra) || !ncol(spectra) || any(!is.finite(spectra)))
    stop("spectra must be a non-empty finite numeric matrix", call. = FALSE)
  offset <- apply(spectra, 1L, min)
  scale <- apply(spectra, 1L, max) - offset
  scale[scale <= 0] <- 1
  list(normalized = sweep(sweep(spectra, 1L, offset, "-"), 1L, scale, "/"),
       scaling = structure(list(offset = offset, scale = scale), class = "radar_scaling"))
}

#' Restore component intensity units
#' @param components Named list of component matrices.
#' @param scaling List containing per-spectrum offset and scale vectors.
#' @return The components in original intensity units. The offset is added to baseline.
#' @export
denormalize_components <- function(components, scaling) {
  if (!is.list(components) || is.null(names(components))) stop("components must be a named list", call. = FALSE)
  out <- lapply(components, function(x) {
    if (is.null(dim(x))) {
      if (length(scaling$scale) != 1L) stop("scaling does not match the components", call. = FALSE)
      x * scaling$scale
    } else {
      if (nrow(x) != length(scaling$scale)) stop("scaling does not match the components", call. = FALSE)
      sweep(x, 1L, scaling$scale, "*")
    }
  })
  if (!is.null(out$baseline)) out$baseline <- if (is.null(dim(out$baseline)))
    out$baseline + scaling$offset else sweep(out$baseline, 1L, scaling$offset, "+")
  out
}
