# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjoberg; Copyright (c) 2026 Bogdan Oancea

.metric_pair <- function(y_true, y_pred) {
  if (!is.numeric(y_true) || !is.numeric(y_pred) || !identical(dim(y_true), dim(y_pred)) || length(y_true) != length(y_pred))
    stop("truth and prediction must have identical numeric shapes", call. = FALSE)
  if (!length(y_true) || any(!is.finite(y_true)) || any(!is.finite(y_pred)))
    stop("truth and prediction must be finite and non-empty", call. = FALSE)
  invisible(NULL)
}
.metric_reduce <- function(x, axis, fun) {
  if (is.null(dim(x)) || is.null(axis)) return(fun(as.numeric(x)))
  nd <- length(dim(x)); if (axis == -1L) axis <- nd
  if (length(axis) != 1L || !(axis %in% seq_len(nd))) stop("invalid axis", call. = FALSE)
  if (nd == 1L) return(fun(x))
  apply(x, setdiff(seq_len(nd), axis), fun)
}

#' Errors for spectra and their components
#' @param y_true,y_pred Numeric truth and predicted vectors or arrays of identical shape.
#' @param axis Axis to reduce; -1 means the final axis, NULL pools all values.
#' @param eps Threshold below which a constant target has undefined R squared.
#' @return Numeric error values, one per spectrum by default.
#' @export
#' @examples
#' truth <- matrix(c(1, 2, 3, 4), nrow = 2)
#' estimate <- truth + 0.1
#' rmse(truth, estimate)
rmse <- function(y_true, y_pred, axis = -1L) {
  .metric_pair(y_true, y_pred)
  sqrt(.metric_reduce((y_pred-y_true)^2, axis, mean))
}
#' @rdname rmse
#' @export
mae <- function(y_true, y_pred, axis = -1L) {
  .metric_pair(y_true, y_pred)
  .metric_reduce(abs(y_pred-y_true), axis, mean)
}
#' @rdname rmse
#' @export
max_error <- function(y_true, y_pred, axis = -1L) {
  .metric_pair(y_true, y_pred)
  .metric_reduce(abs(y_pred-y_true), axis, max)
}
#' @rdname rmse
#' @export
r2 <- function(y_true, y_pred, axis = -1L, eps = 1e-12) {
  .metric_pair(y_true, y_pred)
  ss_res <- .metric_reduce((y_true-y_pred)^2, axis, sum)
  ss_tot <- .metric_reduce(y_true, axis, function(z) sum((z-mean(z))^2))
  out <- 1-ss_res/ss_tot; out[ss_tot <= eps] <- NaN; out
}

#' Per-spectrum errors for each component
#' @param truth,prediction Lists or RADAR objects containing all four components.
#' @return A nested named list of errors or summary statistics.
#' @export
#' @examples
#' dataset <- generate(2, length = 128, seed = 7)
#' parts <- decompose(dataset$X, model = "standard")
#'
#' # One error per spectrum, for each component.
#' per_spectrum_errors(dataset, parts)
#'
#' # evaluate_components() is the same function under its other name.
#' evaluate_components(dataset, parts)
per_spectrum_errors <- function(truth, prediction) {
  out <- lapply(COMPONENTS, function(c) {
    t <- .radar_matrix(truth[[c]], c); p <- .radar_matrix(prediction[[c]], c)
    list(rmse=rmse(t,p), mae=mae(t,p), max_error=max_error(t,p))
  }); stats::setNames(out, COMPONENTS)
}
#' @rdname per_spectrum_errors
#' @export
error_statistics <- function(truth, prediction) {
  lapply(per_spectrum_errors(truth, prediction), function(errs) {
    out <- list()
    for (name in c("rmse","max_error")) {
      v <- errs[[name]]
      out[[paste0(name,"_mean")]] <- mean(v)
      out[[paste0(name,"_median")]] <- stats::median(v)
      out[[paste0(name,"_variance")]] <- mean((v-mean(v))^2)
    }
    out
  })
}
#' @rdname per_spectrum_errors
#' @export
evaluate_components <- function(truth, prediction) {
  errs <- per_spectrum_errors(truth, prediction)
  out <- lapply(COMPONENTS, function(c) {
    t <- as.numeric(truth[[c]]); p <- as.numeric(prediction[[c]])
    c(lapply(errs[[c]], mean), list(r2=r2(t,p)))
  }); stats::setNames(out, COMPONENTS)
}
#' Spectral signal-to-noise ratio
#' @param signal Intensities in a band region.
#' @param noise Intensities in a detrended band-free region.
#' @param predicted_noise Estimated noise to subtract.
#' @param eps Positive stabilizing constant.
#' @return Mean signal divided by population noise standard deviation.
#' @export
#' @examples
#' set.seed(1)
#'
#' signal <- matrix(sin(seq(0, 6, length.out = 64)), nrow = 1)
#' noise <- matrix(rnorm(64, sd = 0.05), nrow = 1)
#' snr(signal, noise)
#'
#' # Supplying a noise estimate scores a denoising result.
#' snr(signal, noise, predicted_noise = noise * 0.9)
snr <- function(signal, noise, predicted_noise = 0, eps = 1e-4) {
  residual <- as.numeric(noise-predicted_noise)
  mean(signal)/(sqrt(mean((residual-mean(residual))^2))+eps)
}
