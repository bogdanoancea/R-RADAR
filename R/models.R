# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjoberg (original architectures and weights)
# Copyright (c) 2026 Bogdan Oancea (RADAR R package)

.radar_architectures <- c("standard", "ensemble")
.radar_components <- c("baseline", "cosmic_rays", "noise", "peaks")

# All tensors use batch-first dimensions; the final dimension is channels.
# Matrix multiplication performs cross-correlation using native R arrays.
.radar_node <- function(ctx, value, parents = NULL, backward = NULL) {
  z <- new.env(parent = emptyenv())
  z$value <- value
  if (ctx$gradient) {
    z$parents <- parents
    z$backward <- backward
    z$grad <- NULL
    ctx$nodes[[length(ctx$nodes) + 1L]] <- z
  }
  z
}
.radar_addgrad <- function(node, value) {
  if (is.null(node$grad)) node$grad <- value else node$grad <- node$grad + value
}
.radar_paramgrad <- function(ctx, layer, key, value) {
  if (is.null(ctx$grads[[layer]])) ctx$grads[[layer]] <- list()
  if (is.null(ctx$grads[[layer]][[key]])) ctx$grads[[layer]][[key]] <- value
  else ctx$grads[[layer]][[key]] <- ctx$grads[[layer]][[key]] + value
}
.radar_patches <- function(x, kernel, stride = 1L, edge = TRUE) {
  d <- dim(x); n <- d[1L]; nc <- tail(d, 1L)
  if (length(d) == 3L) {
    h <- d[2L]; kh <- kernel[1L]; oh <- h
    idx <- matrix(0L, n * h, kh)
    for (j in seq_len(kh)) {
      p <- pmax(1L, pmin(h, seq_len(h) + j - 1L - kh %/% 2L))
      idx[, j] <- as.vector(outer(seq_len(n), (p - 1L) * n, `+`))
    }
    outdim <- c(n, h)
  } else {
    h <- d[2L]; w <- d[3L]; kh <- kernel[1L]; kw <- kernel[2L]
    sh <- stride[1L]; sw <- stride[2L]; oh <- ceiling(h / sh); ow <- ceiling(w / sw)
    ph <- floor(max(0L, (oh - 1L) * sh + kh - h) / 2L)
    pw <- floor(max(0L, (ow - 1L) * sw + kw - w) / 2L)
    grid <- expand.grid(n = seq_len(n), h = seq_len(oh), w = seq_len(ow))
    idx <- matrix(0L, n * oh * ow, kh * kw)
    for (j in seq_len(kw)) for (i in seq_len(kh)) {
      hp <- (grid$h - 1L) * sh + i - ph
      wp <- (grid$w - 1L) * sw + j - pw
      valid <- hp >= 1L & hp <= h & wp >= 1L & wp <= w
      ii <- grid$n + (hp - 1L) * n + (wp - 1L) * n * h
      ii[!valid] <- 0L
      idx[, i + (j - 1L) * kh] <- ii
    }
    outdim <- c(n, oh, ow)
  }
  xm <- matrix(x, ncol = nc)
  p <- ncol(idx); patches <- matrix(0, nrow(idx), p * nc)
  for (j in seq_len(nc)) {
    v <- c(0, xm[, j])
    patches[, (j - 1L) * p + seq_len(p)] <- matrix(v[idx + 1L], nrow(idx))
  }
  list(x = patches, idx = idx, outdim = outdim, dim = d)
}
.radar_unpatch <- function(patches, cache) {
  d <- cache$dim; nc <- tail(d, 1L); p <- ncol(cache$idx)
  dx <- matrix(0, prod(head(d, -1L)), nc)
  ii <- as.vector(cache$idx); valid <- ii > 0L
  for (j in seq_len(nc)) {
    values <- as.vector(patches[, (j - 1L) * p + seq_len(p), drop = FALSE])
    sums <- rowsum(values[valid], ii[valid], reorder = FALSE)
    dx[as.integer(rownames(sums)), j] <- sums[, 1L]
  }
  array(dx, d)
}
.radar_conv <- function(ctx, node, name, stride = 1L) {
  pars <- ctx$weights[[name]]; w <- pars$kernel; b <- pars$bias
  kw <- head(dim(w), -2L)
  cache <- .radar_patches(node$value, kw, stride)
  wm <- matrix(w, ncol = length(b))
  y <- sweep(cache$x %*% wm, 2L, b, `+`)
  value <- array(y, c(cache$outdim, length(b)))
  l1 <- if (grepl("^(baseline|cosmic_ray|peak|noise)_conv$", name)) 1e-5 else ctx$l1
  l2 <- ctx$l2
  if (ctx$loss) {
    ctx$reg_loss <- ctx$reg_loss + l1 * sum(abs(w)) + l2 * (sum(w * w) + sum(b * b)) +
      l2 * sum(value * value) / dim(node$value)[1L]
  }
  .radar_node(ctx, value, list(node), function(g) {
    gm <- matrix(g + 2 * l2 * value / dim(node$value)[1L], ncol = length(b))
    dw <- crossprod(cache$x, gm)
    db <- colSums(gm)
    .radar_paramgrad(ctx, name, "kernel", array(dw, dim(w)) + l1 * sign(w) + 2 * l2 * w)
    .radar_paramgrad(ctx, name, "bias", db + 2 * l2 * b)
    .radar_addgrad(node, .radar_unpatch(gm %*% t(wm), cache))
  })
}
.radar_bn <- function(ctx, node, name) {
  p <- ctx$weights[[name]]; d <- dim(node$value); c <- tail(d, 1L)
  x <- matrix(node$value, ncol = c)
  if (ctx$training) {
    mu <- colMeans(x); centered <- sweep(x, 2L, mu, `-`)
    variance <- colMeans(centered * centered)
    ctx$weights[[name]]$moving_mean <- 0.99 * p$moving_mean + 0.01 * mu
    ctx$weights[[name]]$moving_variance <- 0.99 * p$moving_variance + 0.01 * variance
  } else {
    mu <- p$moving_mean; variance <- p$moving_variance
    centered <- sweep(x, 2L, mu, `-`)
  }
  inv <- 1 / sqrt(variance + 0.001)
  norm <- sweep(centered, 2L, inv, `*`)
  y <- sweep(sweep(norm, 2L, p$gamma, `*`), 2L, p$beta, `+`)
  .radar_node(ctx, array(y, d), list(node), function(g) {
    gm <- matrix(g, ncol = c)
    .radar_paramgrad(ctx, name, "gamma", colSums(gm * norm))
    .radar_paramgrad(ctx, name, "beta", colSums(gm))
    dn <- sweep(gm, 2L, p$gamma, `*`)
    if (ctx$training) {
      m <- nrow(x)
      dx <- sweep(dn, 2L, colMeans(dn), `-`) - sweep(norm, 2L, colMeans(dn * norm), `*`)
      dx <- sweep(dx, 2L, inv, `*`)
    } else dx <- sweep(dn, 2L, inv, `*`)
    .radar_addgrad(node, array(dx, d))
  })
}
.radar_act <- function(ctx, node, kind = "leaky") {
  x <- node$value
  if (kind == "leaky") {
    y <- pmax(x, 0) + 0.3 * pmin(x, 0); deriv <- ifelse(x >= 0, 1, 0.3)
  } else if (kind == "relu") {
    y <- pmax(x, 0); deriv <- as.numeric(x > 0)
  } else { y <- tanh(x); deriv <- 1 - y * y }
  dim(y) <- dim(x); dim(deriv) <- dim(x)
  .radar_node(ctx, y, list(node), function(g) .radar_addgrad(node, g * deriv))
}
.radar_bn_act <- function(ctx, node, prefix) .radar_act(ctx, .radar_bn(ctx, node, paste0(prefix, "_bn")))
.radar_reshape <- function(ctx, node, d) {
  .radar_node(ctx, array(node$value, d), list(node), function(g) .radar_addgrad(node, array(g, dim(node$value))))
}
.radar_add <- function(ctx, nodes, signs = rep(1, length(nodes))) {
  value <- nodes[[1L]]$value * signs[1L]
  for (j in seq_along(nodes)[-1L]) value <- value + signs[j] * nodes[[j]]$value
  .radar_node(ctx, value, nodes, function(g) {
    for (j in seq_along(nodes)) .radar_addgrad(nodes[[j]], signs[j] * g)
  })
}
.radar_pool <- function(ctx, node, size) {
  cache <- .radar_patches(node$value, size)
  p <- ncol(cache$idx); c <- tail(dim(node$value), 1L)
  y <- matrix(0, nrow(cache$x), c)
  for (j in seq_len(c)) y[, j] <- rowMeans(cache$x[, (j - 1L) * p + seq_len(p), drop = FALSE])
  .radar_node(ctx, array(y, dim(node$value)), list(node), function(g) {
    gm <- matrix(g, ncol = c)
    gp <- matrix(0, nrow(gm), p * c)
    for (j in seq_len(c)) gp[, (j - 1L) * p + seq_len(p)] <- gm[, rep(j, p), drop = FALSE] / p
    .radar_addgrad(node, .radar_unpatch(gp, cache))
  })
}
.radar_maxpool <- function(ctx, node) {
  d <- dim(node$value); n <- d[1L]; len <- d[2L]; c <- d[3L]
  y <- matrix(0, n, c); where <- vector("list", n * c)
  for (j in seq_len(c)) for (i in seq_len(n)) {
    v <- node$value[i, , j]; maximum <- max(v)
    k <- which(v == maximum)
    y[i, j] <- maximum
    where[[i + (j - 1L) * n]] <- i + (k - 1L) * n + (j - 1L) * n * len
  }
  .radar_node(ctx, y, list(node), function(g) {
    dx <- array(0, d)
    for (j in seq_along(where)) dx[where[[j]]] <- g[j] / length(where[[j]])
    .radar_addgrad(node, dx)
  })
}
.radar_outer <- function(ctx, signal, key) {
  d <- dim(signal$value); n <- d[1L]; len <- d[2L]; k <- ncol(key$value)
  x <- matrix(signal$value, n, len); value <- array(0, c(n, len, k, 1L))
  for (j in seq_len(k)) value[, , j, 1L] <- x * key$value[, j]
  .radar_node(ctx, value, list(signal, key), function(g) {
    dx <- matrix(0, n, len); dk <- matrix(0, n, k)
    for (j in seq_len(k)) {
      gj <- matrix(g[, , j, 1L], n, len)
      dx <- dx + gj * key$value[, j]
      dk[, j] <- rowSums(gj * x)
    }
    .radar_addgrad(signal, array(dx, d)); .radar_addgrad(key, dk)
  })
}
.radar_attention <- function(ctx, signal, features, prefix, scaler, dim_red_size) {
  key <- .radar_maxpool(ctx, .radar_conv(ctx, features, paste0(prefix, "key_conv")))
  t <- .radar_outer(ctx, signal, key)
  stride <- c(1L, as.integer(sqrt(dim_red_size)))
  t <- .radar_conv(ctx, t, paste0(prefix, "conv2d_1"), stride)
  t <- .radar_bn_act(ctx, t, paste0(prefix, "conv2d_1"))
  .radar_conv(ctx, t, paste0(prefix, "conv2d_2"), stride)
}
.radar_forward <- function(model, x, training = FALSE, gradient = FALSE, loss = gradient) {
  if (is.null(dim(x))) x <- matrix(x, 1L)
  if (length(dim(x)) == 3L) {
    if (dim(x)[3L] != 1L) stop("The input must have one channel.", call. = FALSE)
    x <- matrix(x, dim(x)[1L], dim(x)[2L])
  }
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x) || any(!is.finite(x)))
    stop("newdata must be a finite numeric spectrum or a nonempty spectra matrix.", call. = FALSE)
  ctx <- new.env(parent = emptyenv()); ctx$gradient <- gradient; ctx$training <- training; ctx$loss <- loss
  ctx$weights <- model$weights; ctx$nodes <- list(); ctx$grads <- list()
  ctx$l1 <- model$hyperparameters$l1_param; ctx$l2 <- model$hyperparameters$l2_param
  ctx$reg_loss <- 0
  n <- nrow(x); len <- ncol(x); h <- model$hyperparameters
  inp <- .radar_node(ctx, array(x, c(n, len, 1L)))
  t <- .radar_bn_act(ctx, .radar_conv(ctx, inp, "enc_conv1"), "enc_conv1")
  residual <- .radar_attention(ctx, inp, t, "", h$scaler, h$dim_red_size)
  t <- .radar_bn_act(ctx, residual, "conv2d_2")
  t <- .radar_reshape(ctx, t, c(n, len, h$scaler))
  t <- .radar_bn_act(ctx, .radar_conv(ctx, t, "dec_conv1"), "dec_conv1")
  t <- .radar_bn_act(ctx, .radar_conv(ctx, t, "dec_conv2"), "dec_conv2")
  bl <- .radar_conv(ctx, t, "baseline_conv")
  if (model$architecture == "standard") bl <- .radar_act(ctx, bl, "relu")
  bl <- .radar_act(ctx, .radar_pool(ctx, bl, 33L), "relu")
  cr <- .radar_act(ctx, .radar_conv(ctx, t, "cosmic_ray_conv"), "relu")
  if (model$architecture == "ensemble") {
    reduced <- .radar_add(ctx, list(inp, bl, cr), c(1, -1, -1))
    t <- .radar_bn_act(ctx, .radar_conv(ctx, reduced, "s2_enc_conv1"), "s2_enc_conv1")
    t <- .radar_attention(ctx, reduced, t, "s2_", h$scaler, h$dim_red_size)
    t <- .radar_bn_act(ctx, .radar_add(ctx, list(t, residual)), "s2_conv2d_2")
    t <- .radar_reshape(ctx, t, c(n, len, h$scaler))
    t <- .radar_bn_act(ctx, .radar_conv(ctx, t, "s2_dec_conv1"), "s2_dec_conv1")
    t <- .radar_bn_act(ctx, .radar_conv(ctx, t, "s2_dec_conv2"), "s2_dec_conv2")
  }
  pk <- .radar_conv(ctx, t, "peak_conv")
  if (model$architecture == "standard") pk <- .radar_act(ctx, pk, "relu")
  pk <- .radar_act(ctx, .radar_pool(ctx, pk, 4L), "relu")
  nz <- if (model$architecture == "standard") .radar_add(ctx, list(inp, bl, cr, pk), c(1, -1, -1, -1))
    else .radar_act(ctx, .radar_conv(ctx, t, "noise_conv"), "tanh")
  outputs <- list(baseline = bl, cosmic_rays = cr, noise = nz, peaks = pk)
  list(outputs = outputs, context = ctx)
}

.radar_initial_weights <- function(architecture, scaler, dim_red_size) {
  conv <- function(shape) {
    fan_in <- prod(head(shape, -1L)); fan_out <- tail(shape, 1L) * prod(head(shape, -2L))
    limit <- sqrt(6 / (fan_in + fan_out))
    list(kernel = array(stats::runif(prod(shape), -limit, limit), shape), bias = rep(0, tail(shape, 1L)))
  }
  bn <- function() list(gamma = rep(1, scaler), beta = rep(0, scaler), moving_mean = rep(0, scaler), moving_variance = rep(1, scaler))
  stage <- function() list(enc_conv1 = conv(c(128L, 1L, scaler)), enc_conv1_bn = bn(),
    key_conv = conv(c(64L, scaler, dim_red_size)),
    conv2d_1 = conv(c(32L, floor(dim_red_size / 10), 1L, scaler)), conv2d_1_bn = bn(),
    conv2d_2 = conv(c(16L, floor(dim_red_size / 10), scaler, scaler)), conv2d_2_bn = bn(),
    dec_conv1 = conv(c(128L, scaler, scaler)), dec_conv1_bn = bn(),
    dec_conv2 = conv(c(64L, scaler, scaler)), dec_conv2_bn = bn())
  w <- stage()
  w$baseline_conv <- conv(c(64L, scaler, 1L)); w$cosmic_ray_conv <- conv(c(3L, scaler, 1L))
  if (architecture == "ensemble") {
    s2 <- stage(); names(s2) <- paste0("s2_", names(s2)); w <- c(w, s2)
    w$noise_conv <- conv(c(1L, scaler, 1L))
  }
  w$peak_conv <- conv(c(8L, scaler, 1L)); w
}

#' Construct a native R RADAR convolutional network
#'
#' Networks have the same layers, padding, activations and output order as the
#' original models. New networks are randomly initialized. Use [load_model()]
#' for the pretrained models. All inference and differentiation use R arrays.
#' @param architecture Either `"standard"` or `"ensemble"`.
#' @param scaler Number of hidden channels.
#' @param dim_red_size Key-vector dimension, a perfect square of at least 16.
#' @param l1_param,l2_param Kernel and activity regularization coefficients.
#' @param seed Random seed, or `NULL` to use the current random state.
#' @param name Model name.
#' @return A `radar_model` object.
#' @export
build_model <- function(architecture = "ensemble", scaler = NULL, dim_red_size = NULL,
                        l1_param = 1e-6, l2_param = 1e-6, seed = 0, name = paste0("radar_", architecture)) {
  architecture <- match.arg(architecture, .radar_architectures)
  if (is.null(scaler)) scaler <- if (architecture == "standard") 20L else 10L
  if (is.null(dim_red_size)) dim_red_size <- if (architecture == "standard") 256L else 144L
  if (length(scaler) != 1L || !is.finite(scaler) || scaler < 1 || scaler != as.integer(scaler))
    stop("scaler must be a positive integer.", call. = FALSE)
  if (length(dim_red_size) != 1L || !is.finite(dim_red_size) || dim_red_size < 16L ||
      sqrt(dim_red_size) != floor(sqrt(dim_red_size)))
    stop("dim_red_size must be a perfect square of at least 16.", call. = FALSE)
  if (any(!is.finite(c(l1_param, l2_param))) || any(c(l1_param, l2_param) < 0))
    stop("Regularization parameters must be finite and nonnegative.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  w <- .radar_initial_weights(architecture, as.integer(scaler), as.integer(dim_red_size))
  structure(list(architecture = architecture, name = name,
    hyperparameters = list(scaler = as.integer(scaler), dim_red_size = as.integer(dim_red_size),
      l1_param = l1_param, l2_param = l2_param), weights = w,
    n_parameters = sum(vapply(w, function(p) sum(lengths(p)), numeric(1))),
    origin = "random initialization in R", training = NULL, history = NULL,
    format_version = 1L, backend = "pure R", learning_rate = 1e-5), class = "radar_model")
}
#' @rdname build_model
#' @export
build_standard_model <- function(scaler = 20L, dim_red_size = 256L, l1_param = 1e-6,
                                 l2_param = 1e-6, seed = 0, name = "radar_standard")
  build_model("standard", scaler, dim_red_size, l1_param, l2_param, seed, name)
#' @rdname build_model
#' @export
build_ensemble_model <- function(scaler = 10L, dim_red_size = 144L, l1_param = 1e-6,
                                 l2_param = 1e-6, seed = 0, name = "radar_ensemble")
  build_model("ensemble", scaler, dim_red_size, l1_param, l2_param, seed, name)
#' @rdname build_model
#' @param ... Passed to `build_model()`.
#' @export
new_model <- function(architecture = "ensemble", ...) build_model(architecture, ...)

#' Configure an optimizer for a RADAR model
#' @param model A `radar_model`.
#' @param architecture Architecture name (must agree with `model`).
#' @param learning_rate Positive number, or a function of the zero-based step.
#' @return The configured model.
#' @export
compile_model <- function(model, architecture = model$architecture, learning_rate = 1e-5) {
  if (!inherits(model, "radar_model")) stop("Expected a radar_model.", call. = FALSE)
  if (!identical(architecture, model$architecture)) stop("Architecture does not match the model.", call. = FALSE)
  if (!is.function(learning_rate) && (length(learning_rate) != 1L || !is.finite(learning_rate) || learning_rate <= 0))
    stop("learning_rate must be positive or a schedule function.", call. = FALSE)
  model$learning_rate <- learning_rate; model
}

#' Predict normalized spectral components
#' @param object A `radar_model`.
#' @param newdata Numeric vector, spectra-by-points matrix, or batch-by-points-by-1 array.
#' @param batch_size Maximum spectra per prediction group. Each spectrum is
#'   evaluated separately to bound the memory needed by the attention layer.
#' @param verbose Print progress when positive.
#' @param ... Reserved for future use.
#' @return Named list of four spectra-by-points matrices in network units.
#' @export
#' @examples
#' model <- load_model("standard")
#' dataset <- generate(2, length = 128, seed = 7)
#'
#' # Inputs are expected on the normalized scale here; decompose() is the
#' # higher-level entry point that normalizes and restores units for you.
#' parts <- predict(model, dataset$X)
#' names(parts)
#' str(parts$peaks)
predict.radar_model <- function(object, newdata, batch_size = 32L, verbose = 0L, ...) {
  if (is.null(dim(newdata))) newdata <- matrix(newdata, 1L)
  if (length(dim(newdata)) == 3L && dim(newdata)[3L] == 1L)
    newdata <- matrix(newdata, dim(newdata)[1L], dim(newdata)[2L])
  if (!is.matrix(newdata) || !is.numeric(newdata) || !nrow(newdata) || !ncol(newdata) || any(!is.finite(newdata)))
    stop("newdata must be a nonempty finite numeric spectra matrix.", call. = FALSE)
  if (length(batch_size) != 1L || !is.finite(batch_size) || batch_size < 1L)
    stop("batch_size must be positive.", call. = FALSE)
  out <- setNames(lapply(.radar_components, function(k) matrix(0, nrow(newdata), ncol(newdata))), .radar_components)
  for (i in seq_len(nrow(newdata))) {
    ans <- .radar_forward(object, newdata[i, , drop = FALSE])$outputs
    for (k in .radar_components) out[[k]][i, ] <- as.numeric(ans[[k]]$value)
    if (verbose > 0L) message("Spectrum ", i, "/", nrow(newdata))
  }
  out
}
#' @export
print.radar_model <- function(x, ...) {
  cat("RADAR ", x$architecture, " model (pure R), ", format(x$n_parameters, big.mark = ","), " parameters\n", sep = "")
  cat("Origin: ", x$origin, "\n", sep = ""); invisible(x)
}

.radar_weights_path <- function(architecture) {
  p <- system.file("weights", paste0(architecture, ".rds"), package = "RADAR")
  if (!nzchar(p)) stop("Bundled R weights are missing; install the complete RADAR package.", call. = FALSE)
  p
}


#' Load pretrained or saved native R model weights
#' @param model `"ensemble"`, `"standard"`, an RDS model path, or an existing model.
#' @param architecture Optional architecture name used to validate the model.
#' @details Pretrained models and weights are bundled as native R objects.
#' @param compile Whether to configure the default optimizer when loading.
#' @return A native `radar_model`.
#' @export
#' @examples
#' # A bundled architecture by name; load_pretrained() is the same call
#' # restricted to the bundled names.
#' model <- load_model("standard")
#' model
#' model_info(model)
load_model <- function(model = "ensemble", architecture = NULL, compile = TRUE) {
  if (inherits(model, "radar_model")) {
    if (!is.null(architecture) && !identical(architecture, model$architecture))
      stop("The model architecture does not match the requested architecture.", call. = FALSE)
    return(model)
  }
  if (!is.character(model) || length(model) != 1L) stop("model must be a model name, path, or radar_model.", call. = FALSE)
  if (model %in% .radar_architectures && !file.exists(model)) {
    ans <- readRDS(.radar_weights_path(model))
  } else {
    ans <- readRDS(model)
  }
  if (!inherits(ans, "radar_model") || !identical(ans$format_version, 1L) || !is.list(ans$weights))
    stop("This file is not a supported native RADAR R model.", call. = FALSE)
  if (!is.null(architecture) && !identical(architecture, ans$architecture))
    stop("The model architecture does not match the requested architecture.", call. = FALSE)
  if (compile) ans <- compile_model(ans, learning_rate = ans$learning_rate)
  ans
}
#' @rdname load_model
#' @export
load_pretrained <- function(architecture = "ensemble") load_model(match.arg(architecture, .radar_architectures))

#' Save a native R RADAR model
#' @param model A `radar_model`.
#' @param path Destination RDS file.
#' @param overwrite Replace an existing file.
#' @return `path`, invisibly.
#' @export
#' @examples
#' model <- load_model("standard")
#'
#' path <- file.path(tempdir(), "radar-model.rds")
#' save_model(model, path)
#' restored <- load_model(path)
#' identical(class(restored), class(model))
#' unlink(path)
save_model <- function(model, path, overwrite = TRUE) {
  if (!inherits(model, "radar_model")) stop("Expected a radar_model.", call. = FALSE)
  if (!grepl("\\.rds$", path, ignore.case = TRUE))
    stop("Native R models use RDS serialization; choose a .rds path.", call. = FALSE)
  if (file.exists(path) && !overwrite) stop("The destination exists and overwrite is FALSE.", call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  model$saved_with <- list(RADAR = "0.1.0", R = as.character(getRversion()), platform = R.version$platform)
  model$date_saved <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  saveRDS(model, path, compress = "xz"); invisible(path)
}
#' Inspect a model's architecture and provenance
#' @param model Model name, path, or `radar_model`.
#' @param show Print a readable summary.
#' @return A list describing the model and its training.
#' @export
#' @examples
#' # Accepts a bundled name or a radar_model object.
#' info <- model_info("standard")
#' str(info, max.level = 1)
#'
#' # Printed summary instead of a returned list.
#' model_info("ensemble", show = TRUE)
model_info <- function(model = "ensemble", show = FALSE) {
  m <- load_model(model, compile = FALSE)
  ans <- m[c("architecture", "hyperparameters", "n_parameters", "origin", "training", "saved_with", "date_saved", "backend", "format_version")]
  ans$compatibility <- character()
  if (show) { print(m); if (!is.null(m$training)) print(m$training) }
  ans
}
