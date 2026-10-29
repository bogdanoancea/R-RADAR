# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjöberg (original losses and training schedule)
# Copyright (c) 2026 Bogdan Oancea (RADAR R package)

#' Original RADAR learning-rate levels
#' @export
ORIGINAL_LR_VALUES <- 1.5 * c(7e-4, 5e-4, 3e-4, 1e-4, 7e-5, 5e-5, 3e-5, 1e-5, 7e-6, 5e-6, 3e-6, 1e-6)

#' Original twelve-level learning-rate schedule
#' @param total_steps Positive total number of optimizer steps.
#' @return A function accepting a zero-based optimizer step.
#' @export
original_learning_rate_schedule <- function(total_steps) {
  if (length(total_steps) != 1L || !is.finite(total_steps) || total_steps < 1)
    stop("total_steps must be positive.", call. = FALSE)
  boundaries <- if (total_steps < 12) numeric() else round(total_steps * seq_len(11L) / 12)
  function(step) ORIGINAL_LR_VALUES[1L + sum(step > boundaries)]
}

.radar_loss_inputs <- function(y_true, y_pred) {
  if (!is.numeric(y_true) || !is.numeric(y_pred) || length(y_true) != length(y_pred) ||
      !length(y_true) || any(!is.finite(y_true)) || any(!is.finite(y_pred)) ||
      !identical(dim(y_true), dim(y_pred)))
    stop("Targets and predictions must be finite numeric arrays of equal shape.", call. = FALSE)
}
#' Neural-network losses
#' @param y_true,y_pred Finite numeric arrays of the same shape.
#' @return A scalar loss or coefficient of determination.
#' @export
root_mean_squared_error <- function(y_true, y_pred) {
  .radar_loss_inputs(y_true, y_pred); sqrt(mean((y_pred - y_true)^2))
}
#' @rdname root_mean_squared_error
#' @export
r_square <- function(y_true, y_pred) {
  .radar_loss_inputs(y_true, y_pred)
  1 - sum((y_true - y_pred)^2) / (sum((y_true - mean(y_true))^2) + 1e-7)
}
#' @rdname root_mean_squared_error
#' @export
r_square_loss <- function(y_true, y_pred) {
  .radar_loss_inputs(y_true, y_pred)
  sum((y_true - y_pred)^2) / (sum((y_true - mean(y_true))^2) + 1e-7)
}
#' @rdname root_mean_squared_error
#' @export
joined_loss <- function(y_true, y_pred) root_mean_squared_error(y_true, y_pred) + r_square_loss(y_true, y_pred)

.radar_loss_and_grad <- function(model, x, targets, gradient = TRUE, training = gradient) {
  run <- .radar_forward(model, x, training = training, gradient = gradient, loss = TRUE)
  ctx <- run$context; values <- numeric(4L); names(values) <- .radar_components
  for (name in .radar_components) {
    node <- run$outputs[[name]]; pred <- matrix(node$value, nrow(x), ncol(x))
    truth <- targets[[name]]; err <- pred - truth
    rmse <- sqrt(mean(err * err)); values[name] <- rmse
    g <- if (rmse == 0) err * 0 else err / (length(err) * rmse)
    if (model$architecture == "ensemble") {
      denom <- sum((truth - mean(truth))^2) + 1e-7
      values[name] <- values[name] + sum(err * err) / denom
      g <- g + 2 * err / denom
    }
    if (gradient) .radar_addgrad(node, array(g, dim(node$value)))
  }
  loss <- sum(values) + ctx$reg_loss
  if (gradient) {
    for (j in rev(seq_along(ctx$nodes))) {
      node <- ctx$nodes[[j]]
      if (!is.null(node$grad) && !is.null(node$backward)) node$backward(node$grad)
    }
  }
  list(loss = loss, component_loss = values, gradients = ctx$grads, weights = ctx$weights,
       predictions = lapply(run$outputs, function(node) matrix(node$value, nrow(x), ncol(x))))
}

.radar_training_data <- function(data) {
  if (is.character(data) && length(data) == 1L) data <- load_dataset(data)
  if (is.list(data) && !is.null(data$X)) {
    x <- data$X
    targets <- if (!is.null(data$targets)) data$targets else unclass(data)[.radar_components]
  } else if (is.list(data) && length(data) == 2L) {
    x <- data[[1L]]; targets <- data[[2L]]
  } else stop("data must be a synthetic data set, an (X, targets) list, or a data-set path.", call. = FALSE)
  if (length(dim(x)) == 3L && dim(x)[3L] == 1L) x <- matrix(x, dim(x)[1L], dim(x)[2L])
  if (is.null(dim(x))) x <- matrix(x, 1L)
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x) || any(!is.finite(x)))
    stop("Training X must be a nonempty finite numeric matrix.", call. = FALSE)
  if (!is.list(targets) || length(targets) != 4L) stop("Four component targets are required.", call. = FALSE)
  if (is.null(names(targets))) names(targets) <- .radar_components
  if (!all(.radar_components %in% names(targets))) stop("Targets must contain baseline, cosmic_rays, noise and peaks.", call. = FALSE)
  targets <- lapply(targets[.radar_components], function(y) {
    if (is.null(dim(y)) && nrow(x) == 1L) y <- matrix(y, 1L)
    if (!is.numeric(y) || !identical(dim(y), dim(x)) || any(!is.finite(y)))
      stop("Each target must be a finite numeric matrix with the same dimensions as X.", call. = FALSE)
    y
  })
  list(X = x, targets = targets)
}

.radar_adam <- function(model, gradients, state, rate) {
  state$step <- state$step + 1L
  for (layer in names(gradients)) for (key in names(gradients[[layer]])) {
    g <- gradients[[layer]][[key]]
    if (any(!is.finite(g))) stop("Non-finite gradients; lower the learning rate or inspect the data.", call. = FALSE)
    if (is.null(state$m[[layer]][[key]])) {
      state$m[[layer]][[key]] <- g * 0; state$v[[layer]][[key]] <- g * 0
    }
    state$m[[layer]][[key]] <- 0.9 * state$m[[layer]][[key]] + 0.1 * g
    state$v[[layer]][[key]] <- 0.999 * state$v[[layer]][[key]] + 0.001 * g * g
    alpha <- rate * sqrt(1 - 0.999^state$step) / (1 - 0.9^state$step)
    model$weights[[layer]][[key]] <- model$weights[[layer]][[key]] -
      alpha * state$m[[layer]][[key]] / (sqrt(state$v[[layer]][[key]]) + 1e-7)
  }
  list(model = model, state = state)
}

#' Train or fine-tune a RADAR network entirely in R
#'
#' Reverse-mode differentiation covers every convolution, learned key, outer
#' product, batch-normalization layer, residual connection and component head.
#' Training uses Adam, the original component losses and L1/L2 regularization.
#' Batch-normalization uses batch/spatial statistics and updates moving moments.
#' This CPU implementation favors portability over the speed of a GPU runtime;
#' a small batch size is recommended because the attention layer materializes
#' convolution patches. The default batch size is one to bound memory use.
#'
#' @param data A synthetic data set, `list(X, targets)`, a saved data-set path,
#'   or `NULL` for freshly generated synthetic batches.
#' @param model An architecture name or an existing `radar_model`.
#' @param init `NULL` for random weights, `"pretrained"`, or a saved model path.
#' @param save_to Optional destination `.rds` model file.
#' @param epochs Number of passes over the training data.
#' @param steps_per_epoch Synthetic batches per epoch; ignored for fixed data.
#' @param batch_size Spectra per optimizer step.
#' @param learning_rate A positive number, a function of step, or `"original"`.
#' @param length,min_width,max_width,max_num_peaks,increase_noise,legacy_scaling
#'   Synthetic-data options passed to [generate_dataset()].
#' @param validation_data `NULL`, `"synthetic"`, or a validation data set.
#' @param validation_size,validation_length Size and length of generated validation data.
#' @param seed Random seed, or `NULL` to retain the current random state.
#' @param callbacks List of functions called with `(epoch, model, logs)` after
#'   each epoch. A returned `radar_model` replaces the current model.
#' @param checkpoint_path Optional `.rds` path for the model with the best loss.
#' @param verbose Print epoch losses when positive.
#' @return The trained `radar_model`; its `history` contains epoch losses.
#' @export
train <- function(data = NULL, model = "ensemble", init = NULL, save_to = NULL,
                  epochs = 10L, steps_per_epoch = 1000L, batch_size = 1L,
                  learning_rate = "original", length = 750L, min_width = 5L,
                  max_width = 256L, max_num_peaks = 5L, increase_noise = FALSE,
                  legacy_scaling = FALSE, validation_data = NULL, validation_size = 100L,
                  validation_length = 2000L, seed = 0, callbacks = NULL,
                  checkpoint_path = NULL, verbose = 1L) {
  for (z in c("epochs", "steps_per_epoch", "batch_size")) {
    v <- get(z)
    if (length(v) != 1L || !is.finite(v) || v < 1L || v != as.integer(v)) stop(z, " must be a positive integer.", call. = FALSE)
  }
  if (inherits(model, "radar_model")) m <- model
  else if (is.null(init)) m <- build_model(model, seed = seed)
  else if (identical(init, "pretrained")) m <- load_model(model)
  else m <- load_model(init, architecture = model)
  if (!is.null(seed)) set.seed(seed)
  fixed <- if (is.null(data)) NULL else .radar_training_data(data)
  if (!is.null(fixed)) steps_per_epoch <- ceiling(nrow(fixed$X) / batch_size)
  total_steps <- epochs * steps_per_epoch
  lr_record <- if (is.function(learning_rate)) "R schedule function" else learning_rate
  if (is.character(learning_rate)) {
    if (!identical(learning_rate, "original")) stop("learning_rate must be 'original', positive, or a schedule function.", call. = FALSE)
    learning_rate <- original_learning_rate_schedule(total_steps)
  }
  m <- compile_model(m, learning_rate = learning_rate)
  if (identical(validation_data, "synthetic")) {
    validation_data <- generate_dataset(validation_size, length = validation_length,
      min_width = 5L, max_width = 300L, max_num_peaks = 5L, increase_noise = increase_noise,
      legacy_scaling = legacy_scaling, seed = if (is.null(seed)) NULL else seed + 12345L)
  }
  validation <- if (is.null(validation_data)) NULL else .radar_training_data(validation_data)
  if (!is.null(callbacks) && (!is.list(callbacks) || !all(vapply(callbacks, is.function, logical(1)))))
    stop("callbacks must be a list of R functions.", call. = FALSE)
  state <- list(step = 0L, m = list(), v = list())
  history <- list(loss = numeric(epochs), val_loss = if (!is.null(validation)) numeric(epochs) else NULL,
                  component_loss = matrix(0, epochs, 4L, dimnames = list(NULL, .radar_components)))
  best <- Inf; spectra_seen <- 0L
  prior <- m$training
  m$origin <- if (is.null(prior)) "trained in R" else paste0(m$origin, "; fine-tuned in R")
  for (epoch in seq_len(epochs)) {
    losses <- numeric(steps_per_epoch); part <- matrix(0, steps_per_epoch, 4L); sizes <- integer(steps_per_epoch)
    if (!is.null(fixed)) order <- sample.int(nrow(fixed$X))
    for (step in seq_len(steps_per_epoch)) {
      if (is.null(fixed)) {
        batch <- .radar_training_data(generate_dataset(batch_size, length = length,
          min_width = min_width, max_width = max_width, max_num_peaks = max_num_peaks,
          increase_noise = increase_noise, legacy_scaling = legacy_scaling, seed = NULL))
      } else {
        idx <- order[seq.int((step - 1L) * batch_size + 1L, min(step * batch_size, nrow(fixed$X)))]
        batch <- list(X = fixed$X[idx, , drop = FALSE], targets = lapply(fixed$targets, function(y) y[idx, , drop = FALSE]))
      }
      run <- .radar_loss_and_grad(m, batch$X, batch$targets)
      if (!is.finite(run$loss)) stop("Non-finite training loss; inspect the data and learning rate.", call. = FALSE)
      m$weights <- run$weights
      rate <- if (is.function(learning_rate)) learning_rate(state$step) else learning_rate
      if (length(rate) != 1L || !is.finite(rate) || rate <= 0) stop("The learning-rate schedule returned an invalid value.", call. = FALSE)
      update <- .radar_adam(m, run$gradients, state, rate); m <- update$model; state <- update$state
      losses[step] <- run$loss; part[step, ] <- run$component_loss; sizes[step] <- nrow(batch$X)
      spectra_seen <- spectra_seen + nrow(batch$X)
    }
    history$loss[epoch] <- stats::weighted.mean(losses, sizes)
    history$component_loss[epoch, ] <- colSums(part * sizes) / sum(sizes)
    if (!is.null(validation)) {
      # Evaluate the complete loss (including kernel, bias and activity
      vb <- split(seq_len(nrow(validation$X)), ceiling(seq_len(nrow(validation$X)) / batch_size))
      vl <- vapply(vb, function(idx) {
        run <- .radar_loss_and_grad(m, validation$X[idx, , drop = FALSE],
          lapply(validation$targets, function(y) y[idx, , drop = FALSE]), gradient = FALSE, training = FALSE)
        run$loss
      }, numeric(1))
      history$val_loss[epoch] <- stats::weighted.mean(vl, lengths(vb))
    }
    logs <- list(loss = history$loss[epoch])
    if (!is.null(validation)) logs$val_loss <- history$val_loss[epoch]
    m$history <- lapply(history, function(z) if (is.matrix(z)) z[seq_len(epoch), , drop = FALSE] else if (is.null(z)) NULL else z[seq_len(epoch)])
    # Record provenance before writing a checkpoint: it describes exactly the
    # weights saved at this epoch, with one stable pre-training ancestor.
    m$training <- list(epochs = epoch, steps_per_epoch = steps_per_epoch, batch_size = batch_size,
      learning_rate = lr_record, spectra_seen = spectra_seen,
      data = if (is.null(fixed)) paste0("synthetic stream; length ", length) else paste0("fixed data set of ", nrow(fixed$X), " spectra"),
      final_loss = history$loss[epoch], backend = "R reverse-mode differentiation and Adam")
    if (!is.null(prior)) m$training$started_from <- prior
    monitor <- if (!is.null(validation)) logs$val_loss else logs$loss
    if (!is.null(checkpoint_path) && monitor < best) { save_model(m, checkpoint_path); best <- monitor }
    if (!is.null(callbacks)) for (cb in callbacks) {
      replacement <- cb(epoch, m, logs)
      if (inherits(replacement, "radar_model")) m <- replacement
    }
    if (verbose > 0L) message("Epoch ", epoch, "/", epochs, ": loss = ", signif(logs$loss, 6L),
      if (!is.null(logs$val_loss)) paste0(", validation loss = ", signif(logs$val_loss, 6L)) else "")
  }
  m$history <- history
  if (!is.null(save_to)) save_model(m, save_to)
  m
}
