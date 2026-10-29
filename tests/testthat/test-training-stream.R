test_that("training stream checks pass", {
  # A real pretrained network can be fine-tuned using the native synthetic stream.
  # Checkpoint provenance must describe the checkpoint's weights, not the ancestor.
  initial <- load_model("standard")
  checkpoint <- tempfile(fileext = ".rds")
  called <- 0L
  callback <- function(epoch, model, logs) {
    expect_true(all(epoch == 1L, is.finite(logs$loss), is.finite(logs$val_loss),
              model$training$epochs == 1L, model$training$spectra_seen == 1L))
    called <<- called + 1L
    invisible(NULL)
  }
  trained <- train(model = initial, epochs = 1L, steps_per_epoch = 1L, batch_size = 1L,
    length = 32L, min_width = 2L, max_width = 16L, max_num_peaks = 2L,
    validation_data = "synthetic", validation_size = 1L, validation_length = 32L,
    learning_rate = 1e-6, callbacks = list(callback), checkpoint_path = checkpoint,
    seed = 47L, verbose = 0L)
  expect_true(all(called == 1L, is.finite(trained$history$val_loss), file.exists(checkpoint),
            trained$training$spectra_seen == 1L,
            trained$training$started_from$spectra_seen == 240000000L,
            any(trained$weights$enc_conv1$kernel != initial$weights$enc_conv1$kernel)))
  back <- load_model(checkpoint); unlink(checkpoint)
  expect_true(all(identical(back$training, trained$training),
            identical(back$weights, trained$weights), identical(back$history, trained$history)))
  # Chained fine-tuning preserves exactly one ancestor per training invocation.
  continuation <- train(generate_dataset(1L, length = 32L, seed = 9L), model = trained,
    epochs = 1L, batch_size = 1L, learning_rate = 1e-7, verbose = 0L)
  expect_true(all(identical(continuation$training$started_from, trained$training)))
})

