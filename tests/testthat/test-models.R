test_that("models checks pass", {
  # Check architecture sizes, inference length preservation and original weights.
  expect_true(all(build_standard_model()$n_parameters == 585299,
            build_ensemble_model()$n_parameters == 280592))
  for (architecture in c("standard", "ensemble")) {
    model <- load_model(architecture)
    input <- matrix(seq(0, 1, length.out = 16L), 1L)
    pred <- predict(model, input)
    expect_true(all(identical(names(pred), c("baseline", "cosmic_rays", "noise", "peaks")),
              all(vapply(pred, function(z) identical(dim(z), c(1L, 16L)) && all(is.finite(z)), logical(1)))))
    if (architecture == "standard") expect_true(all(max(abs(Reduce(`+`, pred) - input)) < 1e-12))
  }
  # Finite differences independently verify gradients through all parameter groups,
  # including 2D kernels, max-pooled keys and the two-stage residual connection.
  loss_and_grad <- getFromNamespace(".radar_loss_and_grad", "radar.raman")
  components <- c("baseline", "cosmic_rays", "noise", "peaks")
  for (architecture in c("standard", "ensemble")) {
    model <- build_model(architecture, scaler = 2L, dim_red_size = 16L, seed = 16L)
    set.seed(5L); x <- matrix(runif(16L), 2L, 8L)
    targets <- setNames(replicate(4L, matrix(runif(16L), 2L, 8L), simplify = FALSE), components)
    analytic <- loss_and_grad(model, x, targets)
    for (layer in names(model$weights)) {
      for (key in intersect(names(model$weights[[layer]]), c("kernel", "bias", "gamma", "beta"))) {
        gradient <- analytic$gradients[[layer]][[key]]
        idx <- which.max(abs(gradient)); epsilon <- 1e-6
        plus <- minus <- model
        plus$weights[[layer]][[key]][idx] <- plus$weights[[layer]][[key]][idx] + epsilon
        minus$weights[[layer]][[key]][idx] <- minus$weights[[layer]][[key]][idx] - epsilon
        numerical <- (loss_and_grad(plus, x, targets)$loss - loss_and_grad(minus, x, targets)$loss) / (2 * epsilon)
        expect_true(all(abs(numerical - gradient[idx]) < 1e-4))
      }
    }
    trained <- train(list(x, targets), model = model, epochs = 3L, batch_size = 2L,
                     learning_rate = 1e-6, validation_data = NULL, verbose = 0L)
    expect_true(all(trained$history$loss[3L] < trained$history$loss[1L],
              any(trained$weights$enc_conv1$kernel != model$weights$enc_conv1$kernel),
              any(trained$weights$key_conv$kernel != model$weights$key_conv$kernel)))
    tmp <- tempfile(fileext = ".rds"); save_model(trained, tmp)
    reloaded <- load_model(tmp); unlink(tmp)
    expect_true(all(identical(predict(trained, x), predict(reloaded, x))))
  }
  new_node <- getFromNamespace(".radar_node", "radar.raman")
  max_pool <- getFromNamespace(".radar_maxpool", "radar.raman")
  ctx <- new.env(parent = emptyenv()); ctx$gradient <- TRUE; ctx$nodes <- list()
  node <- new_node(ctx, array(c(2, 2, 1, 3, 0, 3), c(1L, 3L, 2L)))
  pooled <- max_pool(ctx, node); pooled$backward(matrix(c(4, 6), 1L, 2L))
  expect_true(all(identical(as.vector(node$grad), c(2, 2, 0, 3, 0, 3))))
})

