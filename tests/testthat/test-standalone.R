test_that("standalone checks pass", {
  # Both complete pretrained models must load from installed package files alone.
  for (architecture in c("standard", "ensemble")) {
    bundled <- system.file("weights", paste0(architecture, ".rds"),
                           package = "RADAR", mustWork = TRUE)
    stored <- readRDS(bundled)
    model <- load_model(architecture, compile = FALSE)
    expect_true(all(identical(stored$weights, model$weights),
              identical(model$backend, "pure R"),
              identical(model$architecture, architecture)))
    path <- tempfile(fileext = ".rds")
    save_model(model, path)
    restored <- load_model(path, compile = FALSE)
    expect_true(all(identical(restored$weights, model$weights)))
    unlink(path)
  }
  # Model persistence is restricted to native R serialization.
  bad <- tempfile(fileext = ".txt")
  expect_true(all(inherits(try(save_model(model, bad), silent = TRUE), "try-error"),
            !file.exists(bad)))
})

