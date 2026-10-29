test_that("reference checks pass", {
  
  reference <- readRDS(system.file("extdata", "reference_outputs.rds", package = "RADAR"))
  for (arch in c("standard", "ensemble")) {
    net <- load_model(arch)
    for (key in c("X1024", "X1738")) {
      predicted <- stats::predict(net, reference[[key]])
      for (component in COMPONENTS) {
        expected <- reference[[paste(arch, key, component, sep = "_")]]
        error <- max(abs(predicted[[component]] - expected))
        cat(sprintf("%s %s %s: max absolute error %.9g\n", arch, key, component, error))
        expect_true(all(is.finite(error), error < 1e-5))
      }
    }
  }
})

