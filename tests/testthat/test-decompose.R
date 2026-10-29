test_that("decompose checks pass", {
  
  fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
  for (arch in c("standard", "ensemble")) {
    net <- load_model(arch)
    x <- matrix(seq(0, 1, length.out = 128), 2, 64, byrow = TRUE)
    p <- radar.raman::decompose(x, model = net)
    expect_true(all(inherits(p, "radar_components"), identical(dim(p$peaks), dim(x)),
              all(p$peaks >= 0), all(p$cosmic_rays >= 0),
              identical(p$model, arch), p$units == "original"))
    one <- radar.raman::decompose(x[1, ], model = net)
    expect_true(all(is.null(dim(one$peaks)), length(one$peaks) == 64))
    ragged <- radar.raman::decompose(list(x[1, ], x[2, 1:40]), model = net)
    expect_true(all(length(ragged) == 2, ragged[[2]]$length == 40))
    file <- tempfile(fileext = ".rds")
    save_model(net, file)
    same <- stats::predict(load_model(file), x)
    raw <- stats::predict(net, x)
    expect_true(all(isTRUE(all.equal(raw, same, tolerance = 1e-12))))
    unlink(file)
  }
  
  x <- matrix(seq(-1, 4, length.out = 128), 2, 64)
  p <- radar.raman::decompose(x, model = "standard")
  q <- radar.raman::decompose(x * 1234.5 + 300, model = "standard")
  expect_true(all(max(abs(p$reconstruction - x)) < 1e-10,
            max(abs(q$peaks - p$peaks * 1234.5)) < 1e-6,
            max(abs(q$baseline - p$baseline * 1234.5 - 300)) < 1e-6))
  flat <- radar.raman::decompose(rep(5, 64), model = "standard")
  expect_true(all(all(is.finite(flat$peaks)), max(abs(flat$reconstruction - 5)) < 1e-10))
  expect_true(all(fails(radar.raman::decompose(c(1, NA), model = "standard")),
            fails(radar.raman::decompose(1, model = "standard")),
            fails(radar.raman::decompose(1:20, wavenumbers = 1:10)),
            fails(radar.raman::decompose(matrix(numeric(), 0, 20))),
            fails(radar.raman::decompose(1:20, units = "invalid"))))
  
  cube <- array(seq(0, 1, length.out = 2 * 3 * 64), c(2, 3, 64))
  axis <- seq(400, 1800, length.out = 64)
  map <- radar.raman::decompose(cube, model = "standard", wavenumbers = axis)
  expect_true(all(identical(dim(as_map(map)), c(2L, 3L, 64L)),
            identical(dim(band_image(map, 900, 1100)), c(2L, 3L)),
            max(abs(array(map$spectra, c(2, 3, 64)) - cube)) == 0))
  
  # Geometry of source-library maps also survives the high-level decomposition.
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    file <- system.file("extdata", "interoperability", "fixture-map-no-positions.json",
                        package = "radar.raman", mustWork = TRUE)
    source_map <- radar.raman::read(file)
    parts <- radar.raman::decompose(source_map, model = "standard")
    expect_true(all(max(abs(as_map(parts, "spectra")[2, 1, ] - source_map$spectra[4, ])) < 1e-12))
  }
  
  ds <- generate(2, length = 64, seed = 12)
  scores <- evaluate("standard", ds)
  expect_true(all(all(c("rmse", "mae", "max_error", "r2") %in% names(scores$peaks))))
  stats <- evaluate("standard", ds, statistics = TRUE)
  expect_true(all("max_error_mean" %in% names(stats$peaks)))
  
  path <- tempfile(fileext = ".rds")
  invisible(capture.output(radar_cli(c("generate", "-n", "2", "--length", "64", "-o", path))))
  expect_true(all(file.exists(path), inherits(load_dataset(path), "synthetic_dataset")))
  out <- tempfile(fileext = ".rds")
  invisible(capture.output(radar_cli(c("decompose", path, "-m", "standard", "-o", out))))
  expect_true(all(inherits(load_components(out), "radar_components")))
  unlink(c(path, out))
})

