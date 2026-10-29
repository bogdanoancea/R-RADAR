# radar.raman

An R implementation of **RADAR** - **R**aman spectral **A**nalysis using
**D**eep learning for **A**rtifact **R**emoval - and a port of the Python
`radar-raman` package.

RADAR decomposes Raman spectra into baseline, cosmic rays, noise and Raman
peaks. Both network architectures and pretrained weight sets are included.

Inference, differentiation, Adam training, synthetic data generation and
evaluation are implemented in R. Core functionality requires only R; optional R
packages support particular file formats.

## Installation

From R, directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("bogdanoancea/R-RADAR")

# With the vignettes built; needs knitr, rmarkdown and Pandoc.
remotes::install_github("bogdanoancea/R-RADAR", build_vignettes = TRUE)
```

From a local clone, in the terminal:

```sh
R CMD INSTALL .
```

Or install a source archive from R:

```r
install.packages("radar.raman_0.1.0.tar.gz", repos = NULL, type = "source")
```

R 4.1.0 or later is required. No compiler is needed.

## Complete guide

The comprehensive vignette covers the R package, including file formats, models, training, synthetic data, evaluation,
the API reference, command-line tools, troubleshooting, history and credits.

```r
vignette("RADAR", package = "radar.raman")
browseVignettes("radar.raman")
```

Its source is [vignettes/RADAR.Rmd](vignettes/RADAR.Rmd). Building the
vignette requires knitr, rmarkdown and Pandoc; the built source archive
includes the rendered HTML guide.

## Network architectures

The [architecture guide](vignettes/architectures.Rmd) describes Joel Sjöberg's
standard and ensemble networks, with diagrams, layer shapes, output heads,
residual connections, parameter counts and training behavior.

```r
vignette("architectures", package = "radar.raman")
```

## Decompose spectra

```r
library(radar.raman)

# A vector is one spectrum; each row of a matrix is a spectrum.
test <- generate(4, length = 256, seed = 7)
parts <- radar.raman::decompose(test$X, model = "standard")
parts$baseline
parts$cosmic_rays
parts$noise
parts$peaks                     # corrected, denoised Raman bands
plot(parts, index = 1)

# Derived spectra and reconstruction are computed on access.
parts$denoised
parts$despiked
parts$baseline_corrected
parts$reconstruction
parts$residual

# Files can be passed directly; their axes and map metadata are retained.
# parts <- radar.raman::decompose("measurement.spc", output = "components.rds")
radar.raman::write(parts, "components.rds")
back <- load_components("components.rds")
```

The standard model computes noise as the residual, so its four components sum
to the input, up to numerical precision. The ensemble model predicts all four
components independently; its reconstruction can differ from the input.

Inputs are normalized per spectrum before inference, and components are returned
in original intensity units by default. The offset is restored to baseline.
Use `units = "normalized"` for scaled components, or `normalize = FALSE` for
inputs already in the training scale. Unequal-length spectra can be supplied
as a list; the result is then a list of decompositions.

## Spectral files and maps

```r
formats(show = TRUE)
path <- system.file("extdata", "synthetic_spectra.csv", package = "radar.raman")
data <- radar.raman::read(path)
radar.raman::write(data, "converted.jdx")

axis <- seq(400, 1800, length.out = 64)
cube <- array(runif(2 * 3 * 64), dim = c(2, 3, 64))
map <- radar.raman::decompose(cube, model = "standard", wavenumbers = axis)
dim(as_map(map, "peaks"))       # 2 3 64
band_image(map, 900, 1100)      # integrated Raman band, 2 x 3
```

The last array dimension is spectral. Spatial axes follow R's usual
column-major array order. Indices are one-based.

Core formats include CSV/TSV/text, RDS, SPC, JCAMP-DX, MATLAB and NumPy numeric
files. Optional packages add JSON (`jsonlite`), HDF5/NetCDF-4 (`hdf5r`), SQLite
(`DBI`, `RSQLite`), Excel (`readxl`, `openxlsx`), Parquet/Feather (`arrow`) and
vendor XML (`xml2`). Renishaw WDF is read natively. Use `formats()` for the
capabilities and dependencies of each adapter. The vignette describes format limitations.

## Models, training and evaluation

```r
model <- load_model("ensemble")
model_info(model, show = TRUE)
save_model(model, "model.rds")
model <- load_model("model.rds")

scores <- evaluate("standard", test)
scores$peaks
evaluate("standard", test, statistics = TRUE)

# A small training run. Every network layer is differentiated in R.
fresh <- new_model("standard", scaler = 2, dim_red_size = 16, seed = 3)
fitted <- train(test, model = fresh, epochs = 1, batch_size = 2,
                validation_data = NULL, learning_rate = 1e-5, verbose = 0)
fitted$history
save_model(fitted, "trained.rds")

# Fine-tune a bundled network, for example:
# fitted <- train(test, model = "standard", init = "pretrained",
#                 epochs = 1, learning_rate = 1e-6, validation_data = NULL)
```

Omitting training data creates fresh synthetic batches. The original
12-level learning-rate schedule is available through
`original_learning_rate_schedule()`. Model files use RDS and retain weights,
architecture, provenance and training history. Training is CPU based; large
runs require substantial computing time. Synthetic data generation supports
reproducible seeds and optional legacy scaling.

## Command line

The installed package includes an Rscript entry point:

```sh
Rscript /path/to/installed/radar.raman/exec/radar.R --help
Rscript /path/to/installed/radar.raman/exec/radar.R generate -n 4 --length 128 -o test.rds
Rscript /path/to/installed/radar.raman/exec/radar.R decompose measurement.spc -m standard -o result.rds
```

Find the entry point with `system.file("exec", "radar.R", package = "radar.raman")`,
or call `radar_cli(c("formats"))` from R.

## Validation and attribution

Tests cover fixed numerical references, both architectures, gradients, training,
normalization, file round trips and maps. Run `R CMD check` on the source archive.
Run the test suite in RStudio with **Build > Test Package**, or from R:

```r
devtools::test()
```

The suite uses `testthat` and also runs during `R CMD check`.
See AUTHORS.md and LICENSE for credits and `citation("radar.raman")` for the method.

### Provenance

This package is a port of the Python
[radar-raman](https://github.com/bogdanoancea/RADAR) package
([PyPI](https://pypi.org/project/radar-raman/),
[doi:10.5281/zenodo.23760379](https://doi.org/10.5281/zenodo.23760379)), which
in turn builds on Joel Sjoberg's [original
implementation](https://github.com/JoelSjoberg/Comprehensive-Raman-Spectral-Analysis-Using-Deep-Learning).
The network architectures, the synthetic data generator and the pretrained
weights all originate there.

The expected outputs in `inst/extdata/reference_outputs.rds` were produced by
that original code running under legacy Keras 2, so the regression test in
`tests/testthat/test-reference.R` checks this port against the original
implementation rather than against itself. See `inst/extdata/README.md` for the
full chain.
