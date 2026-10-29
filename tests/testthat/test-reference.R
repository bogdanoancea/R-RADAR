# Numerical regression against the ORIGINAL implementation.
#
# Provenance of inst/extdata/reference_outputs.rds
# ------------------------------------------------
# Joel Sjoberg's original repository code (Scripts/essentials.py,
# Scripts/generator.py, Models/*.h5) was run under legacy Keras 2
# (tf-keras 2.21, TF_USE_LEGACY_KERAS=1) to produce reference_outputs.npz,
# shipped in the Python radar-raman package as tests/data/reference_outputs.npz.
# This .rds is a direct conversion of that file: all 18 arrays agree with the
# .npz to 2.2e-16, i.e. float64 epsilon. The expected values are therefore
# outputs of the original implementation and not of this R port, so the checks
# below are genuine port-fidelity tests rather than self-comparisons.
#
#   X1024  6 synthetic spectra of length 1024 from the original
#          train_generator(length = 1024, batch_size = 6, min_width = 10,
#          max_width = 300, increase_noise = True) under numpy seed 3.
#   X1738  3 normalized random-walk curves of length 1738, the length of the
#          live-cell spectra, covering a second input length.
#
# Tolerance
# ---------
# 1e-5 absolute, deliberately identical to the atol used by the Python
# package's test_matches_original_implementation, so both ports are held to
# the same bar. The original stored float32 predictions, so agreement is
# bounded by float32 rounding accumulated through the network. Observed worst
# case is ~6.2e-6 on ensemble X1024 peaks; every other component is at or
# below ~2.3e-6. The tolerance is not tightened because the float32 reference
# does not support a finer claim.
#
# Degenerate references
# ---------------------
# Three of the sixteen component references are identically zero: X1738 is a
# random walk, so it contains no cosmic-ray spikes, and the ensemble model
# resolves no Raman bands in it. Comparing zero against zero cannot detect a
# regression, so those three are named explicitly and every other component is
# additionally required to carry signal on both sides. A head that silently
# begins returning zeros is then a failure instead of a vacuous pass.

DEGENERATE_REFERENCES <- c("standard_X1738_cosmic_rays",
                           "ensemble_X1738_cosmic_rays",
                           "ensemble_X1738_peaks")

test_that("predictions reproduce the original implementation", {

  reference <- readRDS(system.file("extdata", "reference_outputs.rds",
                                   package = "radar.raman", mustWork = TRUE))

  for (arch in c("standard", "ensemble")) {
    net <- load_model(arch)
    for (key in c("X1024", "X1738")) {
      predicted <- stats::predict(net, reference[[key]])
      for (component in COMPONENTS) {

        label <- paste(arch, key, component, sep = "_")
        expected <- reference[[label]]
        got <- predicted[[component]]

        expect_true(all(is.finite(got)), info = label)
        expect_lt(max(abs(got - expected)), 1e-5, label = label)

        if (label %in% DEGENERATE_REFERENCES) {
          # Assert the documented degeneracy still holds, so the list above
          # cannot go stale without the suite saying so.
          expect_equal(max(abs(expected)), 0,
                       label = paste0(label, " (reference magnitude)"))
        } else {
          # Require signal on both sides, so the comparison above cannot pass
          # merely because reference and prediction are both zero.
          expect_gt(max(abs(expected)), 0,
                    label = paste0(label, " (reference magnitude)"))
          expect_gt(max(abs(got)), 0,
                    label = paste0(label, " (prediction magnitude)"))
        }
      }
    }
  }
})
