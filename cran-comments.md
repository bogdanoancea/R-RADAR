## Test environments

* local: macOS 26.6.2 (aarch64-apple-darwin23), R 4.6.1
* win-builder: Windows Server 2022 x64, R Under development (unstable)
  (2026-11-12 r90881 ucrt)
* GitHub Actions: ubuntu-latest, macos-latest, windows-latest, R release

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes

The NOTE is "New submission", which is expected for a first submission.

On win-builder the same NOTE also lists possibly misspelled words in
DESCRIPTION:

    Raman (3:8, 12:25, 13:5)
    convolutional (13:68)
    pretrained (13:35)

All three are spelled correctly. "Raman" is the surname of C. V. Raman, after
whom Raman spectroscopy and the Raman effect are named. "Convolutional" and
"pretrained" are standard terms in the neural network literature.

## Other comments

This package is an R port of the Python package radar-raman. It contains no
compiled code.

`inst/weights` holds the two pretrained convolutional networks as native R
objects, about 3 MB in total. They are bundled deliberately: the package
performs inference out of the box, with no download step and no conversion
from another framework's format.

The regression test in `tests/testthat/test-reference.R` compares this
implementation's predictions against fixed outputs of the original
implementation, shipped in `inst/extdata/reference_outputs.rds`. Agreement is
at float32 rounding level, and the 1e-5 tolerance matches the one used by the
Python package, so both ports are held to the same bar.

Examples and tests depend only on the hard dependencies. The optional packages
in Suggests enable additional spectral file formats, and both the examples and
the test suite narrow their scope when those packages are absent rather than
failing or skipping wholesale.
