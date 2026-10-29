# Example and validation data

synthetic_spectra.csv contains five spectra with 1,024 points and an axis from
400 to 1,800 inverse centimetres. synthetic_spectra_truth.rds contains their
known components on a normalized scale. CSV intensities equal 4000 * X + 600.

reference_outputs.rds contains fixed numerical regression inputs and expected
outputs for both architectures at lengths 1,024 and 1,738. Tests use an
absolute tolerance of 1e-5. All fixtures are included in this package.
