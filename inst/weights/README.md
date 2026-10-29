# Bundled pretrained models

standard.rds and ensemble.rds contain complete native radar_model objects,
including architectures, hyperparameters and trained numerical weights.
load_model() reads these installed files directly. No download or conversion
is needed for installation or use. Standard: 585,299 parameters; ensemble:
280,592 parameters. Training used 240 million synthetic spectra of length 750.

Architectures and pretrained weights: Copyright (c) 2025 Joel Sjöberg.
R implementation: Copyright (c) 2026 Bogdan Oancea. BSD 3-Clause license.
