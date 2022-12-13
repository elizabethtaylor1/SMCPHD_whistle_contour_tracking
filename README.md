
# SMC-PHD for whistle contour tracking
[![DOI](https://zenodo.org/badge/284840553.svg)](https://zenodo.org/badge/latestdoi/284840553)

A Matlab package for tracking dolphin whistle contours.

Copyright (c) 2020, Pina Gruden

This package provides an implementation of the Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD) filter that uses radial basis function network motion model to track dolphin whistle contours.

See reference paper: Gruden, P. and White, P. (2020) Automated extraction of dolphin whistles - a Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD) approach; Accepted for publication in the Journal of the Acoustical Society of America (28 September 2020).

Steps:
1. Train RBF network (RBF_netork_training folder). Run prepare input/output first, then RBFnet_TrainTest_v2.
2. Develop GMM (build_gmm.m)
3. Run train_parameters to get birthpdf.m, requires picking to have informed or not informed priors
4. Run SMCPHD (Run_SMCPHD.m).