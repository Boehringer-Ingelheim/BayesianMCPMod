# Changelog

## BayesianMCPMod 1.3.0 (XX-Feb-2026)

CRAN release: 2026-02-23

- Fixed a bug that would occur when predicting from the beta model shape
  outside of the original dose range.
- Fixed a bug in which the MED assessment could not be performed when
  specifying a negative direction of beneficial effect and an evidence
  level other than 0.5.
- Added functions and vignettes for the binary endpoint case.
- Added functionality to
  [`assessDesign()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/assessDesign.md)
  to provide custom simulated data and custom model estimates enabling
  complex data simulation and analysis methods.
- Added argument to
  [`assessDesign()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/assessDesign.md)
  for number of bootstrap samples in case `evidence_level` is provided.
- Added functionality to
  [`plot.modelFits()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/plot.modelFits.md)
  to plot effect sizes.
- Added calls to [`set.seed()`](https://rdrr.io/r/base/Random.html) in
  vignette’s code blocks to facilitate individual code block
  reproducibility.

## BayesianMCPMod 1.2.0 (28-Aug-2025)

CRAN release: 2025-08-28

- Fixed a bug in
  [`performBayesianMCPMod()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/performBayesianMCPMod.md)
  where the model significance status from the MCP step was sometimes
  not correctly assigned to the fitted model in the Mod step.
- Fixed a bug in `print.modelFit()` where sometimes the coefficients for
  the fitted model shapes were not printed correctly.
- Fixed a bug in
  [`getMED()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getMED.md)
  where quantile and evidence level could sometimes not be matched due
  to floating-point precision issues when using bootstrapped quantiles.
- Changed functions
  [`getPosterior()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getPosterior.md),
  [`getCritProb()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getCritProb.md),
  and
  [`getContr()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getContr.md)
  to accept a covariance matrix instead of a standard deviation vector
  as argument.
- Added support for none-zero off-diagonal covariance matrices in the
  MCP step.
- Added bootstrapped differences to
  [`getBootstrapSamples()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getBootstrapSamples.md).
- Added average MED identification rate as attribute to
  [`assessDesign()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/assessDesign.md)
  output.
- Made the `future.apply` package optional.
- Re-worked vignettes and improved the output of print functions.

## BayesianMCPMod 1.1.0 (07-Mar-2025)

CRAN release: 2025-03-07

- Fixed a bug in
  [`plot.modelFits()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/plot.modelFits.md)
  that would plot credible bands based on incorrectly selected
  bootstrapped quantiles.
- Added
  [`getMED()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getMED.md),
  a function to assess the minimally efficacious dose (MED) and
  integrated
  [`getMED()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getMED.md)
  into
  [`assessDesign()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/assessDesign.md)
  and
  [`performBayesianMCPMod()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/performBayesianMCPMod.md).
- Added parallel processing using the future framework.
- Modified the handling of the fit of an average model: Now,
  [`getModelFits()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getModelFits.md)
  has an argument to fit an average model and this will be carried
  forward for all subsequent functions.
- Re-introduced
  [`getBootstrapSamples()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getBootstrapSamples.md),
  a separate function for bootstrapping samples from the posterior
  distributions of the dose levels.
- Adapted the vignettes to new features.

## BayesianMCPMod 1.0.2 (06-Feb-2025)

CRAN release: 2025-02-06

- Addition of new vignette comparing frequentist and Bayesian MCPMod
  using vague priors.
- Extension of
  [`getPosterior()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getPosterior.md)
  to allow the input of a fully populated variance-covariance matrix.
- Added the non-monotonic model shapes beta and quadratic.
- New argument in
  [`assessDesign()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/assessDesign.md)
  to optionally skip the Mod part of MCPMod.
- Additional tests.

## BayesianMCPMod 1.0.1 (03-Apr-2024)

CRAN release: 2024-04-05

- Re-submission of the `BayesianMCPMod` package.
- Removed a test that occasionally failed on the fedora CRAN test
  system.
- Fixed a bug in
  [`getBootstrapQuantiles()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getBootstrapQuantiles.md)
  that would return wrong bootstrapped quantiles.
- Added
  [`getBootstrapSamples()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getBootstrapSamples.md),
  a separate function for bootstrapping samples.

## BayesianMCPMod 1.0.0 (31-Dec-2023)

CRAN release: 2024-01-08

- Initial release of the `BayesianMCPMod` package.
- Special thanks to Jana Gierse, Bjoern Bornkamp, Chen Yao, Marius
  Thomas & Mitchell Thomann for their review and valuable comments.
- Thanks to Kevin Kunzmann for R infrastructure support and to Frank
  Fleischer for methodological support.
