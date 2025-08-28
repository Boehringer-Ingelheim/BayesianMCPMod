# BayesianMCPMod CRAN Comments

Maintainer: <stephan.wojciekowski@boehringer-ingelheim.com>

## Test environments & R CMD check results
- Local aarch64-apple-darwin20, R version 4.5.1 (2025-06-13)
  - 0 errors √ | 0 warnings √ | 0 notes √
- Winbuilder x86_64-w64-mingw32, Windows Server, R under development (unstable) (2025-08-27 r88724 ucrt)
  - Status: OK
- Macbuilder aarch64-apple-darwin20, macOS Ventura 13.3.1, R version 4.4.2 (2024-10-31)
  - Status: OK
- GitHub Action Linux: R version 4.5.1 (2025-06-13), x86_64-pc-linux-gnu
  - Status: OK
- GitHub Action Mac: R version 4.5.1 (2025-06-13), aarch64-apple-darwin20
  - Status: OK
- GitHub Action Windows: R version 4.5.1 (2025-06-13 ucrt), x86_64-w64-mingw32
  - Status: 1 NOTE
    - checking for detritus in the temp directory ... NOTE
    - This seems to be a github related note that I could not reproduce elsewhere.

### From NEWS.md: BayesianMCPMod 1.2.0 (28-Aug-2025)

* Fixed a bug in performBayesianMCPMod() where the model significance status from the MCP step was sometimes not correctly assigned to the fitted model in the Mod step
* Fixed a bug in print.modelFit() where sometimes the coefficients for the fitted model shapes were not printed correctly
* Fixed a bug in getMED() where quantile and evidence level could sometimes not be matched due to floating-point precision issues when using bootstrapped quantiles
* Changed functions getPosterior(), getCritProb(), and getContr() to accept a covariance matrix instead of a standard deviation vector as argument
* Added support for none-zero off-diagonal covariance matrices in the MCP step
* Added bootstrapped differences to getBootstrapSamples()
* Added average MED identification rate as attribute to assessDesign() output
* Made the future.apply package optional
* Re-worked vignettes and improved the output of print functions