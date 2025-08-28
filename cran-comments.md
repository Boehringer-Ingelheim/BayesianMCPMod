# BayesianMCPMod CRAN Comments

## Maintainer

'Stephan Wojciekowski <stephan.wojciekowski@boehringer-ingelheim.com>'

## Test environments
- Local aarch64-apple-darwin20, R 4.4.2

- Winbuilder x86_64-w64-mingw32, Windows Server, R under development (unstable) (2025-02-05 r87692 ucrt)
- Macbuilder aarch64-apple-darwin20, macOS Ventura 13.3.1, R version 4.4.2 (2024-10-31)

- Rhub Linux
- Rhub Mac
- Rhub windows

## R CMD check results

### Local aarch64-apple-darwin20, R 4.4.2
0 errors √ | 0 warnings √ | 0 notes √

### Winbuilder x86_64-w64-mingw32, macOS Ventura 13.3.1, R under development (unstable) (2025-02-05 r87692 ucrt)

* DONE
Status: OK

### Macbuilder aarch64-apple-darwin20, macOS Ventura 13.3.1, R version 4.4.2 (2024-10-31)

* DONE
Status: OK

### Rhub Linux / ubuntu-next

* DONE
Status: OK

### Rhub Mac / macos-arm64

* DONE
Status: OK

### Rhub Windows / windows

* DONE
Status: OK
   
## From NEWS.md

### BayesianMCPMod 1.2.0 (28-Aug-2025)

* Fixed a bug in performBayesianMCPMod() where the model significance status from the MCP step was sometimes not correctly assigned to the fitted model in the Mod step
* Fixed a bug in print.modelFit() where sometimes the coefficients for the fitted model shapes were not printed correctly
* Fixed a bug in getMED() where quantile and evidence level could sometimes not be matched due to floating-point precision issues when using bootstrapped quantiles
* Changed functions getPosterior(), getCritProb(), and getContr() to accept a covariance matrix instead of a standard deviation vector as argument
* Added support for none-zero off-diagonal covariance matrices in the MCP step
* Added bootstrapped differences to getBootstrapSamples()
* Added average MED identification rate as attribute to assessDesign() output
* Made the future.apply package optional
* Re-worked vignettes and improved the output of print functions