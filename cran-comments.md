# BayesianMCPMod CRAN Comments

## Maintainer

'Stephan Wojciekowski <stephan.wojciekowski@boehringer-ingelheim.com>'

## Test environments
- Local aarch64-apple-darwin20, R 4.4.2

- Winbuilder x86_64-w64-mingw32, Windows Server, R under development (unstable) (2025-02-05 r87692 ucrt)
- Macbuilder aarch64-apple-darwin20, macOS Ventura 13.3.1, R version 4.4.2 (2024-10-31)

- Github Linux / ubuntu-latest (release)
- Github Mac / macos-latest (release)
- Github windows, all R versions on GitHub Actions windows-latest

## R CMD check results

### Local aarch64-apple-darwin20, Windows Server, R 4.4.2
0 errors √ | 0 warnings √ | 0 notes √

### Winbuilder x86_64-w64-mingw32, macOS Ventura 13.3.1, R under development (unstable) (2025-02-05 r87692 ucrt)

* DONE
Status: OK

### Macbuilder aarch64-apple-darwin20, macOS Ventura 13.3.1, R version 4.4.2 (2024-10-31)

* DONE
Status: OK

### Github Linux / ubuntu-latest (release)

* DONE
Status: OK

### Github Mac / macos-latest (release)

* DONE
Status: OK

### Windows / windows-latest (release)

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'Rscript420e8038' 'Rscripta4ce8028'
* DONE
Status: 1 NOTE

-> This note seems is related with the parallelization on the github server and does not occur on the Winbuilder server.
   
## From NEWS.md

## BayesianMCPMod 1.1.0 (07-Mar-2025)

- Fixed a bug in plot.modelFits() that would plot credible bands based on incorrectly selected bootstrapped quantiles
- Added getMED(), a function to assess the minimally efficacious dose (MED) and integrated getMED() into assessDesign() and performBayesianMCPMod
- Added parallel processing using the future framework
- Modified the handling of the fit of an average model: Now, getModelFits() has an argument to fit an average model and this will be carried forward for all subsequent functions
- Re-introduced getBootstrapSamples(), a separate function for bootstrapping samples from the posterior distributions of the dose levels
- Adapted the vignettes to new features