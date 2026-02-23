# BayesianMCPMod CRAN Comments

Maintainer: <stephan.wojciekowski@boehringer-ingelheim.com>

## Test environments & R CMD check results
- Local aarch64-apple-darwin20, R version 4.5.1 (2025-06-13)
  - 0 errors √ | 0 warnings √ | 0 notes √
- Local x86_64-pc-linux-gnu, R version 4.5.0 (2025-04-11)
  - 0 errors √ | 0 warnings √ | 0 notes √
- Winbuilder x86_64-w64-mingw32, Windows Server, R under development (unstable) (2025-08-27 r88724 ucrt)
  - Status: OK
- Macbuilder aarch64-apple-darwin20, macOS Ventura 13.3.1, R version 4.4.2 (2024-10-31)
  - Status: OK
- GitHub Action Linux: using R version 4.5.2 (2025-10-31), x86_64-pc-linux-gnu
  - Status: OK
- GitHub Action Mac: R version 4.5.2 (2025-10-31), aarch64-apple-darwin20
  - Status: OK
- GitHub Action Windows: R version 4.5.2 (2025-10-31 ucrt), x86_64-w64-mingw32
  - Status: 1 NOTE
    - checking for detritus in the temp directory ... NOTE
    - Found the following files/directories: 'Rscript1230a75fa' 'Rscripta0ca75fa'
    - This seems to be a Windows github action related NOTE that I could not reproduce elsewhere.

## BayesianMCPMod 1.3.0 (XX-Feb-2026)

* Fixed a bug that would occur when predicting from the beta model shape outside of the original dose range.
* Fixed a bug in which the MED assessment could not be performed when specifying a negative direction of beneficial effect and an evidence level other than 0.5.
* Added functions and vignettes for the binary endpoint case.
* Added functionality to `assessDesign()` to provide custom simulated data and custom model estimates enabling complex data simulation and analysis methods.
* Added argument to `assessDesign()` for number of bootstrap samples in case `evidence_level` is provided.
* Added functionality to `plot.modelFits()` to plot effect sizes.
* Added calls to `set.seed()` in vignette's code blocks to facilitate individual code block reproducibility.