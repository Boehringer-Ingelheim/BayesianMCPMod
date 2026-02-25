# BayesianMCPMod CRAN Comments

Maintainer: <stephan.wojciekowski@boehringer-ingelheim.com>

## Test environments & R CMD check results
- Local aarch64-apple-darwin20, R version 4.5.1 (2025-06-13)
  - 1 NOTE:
    - Days since last update: 2
      - This is due to the hotfix update, please see below
- Local x86_64-pc-linux-gnu, R version 4.5.0 (2025-04-11)
  - 1 NOTE:
    - Days since last update: 2
      - This is due to the hotfix update, please see below
- Winbuilder x86_64-w64-mingw32, Windows Server 2022 x64 (build 20348), R Under development (unstable) (2026-02-24 r89478 ucrt)
  - 2 NOTEs:
    - Days since last update: 2
      - This is due to the hotfix update, please see below
    - Examples with CPU (user + system) or elapsed time > 10s
                   user system elapsed
      assessDesign 9.66.  0.49   10.14
      - This was not reproducible in other environments, assessDesign examples should run in < 1s
- Macbuilder aarch64-apple-darwin23, macOS 26.2 (25C56), R Under development (unstable) (2026-02-14 r89420)
  - Status: OK
- GitHub Action Linux: R version 4.5.2 (2025-10-31), x86_64-pc-linux-gnu
  - Status: OK
- GitHub Action Mac: R version 4.5.2 (2025-10-31), aarch64-apple-darwin20
  - Status: OK
- GitHub Action Windows: R version 4.5.2 (2025-10-31 ucrt), x86_64-w64-mingw32
  - Status: 1 NOTE
    - checking for detritus in the temp directory ... NOTE
    - Found the following files/directories: 'Rscript1230a75fa' 'Rscripta0ca75fa'
    - This seems to be a Windows github action related NOTE that I could not reproduce elsewhere.

## BayesianMCPMod 1.3.1 (25-Feb-2026)

* Fixed a newly introduced bug that would occur if the R package `future.apply` was not installed. 
* Added flexibility to bootstrapped credible bands in `plot.modelFits()`.