# BayesianMCPMod CRAN Comments

Maintainer: <stephan.wojciekowski@boehringer-ingelheim.com>

## Test environments & R CMD check results
- Local aarch64-apple-darwin20, R version 4.5.1 (2025-06-13)
  - 0 errors ✔ | 0 warnings ✔ | 0 notes ✔
- Winbuilder x86_64-w64-mingw32,  Under development (unstable) (2026-05-12 r90049 ucrt)
  - Status: OK
- Macbuilder aarch64-apple-darwin23, 4.6.0 Patched (2026-04-24 r89963)
  - Status: OK
- GitHub Action Linux:
  - Status: OK
- GitHub Action Mac:
  - Status: OK
- GitHub Action Windows:
  - Status: OK
  
## BayesianMCPMod 1.3.2 (23-Mar-2026)

* Included Firth's penalized regression model for binary endpoints in case of separation.