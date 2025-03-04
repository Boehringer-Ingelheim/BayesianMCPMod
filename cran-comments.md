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

### BayesianMCPMod 1.0.2 (06-Feb-2025)

- Addition of new vignette comparing frequentist and Bayesian MCPMod using vague priors
- Extension of getPosterior to allow the input of a fully populated variance-covariance matrix
- Added the non-monotonic model shapes beta and quadratic
- New argument in assessDesign() to skip the Mod part of Bayesian MCPMod
- Additional tests 

### BayesianMCPMod 1.0.1 (03-Apr-2024)

- Re-submission of the 'BayesianMCPMod' package
- Removed a test that occasionally failed on the fedora CRAN test system
- Fixed a bug that would return wrong bootstrapped quantiles in getBootstrappedQuantiles()
- Added getBootstrapSamples(), a separate function for bootstrapping samples

### BayesianMCPMod 1.0.0 (31-Dec-2023)

- Initial release of the 'BayesianMCPMod' package
- Special thanks to Jana Gierse, Bjoern Bornkamp, Chen Yao, Marius Thoma & Mitchell Thomann for their review and valuable comments
- Thanks to Kevin Kunzmann for R infrastructure support and to Frank Fleischer for methodological support
