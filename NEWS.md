## BayesianMCPMod 1.0.2 (xx-Aug-2024)

- addition of new vignette comparing frequentist and bayesian MCPMod (using non-informative prior) 
- extension of posterior function to allow also input of fully populated variance-covariance matrices 
- added non-monotonic models (beta and quadratic)  
- more tests 


## BayesianMCPMod 1.0.1 (03-Apr-2024)

- Re-submission of the 'BayesianMCPMod' package
- Removed a test that occasionally failed on the fedora CRAN test system
- Fixed a bug that would return wrong bootstrapped quantiles in getBootstrapQuantiles()
- Added getBootstrapSamples(), a separate function for bootstrapping samples

## BayesianMCPMod 1.0.0 (31-Dec-2023)

- Initial release of the 'BayesianMCPMod' package
- Special thanks to Jana Gierse, Bjoern Bornkamp, Chen Yao, Marius Thomas & Mitchell Thomann for their review and valuable comments
- Thanks to Kevin Kunzmann for R infrastructure support and to Frank Fleischer for methodological support