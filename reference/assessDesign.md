# assessDesign

This function performs simulation based trial design evaluations for a
set of specified dose-response models

## Usage

``` r
assessDesign(
  n_patients,
  mods,
  prior_list,
  sd = NULL,
  contr = NULL,
  dr_means = NULL,
  data_sim = NULL,
  estimates_sim = NULL,
  n_sim = 1000,
  alpha_crit_val = 0.05,
  modeling = FALSE,
  simple = TRUE,
  avg_fit = TRUE,
  reestimate = FALSE,
  delta = NULL,
  evidence_level = NULL,
  n_bs_samples = 1000,
  med_selection = c("avgFit", "bestFit"),
  probability_scale = FALSE
)
```

## Arguments

- n_patients:

  Vector specifying the planned number of patients per dose group. A
  minimum of 2 patients are required in each group.

- mods:

  An object of class `Mods` as specified in the `DoseFinding` package.

- prior_list:

  A prior_list object specifying the utilized prior for the different
  dose groups

- sd:

  A positive value, specification of assumed sd. Not required if either
  `data_sim` or `estimates_sim` is provided. Also not required in case
  of binary endpoint. Default NULL

- contr:

  An object of class `optContr` as created by the
  `DoseFinding::getContr` function. Allows specification of a fixed
  contrasts matrix. Default NULL.

- dr_means:

  A vector, allows specification of individual (not model based) assumed
  effects per dose group. Default NULL.

- data_sim:

  An optional data frame for custom simulated data. Must follow the data
  structure as provided by
  [`simulateData()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/simulateData.md).
  Default NULL.

- estimates_sim:

  An optional named list of 1) list of vectors for the estimated means
  per dose group (`estimates_sim$mu_hats`) and 2) of list of matrices
  for the covariance matrices specifying the (estimated) variabilities
  (`estimates_sim$S_hats`). Dimensions of entries must match the number
  of dose levels. Default NULL.

- n_sim:

  Number of simulations to be performed

- alpha_crit_val:

  (Un-adjusted) Critical value to be used for the MCP testing step.
  Passed to the `getCritProb` function for the calculation of adjusted
  critical values (on the probability scale). Default 0.05.

- modeling:

  Boolean variable defining whether the Mod part of Bayesian MCP-Mod
  will be performed in the assessment. More heavy on resources. Default
  FALSE.

- simple:

  Boolean variable defining whether simplified fit will be applied, see
  [`?getModelFits`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getModelFits.md).
  Set automatically to TRUE if argument `delta` is provided. Passed to
  `getModelFits` Default TRUE.

- avg_fit:

  Boolean variable, defining whether an average fit (based on
  generalized AIC weights) should be performed in addition to the
  individual models. Default TRUE.

- reestimate:

  Boolean variable defining whether critical value should be calculated
  with re-estimated contrasts (see `getCritProb` function for more
  details). Default FALSE.

- delta:

  A numeric value for the threshold Delta for the MED assessment. If
  NULL, no MED assessment is performed. Default NULL.

- evidence_level:

  A numeric value between 0 and 1 for the evidence level gamma for the
  MED assessment. Only required for Bayesian MED assessment, see
  [`?getMED`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/getMED.md)
  for details. If NULL, MED assessment will be performed on the fitted
  model according to the argument `med_selection`. Default NULL.

- n_bs_samples:

  Number of bootstrap samples for the MED assessment if `evidence_level`
  is provided. Note that more extreme quantiles (i.e., quantiles closer
  to 0 or 1) tend to require more bootstrap samples to maintain
  precision. Default 1000.

- med_selection:

  A string, either `"avgFit"` or `"bestFit"`, for the method of MED
  selection. Default `"avgFit"`.

- probability_scale:

  A boolean to specify if the trial has a continuous or a binary
  outcome. Setting to TRUE will transform calculations from the logit
  scale to the probability scale, which can be desirable for a binary
  outcome. Default FALSE.

## Value

Returns success probabilities for the different assumed dose-response
shapes, attributes also includes information around average success rate
(across all assumed models) and prior Effective sample size.

## Examples

``` r
mods <- DoseFinding::Mods(linear      = NULL,
                          emax        = c(0.5, 1.2),
                          exponential = 2,
                          betaMod     = c(1, 1),
                          doses       = c(0, 0.5, 2,4, 8),
                          maxEff      = 6)
                          
sd <- 12
prior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 12), sigma = 2),
                   DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
                   DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2) ,
                   DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2) ,
                   DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2))
n_patients <- c(40, 60, 60, 60, 60)
dose_levels  <- c(0, 0.5, 2, 4, 8)

success_probabilities <- assessDesign(
  n_patients  = n_patients,
  mods        = mods,
  prior_list  = prior_list,
  sd          = sd,
  n_sim       = 1e2) # speed up example run time

success_probabilities
#> $linear
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.86 
#>   N Simulations:          100 
#>    Model Shape:         lin emax1 emax2   exp betaM 
#>    Significance Freq:  0.81  0.60  0.75  0.74  0.31 
#> 
#> $emax1
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.88 
#>   N Simulations:          100 
#>    Model Shape:         lin emax1 emax2   exp betaM 
#>    Significance Freq:  0.57  0.81  0.78  0.24  0.74 
#> 
#> $emax2
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.9 
#>   N Simulations:          100 
#>    Model Shape:         lin emax1 emax2   exp betaM 
#>    Significance Freq:  0.70  0.81  0.83  0.43  0.79 
#> 
#> $exponential
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.85 
#>   N Simulations:          100 
#>    Model Shape:         lin emax1 emax2   exp betaM 
#>    Significance Freq:  0.81  0.33  0.54  0.83  0.04 
#> 
#> $betaMod
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.86 
#>   N Simulations:          100 
#>    Model Shape:         lin emax1 emax2   exp betaM 
#>    Significance Freq:  0.28  0.61  0.67  0.04  0.85 
#> 
#> attr(,"avgSuccessRate")
#> [1] 0.87
#> attr(,"placEff")
#> [1] 0
#> attr(,"maxEff")
#> [1] 6
#> attr(,"sampleSize")
#> [1] 40 60 60 60 60
#> attr(,"priorESS")
#> Ctrl DG_1 DG_2 DG_3 DG_4 
#>    0    0    0    0    0 

if (interactive()) { # showcasing further functionality

## Analysis with custom data
data_sim <- simulateData(
  n_patients        = n_patients,
  dose_levels       = dose_levels,
  sd                = sd,
  mods              = mods,
  n_sim             = 10)

success_probabilities_cd <- assessDesign(
  n_patients  = n_patients,
  mods        = mods,
  prior_list  = prior_list,
  data_sim    = data_sim,
  sd          = sd,
  n_sim       = 1e2) # speed up example run time

success_probabilities_cd

## Analysis with custom dose response relationship
custom_dr_means <- c(1, 2, 3, 4, 5)

success_probs_custom_dr <- assessDesign(
  n_patients  = n_patients,
  mods        = mods,
  prior_list  = prior_list,
  dr_means    = custom_dr_means,
  sd          = sd,
  n_sim       = 1e2) # speed up example run time

success_probs_custom_dr

## Analysis with custom estimates for means and variabilies
## No simulated data, only simulated model estimates
estimates_sim <- list(mu_hats = replicate(100, list(c(1, 2, 3, 4, 5) + rnorm(5, 0, 1))),
                      S_hats  = list(diag(1, 5)))

success_probs_custom_est <- assessDesign(
  n_patients    = n_patients,
  mods          = mods,
  prior_list    = prior_list,
  estimates_sim = estimates_sim)

success_probs_custom_est

}

if (interactive()) { # takes typically > 5 seconds

# with MED estimation without bootstrapping
# see ?getMED for details

success_probabilities <- assessDesign(
  n_patients     = n_patients,
  mods           = mods,
  prior_list     = prior_list,
  sd             = sd,
  modeling       = TRUE,
  n_sim          = 10, # speed up example run time
  delta          = 7)

  success_probabilities

# with MED estimation with bootstrapping

success_probabilities <- assessDesign(
  n_patients     = n_patients,
  mods           = mods,
  prior_list     = prior_list,
  sd             = sd,
  modeling       = TRUE,
  n_sim          = 10, # speed up example run time
  delta          = 7,
  evidence_level = 0.8)

  success_probabilities

}
```
