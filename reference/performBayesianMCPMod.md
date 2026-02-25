# performBayesianMCPMod

Performs Bayesian MCP Test step and modeling in a combined fashion. See
performBayesianMCP() function for MCP Test step and getModelFits() for
the modeling step

## Usage

``` r
performBayesianMCPMod(
  posterior_list,
  contr,
  crit_prob_adj,
  simple = FALSE,
  avg_fit = TRUE,
  delta = NULL,
  evidence_level = NULL,
  med_selection = c("avgFit", "bestFit"),
  n_samples = 1000,
  probability_scale = FALSE
)
```

## Arguments

- posterior_list:

  An object of class 'postList' or a list of 'postList' objects as
  created by getPosterior() containing information about the (mixture)
  posterior distribution per dose group

- contr:

  An object of class 'optContr' as created by the getContr() function.
  It contains the contrast matrix to be used for the testing step.

- crit_prob_adj:

  A getCritProb object, specifying the critical value to be used for the
  testing (on the probability scale).

- simple:

  Boolean variable, defining whether simplified fit will be applied.
  Passed to the getModelFits() function. Default FALSE.

- avg_fit:

  Boolean variable, defining whether an average fit (based on
  generalized AIC weights) should be performed in addition to the
  individual models. Default TRUE.

- delta:

  A numeric value for the threshold Delta for the MED assessment. If
  NULL, no MED assessment is performed. Default NULL.

- evidence_level:

  A numeric value between 0 and 1 for the evidence level gamma for the
  MED assessment. Only required for Bayesian MED assessment, see ?getMED
  for details. Default NULL.

- med_selection:

  A string, either "avgFit" or "bestFit" based on the lowest gAIC, for
  the method of MED selection. Default "avgFit".

- n_samples:

  A numerical for the number of bootstrapped samples in case the
  Bayesian MED assessment is performed. Default 1e3.

- probability_scale:

  A boolean to specify if the trial has a continuous or a binary
  outcome. Setting to TRUE will transform calculations from the logit
  scale to the probability scale, which can be desirable for a binary
  outcome. Default FALSE.

## Value

Bayesian MCP test result as well as modeling result.

## Examples

``` r
mods <- DoseFinding::Mods(linear      = NULL,
                          emax        = c(0.5, 1.2),
                          exponential = 2,
                          doses       = c(0, 0.5, 2,4, 8))
dose_levels  <- c(0, 0.5, 2, 4, 8)
sd_posterior <- c(2.8, 3, 2.5, 3.5, 4)
contr_mat <- getContr(
  mods          = mods,
  dose_levels   = dose_levels,
  cov_posterior = diag(sd_posterior)^2)
critVal <- getCritProb(
  mods           = mods,
  dose_weights   = c(50, 50, 50, 50, 50), #reflecting the planned sample size
  dose_levels    = dose_levels,
  alpha_crit_val = 0.6) # unreasonable alpha chosen for this example, rather choose 0.05
prior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 5), sigma = 2),
                   DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
                   DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2) ,
                   DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2) ,
                   DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2))
mu <- c(0, 1, 1.5, 2, 2.5)
S_hat <- diag(c(5, 4, 6, 7, 8)^2)
posterior_list <- getPosterior(
  prior_list = prior_list,
  mu_hat     = mu,
  S_hat      = S_hat,
  calc_ess   = TRUE)

performBayesianMCPMod(posterior_list = posterior_list,
                      contr          = contr_mat,
                      crit_prob_adj  = critVal,
                      simple         = FALSE,
                      delta          = 1.1)
#> Bayesian Multiple Comparison Procedure
#>   Significant:                   1 
#>   Critical Probability:          0.5706743 
#>   Maximum Posterior Probability: 0.6412959 
#> Posterior Probabilities for Model Shapes
#>                        lin     emax1     emax2       exp
#>   Posterior Prob 0.6191201 0.6412959 0.6410572 0.5839261 
#>   Significant            1         1         1         1 
#> Average Posterior ESS
#>   Dose Level:   Ctrl DG_1 DG_2 DG_3 DG_4 
#>   Avg Post ESS:  0.3  0.3  0.1  0.1  0.1 
#> MED Assessment
#>   Selection Method:    avgFit 
#>   Identification Rate: 1 
#>    Dose Level: 0.5 2.0 4.0 8.0 
#>    MED Freq:     0   0   1   0 
#>   MED not reached Freq:        0 
#>   No success in MCP step Freq: 0 
#> Model Coefficients
#>   emax   e0 = 0, eMax = 2.2, ed50 = 0.8 
#>   exp    e0 = 0.5, e1 = 3.4, delta = 16 
#>   lin    e0 = 0.5, delta = 0.3 
#> Dose Levels
#>   Ctrl = 0, DG_1 = 0.5, DG_2 = 2, DG_3 = 4, DG_4 = 8 
#> Predictions, Maximum Effect, gAIC, avgFit Model Weights & Significance
#>          Ctrl DG_1 DG_2 DG_3 DG_4 mEff gAIC    w Sign
#>   avgFit  0.4  0.7  1.2  1.6  2.6  2.2   NA   NA   NA 
#>   emax    0.0  0.9  1.6  1.9  2.1  2.0  6.0  0.2  1.0 
#>   exp     0.5  0.6  1.0  1.5  2.7  2.2  6.0  0.2  1.0 
#>   lin     0.5  0.6  1.0  1.6  2.7  2.3  4.0  0.6  1.0 
```
