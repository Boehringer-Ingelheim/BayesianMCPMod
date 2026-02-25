# performBayesianMCP

Performs Bayesian MCP Test step, as described in Fleischer et al.
(2022). Tests for a dose-response effect using a model-based multiple
contrast test based on the (provided) posterior distribution. In
particular for every dose-response candidate the posterior probability
is calculated that the contrast is bigger than 0 (based on the posterior
distribution of the dose groups). In order to obtain significant test
decision we consider the maximum of the posterior probabilities across
the different models. This maximum is compared with a (multiplicity
adjusted) critical value (on the probability scale).

## Usage

``` r
performBayesianMCP(posterior_list, contr, crit_prob_adj)
```

## Arguments

- posterior_list:

  An object derived with getPosterior with information about the
  (mixture) posterior distribution per dose group

- contr:

  An object of class 'optContr' as created by the getContr() function.
  It contains the contrast matrix to be used for the testing step.

- crit_prob_adj:

  A getCritProb object, specifying the critical value to be used for the
  testing (on the probability scale)

## Value

Bayesian MCP test result, with information about p-values for the
individual dose-response shapes and overall significance

## References

Fleischer F, Bossert S, Deng Q, Loley C, Gierse J. 2022. Bayesian
MCPMod. Pharmaceutical Statistics. 21(3): 654-670. doi:10.1002/pst.2193

## Examples

``` r
mods <- DoseFinding::Mods(linear      = NULL,
                          emax        = c(0.5, 1.2),
                          exponential = 2,
                          doses       = c(0, 0.5, 2,4, 8))
dose_levels  <- c(0, 0.5, 2, 4, 8)
sd_posterior <- c(2.8,3,2.5,3.5,4)
contr_mat <- getContr(
  mods          = mods,
  dose_levels   = dose_levels,
  cov_posterior = diag(sd_posterior)^2)
critVal <- getCritProb(
  mods           = mods,
  dose_weights   = c(50, 50, 50, 50, 50), #reflecting the planned sample size
  dose_levels    = dose_levels,
  alpha_crit_val = 0.05)
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

performBayesianMCP(posterior_list = posterior_list,
                   contr          = contr_mat,
                   crit_prob_adj  = critVal)
#> Bayesian Multiple Comparison Procedure
#>   Significant:                   0 
#>   Critical Probability:          0.9751629 
#>   Maximum Posterior Probability: 0.6412959 
#> Posterior Probabilities for Model Shapes
#>                        lin     emax1     emax2       exp
#>   Posterior Prob 0.6191201 0.6412959 0.6410572 0.5839261 
#>   Significant            0         0         0         0 
#> Average Posterior ESS
#>   Dose Level:   Ctrl DG_1 DG_2 DG_3 DG_4 
#>   Avg Post ESS:  0.3  0.3  0.1  0.1  0.1 
```
