# getPosterior

Either the patient level data or both mu_hat as well as S_hat must to be
provided. If patient level data is provided mu_hat and S_hat are
calculated within the function using a linear model. This function
calculates the posterior distribution. Depending on the input for S_hat
this step is either performed for every dose group independently via the
RBesT function postmix() or the mvpostmix() function of the DoseFinding
package is utilized. In the latter case conjugate posterior mixture of
multivariate normals are calculated (DeGroot 1970, Bernardo and Smith
1994)

## Usage

``` r
getPosterior(
  prior_list,
  data = NULL,
  mu_hat = NULL,
  S_hat = NULL,
  calc_ess = FALSE,
  probability_scale = attr(data, "probability_scale")
)
```

## Arguments

- prior_list:

  a prior list with information about the prior to be used for every
  dose group

- data:

  dataframe containing the information of dose and response. Also a
  simulateData object can be provided. Default NULL.

- mu_hat:

  vector of estimated mean values (per dose group). Default NULL.

- S_hat:

  covariance matrix specifying the (estimated) variability. The
  variance-covariance matrix should be provided and the dimension of the
  matrix needs to match the number of dose groups. Default NULL.

- calc_ess:

  boolean variable, indicating whether effective sample size should be
  calculated. Default FALSE.

- probability_scale:

  A boolean to specify if the trial has a continuous or a binary
  outcome. Setting to TRUE will transform calculations from the logit
  scale to the probability scale, which can be desirable for a binary
  outcome. Default `attr(data, "probability_scale")`.

## Value

posterior_list, a posterior list object is returned with information
about (mixture) posterior distribution per dose group (more detailed
information about the conjugate posterior in case of covariance input
for S_hat is provided in the attributes)

## Details

Kindly note that one can sample from the `posterior_list` with
`lapply(posterior_list, RBesT::rmix, n = 10)`.

## References

BERNARDO, Jl. M., and Smith, AFM (1994). Bayesian Theory. 81.

## Examples

``` r
prior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 5), sigma = 2),
                   DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
                   DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2) ,
                   DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2) ,
                   DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2))
                   
mu_hat <- c(0, 1, 1.5, 2, 2.5)
S_hat  <- diag(c(5, 4, 6, 7, 8)^2)

posterior_list <- getPosterior(
   prior_list = prior_list,
   mu_hat     = mu_hat,
   S_hat      = S_hat)

summary(posterior_list)
#>          mean       sd       2.5%    50.0%     97.5%
#> Ctrl 0.000000 3.535534  -6.929519 0.000000  6.929519
#> DG_1 1.000000 3.794733  -6.437540 1.000000  8.437540
#> DG_2 1.431210 5.267373  -8.892652 1.431210 11.755072
#> DG_3 1.798235 5.905630  -9.776588 1.798235 13.373058
#> DG_4 2.362661 6.813267 -10.991096 2.362661 15.716418
```
