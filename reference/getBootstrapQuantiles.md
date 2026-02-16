# getBootstrapQuantiles

A function for the calculation of bootstrapped model predictions.
Samples from the posterior distribution are drawn (via the RBesT
function rmix()) and for every sample the simplified fitting step (see
getModelFits() function) and a prediction is performed. These fits are
then used to identify the specified quantiles. This approach can be
considered as the Bayesian equivalent of the frequentist bootstrap
approach described in O'Quigley et al. (2017). Instead of drawing n
bootstrap samples from the sampling distribution of the trial
dose-response estimates, here the samples are directly taken from the
posterior distribution.

## Usage

``` r
getBootstrapQuantiles(
  model_fits,
  quantiles,
  n_samples = 1000,
  doses = NULL,
  probability_scale = attr(model_fits, "probability_scale")
)
```

## Arguments

- model_fits:

  An object of class modelFits, i.e. information about fitted models &
  corresponding model coefficients as well as the posterior distribution
  that was the basis for the model fitting

- quantiles:

  A vector of quantiles that should be evaluated

- n_samples:

  Number of samples that should be drawn as basis for the bootstrapped
  quantiles

- doses:

  A vector of doses for which a prediction should be performed. If NULL,
  the dose levels of the model_fits will be used. Default NULL.

- probability_scale:

  A boolean variable to specify if the trial has a continuous or a
  binary outcome. Setting to TRUE will transform predictions from the
  logit scale to the probability scale, which can be desirable for a
  binary outcome. Default `attr(model_fits, "probability_scale")`.

## Value

A tibble with columns for model, dose, and bootstrapped samples

## References

O'Quigley J, Iasonos A, Bornkamp B. 2017. Handbook of Methods for
Designing, Monitoring, and Analyzing Dose-Finding Trials (1st ed.).
Chapman and Hall/CRC. doi:10.1201/9781315151984

## Examples

``` r
posterior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 1), sigma = 2),
                       DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 3, s = 1.2), sigma = 2),
                       DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 4, s = 1.5), sigma = 2) ,
                       DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 6, s = 1.2), sigma = 2) ,
                       DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 6.5, s = 1.1), sigma = 2))
models         <- c("exponential", "linear")
dose_levels    <- c(0, 1, 2, 4, 8)
model_fits     <- getModelFits(models      = models,
                               posterior   = posterior_list,
                               dose_levels = dose_levels,
                               simple      = TRUE)

bs_quantiles <- getBootstrapQuantiles(model_fits = model_fits,
                                      quantiles  = c(0.025, 0.5, 0.8, 0.975),
                                      n_samples  = 10, # speeding up example run time
                                      doses      = c(0, 6, 8))
                      
bs_quantiles
#> # A tibble: 72 × 5
#>    model   dose sample_type q_prob   q_val
#>    <chr>  <dbl> <chr>        <dbl>   <dbl>
#>  1 avgFit     0 abs          0.025 -0.0692
#>  2 avgFit     0 abs          0.5    1.68  
#>  3 avgFit     0 abs          0.8    2.18  
#>  4 avgFit     0 abs          0.975  3.73  
#>  5 avgFit     0 diff         0.025  0     
#>  6 avgFit     0 diff         0.5    0     
#>  7 avgFit     0 diff         0.8    0     
#>  8 avgFit     0 diff         0.975  0     
#>  9 avgFit     6 abs          0.025  5.03  
#> 10 avgFit     6 abs          0.5    5.82  
#> # ℹ 62 more rows
```
