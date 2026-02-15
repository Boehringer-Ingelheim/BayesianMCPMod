# getBootstrapSamples

A function to return bootstrap samples from the fitted dose-response
models. Samples from the posterior distribution are drawn (via the RBesT
function rmix()) and for every sample the simplified fitting step (see
getModelFits() function) and a prediction is performed. These samples
are returned by this function. This approach can be considered as the
Bayesian equivalent of the frequentist bootstrap approach described in
O'Quigley et al. (2017). Instead of drawing n bootstrap samples from the
sampling distribution of the trial dose-response estimates, here the
samples are directly taken from the posterior distribution.

## Usage

``` r
getBootstrapSamples(
  model_fits,
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

- n_samples:

  Number of samples that should be drawn

- doses:

  A vector of doses for which a prediction should be performed

- probability_scale:

  A boolean variable to specify if the trial has a continuous or a
  binary outcome. Setting to TRUE will transform predictions from the
  logit scale to the probability scale, which can be desirable for a
  binary outcome. Default `attr(model_fits, "probability_scale")`.

## Value

A tibble with columns for sample_id, model, dose, sample, and
sample_diff

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
                           
bs_samples <- getBootstrapSamples(model_fits = model_fits,
                                  n_samples  = 10, # speeding up example run time
                                  doses      = c(0, 6, 8))
                      
bs_samples
#> # A tibble: 90 × 5
#>    sample_id model        dose   abs  diff
#>        <int> <chr>       <dbl> <dbl> <dbl>
#>  1         1 avgFit          0  1.07  0   
#>  2         1 exponential     0  1.36  0   
#>  3         1 linear          0  1.07  0   
#>  4         1 avgFit          6  6.82  5.75
#>  5         1 exponential     6  6.50  5.14
#>  6         1 linear          6  6.82  5.75
#>  7         1 avgFit          8  8.73  7.66
#>  8         1 exponential     8  8.69  7.33
#>  9         1 linear          8  8.73  7.66
#> 10         2 avgFit          0  2.01  0   
#> # ℹ 80 more rows
```
