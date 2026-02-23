# predict.modelFits

This function performs model predictions based on the provided model and
dose specifications

## Usage

``` r
# S3 method for class 'modelFits'
predict(
  object,
  doses = NULL,
  probability_scale = attr(object, "probability_scale"),
  ...
)
```

## Arguments

- object:

  A modelFits object containing information about the fitted model
  coefficients

- doses:

  A vector specifying the doses for which a prediction should be done

- probability_scale:

  A boolean variable to specify if the trial has a continuous or a
  binary outcome. Setting to TRUE will transform predictions from the
  logit scale to the probability scale, which can be desirable for a
  binary outcome. Default FALSE.

- ...:

  Currently without function

## Value

a list with the model predictions for the specified models and doses

## Examples

``` r
posterior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 1), sigma = 2),
                       DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 3, s = 1.2), sigma = 2),
                       DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 4, s = 1.5), sigma = 2) ,
                       DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 6, s = 1.2), sigma = 2) ,
                       DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 6.5, s = 1.1), sigma = 2))
models         <- c("emax", "exponential", "sigEmax", "linear", "betaMod")
dose_levels    <- c(0, 1, 2, 4, 8)
fit            <- getModelFits(models      = models,
                               posterior   = posterior_list,
                               dose_levels = dose_levels)

predict(fit, doses = c(0, 1, 3, 4, 6, 8))
#> $avgFit
#> [1] 0.1666128 2.8524723 4.9940076 5.5700609 6.3209477 6.6915166
#> 
#> $betaMod
#> [1] 0.006750215 2.896634304 5.209887319 5.930220008 6.745178439 6.507408218
#> 
#> $emax
#> [1] -0.008147769  2.977266571  5.136083680  5.647278695  6.271281086
#> [6]  6.637935285
#> 
#> $exponential
#> [1] 1.636554 2.207315 3.461645 4.150114 5.663125 7.377590
#> 
#> $linear
#> [1] 1.419597 2.167297 3.662695 4.410394 5.905792 7.401191
#> 
#> $sigEmax
#> [1] 0.00731386 2.89049045 5.22947029 5.72787171 6.28978958 6.59200715
#> 
#> attr(,"doses")
#> [1] 0 1 3 4 6 8
```
