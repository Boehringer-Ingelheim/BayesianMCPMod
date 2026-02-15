# getMED

This function provides information on the minimally efficacious dose
(MED). The MED evaluation can either be based on the fitted model shapes
(model_fits) or on bootstrapped quantiles (bs_quantiles).

## Usage

``` r
getMED(
  delta,
  evidence_level = 0.5,
  dose_levels = NULL,
  model_fits = NULL,
  bs_quantiles = NULL,
  probability_scale = attr(model_fits, "probability_scale")
)
```

## Arguments

- delta:

  A numeric value for the threshold Delta.

- evidence_level:

  A numeric value between 0 and 1 for the evidence level gamma. Used for
  the bs_quantiles-based evaluation and not used for the
  model_fits-based evaluation. Default 0.5.

- dose_levels:

  A vector of numerics containing the different dosage levels. Default
  NULL.

- model_fits:

  An object of class modelFits as created with getModelFits(). Default
  NULL.

- bs_quantiles:

  A dataframe created with getBootstrapQuantiles(). Default NULL.

- probability_scale:

  A boolean variable to specify if the trial has a continuous or a
  binary outcome. Setting to TRUE will transform predictions from the
  logit scale to the probability scale, which can be desirable for a
  binary outcome. Default `attr(model_fits, "probability_scale")`.

## Value

A matrix with rows for MED reached, MED, and MED index in the vector of
dose levels and columns for the dose-response shapes.

## Details

The function assumes that the 1st dose group is the control dose group.

The bootstrap approach allows for an MED based on decision rules of the
form \$\$\widehat{\text{MED}} = \text{arg min}\_{d\in\\d_1, \dots,
d_k\\} \left\\ \text{Pr}\left(f(d, \hat\theta) - f(d_1, \hat\theta) \>
\Delta\right) \> \gamma \right\\ .\$\$ The model-shape approach takes
the point estimate of the model into account.

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

# MED based on the model_fit:
getMED(delta = 5, model_fits = model_fits)
#>             avgFit exponential linear
#> med_reached      1           1      1
#> med              8           8      8
                               
# MED based on bootstrapped quantiles
bs_quantiles <- getBootstrapQuantiles(model_fits = model_fits,
                                      quantiles  = c(0.025, 0.2, 0.5, 0.8),
                                      n_samples  = 100) # speeding up example run time
                                      
getMED(delta          = 5,
       evidence_level = 0.8,
       bs_quantiles   = bs_quantiles)
#>             avgFit exponential linear
#> med_reached      0           0      0
#> med             NA          NA     NA
       
# MED on the probability scale
getMED(delta = 0.1, model_fits = model_fits, probability_scale = TRUE)
#>             avgFit exponential linear
#> med_reached      1           1      1
#> med              2           2      2
```
