# getModelFits

Fits dose-response curves for the specified dose-response models, based
on the posterior distributions. For the simplified fit, multivariate
normal distributions will be approximated and reduced by one-dimensional
normal distributions. For the default case, the Nelder-Mead algorithm is
used. In detail, for both approaches the mean vector \\\theta^{Y}\\ and
the covariance \\\Sigma\\ of the (mixture) posterior distributions and
the corresponding posterior weights \\\tilde{\omega}\_{l}\\ for \\l \in
{1,...,L}\\ are used as basis For the full fit a GLS estimator is used
to minimize the following expression for the respective dose-response
models \\m\\ \$\$ \hat{\theta}\_{m}=\text{arg min}\_{\theta\_{m}}
\sum\_{l=1}^{L}
\tilde{\omega}\_{l}(\theta\_{l\_{i}}^{Y}-f(dose\_{i},\hat{\theta}\_{m}))'\Sigma\_{l}^{-1}(\theta\_{l\_{i}}^{Y}-f(dose\_{i},\hat{\theta}\_{m}))\$\$
Therefore the function nloptr of the nloptr package is utilized. In the
simplified case \\L=1\\, as the dimension of the posterior is reduced to
1 first. The generalized AIC values are calculated via the formula
\$\$gAIC\_{m} = \sum\_{l=1}^{L} \tilde{\omega}\_{l} \sum\_{i=0}^{K}
\frac{1}{\Sigma\_{l\_{i,i}}} (\theta\_{l_i}^Y -
f(dose\_{i},\hat{\theta}\_{m}))^2 + 2p \$\$ where \\p\\ denotes the
number of estimated parameters and \\K\\ the number of active dose
levels. Here as well for the simplified case the formula reduces to one
summand as \\L=1\\. Corresponding gAIC based weights for model \\M\\ are
calculated as outlined in Schorning et al. (2016) \$\$ \Omega_I (M) =
\frac{\exp(-0.5 gAIC\_{M})}{\sum\_{m=1}^{Q} \exp(-0.5 gAIC\_{m})} \$\$
where \\Q\\ denotes the number of models included in the averaging
procedure.

## Usage

``` r
getModelFits(
  models,
  dose_levels,
  posterior,
  avg_fit = TRUE,
  simple = FALSE,
  probability_scale = FALSE
)
```

## Arguments

- models:

  A Mods object as created with
  [`DoseFinding::Mods()`](https://openpharma.github.io/DoseFinding/reference/Mods.html)
  or a vector of model names for which a fit will be performed.
  Implemented model shapes are `"linear"`, `"exponential"`,
  `"logistic"`, `"emax"`, `"sigEmax"`, `"quadratic"`, and `"betaMod"`.

- dose_levels:

  A vector containing the different dosage levels.

- posterior:

  A getPosterior object, containing the (multivariate) posterior
  distribution per dosage level.

- avg_fit:

  Boolean variable, defining whether an average fit (based on
  generalized AIC weights) should be performed in addition to the
  individual models. Default TRUE.

- simple:

  Boolean variable, defining whether simplified fit will be applied.
  Default FALSE.

- probability_scale:

  A boolean variable to specify if the predicted dose-response should be
  on the logit scale or the probability scale. Setting to TRUE will
  transform predictions from the logit scale to the probability scale,
  which can be desirable for a binary outcome. Default FALSE.

## Value

An object of class modelFits. A list containing information about the
fitted model coefficients, the prediction per dose group as well as
maximum effect and generalized AIC (and corresponding weight) per model.

## References

Schorning K, Bornkamp B, Bretz F, Dette H. 2016. Model selection versus
model averaging in dose finding studies. Stat Med; 35; 4021-4040.

## Examples

``` r
posterior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 1), sigma = 2),
                       DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 3, s = 1.2), sigma = 2),
                       DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 4, s = 1.5), sigma = 2),
                       DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 6, s = 1.2), sigma = 2),
                       DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 6.5, s = 1.1), sigma = 2))
models         <- c("emax", "exponential", "sigEmax", "linear")
dose_levels    <- c(0, 1, 2, 4, 8)

fit        <- getModelFits(models      = models,
                           posterior   = posterior_list,
                           dose_levels = dose_levels)
                           
fit
#> Model Coefficients
#>   emax   e0 = 0, eMax = 8.1, ed50 = 1.7 
#>   exp    e0 = 1.6, e1 = 8.8, delta = 16 
#>   lin    e0 = 1.4, delta = 0.7 
#>   sigE   e0 = 0, eMax = 7.5, ed50 = 1.5, h = 1.2 
#> Dose Levels
#>   Ctrl = 0, DG_1 = 1, DG_2 = 2, DG_3 = 4, DG_4 = 8 
#> Predictions, Maximum Effect, gAIC & avgFit Model Weights
#>          Ctrl DG_1 DG_2 DG_3 DG_4 mEff gAIC    w
#>   avgFit  0.2  2.8  4.2  5.5  6.7  6.5   NA   NA 
#>   emax    0.0  3.0  4.3  5.6  6.6  6.6  6.2  0.6 
#>   exp     1.6  2.2  2.8  4.2  7.4  5.7 12.8  0.0 
#>   lin     1.4  2.2  2.9  4.4  7.4  6.0  9.4  0.1 
#>   sigE    0.0  2.9  4.4  5.7  6.6  6.6  8.1  0.2 
                           
fit_simple <- getModelFits(models      = models,
                           posterior   = posterior_list,
                           dose_levels = dose_levels,
                           simple      = TRUE)
                           
fit_simple
#> Model Coefficients
#>   emax   e0 = 0, eMax = 8.1, ed50 = 1.7 
#>   exp    e0 = 1.6, e1 = 8.8, delta = 16 
#>   lin    e0 = 1.4, delta = 0.7 
#>   sigE   e0 = 0, eMax = 7.5, ed50 = 1.5, h = 1.2 
#> Dose Levels
#>   Ctrl = 0, DG_1 = 1, DG_2 = 2, DG_3 = 4, DG_4 = 8 
#> Predictions, Maximum Effect, gAIC & avgFit Model Weights
#>          Ctrl DG_1 DG_2 DG_3 DG_4 mEff gAIC    w
#>   avgFit  0.2  2.8  4.2  5.5  6.7  6.5   NA   NA 
#>   emax    0.0  3.0  4.3  5.6  6.6  6.6  6.2  0.6 
#>   exp     1.6  2.2  2.8  4.2  7.4  5.7 12.8  0.0 
#>   lin     1.4  2.2  2.9  4.4  7.4  6.0  9.4  0.1 
#>   sigE    0.0  2.9  4.4  5.7  6.6  6.6  8.1  0.2 
```
