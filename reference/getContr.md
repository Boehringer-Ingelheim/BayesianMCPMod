# getContr

This function calculates contrast vectors that are optimal for detecting
certain alternatives via applying the function optContr() of the
DoseFinding package. Hereby, 4 different options can be distinguished
that are automatically executed based on the input that is provided

1.  Bayesian approach: If dose_weights and a prior_list are provided an
    optimized contrasts for the posterior sample size is calculated. In
    detail, in a first step the dose_weights (typically the number of
    patients per dose group) and the prior information is combined by
    calculating for each dose group a posterior effective sample. Based
    on this posterior effective sample sizes the allocation ratio is
    derived, which allows for a calculation on pseudo-optimal contrasts
    via regular MCPMod are calculated from the regular MCPMod for these
    specific weights

2.  Frequentist approach: If only dose_weights are provided optimal
    contrast vectors are calculated from the regular MCPMod for these
    specific weights

3.  Bayesian approach + re-estimation: If only a cov_posterior (i.e.
    variability of the posterior distribution) is provided,
    pseudo-optimal contrasts based on these posterior weights will be
    calculated

4.  Frequentist approach+re-estimation: If only a cov_new_trial (i.e.
    the estimated variability of a new trial) is provided, optimal
    contrast vectors are calculated from the regular MCPMod for this
    specific covariance matrix.

## Usage

``` r
getContr(
  mods,
  dose_levels,
  dose_weights = NULL,
  prior_list = NULL,
  cov_posterior = NULL,
  cov_new_trial = NULL
)
```

## Arguments

- mods:

  An object of class 'Mods' as created by the function
  'DoseFinding::Mods()'

- dose_levels:

  Vector containing the different dosage levels.

- dose_weights:

  Vector specifying weights for the different doses. Please note that in
  case this information is provided together with a prior (i.e.
  Option 1) is planned these two inputs should be provided on the same
  scale (e.g. patient numbers). Default NULL

- prior_list:

  A list of objects of class 'normMix' as created with
  'RBesT::mixnorm()'. Only required as input for Option 1. Default NULL

- cov_posterior:

  A covariance matrix with information about the variability of the
  posterior distribution, only required for Option 3. Default NULL

- cov_new_trial:

  A covariance matrix with information about the observed variability,
  only required for Option 4. Default NULL

## Value

An object of class 'optContr' as provided by the function
'DoseFinding::optContr()'.

## Examples

``` r
dose_levels  <- c(0, 0.5, 2, 4, 8)
mods <- DoseFinding::Mods(
  linear      = NULL,
  emax        = c(0.5, 1.2),
  exponential = 2,
  doses       = dose_levels,
  maxEff      = 6)
cov_posterior <- diag(c(2.8, 3, 2.5, 3.5, 4)^2)

contr_mat <- getContr(
  mods         = mods,
  dose_levels  = dose_levels,
  cov_posterior = cov_posterior)
```
