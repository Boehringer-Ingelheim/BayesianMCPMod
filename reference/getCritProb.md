# getCritProb

This function calculates multiplicity adjusted critical values. The
critical values are calculated in such a way that when using
non-informative priors the actual error level for falsely declaring a
significant trial in the Bayesian MCPMod is controlled (by the specified
alpha level). Hereby optimal contrasts of the frequentist MCPMod are
applied and two options can be distinguished

1.  Frequentist approach: If only dose_weights are provided optimal
    contrast vectors are calculated from the regular MCPMod for these
    specific weights and the corresponding critical value for this set
    of contrasts is calculated via the critVal() function of the
    DoseFinding package.

2.  Frequentist approach + re-estimation: If only a cov_new_trial (i.e.
    the covariance matrix of a new trial) is provided, optimal contrast
    vectors are calculated from the regular MCPMod for this specific
    matrix. Here as well the critical value for this set of contrasts is
    calculated via the critVal() function of the DoseFinding package.

## Usage

``` r
getCritProb(
  mods,
  dose_levels,
  dose_weights = NULL,
  cov_new_trial = NULL,
  alpha_crit_val = 0.025
)
```

## Arguments

- mods:

  An object of class "Mods" as specified in the DoseFinding package.

- dose_levels:

  Vector containing the different dosage levels.

- dose_weights:

  Vector specifying weights for the different doses, only required for
  Option i). Default NULL

- cov_new_trial:

  A covariance matrix, only required for Option ii). Default NULL

- alpha_crit_val:

  Significance level. Default set to 0.025.

## Value

Multiplicity adjusted critical value on the probability scale.

## Examples

``` r
mods <- DoseFinding::Mods(linear      = NULL,
                          emax        = c(0.5, 1.2),
                          exponential = 2,
                          doses       = c(0, 0.5, 2,4, 8))
dose_levels <- c(0, 0.5, 2, 4, 8)
critVal <- getCritProb(
  mods           = mods,
  dose_weights   = c(50,50,50,50,50), #reflecting the planned sample size
  dose_levels    = dose_levels,
  alpha_crit_val = 0.05)
```
