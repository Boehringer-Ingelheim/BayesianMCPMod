# simulateData

Function to simulate patient level data for a normally distributed
endpoint

## Usage

``` r
simulateData(
  n_patients,
  dose_levels,
  sd = NULL,
  mods = NULL,
  n_sim = 1000,
  true_model = NULL,
  dr_means = NULL,
  probability_scale = FALSE
)
```

## Arguments

- n_patients:

  Vector containing number of patients as a numerical value per
  dose-group.

- dose_levels:

  Vector containing the different dosage levels.

- sd:

  Standard deviation on patient level. Can be NULL if
  `probability_scale` is TRUE. Default NULL.

- mods:

  An object of class "Mods" as specified in the DoseFinding package. Can
  be NULL if ´dr_means´ is not NULL. Default NULL.

- n_sim:

  Number of simulations to be performed, Default is 1000

- true_model:

  A character for model name, e.g. "emax". Assumed true underlying
  model. If NULL, all dose-response models included in the mods input
  parameter will be used. Default NULL.

- dr_means:

  an optional vector, with information about assumed effects per dose
  group. Default NULL.

- probability_scale:

  A boolean to specify if the trial has a continuous or a binary
  outcome. Setting to TRUE will transform calculations from the logit
  scale to the probability scale, which can be desirable for a binary
  outcome. Default FALSE.

## Value

A list object, containing patient level simulated data for all assumed
true models. Also providing information about simulation iteration,
patient number as well as dosage levels.

## Examples

``` r
models <- DoseFinding::Mods(linear      = NULL,
                            linlog      = NULL,
                            emax        = c(0.5, 1.2),
                            exponential = 2, 
                            doses       = c(0, 0.5, 2,4, 8),
                            maxEff      = 6)
dose_levels <- c(0, 0.5, 2, 4, 8)
sd          <- 12
n_patients  <- c(40, 60, 60, 60, 60)

sim_data <- simulateData(n_patients  = n_patients,
                         dose_levels = dose_levels,
                         sd          = sd,
                         mods        = models)

head(sim_data)
#>   simulation ptno dose     linear     linlog      emax1      emax2 exponential
#> 1          1    1    0  13.081588  13.081588  13.081588  13.081588   13.081588
#> 2          1    2    0   4.215108   4.215108   4.215108   4.215108    4.215108
#> 3          1    3    0 -13.124812 -13.124812 -13.124812 -13.124812  -13.124812
#> 4          1    4    0  -7.065638  -7.065638  -7.065638  -7.065638   -7.065638
#> 5          1    5    0  -6.195623  -6.195623  -6.195623  -6.195623   -6.195623
#> 6          1    6    0  17.412528  17.412528  17.412528  17.412528   17.412528

# custom response "model" shape
custom_dose_response <- c(1, 2, 3, 4, 5)
sim_data_custom_dr   <- simulateData(n_patients  = n_patients,
                                     dose_levels = dose_levels,
                                     sd          = sd,
                                     dr_means    = custom_dose_response)

head(sim_data_custom_dr)
#>   simulation ptno dose dr_response
#> 1          1    1    0   20.716093
#> 2          1    2    0    4.075376
#> 3          1    3    0    1.454067
#> 4          1    4    0   13.761368
#> 5          1    5    0   -2.650874
#> 6          1    6    0   10.975242
```
