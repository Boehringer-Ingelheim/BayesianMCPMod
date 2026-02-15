# Simulation Example of Bayesian MCPMod for Continuous Data

``` r
library(BayesianMCPMod)
library(clinDR)
library(dplyr)

set.seed(7015)
```

## Background and Data

In this vignette, we will show the use of the `BayesianMCPMod` R package
for trial planning for continuously distributed data. As in the
[analysis example
vignette](https://boehringer-ingelheim.github.io/BayesianMCPMod/articles/analysis_normal.md),
we focus on the indication MDD and make use of historical data that is
included in the clinDR package. More specifically, trial results for
BRINTELLIX will be utilized to establish an informative prior for the
control group.

A more general overview of the R package was provided with a
[poster](https://github.com/Boehringer-Ingelheim/BayesianMCPMod/blob/36763ee5325955ca5d76a6c140ea3881c74e5fda/inst/PSI_Poster_2024.pdf)
presented during the PSI 2024 Conference.

This package makes use of the
[future](https://cran.r-project.org/package=future) framework for
parallel processing, which can be set up as follows:

``` r
future::plan(future::multisession)
```

## Calculation of a MAP Prior

In a first step, a meta analytic predictive prior will be calculated
using historical data from 5 trials with main endpoint Change from
baseline in MADRS score after 8 weeks. Please note that only information
from the control group will be integrated leading to an informative
mixture prior for the control group, while for the active groups a
non-informative prior will be specified.

``` r
data("metaData")
testdata    <- as.data.frame(metaData)
dataset     <- filter(testdata, bname == "BRINTELLIX")
histcontrol <- filter(dataset, dose == 0, primtime == 8, indication == "MAJOR DEPRESSIVE DISORDER")

hist_data <- data.frame(
  trial = histcontrol$nctno,
  est   = histcontrol$rslt,
  se    = histcontrol$se,
  sd    = histcontrol$sd,
  n     = histcontrol$sampsize)

sd_tot <- with(hist_data, sum(sd * n) / sum(n))
```

We will make use of the same `getPriorList()` function as in the
[analysis example
vignette](https://boehringer-ingelheim.github.io/BayesianMCPMod/articles/analysis_normal.md)
to create a MAP prior.

``` r
dose_levels <- c(0, 2.5, 5, 10, 20)

prior_list  <- getPriorList(
  hist_data     = hist_data,
  dose_levels   = dose_levels,
  robust_weight = 0.3)
```

Kindly note that a vague prior could be implemented via

``` r
prior_list_vague <- rep(list(RBesT::mixnorm(comp1 = c(w = 1, m = 0, n = 1),
                                            sigma = sd_sim, param = "mn")),
                        times = length(dose_levels))
names(prior_list_vague) <- c("Ctrl", "DG_1", "DG_2", "DG_3", "DG_4")
```

## Specification of the New Trial Design

For the hypothetical new trial, we plan with 4 active dose levels and we
specify a broad set of potential dose-response relationships, including
a linear, an exponential, an emax and 2 sigEMAX models.  
Furthermore, we assume a maximum effect of -3 on top of control
(i.e. assuming that active treatment can reduce the MADRS score after 8
weeks by up to 15.8) and plan a trial with 80 patients for all active
groups and 60 patients for control.

``` r
exp     <- DoseFinding::guesst(
  d     = 5,
  p     = c(0.2),
  model = "exponential",
  Maxd  = max(dose_levels))

emax    <- DoseFinding::guesst(
  d     = 2.5,
  p     = c(0.9),
  model = "emax")

sigemax <- DoseFinding::guesst(
  d     = c(2.5, 5),
  p     = c(0.1, 0.6),
  model = "sigEmax")

sigemax2 <- DoseFinding::guesst(
  d     = c(2, 4),
  p     = c(0.3, 0.8),
  model = "sigEmax")

mods <- DoseFinding::Mods(
  linear      = NULL,
  emax        = emax,
  exponential = exp,
  sigEmax     = rbind(sigemax, sigemax2),
  doses       = dose_levels,
  maxEff      = -3,
  placEff     = -12.8)

n_patients <- c(60, 80, 80, 80, 80)
```

## Calculation of the Success Probabilities

To calculate success probabilities for the different assumed
dose-response models and the specified trial design we will apply the
assessDesign function. For illustration purposes, the number of
simulated trial results is reduced to 100 in this example.

``` r
set.seed(7015) # re-sets seed only for this example; remove in your analysis script
success_probabilities <- assessDesign(
  n_patients  = n_patients,
  mods        = mods,
  prior_list  = prior_list,
  sd          = sd_tot,
  n_sim       = 100) # speed up example run-time

success_probabilities
#> $linear
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.68 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.59  0.21  0.59  0.48  0.36 
#> 
#> $emax
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.82 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.26  0.78  0.23  0.36  0.65 
#> 
#> $exponential
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.63 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.61  0.16  0.60  0.42  0.30 
#> 
#> $sigEmax1
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.82 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.63  0.52  0.58  0.73  0.76 
#> 
#> $sigEmax2
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.84 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.43  0.73  0.37  0.61  0.77 
#> 
#> attr(,"avgSuccessRate")
#> [1] 0.758
#> attr(,"placEff")
#> [1] -12.8
#> attr(,"maxEff")
#> [1] -3
#> attr(,"sampleSize")
#> [1] 60 80 80 80 80
#> attr(,"priorESS")
#>  Ctr DG_1 DG_2 DG_3 DG_4 
#> 20.6  1.0  1.0  1.0  1.0
```

As an alternative, we will evaluate a design with the same overall
sample size but allocating more patients on the highest dose group and
control.

``` r
set.seed(7015) # re-sets seed only for this example; remove in your analysis script
success_probabilities_uneq <- assessDesign(
  n_patients  = c(80, 60, 60, 60, 120),
  mods        = mods,
  prior_list  = prior_list,
  sd          = sd_tot,
  n_sim       = 100) # speed up example run-time
success_probabilities_uneq
#> $linear
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.8 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.76  0.33  0.75  0.64  0.53 
#> 
#> $emax
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.88 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.45  0.82  0.39  0.56  0.73 
#> 
#> $exponential
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.79 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.76  0.30  0.76  0.64  0.51 
#> 
#> $sigEmax1
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.86 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.73  0.51  0.69  0.82  0.82 
#> 
#> $sigEmax2
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.89 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.57  0.77  0.56  0.68  0.82 
#> 
#> attr(,"avgSuccessRate")
#> [1] 0.844
#> attr(,"placEff")
#> [1] -12.8
#> attr(,"maxEff")
#> [1] -3
#> attr(,"sampleSize")
#> [1]  80  60  60  60 120
#> attr(,"priorESS")
#>  Ctr DG_1 DG_2 DG_3 DG_4 
#> 20.6  1.0  1.0  1.0  1.0
```

For this specific trial setting the adapted allocation ratio leads to
increased success probabilities under all assumed dose response
relationships.

Instead of specifying the assumed effects via the models, it is also
possible to directly specify the effects for the individual dose levels
via the dr_means input. This allows e.g. also the simulation of
scenarios with a prior-data conflict.

``` r
set.seed(7015) # re-sets seed only for this example; remove in your analysis script
success_probabilities_dr <- assessDesign(
  n_patients  = c(60, 80, 80, 80, 80),
  mods        = mods,
  prior_list  = prior_list,
  sd          = sd_tot,
  dr_means    = c(-12, -14, -15, -16, -17),
  n_sim       = 100) # speed up example run-time
success_probabilities_dr
#> $dr_response
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.99 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.87  0.90  0.86  0.94  0.97 
#> 
#> attr(,"avgSuccessRate")
#> [1] 0.99
#> attr(,"placEff")
#> [1] -12
#> attr(,"maxEff")
#> [1] 5
#> attr(,"sampleSize")
#> [1] 60 80 80 80 80
#> attr(,"priorESS")
#>  Ctr DG_1 DG_2 DG_3 DG_4 
#> 20.6  1.0  1.0  1.0  1.0
```

## Assessment of the Minimally Efficacious Dose

The assessment of the Minimally Efficacious Dose (MED) is integrated in
the
[`assessDesign()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/assessDesign.md)
function via the arguments `delta` and `evidence_level`. If only the
argument `delta` is provided, the estimated model shapes will be used to
assess the MED. If both the arguments `delta` and `evidence_level` are
provided, a Bayesian decision rule of the form
$$\widehat{\text{MED}} = \text{arg min}_{d \in \{ d_{1},\ldots,d_{k}\}}\left\{ \text{Pr}\left( f\left( d,\widehat{\theta} \right) - f\left( d_{1},\widehat{\theta} \right) > \Delta \right) > \gamma \right\}$$
will be applied, see also `?getMED()`. The computational cost for the
Bayesian decision rule within the
[`assessDesign()`](https://boehringer-ingelheim.github.io/BayesianMCPMod/reference/assessDesign.md)
function is rather high depending on the number of simulated trial
outcomes, as the required quantiles need to be bootstrapped for each
true model shape and each simulated outcome. Thus, using the Bayesian
decision rule when assessing the trial’s design is only recommended when
using parallel computing.

``` r
set.seed(7015) # re-sets seed only for this example; remove in your analysis script
success_probabilities_med <- assessDesign(
  n_patients  = c(60, 80, 80, 80, 80),
  mods        = mods,
  prior_list  = prior_list,
  sd          = sd_tot,
  delta       = 2,
  n_sim       = 100) # speed up example run-time
success_probabilities_med
#> $linear
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.68 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.59  0.21  0.59  0.48  0.36 
#> MED Assessment
#>   Selection Method:    avgFit 
#>   Identification Rate: 0.68 
#>    Dose Level:  2.5  5.0 10.0 20.0 
#>    MED Freq:   0.00 0.03 0.19 0.46 
#>   MED not reached Freq:        0 
#>   No success in MCP step Freq: 0.32 
#> 
#> $emax
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.82 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.26  0.78  0.23  0.36  0.65 
#> MED Assessment
#>   Selection Method:    avgFit 
#>   Identification Rate: 0.8 
#>    Dose Level:  2.5  5.0 10.0 20.0 
#>    MED Freq:   0.51 0.09 0.14 0.06 
#>   MED not reached Freq:        0.02 
#>   No success in MCP step Freq: 0.18 
#> 
#> $exponential
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.63 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.61  0.16  0.60  0.42  0.30 
#> MED Assessment
#>   Selection Method:    avgFit 
#>   Identification Rate: 0.63 
#>    Dose Level:  2.5  5.0 10.0 20.0 
#>    MED Freq:   0.00 0.02 0.16 0.45 
#>   MED not reached Freq:        0 
#>   No success in MCP step Freq: 0.37 
#> 
#> $sigEmax1
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.82 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.63  0.52  0.58  0.73  0.76 
#> MED Assessment
#>   Selection Method:    avgFit 
#>   Identification Rate: 0.82 
#>    Dose Level:  2.5  5.0 10.0 20.0 
#>    MED Freq:   0.02 0.17 0.37 0.26 
#>   MED not reached Freq:        0 
#>   No success in MCP step Freq: 0.18 
#> 
#> $sigEmax2
#> Bayesian Multiple Comparison Procedure
#>   Estimated Success Rate: 0.84 
#>   N Simulations:          100 
#>    Model Shape:         lin  emax   exp sigE1 sigE2 
#>    Significance Freq:  0.43  0.73  0.37  0.61  0.77 
#> MED Assessment
#>   Selection Method:    avgFit 
#>   Identification Rate: 0.83 
#>    Dose Level:  2.5  5.0 10.0 20.0 
#>    MED Freq:   0.25 0.17 0.22 0.19 
#>   MED not reached Freq:        0.01 
#>   No success in MCP step Freq: 0.16 
#> 
#> attr(,"avgSuccessRate")
#> [1] 0.758
#> attr(,"avgMEDIdentificationRate")
#> [1] 0.752
#> attr(,"placEff")
#> [1] -12.8
#> attr(,"maxEff")
#> [1] -3
#> attr(,"sampleSize")
#> [1] 60 80 80 80 80
#> attr(,"priorESS")
#>  Ctr DG_1 DG_2 DG_3 DG_4 
#> 20.6  1.0  1.0  1.0  1.0
```
