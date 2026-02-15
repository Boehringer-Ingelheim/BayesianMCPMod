# Simulation Example of Bayesian MCPMod for binary Data

## Background and Data

In this vignette, we will show the use of the `BayesianMCPMod` R package
for trial planning for binary data.

For this setting the main specifications need to happen on the logit
scale (see also add ([binary analysis example
vignette](https://boehringer-ingelheim.github.io/BayesianMCPMod/articles/binary_endpoint.md)))

As in the ([analysis example
vignette](https://boehringer-ingelheim.github.io/BayesianMCPMod/articles/analysis_normal.md)),
we make use of historical data that is included in the clinDR package.
More specifically, trial results for XELJANZ will be utilized to
establish an informative prior for the control group.

This package makes use of the
[future](https://cran.r-project.org/package=future) framework for
parallel processing, which can be set up as follows:

``` r
future::plan(future::multisession)
```

``` r
suppressPackageStartupMessages({
library(BayesianMCPMod) 
library(RBesT)
library(clinDR)
library(dplyr)
library(tibble)
library(reactable)
library(DoseFinding)
library(BayesianMCPMod)
})


set.seed(7015)
```

``` r
data("metaData")
testdata    <- as.data.frame(metaData)
dataset     <- filter(testdata, bname == "XELJANZ")
histcontrol <- filter(dataset, dose == 0, primtime == 12, indication == "RHEUMATOID ARTHRITIS")

hist_data <- data.frame(
  trial = histcontrol$nctno,
  est   = histcontrol$rslt,
  se    = histcontrol$se,
  sd    = histcontrol$sd,
  r     = round(histcontrol$sampsize*histcontrol$rslt),
  n     = histcontrol$sampsize)

sd_tot <- with(hist_data, sum(sd * n) / sum(n))
```

## Calculation of a MAP Prior

In a first step, a meta analytic prior will be calculated. The approach
to establish the prior is the same as outlined in the ([binary analysis
example
vignette](https://boehringer-ingelheim.github.io/BayesianMCPMod/articles/binary_endpoint.md))).
Please note that only information from the control group will be
integrated leading to an informative mixture prior for the control
group, while for the active groups a non-informative prior will be
specified.

``` r
dose_levels <- c(0, 2.5, 5, 10,20)

# 1) Establish MAP prior (beta mixture distribution) 
set.seed(7015) # re-sets seed only for this example; remove in your analysis script
map <- gMAP(cbind(hist_data$r, hist_data$n - hist_data$r) ~ 1|histcontrol$nctno, family = binomial, tau.dist = "HalfNormal",
            tau.prior = 0.5, beta.prior = (1/sqrt(0.1*0.9)), warmup = 1000, iter = 10000, chains = 2, thin = 1)
```

    ## Assuming default prior location   for beta: 0

``` r
map
```

    ## Generalized Meta Analytic Predictive Prior Analysis
    ## 
    ## Call:  gMAP(formula = cbind(hist_data$r, hist_data$n - hist_data$r) ~ 
    ##     1 | histcontrol$nctno, family = binomial, tau.dist = "HalfNormal", 
    ##     tau.prior = 0.5, beta.prior = (1/sqrt(0.1 * 0.9)), iter = 10000, 
    ##     warmup = 1000, thin = 1, chains = 2)
    ## 
    ## Exchangeability tau strata: 1 
    ## Prediction tau stratum    : 1 
    ## Maximal Rhat              : 1 
    ## 
    ## Between-trial heterogeneity of tau prediction stratum
    ##   mean     sd   2.5%    50%  97.5% 
    ## 0.3970 0.3010 0.0156 0.3340 1.1200 
    ## 
    ## MAP Prior MCMC sample
    ##   mean     sd   2.5%    50%  97.5% 
    ## 0.1780 0.1130 0.0363 0.1550 0.4900

``` r
prior <- automixfit(map) #fits mixture distribution from MCMC samples from above
ess(prior)
```

    ## [1] 17.48444

``` r
#ess(prior)
p<-summary(prior)[1]

#ii) Robustify prior
prior.rob<-RBesT::robustify(
      priormix = prior,
      mean     = 0.5,
      weight   = 0.2)

#ess(prior.rob)

#iii) Translate prior to logit scale (to approximate via normal mixture model)
r <- rmix(prior.rob, n=1e4)
log.r <- logit(r)
prior.ctr<- automixfit(log.r, type = "norm")
sigma(prior.ctr)<-sqrt(1/(p*(1-p)))
#ess(prior.ctr, sigma = sqrt(1/(p*(1-p))))
#Specification of reference scale (this follows the idea of [@Neuenschwander2016]). 


#Specify a prior list
prior_trt <- RBesT::mixnorm(
    comp1 = c(w = 1, m = logit(summary(prior)[1]), n = 1),
    sigma = sqrt(1/(p*(1-p))),
    param = "mn")
  
  prior_list <- c(list(prior.ctr),
                  rep(x     = list(prior_trt),
                      times = length(dose_levels[-1])))

dose_names <- c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
names(prior_list) <- dose_names
```

Kindly note that a vague prior could be implemented via

``` r
prior_list_vague <- rep(list(RBesT::mixnorm(comp1 = c(w = 1, m = logit(p), n=1),
                                            sigma = sqrt(1/(p*(1-p))), param = "mn")),
                        times = length(dose_levels))
names(prior_list_vague) <- c("Ctrl", "DG_1", "DG_2", "DG_3", "DG_4")
```

## Specification of the New Trial Design

For the hypothetical new trial, we plan with 4 active dose levels and we
specify a broad set of potential dose-response relationships, including
a linear, an exponential, an emax, a logistic and a sigEMAX models.  
Furthermore, we assume a maximum response rate of 48% (resp. 30% on top
of control) and plan a trial with 40 patients for all active groups and
30 patients for control. Please note that the response rates need to be
provided on the **logit scale**.

``` r
n_patients <- c(30, 40, 40, 40, 40)

models <- Mods(
  linear = NULL,
  sigEmax = c(10, 5),
  logistic = c(11, 15),
  exponential = 10,
  emax = 2,
  doses = dose_levels,
  placEff = RBesT::logit(0.18),
  maxEff = (RBesT::logit(0.48) - RBesT::logit(0.18))
)
```

## Calculation of the Success Probabilities

To calculate success probabilities for the different assumed
dose-response models and the specified trial design we will apply the
assessDesign function. We are not only interested in the success
probability for the testing step, but also the assessment of the
Minimally Efficacious Dose (MED). This effect delta is provided on the
probability scale. In our case we would like to see a difference of 20%
(compared to control) to claim efficacy.

For illustration purposes, the number of simulated trial results is
reduced to 100 in this example.

``` r
set.seed(7015) # re-sets seed only for this example; remove in your analysis script
success_probabilities <- assessDesign(
  n_patients  = n_patients,
  mods        = models,
  prior_list  = prior_list,
  n_sim       = 100, probability_scale=TRUE,delta          = 0.2
  ) 
success_probabilities
```

    ## $linear
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.96 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.91 0.90 0.91 0.89 0.73 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.88 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.00 0.00 0.05 0.83 
    ##   MED not reached Freq:        0.08 
    ##   No success in MCP step Freq: 0.04 
    ## 
    ## $sigEmax
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.99 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.98 0.97 0.98 0.97 0.72 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.92 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.00 0.00 0.10 0.82 
    ##   MED not reached Freq:        0.07 
    ##   No success in MCP step Freq: 0.01 
    ## 
    ## $logistic
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.99 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.96 0.94 0.96 0.95 0.77 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.92 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.00 0.02 0.13 0.77 
    ##   MED not reached Freq:        0.07 
    ##   No success in MCP step Freq: 0.01 
    ## 
    ## $exponential
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.97 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.92 0.92 0.93 0.97 0.52 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.88 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.00 0.00 0.00 0.88 
    ##   MED not reached Freq:        0.09 
    ##   No success in MCP step Freq: 0.03 
    ## 
    ## $emax
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.94 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.80 0.66 0.80 0.64 0.92 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.87 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.27 0.28 0.17 0.15 
    ##   MED not reached Freq:        0.07 
    ##   No success in MCP step Freq: 0.06 
    ## 
    ## attr(,"avgSuccessRate")
    ## [1] 0.97
    ## attr(,"avgMEDIdentificationRate")
    ## [1] 0.894
    ## attr(,"placEff")
    ## [1] -1.516347
    ## attr(,"maxEff")
    ## [1] 1.436305
    ## attr(,"sampleSize")
    ## [1] 30 40 40 40 40
    ## attr(,"priorESS")
    ##       Ctr DG_1.mean DG_2.mean DG_3.mean DG_4.mean 
    ##      10.6       1.0       1.0       1.0       1.0

As an alternative, we will evaluate a design with the same overall
sample size but allocating more patients on the highest dose group and
control.

``` r
set.seed(7015) # re-sets seed only for this example; remove in your analysis script
success_probabilities_uneq <- assessDesign(
  n_patients  = c(40, 30, 30, 30, 50),
  mods        = models,
  prior_list  = prior_list,
  n_sim       = 100,probability_scale=TRUE,delta          = 0.2) # speed up example run-time
success_probabilities_uneq
```

    ## $linear
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.98 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.94 0.95 0.94 0.93 0.86 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.88 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.00 0.02 0.09 0.77 
    ##   MED not reached Freq:        0.1 
    ##   No success in MCP step Freq: 0.02 
    ## 
    ## $sigEmax
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.99 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.99 0.97 0.99 0.98 0.79 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.95 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.00 0.00 0.12 0.83 
    ##   MED not reached Freq:        0.04 
    ##   No success in MCP step Freq: 0.01 
    ## 
    ## $logistic
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.98 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.98 0.96 0.98 0.96 0.87 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.89 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.01 0.02 0.10 0.76 
    ##   MED not reached Freq:        0.09 
    ##   No success in MCP step Freq: 0.02 
    ## 
    ## $exponential
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.99 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.99 0.99 0.99 0.99 0.75 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.95 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.00 0.00 0.04 0.91 
    ##   MED not reached Freq:        0.04 
    ##   No success in MCP step Freq: 0.01 
    ## 
    ## $emax
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.99 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.88 0.78 0.88 0.76 0.97 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.91 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.14 0.28 0.23 0.26 
    ##   MED not reached Freq:        0.08 
    ##   No success in MCP step Freq: 0.01 
    ## 
    ## attr(,"avgSuccessRate")
    ## [1] 0.986
    ## attr(,"avgMEDIdentificationRate")
    ## [1] 0.916
    ## attr(,"placEff")
    ## [1] -1.516347
    ## attr(,"maxEff")
    ## [1] 1.436305
    ## attr(,"sampleSize")
    ## [1] 40 30 30 30 50
    ## attr(,"priorESS")
    ##       Ctr DG_1.mean DG_2.mean DG_3.mean DG_4.mean 
    ##      10.6       1.0       1.0       1.0       1.0

As a second alternative we will use the vague prior:

``` r
set.seed(7015) # re-sets seed only for this example; remove in your analysis script
success_probabilities_vague <- assessDesign(
  n_patients  = n_patients,
  mods        = models,
  prior_list  = prior_list_vague,
  n_sim       = 100, probability_scale=TRUE,delta          = 0.2
  ) 
success_probabilities_vague
```

    ## $linear
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.91 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.89 0.85 0.88 0.86 0.56 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.85 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.00 0.01 0.09 0.75 
    ##   MED not reached Freq:        0.06 
    ##   No success in MCP step Freq: 0.09 
    ## 
    ## $sigEmax
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.99 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.97 0.97 0.97 0.96 0.65 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.91 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.00 0.00 0.13 0.78 
    ##   MED not reached Freq:        0.08 
    ##   No success in MCP step Freq: 0.01 
    ## 
    ## $logistic
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.97 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.94 0.93 0.94 0.92 0.69 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.9 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.00 0.05 0.14 0.71 
    ##   MED not reached Freq:        0.07 
    ##   No success in MCP step Freq: 0.03 
    ## 
    ## $exponential
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.96 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.91 0.91 0.92 0.96 0.49 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.86 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.00 0.00 0.01 0.85 
    ##   MED not reached Freq:        0.1 
    ##   No success in MCP step Freq: 0.04 
    ## 
    ## $emax
    ## Bayesian Multiple Comparison Procedure
    ##   Estimated Success Rate: 0.87 
    ##   N Simulations:          100 
    ##    Model Shape:        lin sigE  log  exp emax 
    ##    Significance Freq: 0.68 0.54 0.68 0.53 0.84 
    ## MED Assessment on Probability Scale
    ##   Selection Method:    avgFit 
    ##   Identification Rate: 0.83 
    ##    Dose Level:  2.5  5.0 10.0 20.0 
    ##    MED Freq:   0.23 0.27 0.11 0.22 
    ##   MED not reached Freq:        0.04 
    ##   No success in MCP step Freq: 0.13 
    ## 
    ## attr(,"avgSuccessRate")
    ## [1] 0.94
    ## attr(,"avgMEDIdentificationRate")
    ## [1] 0.87
    ## attr(,"placEff")
    ## [1] -1.516347
    ## attr(,"maxEff")
    ## [1] 1.436305
    ## attr(,"sampleSize")
    ## [1] 30 40 40 40 40
    ## attr(,"priorESS")
    ## Ctrl.mean DG_1.mean DG_2.mean DG_3.mean DG_4.mean 
    ##         1         1         1         1         1
