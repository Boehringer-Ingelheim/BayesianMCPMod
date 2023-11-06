library(clinDR)
library(dplyr)

set.seed(2023)
data("metaData")
testdata    <- as.data.frame(metaData)
dataset     <- filter(testdata, bname == "BRINTELLIX")
histcontrol <- filter(dataset, dose == 0, primtime == 8, indication == "MAJOR DEPRESSIVE DISORDER")

##Create MAP Prior
hist_data <- data.frame(
  trial = histcontrol$nctno,
  est   = histcontrol$rslt,
  se    = histcontrol$se,
  sd    = histcontrol$sd,
  n     = histcontrol$sampsize)

dose_levels <- c(0, 2.5, 5, 10, 20)

prior_list <- getPriorList(
  hist_data   = hist_data,
  dose_levels = dose_levels,
  robustify_weight = 0.5)
#Pre-Specification (B)MCPMod 

## candidate models for MCPMod
# linear function - no guestimates needed
exp     <- DoseFinding::guesst(d     = 5,
                               p     = c(0.2),
                               model = "exponential",
                               Maxd  = max(dose_levels))
emax    <- DoseFinding::guesst(d     = 2.5,
                               p     = c(0.9),
                               model = "emax")
sigemax<-  DoseFinding::guesst(d     = c(2.5,5),
                               p     = c(0.1,0.6),
                               model = "sigEmax")
#beta <- DoseFinding::guesst(d=5, p=0.8, model="betaMod", dMax=1, scal=1.2, Maxd=20)

mods <- DoseFinding::Mods(
  linear      = NULL,
  emax        = emax,
  exponential = exp,
  sigEmax     = sigemax,
  #betaMod     = beta,
  doses       = dose_levels,
  maxEff      = -3,
  placEff     = -12.8)

n_patients = c(60, 80, 80, 80, 80)

test_that("assessDesign works as intented", {
  success_probabilities <- assessDesign(
    n_patients = n_patients,
    mods = mods,
    prior_list = prior_list)
  success_probabilities_uneq <- assessDesign(
    n_patients = c(80, 60, 60, 60, 120),
    mods = mods,
    prior_list = prior_list)
  
  expect_type(success_probabilities, "list")
  expect_type(success_probabilities_uneq, "list")
})