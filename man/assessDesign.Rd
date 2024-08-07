% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BMCPMod.R
\name{assessDesign}
\alias{assessDesign}
\title{assessDesign
.}
\usage{
assessDesign(
  n_patients,
  mods,
  prior_list,
  sd,
  n_sim = 1000,
  alpha_crit_val = 0.05,
  simple = TRUE,
  reestimate = FALSE,
  contr = NULL,
  dr_means = NULL
)
}
\arguments{
\item{n_patients}{Vector specifying the planned number of patients per dose group}

\item{mods}{An object of class "Mods" as specified in the DoseFinding package.}

\item{prior_list}{A prior_list object specifying the utilized prior for the different dose groups}

\item{sd}{A positive value, specification of assumed sd}

\item{n_sim}{Number of simulations to be performed}

\item{alpha_crit_val}{(Un-adjusted) Critical value to be used for the MCP testing step. Passed to the getCritProb() function for the calculation of adjusted critical values (on the probability scale). Default is 0.05.}

\item{simple}{Boolean variable defining whether simplified fit will be applied. Passed to the getModelFits function. Default FALSE.}

\item{reestimate}{Boolean variable defining whether critical value should be calculated with re-estimated contrasts (see getCritProb function for more details). Default FALSE}

\item{contr}{An object of class 'optContr' as created by the getContr() function. Allows specification of a fixed contrasts matrix. Default NULL}

\item{dr_means}{A vector, allows specification of  individual (not model based) assumed effects per dose group. Default NULL}
}
\value{
Returns success probabilities for the different assumed dose-response shapes, attributes also includes information around average success rate (across all assumed models) and prior Effective sample size
}
\description{
This function performs simulation based trial design evaluations for a set of specified dose-response models
}
\examples{
if (interactive()) { # takes typically > 5 seconds

mods <- DoseFinding::Mods(linear      = NULL,
                          linlog      = NULL,
                          emax        = c(0.5, 1.2),
                          exponential = 2,
                          doses       = c(0, 0.5, 2,4, 8),
                          maxEff      = 6)
sd <- 12
prior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 12), sigma = 2),
                   DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 1, s = 12), sigma = 2),
                   DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.2, s = 11), sigma = 2) ,
                   DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 1.3, s = 11), sigma = 2) ,
                   DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 2, s = 13), sigma = 2))
n_patients <- c(40, 60, 60, 60, 60)

success_probabilities <- assessDesign(
  n_patients  = n_patients,
  mods        = mods,
  prior_list  = prior_list,
  sd          = sd,
  n_sim       = 1e2) # speed up exammple run time

success_probabilities

}

}
