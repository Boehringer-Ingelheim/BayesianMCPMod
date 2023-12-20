etPriorList <- function (
    
  hist_data,
  dose_levels,
  dose_names    = NULL,
  robust_weight = 0.5
  
) {
  
  sd_tot <- with(hist_data, sum(sd * n) / sum(n))
  
  gmap <- RBesT::gMAP(
    formula    = cbind(est, se) ~ 1 | trial,
    weights    = hist_data$n,
    data       = hist_data,
    family     = gaussian,
    beta.prior = cbind(0, 100 * sd_tot),
    tau.dist   = "HalfNormal",
    tau.prior  = cbind(0, sd_tot / 4))
  
  prior_ctr <- RBesT::automixfit(gmap)
  
  if (!is.null(robust_weight)) {
    
    prior_ctr <- suppressMessages(RBesT::robustify(
      priormix = prior_ctr,
      weight   = robust_weight,
      sigma    = sd_tot))
    
  }
  
  prior_trt <- RBesT::mixnorm(
    comp1 = c(w = 1, m = summary(prior_ctr)[1], n = 1),
    sigma = sd_tot,
    param = "mn")
  
  prior_list <- c(list(prior_ctr),
                  rep(x     = list(prior_trt),
                      times = length(dose_levels[-1])))
  
  if (is.null(dose_names)) {
    
    dose_names <- c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
    
  }
  
  names(prior_list) <- dose_names
  
  return (prior_list)
  
}
# Create minimal test case
n_hist_trials = 1

hist_data <- data.frame(
  trial = seq(1, n_hist_trials, 1),
  est   = rep(1, n_hist_trials),
  se    = rep(1, n_hist_trials),
  sd    = rep(1, n_hist_trials),
  n     = rep(1, n_hist_trials)
)

n_patients <- c(2, 1)
dose_levels <- c(0, 2.5)
mean <- c(8, 12)
sd <- c(0.5, 0.8)

mods <- DoseFinding::Mods(
  linear = NULL, 
  doses = dose_levels
)

prior_list <- getPriorList(
  hist_data   = hist_data,
  dose_levels = dose_levels,
  robustify_weight = 0.5
)

n_sim = 1
alpha_crit_val = 0.05
simple = TRUE

# eval_design <- assessDesign(
#   n_patients = n_patients, 
#   mods = mods, 
#   prior_list = prior_list,
#   n_sim = n_sim,
#   alpha_crit_val = alpha_crit_val,
#   simple = TRUE
# )