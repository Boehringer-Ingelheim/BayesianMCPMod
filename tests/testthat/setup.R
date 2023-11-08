
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
