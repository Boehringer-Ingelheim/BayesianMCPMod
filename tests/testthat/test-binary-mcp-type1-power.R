test_that("MCP step sanity: low type-I error under null and decent power under monotone alternative (binary)", {
  skip_if_not_installed("RBesT")
  skip_if_not_installed("DoseFinding")
  
  set.seed(202)
  
  dose_levels <- c(0, 1, 2, 4)
  n_per_dose  <- 90
  R_null <- 40
  R_alt  <- 40
  
  mods <- DoseFinding::Mods(
    linear = NULL,
    emax   = c(0.5, 1.2),
    doses  = dose_levels,
    maxEff = 2
  )
  
  # Mildly regularizing prior centered near a reasonable baseline.
  # This keeps separation-ish cases from blowing up while still letting data drive.
  make_prior <- function(p_center) {
    setNames(
      lapply(seq_along(dose_levels), function(i) {
        RBesT::mixnorm(comp1 = c(w = 1, m = qlogis(p_center[i]), s = 2.0), sigma = 2)
      }),
      c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
    )
  }
  
  # Extract "significant" robustly.
  # Prefer structure if present; fall back to parsing printed output.
  extract_significant <- function(bmcp_obj) {
    # 1) Common possibilities
    if (is.list(bmcp_obj)) {
      for (nm in c("Significant", "significant", "sig")) {
        if (!is.null(bmcp_obj[[nm]])) return(as.integer(bmcp_obj[[nm]]))
      }
    }
    if (is.matrix(bmcp_obj) || is.data.frame(bmcp_obj)) {
      rn <- rownames(bmcp_obj)
      if (!is.null(rn) && "Significant" %in% rn) {
        return(as.integer(bmcp_obj["Significant", 1]))
      }
      if ("Significant" %in% colnames(bmcp_obj)) {
        return(as.integer(bmcp_obj[1, "Significant"]))
      }
    }
    
    # 2) Parse printed output (S3 print shows "Significant: <0/1>")
    out <- paste(capture.output(print(bmcp_obj)), collapse = "\n")
    m <- regexec("Significant:\\s*([01])", out)
    reg <- regmatches(out, m)[[1]]
    if (length(reg) >= 2) return(as.integer(reg[2]))
    
    stop("Could not extract 'Significant' indicator from BayesianMCP object")
  }
  
  run_one <- function(p_vec, crit_prob_adj = 0.95) {
    dat <- data.frame(
      simulation = 1L,
      dose       = rep(dose_levels, each = n_per_dose),
      response   = unlist(lapply(p_vec, function(p) rbinom(n_per_dose, 1, p)))
    )
    attr(dat, "probability_scale") <- TRUE
    
    prior_list <- make_prior(p_center = p_vec)
    
    post <- getPosterior(
      prior_list        = prior_list,
      data              = dat,
      probability_scale = TRUE
    )
    
    contr <- getContr(
      mods         = mods,
      dose_levels  = dose_levels,
      dose_weights = rep(n_per_dose, length(dose_levels)),
      prior_list   = prior_list
    )
    
    res <- performBayesianMCPMod(
      posterior_list    = post,
      contr             = contr,
      crit_prob_adj     = crit_prob_adj,
      simple            = TRUE,
      avg_fit           = TRUE,
      probability_scale = TRUE
    )
    
    extract_significant(res$BayesianMCP)
  }
  
  # ---- Null (flat) ----
  p_null <- rep(0.25, length(dose_levels))
  sig_null <- replicate(R_null, run_one(p_null, crit_prob_adj = 0.95))
  
  type1 <- mean(sig_null == 1)
  
  # Type-I error should be "low".
  # With only 40 reps, allow stochastic wiggle, but it should not be huge.
  expect_true(type1 <= 0.20)
  
  # ---- Alternative (monotone increasing) ----
  p_alt <- c(0.10, 0.18, 0.32, 0.55)
  sig_alt <- replicate(R_alt, run_one(p_alt, crit_prob_adj = 0.95))
  
  power <- mean(sig_alt == 1)
  
  # With n_per_dose=90 and a pretty clear signal, we expect decent power.
  expect_true(power >= 0.55)
})
