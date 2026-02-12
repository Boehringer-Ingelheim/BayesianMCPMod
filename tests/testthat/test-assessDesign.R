# Tests for assessDesign --------------------------------------------------



test_that("base case input throws no error and has correct properties", {
  expect_no_error(
    eval_design <- assessDesign(
      n_patients = n_patients,
      mods = mods,
      sd = sd,
      prior_list = prior_list,
      n_sim = n_sim,
      alpha_crit_val = alpha_crit_val,
      simple = TRUE
    )
  )
  
  # assessDesign should give results for each model in mods
  expect_equal(
    names(eval_design), names(mods)
  )
  
  # assessDesign result should have rows = n_sim
  expect_equal(
    attr(eval_design$linear, "dim")[1],
    n_sim
  )
  
  # assessDesign result (in this base case) should have crit_prob = 1 - alpha_crit_val
  expect_equal(
    attr(eval_design$linear, "critProb"),
    1 - alpha_crit_val
  )
  
  contr_mat <- getContr(
    mods = mods,
    dose_levels = dose_levels,
    dose_weights = n_patients,
    prior_list = prior_list
  )
  
  expect_no_error(
    eval_design <- assessDesign(
      n_patients = n_patients,
      mods = mods,
      sd = sd,
      prior_list = prior_list,
      n_sim = n_sim,
      alpha_crit_val = alpha_crit_val,
      simple = TRUE,
      modeling = TRUE
    )
  )
  
  # assessDesign result should have rows = n_sim
  expect_equal(
    attr(eval_design$linear$BayesianMCP, "dim")[1],
    n_sim
  )
  
  # assessDesign result (in this base case) should have crit_prob = 1 - alpha_crit_val
  expect_equal(
    attr(eval_design$linear$BayesianMCP, "critProb"),
    1 - alpha_crit_val
  )
  
  expect_no_error(
    assessDesign(
      n_patients = n_patients,
      mods = mods,
      sd = sd,
      prior_list = prior_list,
      n_sim = n_sim,
      alpha_crit_val = alpha_crit_val,
      simple = TRUE,
      reestimate = TRUE,
      contr = contr_mat
    )
  )
  
  
  sd_tot <- 9.4
  
  dose_levels <- c(0, 2.5, 5, 10, 20)
  
  prior_list <- lapply(dose_levels, function(dose_group) {
    RBesT::mixnorm(weak = c(w = 1, m = 0, s = 200), sigma = 10)
  })
  
  names(prior_list) <- c("Ctr", paste0("DG_", dose_levels[-1]))
  
  exp <- DoseFinding::guesst(
    d     = 5,
    p     = c(0.2),
    model = "exponential",
    Maxd  = max(dose_levels)
  )
  
  emax <- DoseFinding::guesst(
    d     = 2.5,
    p     = c(0.9),
    model = "emax"
  )
  
  sigemax <- DoseFinding::guesst(
    d     = c(2.5, 5),
    p     = c(0.1, 0.6),
    model = "sigEmax"
  )
  
  sigemax2 <- DoseFinding::guesst(
    d     = c(2, 4),
    p     = c(0.3, 0.8),
    model = "sigEmax"
  )
  
  mods <- DoseFinding::Mods(
    linear      = NULL,
    emax        = emax,
    exponential = exp,
    sigEmax     = rbind(sigemax, sigemax2),
    doses       = dose_levels,
    maxEff      = -3,
    placEff     = -12.8
  )
  
  n_patients <- c(60, 80, 80, 80, 80)
  
  expect_no_error(
    assessDesign(
      n_patients = n_patients,
      mods       = mods,
      prior_list = prior_list,
      sd         = sd_tot,
      n_sim      = 10,
      reestimate = TRUE
    )
  )
})


### n_patients param ###

test_that("assessDesign validates n_patients parameter input and give appropriate error messages", {
  # assertions that aren't tested here for sake of brevity
  # n_patients should be a non-NULL numeric vector
  
  expect_error(
    assessDesign(n_patients = n_patients[-1], sd = sd, mods = mods, prior_list = prior_list, n_sim = n_sim)
  )
  
  expect_error(
    assessDesign(n_patients = rep(1, length(n_patients)), sd = sd, mods = mods, prior_list = prior_list, n_sim = n_sim),
  )
})

### mods param ###

test_that("assessDesign validates mods parameter input and give appropriate error messages", {
  # assertions that aren't tested here for sake of brevity
  # mods should be non-NULL object of class "Mods" from {DoseFinding}
  
  
  # checking that DoseFinding didn't change how they named their 'doses' attribute
  expect_true(
    "doses" %in% names(attributes(mods))
  )
  
  mods2 <- mods
  attr(mods2, "doses") <- 0
  expect_error(
    assessDesign(n_patients = n_patients, mods = mods2, sd = sd, prior_list = prior_list, n_sim = n_sim)
  )
  rm(mods2)
})

## prior_list param ###

test_that("assessDesign validates prior_list parameter input and give appropriate error messages", {
  # assertions that aren't tested here for sake of brevity
  # prior_list should be a non-NULL named list with length = number of dose levels
  # length(attr(prior_list, "dose_levels")) == n_patients (see above)
  
  # checking that we didn't change how we named the 'dose_levels' attribute
  expect_true(
    "doses" %in% names(attributes(mods))
  )
})

test_that("assessDesign: input validation branches are triggered", {
  skip_if_not_installed("DoseFinding")
  skip_if_not_installed("RBesT")
  
  dose_levels <- c(0, 1, 2)
  mods <- DoseFinding::Mods(
    linear = NULL,
    doses  = dose_levels,
    maxEff = 1
  )
  
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 2), sigma = 2)
    }),
    c("Ctr", "DG_1", "DG_2")
  )
  
  # n_patients < 2
  expect_error(
    assessDesign(
      n_patients = c(1, 2, 2),
      mods = mods,
      prior_list = prior_list,
      n_sim = 1
    )
  )
  
  # length mismatch
  expect_error(
    assessDesign(
      n_patients = c(10, 10),
      mods = mods,
      prior_list = prior_list,
      n_sim = 1
    )
  )
  
  # conflicting inputs
  expect_error(
    assessDesign(
      n_patients = c(10, 10, 10),
      mods = mods,
      prior_list = prior_list,
      sd = 1,
      data_sim = data.frame(),
      n_sim = 1
    )
  )
  
  expect_error(
    assessDesign(
      n_patients = c(10, 10, 10),
      mods = mods,
      prior_list = prior_list,
      sd = 1,
      estimates_sim = list(),
      n_sim = 1
    )
  )
})


test_that("assessDesign: binary endpoint runs and returns per-true-model results with avgSuccessRate attribute", {
  skip_if_not_installed("DoseFinding")
  skip_if_not_installed("RBesT")
  
  set.seed(401)
  
  dose_levels <- c(0, 1, 2)
  n_patients  <- c(20, 20, 20)
  
  mods <- DoseFinding::Mods(
    linear = NULL,
    doses  = dose_levels,
    maxEff = 1
  )
  
  # Logit-scale priors
  p_true <- c(0.10, 0.20, 0.40)
  prior_list <- setNames(
    lapply(seq_along(p_true), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = qlogis(p_true[i]), s = 2), sigma = 2)
    }),
    c("Ctr", "DG_1", "DG_2")
  )
  
  res <- assessDesign(
    n_patients = n_patients,
    mods = mods,
    prior_list = prior_list,
    n_sim = 2,
    probability_scale = TRUE,
    modeling = FALSE
  )
  
  # Returned object is a list with names == true underlying model names
  expect_true(is.list(res))
  expect_true(length(res) >= 1)
  expect_true(!is.null(names(res)))
  
  # avgSuccessRate is stored as attribute on the returned list
  expect_true(!is.null(attr(res, "avgSuccessRate")))
  expect_true(is.numeric(attr(res, "avgSuccessRate")))
  expect_true(attr(res, "avgSuccessRate") >= 0 && attr(res, "avgSuccessRate") <= 1)
  
  # Each element (when modeling=FALSE) is the BayesianMCP result and should have successRate attribute
  for (nm in names(res)) {
    expect_true(!is.null(attr(res[[nm]], "successRate")))
    sr <- attr(res[[nm]], "successRate")
    expect_true(is.numeric(sr))
    expect_true(sr >= 0 && sr <= 1)
  }
})


test_that("assessDesign: custom data_sim path requires model columns and returns results; contr warning is expected", {
  skip_if_not_installed("DoseFinding")
  skip_if_not_installed("RBesT")
  
  set.seed(402)
  
  dose_levels <- c(0, 1)
  n_patients  <- c(30, 30)
  
  mods <- DoseFinding::Mods(
    linear = NULL,
    doses  = dose_levels,
    maxEff = 1
  )
  
  prior_list <- setNames(
    lapply(seq_along(dose_levels), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 2), sigma = 2)
    }),
    c("Ctr", "DG_1")
  )
  
  # IMPORTANT: data_sim must match simulateData() structure:
  # first 3 columns: simulation, ptno, dose
  # then ONE OR MORE model columns (true model responses), e.g. "lin"
  data_sim <- data.frame(
    simulation = 1L,
    ptno = seq_len(sum(n_patients)),
    dose = rep(dose_levels, times = n_patients),
    lin  = rnorm(sum(n_patients))
  )
  
  # assessDesign will message about providing contr; we allow that message
  expect_message(
    res <- assessDesign(
      n_patients = n_patients,
      mods = mods,
      prior_list = prior_list,
      data_sim = data_sim,
      n_sim = 1
    ),
    regexp = "Consider to provide 'contr'"
  )
  
  expect_true(is.list(res))
  expect_true("lin" %in% names(res))
  expect_true(!is.null(attr(res, "avgSuccessRate")))
})


test_that("assessDesign: modeling + delta attaches MED in a robust form and sets avgMEDIdentificationRate", {
  skip_if_not_installed("DoseFinding")
  skip_if_not_installed("RBesT")
  
  set.seed(403)
  
  dose_levels <- c(0, 1, 2)
  n_patients  <- c(25, 25, 25)
  
  mods <- DoseFinding::Mods(
    linear = NULL,
    doses  = dose_levels,
    maxEff = 1
  )
  
  p_flat <- c(0.25, 0.25, 0.25)
  
  prior_list <- setNames(
    lapply(seq_along(p_flat), function(i) {
      RBesT::mixnorm(comp1 = c(w = 1, m = qlogis(p_flat[i]), s = 1.5), sigma = 2)
    }),
    c("Ctr", "DG_1", "DG_2")
  )
  
  res_avg <- assessDesign(
    n_patients = n_patients,
    mods = mods,
    prior_list = prior_list,
    n_sim = 3,
    probability_scale = TRUE,
    modeling = TRUE,
    delta = 0.2,
    med_selection = "avgFit"
  )
  
  expect_true(is.list(res_avg))
  expect_true(length(res_avg) >= 1)
  
  # Top-level identification rate attribute
  ir <- attr(res_avg, "avgMEDIdentificationRate")
  expect_true(!is.null(ir))
  expect_true(is.numeric(ir))
  expect_true(ir >= 0 && ir <= 1)
  
  # MED is stored as attribute on each element, but its structure may vary.
  # We just need it to contain a med_reached indicator (0/1) in some form.
  extract_med_reached <- function(med) {
    if (is.null(med)) return(NULL)
    
    # case 1: named atomic vector/list
    if (!is.null(names(med)) && "med_reached" %in% names(med)) {
      return(as.numeric(med[["med_reached"]]))
    }
    
    # case 2: matrix/data.frame with rowname
    if ((is.matrix(med) || is.data.frame(med)) && !is.null(rownames(med)) &&
        "med_reached" %in% rownames(med)) {
      return(as.numeric(med["med_reached", ]))
    }
    
    # case 3: matrix/data.frame with column name
    if ((is.matrix(med) || is.data.frame(med)) && "med_reached" %in% colnames(med)) {
      return(as.numeric(med[, "med_reached"]))
    }
    
    # case 4: 2-row matrix without rownames where first row is med_reached
    if (is.matrix(med) && nrow(med) >= 1) {
      # Heuristic: first row should be 0/1-like if it's med_reached
      cand <- as.numeric(med[1, ])
      if (all(cand %in% c(0, 1, NA))) return(cand)
    }
    
    NULL
  }
  
  any_has_med <- FALSE
  any_has_reached <- FALSE
  
  for (nm in names(res_avg)) {
    med <- attr(res_avg[[nm]], "MED")
    if (!is.null(med)) any_has_med <- TRUE
    
    mr <- extract_med_reached(med)
    if (!is.null(mr)) {
      any_has_reached <- TRUE
      mr <- mr[!is.na(mr)]
      expect_true(all(mr %in% c(0, 1)))
    }
  }
  
  # Invariant expectations: at least one element should carry MED,
  # and at least one should expose a med_reached indicator in some form.
  expect_true(any_has_med)
  expect_true(any_has_reached)
})
