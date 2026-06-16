# **************************************************************************** #
# Project: Thesis course
# Project Title: A Unified Statistical Framework for Efficacy-Integrated Dose 
# Optimisation Designs in Early-Phase Oncology Trials
# Do File: Script that generated simulated data for TEPI, STEIN and RWR designs
# Date: 13th March 2026
# Last date editted: 16th June 2026
# **************************************************************************** #

# This script implements a simulation-based comparison of Phase I/II dose-finding
# designs under a unified set of LORDs benchmark scenarios.

# Main features:
#   - Common truth scenarios (LORDs A-D)
#   - Trinomial outcomes:
#       Neutral  = no efficacy, no toxicity
#       Success  = efficacy without toxicity
#       Toxicity = toxicity
#   - Two data-generating mechanisms:
#       base    = logit/logit
#       misspec = probit/probit
#   - Designs compared:
#       Model-assisted: TEPI, STEIN, BOIN12
#   - Final analysis:
#       empirical summary, logit GLM, probit GLM, Firth logistic regression
#   - No early stopping, to keep all designs comparable under N = 60

# Designs:
#   - TEPI
#   - STEIN
#   - RWR (Random Walk Rule active control)
#
# Data-generating mechanisms:
#   - Base case: logit/logit
#   - Misspecification 1: probit/probit
#   - Misspecification 2: quadratic efficacy + logistic toxicity
#
# Outputs:
#   - all_reps_tepi_stein_rwr.csv
#   - all_allocations_tepi_stein_rwr.csv

# 0. Load packages ----

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(logistf)
library(patchwork)
library(here)
library(scales)
library(tibble)

set.seed(2526)

# 1. Settings ----

N_total     <- 60
cohort_size <- 3
startdose   <- 1
nsim        <- 1000

designs <- c("TEPI", "STEIN", "RWR")

design_cols <- c(
  TEPI  = "#56B4E9",
  STEIN = "#009E73",
  RWR   = "#CC79A7"
)

scenario_labels_long <- c(
  A = "A: Narrow therapeutic window",
  B = "B: Wide therapeutic window",
  C = "C: Safe / low toxicity",
  D = "D: Unsafe / high toxicity"
)

scenario_labels_short <- c(
  A = "A\nNarrow",
  B = "B\nWide",
  C = "C\nSafe / low tox",
  D = "D\nUnsafe / high tox"
)

dir.create(here("Results"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("Results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("Results", "figures"), recursive = TRUE, showWarnings = FALSE)

# 2. LORDs dose grid ----

dose_log <- c(-1.20, -0.23, 0.92, 2.02, 3.00, 3.69, 4.38, 5.08, 5.77)
dose_mg  <- c(0.3, 0.8, 2.5, 7.5, 20, 40, 80, 160, 320)

dose_lower <- min(dose_log)
dose_upper <- max(dose_log)

# 3. LORDs scenarios A-D ----

lords_scenarios <- list(
  A = list(
    theta = c(theta1 = 0.855, theta2 = 0.566, theta3 = -5.768, theta4 = 1.000),
    mined_log = 0.92,
    obd_log   = 2.75,
    mtd_log   = 4.38,
    label     = "Scenario A: narrow therapeutic window"
  ),
  B = list(
    theta = c(theta1 = 2.017, theta2 = 2.827, theta3 = -11.537, theta4 = 2.000),
    mined_log = -0.22,
    obd_log   = 2.04,
    mtd_log   = 5.08,
    label     = "Scenario B: wide therapeutic window"
  ),
  C = list(
    theta = c(theta1 = -3.539, theta2 = 1.124, theta3 = -26.618, theta4 = 3.674),
    mined_log = 4.38,
    obd_log   = 6.03,
    mtd_log   = 6.87,
    label     = "Scenario C: safe"
  ),
  D = list(
    theta = c(theta1 = 1.437, theta2 = 0.125, theta3 = -1.525, theta4 = 1.227),
    mined_log = -1.02,
    obd_log   = -1.80,
    mtd_log   = 0.11,
    label     = "Scenario D: unsafe"
  )
)

scenario_names <- names(lords_scenarios)

# 4. Data-generating mechanisms ----
dgp_settings <- list(
  base = list(
    eff_model = "cr_linear",
    eff_link  = "logit",
    tox_model = "linear",
    tox_link  = "logit",
    label     = "Base case: logit/logit"
  ),
  misspec1 = list(
    eff_model = "cr_linear",
    eff_link  = "probit",
    tox_model = "linear",
    tox_link  = "probit",
    label     = "Misspecification 1: probit/probit"
  ),
  misspec2 = list(
    eff_model = "quadratic",
    eff_link  = "logit",
    tox_model = "linear",
    tox_link  = "logit",
    label     = "Misspecification 2: quadratic efficacy + logistic toxicity"
  )
)

dgp_labels_all <- vapply(dgp_settings, `[[`, character(1), "label")

# 5. Helper functions for truth probabilities and patient outcome generation ----

inv_link <- function(eta, link = c("logit", "probit")) {
  link <- match.arg(link)
  if (link == "logit") return(plogis(eta))
  pnorm(eta)
}

pE_cond_cr_linear <- function(d, theta, link = "logit") {
  inv_link(theta["theta1"] + theta["theta2"] * d, link = link)
}

pE_cond_quadratic <- function(d, theta, dU = max(dose_log)) {
  dmax <- (qlogis(0.95) - theta["theta1"]) / theta["theta2"]
  
  p_left <- plogis(theta["theta1"] + theta["theta2"] * d)
  p_right <- ((0.8 - 0.95) * (d - dmax)^2) / ((dU - dmax)^2) + 0.95
  
  out <- ifelse(d <= dmax, p_left, p_right)
  pmin(pmax(out, 0), 1)
}

pT_model <- function(d, theta, tox_model = "linear", tox_link = "logit") {
  inv_link(theta["theta3"] + theta["theta4"] * d, link = tox_link)
}

truth_probs <- function(d, theta, dgp) {
  
  pE_cond <- switch(
    dgp$eff_model,
    cr_linear = pE_cond_cr_linear(d, theta, link = dgp$eff_link),
    quadratic = pE_cond_quadratic(d, theta, dU = max(dose_log))
  )
  
  pT <- pT_model(d, theta, tox_model = dgp$tox_model, tox_link = dgp$tox_link)
  
  tibble(
    dose_log = d,
    pE_cond  = pE_cond,
    pT       = pT,
    neutral  = (1 - pE_cond) * (1 - pT),
    success  = pE_cond * (1 - pT),
    toxicity = pT
  )
}

simulate_one_patient <- function(d, theta, dgp) {
  
  probs <- truth_probs(d = d, theta = theta, dgp = dgp)
  
  pE_cond <- probs$pE_cond
  pT      <- probs$pT
  
  tox <- rbinom(1, 1, pT)
  
  eff <- if (tox == 1) {
    0
  } else {
    rbinom(1, 1, pE_cond)
  }
  
  outcome <- if (tox == 1) {
    "toxicity"
  } else if (eff == 1) {
    "success"
  } else {
    "neutral"
  }
  
  tibble(
    dose_log = d,
    tox = tox,
    eff = eff,
    outcome = outcome
  )
}


# ------------------------------------------------------------------------------
# 6. Truth tables
# ------------------------------------------------------------------------------

make_truth_table <- function(scenario_name, dgp_name) {
  sc  <- lords_scenarios[[scenario_name]]
  dgp <- dgp_settings[[dgp_name]]
  
  truth_probs(
    d = dose_log,
    theta = sc$theta,
    dgp = dgp
  ) %>%
    mutate(
      scenario = scenario_name,
      dgp = dgp$label
    )
}

truth_tables_all <- map_dfr(
  names(dgp_settings),
  function(dgp_name) {
    map_dfr(scenario_names, ~ make_truth_table(.x, dgp_name))
  }
)

write.csv(
  truth_tables_all,
  here("Results", "tables", "truth_tables_all_tepi_stein_rwr.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7. Dose-level summaries and utility helpers
# ------------------------------------------------------------------------------

dose_stats <- function(trial, d_grid) {
  trial %>%
    group_by(dose_log) %>%
    summarise(
      n_treated = n(),
      tox = sum(tox),
      eff = sum(eff),
      .groups = "drop"
    ) %>%
    right_join(tibble(dose_log = d_grid), by = "dose_log") %>%
    arrange(match(dose_log, d_grid)) %>%
    mutate(
      n_treated = replace_na(n_treated, 0L),
      tox = replace_na(tox, 0L),
      eff = replace_na(eff, 0L)
    )
}

posterior_beta_mean <- function(x, n, a = 1, b = 1) {
  (a + x) / (a + b + n)
}

safe_nearest_idx <- function(x, grid) {
  if (length(x) == 0 || all(is.na(x))) return(NA_integer_)
  which.min(abs(grid - x))
}


# ------------------------------------------------------------------------------
# 8. Design decision functions
# ------------------------------------------------------------------------------

# ==============================================================================
# TEPI dose-assignment function
# ==============================================================================

choose_next_dose_tepi <- function(trial, current_idx, d_grid,
                                  tox_intervals = list(
                                    c(0.00, 0.15),
                                    c(0.15, 0.33),
                                    c(0.33, 0.40),
                                    c(0.40, 1.00)
                                  ),
                                  eff_intervals = list(
                                    c(0.00, 0.20),
                                    c(0.20, 0.40),
                                    c(0.40, 0.60),
                                    c(0.60, 1.00)
                                  ),
                                  tox_prior = c(1, 1),
                                  eff_prior = c(1, 1)) {
  
  stats <- dose_stats(trial, d_grid)
  row <- stats[current_idx, ]
  
  n_treated <- row$n_treated
  n_tox     <- row$tox
  n_eff     <- row$eff
  
  a_tox_post <- tox_prior[1] + n_tox
  b_tox_post <- tox_prior[2] + n_treated - n_tox
  
  a_eff_post <- eff_prior[1] + n_eff
  b_eff_post <- eff_prior[2] + n_treated - n_eff
  
  decision_table <- matrix(
    c(
      "E",   "E",   "E",   "E",
      "E",   "E",   "E",   "S",
      "D",   "S",   "S",   "S",
      "D",   "D",   "D",   "D"
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("low", "moderate", "high", "unacceptable"),
      c("low", "moderate", "high", "superb")
    )
  )
  
  tox_labels <- c("low", "moderate", "high", "unacceptable")
  eff_labels <- c("low", "moderate", "high", "superb")
  
  jupm_results <- expand.grid(
    tox_idx = seq_along(tox_intervals),
    eff_idx = seq_along(eff_intervals)
  )
  
  jupm_results$tox_label <- tox_labels[jupm_results$tox_idx]
  jupm_results$eff_label <- eff_labels[jupm_results$eff_idx]
  jupm_results$action_raw <- NA_character_
  jupm_results$cell_prob <- NA_real_
  jupm_results$cell_area <- NA_real_
  jupm_results$jupm <- NA_real_
  
  for (k in seq_len(nrow(jupm_results))) {
    
    tox_bounds <- tox_intervals[[jupm_results$tox_idx[k]]]
    eff_bounds <- eff_intervals[[jupm_results$eff_idx[k]]]
    
    tox_lower <- tox_bounds[1]
    tox_upper <- tox_bounds[2]
    eff_lower <- eff_bounds[1]
    eff_upper <- eff_bounds[2]
    
    p_tox_interval <- pbeta(tox_upper, a_tox_post, b_tox_post) -
      pbeta(tox_lower, a_tox_post, b_tox_post)
    
    p_eff_interval <- pbeta(eff_upper, a_eff_post, b_eff_post) -
      pbeta(eff_lower, a_eff_post, b_eff_post)
    
    p_cell <- p_tox_interval * p_eff_interval
    area_cell <- (tox_upper - tox_lower) * (eff_upper - eff_lower)
    jupm_value <- p_cell / area_cell
    
    jupm_results$action_raw[k] <- decision_table[
      jupm_results$tox_label[k],
      jupm_results$eff_label[k]
    ]
    jupm_results$cell_prob[k] <- p_cell
    jupm_results$cell_area[k] <- area_cell
    jupm_results$jupm[k] <- jupm_value
  }
  
  winner_idx <- which.max(jupm_results$jupm)
  winner <- jupm_results[winner_idx, ]
  
  action_raw <- winner$action_raw
  
  action <- dplyr::case_when(
    action_raw %in% c("E", "EU") ~ "E",
    action_raw %in% c("D", "DUE", "DUT") ~ "D",
    action_raw == "S" ~ "S",
    TRUE ~ "S"
  )
  
  next_idx <- current_idx
  
  if (action == "E" && current_idx < length(d_grid)) {
    next_idx <- current_idx + 1
  }
  
  if (action == "D" && current_idx > 1) {
    next_idx <- current_idx - 1
  }
  
  next_idx <- max(1, min(next_idx, length(d_grid)))
  
  list(
    next_idx = next_idx,
    action = action,
    action_raw = action_raw,
    current_idx = current_idx,
    tox_posterior_mean = a_tox_post / (a_tox_post + b_tox_post),
    eff_posterior_mean = a_eff_post / (a_eff_post + b_eff_post),
    winning_cell = winner,
    jupm_table = jupm_results
  )
}

# ==============================================================================
# STEIN dose-assignment function
# ==============================================================================

choose_next_dose_stein <- function(trial, current_idx, d_grid,
                                   phi0 = 0.20,
                                   phi1 = 0.75 * phi0,
                                   phi2 = 1.25 * phi0,
                                   psi1 = 0.20,
                                   psi2 = 0.60) {
  
  stats <- dose_stats(trial, d_grid)
  row <- stats[current_idx, ]
  
  n_j <- row$n_treated
  x_j <- row$tox
  y_j <- row$eff
  
  if (n_j == 0) {
    return(list(
      next_idx = current_idx,
      action = "S",
      current_idx = current_idx,
      phi_L = NA_real_,
      phi_U = NA_real_,
      psi = NA_real_,
      p_hat = NA_real_,
      q_hat = NA_real_,
      admissible_set = current_idx,
      post_prob_table = NULL
    ))
  }
  
  p_hat_j <- x_j / n_j
  q_hat_j <- y_j / n_j
  
  phi_L <- log((1 - phi1) / (1 - phi0)) /
    log((phi0 * (1 - phi1)) / (phi1 * (1 - phi0)))
  
  phi_U <- log((1 - phi0) / (1 - phi2)) /
    log((phi2 * (1 - phi0)) / (phi0 * (1 - phi2)))
  
  psi <- log((1 - psi1) / (1 - psi2)) /
    log((psi2 * (1 - psi1)) / (psi1 * (1 - psi2)))
  
  if (p_hat_j <= phi_L) {
    A_j <- c(current_idx - 1, current_idx, current_idx + 1)
  } else if (p_hat_j < phi_U) {
    A_j <- c(current_idx - 1, current_idx)
  } else {
    A_j <- c(current_idx - 1)
  }
  
  A_j <- A_j[A_j >= 1 & A_j <= length(d_grid)]
  A_j <- sort(unique(A_j))
  
  if (p_hat_j >= phi_U) {
    next_idx <- max(1, current_idx - 1)
    
    return(list(
      next_idx = next_idx,
      action = "D",
      current_idx = current_idx,
      phi_L = phi_L,
      phi_U = phi_U,
      psi = psi,
      p_hat = p_hat_j,
      q_hat = q_hat_j,
      admissible_set = A_j,
      post_prob_table = NULL
    ))
  }
  
  if (p_hat_j < phi_U && q_hat_j >= psi) {
    return(list(
      next_idx = current_idx,
      action = "S",
      current_idx = current_idx,
      phi_L = phi_L,
      phi_U = phi_U,
      psi = psi,
      p_hat = p_hat_j,
      q_hat = q_hat_j,
      admissible_set = A_j,
      post_prob_table = NULL
    ))
  }
  
  post_prob_table <- data.frame(
    dose_idx = A_j,
    n_treated = NA_integer_,
    n_eff = NA_integer_,
    post_prob_eff_above_psi = NA_real_
  )
  
  for (k in seq_along(A_j)) {
    jprime <- A_j[k]
    row_jprime <- stats[jprime, ]
    
    n_jprime <- row_jprime$n_treated
    y_jprime <- row_jprime$eff
    
    a_post <- 1 + y_jprime
    b_post <- 1 + n_jprime - y_jprime
    post_prob <- 1 - pbeta(psi, a_post, b_post)
    
    post_prob_table$n_treated[k] <- n_jprime
    post_prob_table$n_eff[k] <- y_jprime
    post_prob_table$post_prob_eff_above_psi[k] <- post_prob
  }
  
  max_prob <- max(post_prob_table$post_prob_eff_above_psi)
  candidate_rows <- post_prob_table[
    post_prob_table$post_prob_eff_above_psi == max_prob, ,
    drop = FALSE
  ]
  
  next_idx <- min(candidate_rows$dose_idx)
  
  action <- if (next_idx > current_idx) {
    "E"
  } else if (next_idx < current_idx) {
    "D"
  } else {
    "S"
  }
  
  list(
    next_idx = next_idx,
    action = action,
    current_idx = current_idx,
    phi_L = phi_L,
    phi_U = phi_U,
    psi = psi,
    p_hat = p_hat_j,
    q_hat = q_hat_j,
    admissible_set = A_j,
    post_prob_table = post_prob_table
  )
}

# ==============================================================================
# Random Walk Rule (RWR) based on Ivanova (2003)
# ==============================================================================

choose_next_dose_rwr <- function(trial, current_idx, d_grid,
                                 gamma_tox = 0.20,
                                 cohort_size = 3,
                                 tox_prior = c(1, 1)) {
  
  stats <- dose_stats(trial, d_grid)
  
  current_dose <- d_grid[current_idx]
  
  recent_group <- trial %>%
    filter(dose_log == current_dose) %>%
    tail(cohort_size)
  
  if (nrow(recent_group) == 0) {
    return(list(
      next_idx = current_idx,
      action = "S",
      current_idx = current_idx,
      safe_idx = current_idx,
      recent_group_summary = NULL
    ))
  }
  
  n_tox <- sum(recent_group$tox == 1, na.rm = TRUE)
  n_success <- sum(recent_group$tox == 0 & recent_group$eff == 1, na.rm = TRUE)
  n_neutral <- sum(recent_group$tox == 0 & recent_group$eff == 0, na.rm = TRUE)
  
  proposed_idx <- current_idx
  action <- "S"
  
  if (n_tox >= 2) {
    proposed_idx <- max(1, current_idx - 1)
    action <- "D"
  } else if (n_neutral >= 2) {
    proposed_idx <- min(length(d_grid), current_idx + 1)
    action <- "E"
  } else {
    proposed_idx <- current_idx
    action <- "S"
  }
  
  tox_est <- stats %>%
    mutate(
      ptox_post_mean = (tox_prior[1] + tox) / (tox_prior[1] + tox_prior[2] + n_treated)
    )
  
  safe_doses <- which(tox_est$ptox_post_mean <= gamma_tox)
  
  if (length(safe_doses) == 0) {
    safe_idx <- 1
  } else {
    safe_idx <- max(safe_doses)
  }
  
  next_idx <- min(proposed_idx, safe_idx)
  
  if (next_idx < current_idx) {
    action <- "D"
  } else if (next_idx > current_idx) {
    action <- "E"
  } else {
    action <- "S"
  }
  
  list(
    next_idx = next_idx,
    action = action,
    current_idx = current_idx,
    proposed_idx = proposed_idx,
    safe_idx = safe_idx,
    gamma_tox = gamma_tox,
    recent_group_summary = tibble::tibble(
      n_tox = n_tox,
      n_success = n_success,
      n_neutral = n_neutral
    ),
    toxicity_estimates = tox_est
  )
}


# ------------------------------------------------------------------------------
# 9. Sequential adaptive trial simulation
# ------------------------------------------------------------------------------

simulate_adaptive_trial <- function(design,
                                    scenario_name,
                                    dgp_name = c("base", "misspec1", "misspec2"),
                                    n_total = 60,
                                    cohort_size = 3,
                                    startdose = 1,
                                    seed = NULL) {
  
  dgp_name <- match.arg(dgp_name)
  
  if (!is.null(seed)) set.seed(seed)
  
  sc  <- lords_scenarios[[scenario_name]]
  dgp <- dgp_settings[[dgp_name]]
  
  current_idx <- startdose
  trial <- tibble()
  available <- rep(TRUE, length(dose_log))
  
  while (nrow(trial) < n_total) {
    
    n_this <- min(cohort_size, n_total - nrow(trial))
    
    current_dose_log <- dose_log[current_idx]
    current_dose_mg  <- dose_mg[current_idx]
    patient_ids <- nrow(trial) + seq_len(n_this)
    
    cohort <- map_dfr(seq_len(n_this), ~ simulate_one_patient(
      d = current_dose_log,
      theta = sc$theta,
      dgp = dgp
    )) %>%
      mutate(
        design = design,
        scenario = scenario_name,
        dgp = dgp$label,
        patient_id = patient_ids,
        cohort_id = ceiling(patient_id / cohort_size),
        dose_idx = current_idx,
        dose_log = current_dose_log,
        dose_mg = current_dose_mg
      )
    
    trial <- bind_rows(trial, cohort)
    
    upd <- switch(
      design,
      TEPI = choose_next_dose_tepi(
        trial = trial,
        current_idx = current_idx,
        d_grid = dose_log
      ),
      STEIN = choose_next_dose_stein(
        trial = trial,
        current_idx = current_idx,
        d_grid = dose_log
      ),
      RWR = choose_next_dose_rwr(
        trial = trial,
        current_idx = current_idx,
        d_grid = dose_log,
        cohort_size = cohort_size
      )
    )
    
    current_idx <- upd$next_idx
  }
  
  stats <- dose_stats(trial, dose_log) %>%
    mutate(
      p_success_emp = ifelse(n_treated > 0, eff / n_treated, NA_real_),
      p_tox_emp     = ifelse(n_treated > 0, tox / n_treated, NA_real_)
    )
  
  final_idx <- current_idx
  
  cand <- which(!is.na(stats$p_success_emp) & stats$p_tox_emp <= 0.20)
  if (length(cand) == 0) cand <- which(!is.na(stats$p_success_emp))
  
  if (length(cand) > 0) {
    final_idx <- cand[which.max(stats$p_success_emp[cand])]
  }
  
  list(
    trial = trial,
    final_idx = final_idx,
    final_dose_log = dose_log[final_idx],
    n_treated = stats$n_treated,
    tox_by_dose = stats$tox,
    eff_by_dose = stats$eff
  )
}


# ------------------------------------------------------------------------------
# 10. Final OBD analysis helpers
# ------------------------------------------------------------------------------

empirical_summary <- function(trial, d_grid) {
  trial %>%
    group_by(dose_log) %>%
    summarise(
      p_success_emp = mean(outcome == "success"),
      p_tox_emp     = mean(outcome == "toxicity"),
      n_treated     = n(),
      .groups = "drop"
    ) %>%
    right_join(tibble(dose_log = d_grid), by = "dose_log") %>%
    arrange(match(dose_log, d_grid))
}

fit_pair <- function(trial, link = c("logit", "probit"), method = c("glm", "firth")) {
  
  link <- match.arg(link)
  method <- match.arg(method)
  
  dat_eff <- subset(trial, tox == 0)
  
  if (nrow(dat_eff) < 2 || nrow(trial) < 4) {
    return(list(tox_fit = NULL, eff_fit = NULL))
  }
  
  if (method == "glm") {
    tox_fit <- try(
      glm(tox ~ dose_log, data = trial, family = binomial(link = link)),
      silent = TRUE
    )
    
    eff_fit <- try(
      glm(eff ~ dose_log, data = dat_eff, family = binomial(link = link)),
      silent = TRUE
    )
  } else {
    if (link != "logit") {
      stop("Firth estimation is implemented here for logistic regression only.")
    }
    
    tox_fit <- try(logistf::logistf(tox ~ dose_log, data = trial), silent = TRUE)
    eff_fit <- try(logistf::logistf(eff ~ dose_log, data = dat_eff), silent = TRUE)
  }
  
  if (inherits(tox_fit, "try-error")) tox_fit <- NULL
  if (inherits(eff_fit, "try-error")) eff_fit <- NULL
  
  list(tox_fit = tox_fit, eff_fit = eff_fit)
}

predict_response_safe <- function(fit, newdat) {
  if (is.null(fit)) return(rep(NA_real_, nrow(newdat)))
  
  out <- try(
    predict(fit, newdata = newdat, type = "response"),
    silent = TRUE
  )
  
  if (inherits(out, "try-error")) {
    eta <- as.numeric(cbind(1, newdat$dose_log) %*% coef(fit))
    out <- plogis(eta)
  }
  
  as.numeric(out)
}

predict_success_curve <- function(fits, d_grid, label) {
  newdat <- tibble(dose_log = d_grid)
  
  pT_hat <- predict_response_safe(fits$tox_fit, newdat)
  pE_hat <- predict_response_safe(fits$eff_fit, newdat)
  
  tibble(
    dose_log = d_grid,
    success_hat = (1 - pT_hat) * pE_hat,
    method = label
  )
}

recommend_dose <- function(curve_df) {
  if (all(is.na(curve_df$success_hat))) return(NA_real_)
  curve_df$dose_log[which.max(curve_df$success_hat)]
}


# ------------------------------------------------------------------------------
# 11. Simulate one replicate
# ------------------------------------------------------------------------------

simulate_one_replicate <- function(design,
                                   scenario_name,
                                   dgp_name,
                                   n_total = 60,
                                   cohort_size = 3,
                                   startdose = 1,
                                   seed = NULL) {
  
  sc  <- lords_scenarios[[scenario_name]]
  dgp <- dgp_settings[[dgp_name]]
  
  trial_obj <- simulate_adaptive_trial(
    design = design,
    scenario_name = scenario_name,
    dgp_name = dgp_name,
    n_total = n_total,
    cohort_size = cohort_size,
    startdose = startdose,
    seed = seed
  )
  
  trial <- trial_obj$trial
  
  true_obd_log <- sc$obd_log
  true_obd_idx <- which.min(abs(dose_log - true_obd_log))
  
  selected_idx <- which.min(abs(dose_log - trial_obj$final_dose_log))
  selected_log <- trial_obj$final_dose_log
  
  emp <- empirical_summary(trial, dose_log)
  
  fit_logit  <- fit_pair(trial, link = "logit",  method = "glm")
  fit_probit <- fit_pair(trial, link = "probit", method = "glm")
  fit_firth  <- fit_pair(trial, link = "logit",  method = "firth")
  
  curve_true <- truth_probs(dose_log, sc$theta, dgp) %>%
    transmute(dose_log, success_hat = success, method = "Truth")
  
  curve_emp <- emp %>%
    transmute(dose_log, success_hat = p_success_emp, method = "Empirical")
  
  curve_logit  <- predict_success_curve(fit_logit, dose_log, "Logit GLM")
  curve_probit <- predict_success_curve(fit_probit, dose_log, "Probit GLM")
  curve_firth  <- predict_success_curve(fit_firth, dose_log, "Firth logit")
  
  curves <- bind_rows(
    curve_true,
    curve_emp,
    curve_logit,
    curve_probit,
    curve_firth
  ) %>%
    mutate(
      design = design,
      scenario = scenario_name,
      dgp = dgp$label
    )
  
  emp_rec    <- recommend_dose(curve_emp)
  logit_rec  <- recommend_dose(curve_logit)
  probit_rec <- recommend_dose(curve_probit)
  firth_rec  <- recommend_dose(curve_firth)
  
  summary_out <- tibble(
    design = design,
    scenario = scenario_name,
    dgp = dgp$label,
    
    true_obd_log = true_obd_log,
    true_obd_idx = true_obd_idx,
    
    selected_log = selected_log,
    selected_idx = selected_idx,
    
    bias_log = selected_log - true_obd_log,
    bias_idx = selected_idx - true_obd_idx,
    
    below_true = as.integer(selected_idx < true_obd_idx),
    exact_true = as.integer(selected_idx == true_obd_idx),
    within_one = as.integer(abs(selected_idx - true_obd_idx) <= 1),
    
    n_patients = nrow(trial),
    n_tox = sum(trial$tox),
    n_eff = sum(trial$eff),
    n_neutral = sum(trial$outcome == "neutral"),
    
    emp_rec_idx    = safe_nearest_idx(emp_rec, dose_log),
    logit_rec_idx  = safe_nearest_idx(logit_rec, dose_log),
    probit_rec_idx = safe_nearest_idx(probit_rec, dose_log),
    firth_rec_idx  = safe_nearest_idx(firth_rec, dose_log)
  )
  
  list(
    summary = summary_out,
    curves = curves,
    trial = trial
  )
}


# ------------------------------------------------------------------------------
# 12. Run many replicates
# ------------------------------------------------------------------------------

run_many_replicates <- function(design,
                                scenario_name,
                                dgp_name,
                                n_rep = 1000,
                                n_total = 60,
                                cohort_size = 3,
                                startdose = 1) {
  
  sim_list <- map(seq_len(n_rep), function(r) {
    seed_r <- 1000 + r
    
    sim <- simulate_one_replicate(
      design = design,
      scenario_name = scenario_name,
      dgp_name = dgp_name,
      n_total = n_total,
      cohort_size = cohort_size,
      startdose = startdose,
      seed = seed_r
    )
    
    summary_r <- sim$summary %>% mutate(rep = r)
    allocation_r <- sim$trial %>% mutate(rep = r)
    
    list(
      summary = summary_r,
      allocation = allocation_r
    )
  })
  
  list(
    summary = map_dfr(sim_list, "summary"),
    allocation = map_dfr(sim_list, "allocation")
  )
}


# ------------------------------------------------------------------------------
# 13. Run TEPI, STEIN, and RWR across all scenarios and DGPs
# ------------------------------------------------------------------------------

all_results_list <- list()
counter <- 1

for (des in designs) {
  for (sc in scenario_names) {
    for (dgp_name in names(dgp_settings)) {
      
      message(
        "Running design = ", des,
        ", scenario = ", sc,
        ", DGP = ", dgp_name
      )
      
      all_results_list[[counter]] <- run_many_replicates(
        design = des,
        scenario_name = sc,
        dgp_name = dgp_name,
        n_rep = nsim,
        n_total = N_total,
        cohort_size = cohort_size,
        startdose = startdose
      )
      
      counter <- counter + 1
    }
  }
}

all_reps <- map_dfr(all_results_list, "summary")
all_allocations <- map_dfr(all_results_list, "allocation")


# ------------------------------------------------------------------------------
# 14. Save outputs
# ------------------------------------------------------------------------------

write.csv(
  all_reps,
  here("Results", "tables", "all_reps_tepi_stein_rwr.csv"),
  row.names = FALSE
)

write.csv(
  all_allocations,
  here("Results", "tables", "all_allocations_tepi_stein_rwr.csv"),
  row.names = FALSE
)