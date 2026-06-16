# ==============================================================================
# BOIN12-ONLY SIMULATION SCRIPT
# ==============================================================================
# This script runs BOIN12 only, separately from TEPI, STEIN, and RWR.
# It is adapted from the full simulation framework.
#
# Main adaptation:
#   - if the BOIN12 RDS lookup table does not contain the exact n_treated value,
#     the lookup falls back to the nearest lower available n_treated block.
#
# Outputs:
#   - all_reps_boin12.csv
#   - all_allocations_boin12.csv
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(logistf)
library(patchwork)
library(here)
library(scales)
library(tibble)
library(stringr)

set.seed(2526)


# ------------------------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------------------------

N_total     <- 60
cohort_size <- 3
startdose   <- 1
nsim        <- 1000

designs <- c("BOIN12")

design_cols <- c(
  BOIN12 = "#E69F00"
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


# ------------------------------------------------------------------------------
# 2. LORDs dose grid
# ------------------------------------------------------------------------------

dose_log <- c(-1.20, -0.23, 0.92, 2.02, 3.00, 3.69, 4.38, 5.08, 5.77)
dose_mg  <- c(0.3, 0.8, 2.5, 7.5, 20, 40, 80, 160, 320)

dose_lower <- min(dose_log)
dose_upper <- max(dose_log)


# ------------------------------------------------------------------------------
# 3. LORDs scenarios A-D
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# 4. Data-generating mechanisms
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# 5. Helper functions for truth probabilities and patient outcome generation
# ------------------------------------------------------------------------------

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
  here("Results", "tables", "truth_tables_all_boin12.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7. Dose-level summaries and utility helpers
# ------------------------------------------------------------------------------

dose_stats <- function(trial, d_grid) {
  tibble(dose_log = d_grid) %>%
    left_join(
      trial %>%
        group_by(dose_log) %>%
        summarise(
          n_treated = n(),
          tox = sum(tox, na.rm = TRUE),
          eff = sum(eff, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "dose_log"
    ) %>%
    arrange(match(dose_log, d_grid)) %>%
    mutate(
      n_treated = ifelse(is.na(n_treated), 0L, n_treated),
      tox = ifelse(is.na(tox), 0L, tox),
      eff = ifelse(is.na(eff), 0L, eff)
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
# 8. BOIN escalation/de-escalation boundaries
# ------------------------------------------------------------------------------

boin_boundary_table <- tibble(
  phi_T    = c(0.20, 0.25, 0.30, 0.35, 0.40),
  lambda_e = c(0.157, 0.197, 0.236, 0.276, 0.316),
  lambda_d = c(0.238, 0.298, 0.359, 0.419, 0.480)
)

get_boin_boundaries <- function(phi_T) {
  hit <- which(abs(boin_boundary_table$phi_T - phi_T) < 1e-8)
  
  if (length(hit) != 1) {
    stop("phi_T not found in boundary table. Please supply lambda_e and lambda_d manually.")
  }
  
  list(
    lambda_e = boin_boundary_table$lambda_e[hit],
    lambda_d = boin_boundary_table$lambda_d[hit]
  )
}


# ------------------------------------------------------------------------------
# 9. RDS lookup object
# ------------------------------------------------------------------------------
# This is the object from your current script.
# IMPORTANT:
# - It only covers some n_treated values.
# - The lookup function below includes a fallback to the nearest lower available
#   n_treated value so the simulation can run today.
# ------------------------------------------------------------------------------

rds_lookup_boin12_full <- rio::import(here("BOIN12_v1.4.2.0_EscalationDe-escalation Boundaries.xlsx"))

# ------------------------------------------------------------------------------
# 10. Helper to interpret rule strings
# ------------------------------------------------------------------------------

match_rule_value <- function(x, rule) {
  rule <- trimws(rule)
  
  if (is.na(rule)) return(FALSE)
  if (rule == "Any") return(TRUE)
  
  if (str_detect(rule, "^>=\\s*\\d+$")) {
    k <- as.integer(str_extract(rule, "\\d+"))
    return(x >= k)
  }
  
  if (str_detect(rule, "^<=\\s*\\d+$")) {
    k <- as.integer(str_extract(rule, "\\d+"))
    return(x <= k)
  }
  
  if (str_detect(rule, "^\\d+$")) {
    return(x == as.integer(rule))
  }
  
  FALSE
}


# ------------------------------------------------------------------------------
# 11. Lookup function with FALLBACK
# ------------------------------------------------------------------------------
# If exact n_treated is missing, use the nearest lower available n_treated block.
# This is the key workaround so BOIN12 can run today.
# ------------------------------------------------------------------------------

lookup_boin12_rds_export <- function(n_treated, n_tox_obs, n_eff_obs, rds_lookup) {
  
  sub <- rds_lookup %>%
    filter(n_treated == !!n_treated)
  
  if (nrow(sub) == 0) {
    available_n <- sort(unique(rds_lookup$n_treated))
    fallback_n <- max(available_n[available_n <= n_treated], na.rm = TRUE)
    
    if (!is.finite(fallback_n)) {
      stop(
        paste0(
          "No RDS rows found for n_treated = ", n_treated,
          " and no smaller fallback value is available."
        )
      )
    }
    
    message(
      "BOIN12 lookup fallback: using n_treated = ",
      fallback_n, " instead of ", n_treated
    )
    
    sub <- rds_lookup %>%
      filter(n_treated == fallback_n)
  }
  
  sub <- sub %>%
    rowwise() %>%
    mutate(
      tox_match = match_rule_value(n_tox_obs, n_tox),
      eff_match = match_rule_value(n_eff_obs, n_eff)
    ) %>%
    ungroup() %>%
    filter(tox_match, eff_match)
  
  if (nrow(sub) == 0) {
    stop(
      paste0(
        "RDS lookup failed for n_treated=", n_treated,
        ", n_tox=", n_tox_obs,
        ", n_eff=", n_eff_obs
      )
    )
  }
  
  specificity_score <- function(rule) {
    if (rule == "Any") return(0)
    if (str_detect(rule, "^\\d+$")) return(2)
    return(1)
  }
  
  sub <- sub %>%
    mutate(
      tox_spec = vapply(n_tox, specificity_score, numeric(1)),
      eff_spec = vapply(n_eff, specificity_score, numeric(1)),
      total_spec = tox_spec + eff_spec
    ) %>%
    arrange(desc(total_spec))
  
  best <- sub[1, ]
  
  score_num <- ifelse(best$desirability_score == "E", -Inf, as.numeric(best$desirability_score))
  
  list(
    score = score_num,
    matched_tox_rule = best$n_tox,
    matched_eff_rule = best$n_eff
  )
}


# ------------------------------------------------------------------------------
# 12. BOIN12 allocation function
# ------------------------------------------------------------------------------

choose_next_dose_boin12 <- function(trial, current_idx, d_grid,
                                    phi_T = 0.35,
                                    phi_E = 0.25,
                                    lambda_e = NULL,
                                    lambda_d = NULL,
                                    rds_lookup = rds_lookup_boin12_full) {
  
  if (is.null(lambda_e) || is.null(lambda_d)) {
    bounds <- get_boin_boundaries(phi_T)
    lambda_e <- bounds$lambda_e
    lambda_d <- bounds$lambda_d
  }
  
  stats <- dose_stats(trial, d_grid)
  row <- stats[current_idx, ]
  
  n_cur   <- row$n_treated
  tox_cur <- row$tox
  
  p_hat_T_cur <- if (n_cur > 0) tox_cur / n_cur else 0
  
  if (n_cur == 0) {
    return(list(
      next_idx = current_idx,
      action = "S",
      current_idx = current_idx,
      phi_T = phi_T,
      phi_E = phi_E,
      lambda_e = lambda_e,
      lambda_d = lambda_d,
      p_hat_T = NA_real_,
      compared_set = current_idx,
      desirability_table = NULL
    ))
  }
  
  if (p_hat_T_cur >= lambda_d) {
    next_idx <- max(1, current_idx - 1)
    
    return(list(
      next_idx = next_idx,
      action = "D",
      current_idx = current_idx,
      phi_T = phi_T,
      phi_E = phi_E,
      lambda_e = lambda_e,
      lambda_d = lambda_d,
      p_hat_T = p_hat_T_cur,
      compared_set = c(max(1, current_idx - 1)),
      desirability_table = NULL
    ))
  }
  
  if (p_hat_T_cur > lambda_e && p_hat_T_cur < lambda_d && n_cur >= 6) {
    candidate_set <- c(current_idx - 1, current_idx)
  } else {
    candidate_set <- c(current_idx - 1, current_idx, current_idx + 1)
  }
  
  candidate_set <- candidate_set[candidate_set >= 1 & candidate_set <= length(d_grid)]
  candidate_set <- sort(unique(candidate_set))
  
  desirability_table <- tibble(
    dose_idx = candidate_set,
    n_treated = NA_integer_,
    n_tox_obs = NA_integer_,
    n_eff_obs = NA_integer_,
    desirability_score = NA_real_,
    matched_tox_rule = NA_character_,
    matched_eff_rule = NA_character_
  )
  
  for (k in seq_along(candidate_set)) {
    j <- candidate_set[k]
    row_j <- stats[j, ]
    
    n_j   <- as.integer(row_j$n_treated)
    tox_j <- as.integer(row_j$tox)
    eff_j <- as.integer(row_j$eff)
    
    desirability_table$n_treated[k] <- n_j
    desirability_table$n_tox_obs[k] <- tox_j
    desirability_table$n_eff_obs[k] <- eff_j
    
    lookup_res <- lookup_boin12_rds_export(
      n_treated = n_j,
      n_tox_obs = tox_j,
      n_eff_obs = eff_j,
      rds_lookup = rds_lookup
    )
    
    desirability_table$desirability_score[k] <- lookup_res$score
    desirability_table$matched_tox_rule[k] <- lookup_res$matched_tox_rule
    desirability_table$matched_eff_rule[k] <- lookup_res$matched_eff_rule
  }
  
  max_score <- max(desirability_table$desirability_score)
  winners <- desirability_table %>%
    filter(desirability_score == max_score)
  
  next_idx <- min(winners$dose_idx)
  
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
    phi_T = phi_T,
    phi_E = phi_E,
    lambda_e = lambda_e,
    lambda_d = lambda_d,
    p_hat_T = p_hat_T_cur,
    compared_set = candidate_set,
    desirability_table = desirability_table
  )
}


# ------------------------------------------------------------------------------
# 13. Sequential adaptive trial simulation
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
    
    upd <- choose_next_dose_boin12(
      trial = trial,
      current_idx = current_idx,
      d_grid = dose_log
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
# 14. Final OBD analysis helpers
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
# 15. Simulate one replicate
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
# 16. Run many replicates
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
# 17. Run BOIN12 across all scenarios and DGPs
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

all_reps_b <- map_dfr(all_results_list, "summary")
all_allocations_b <- map_dfr(all_results_list, "allocation")


# ------------------------------------------------------------------------------
# 18. Save outputs
# ------------------------------------------------------------------------------

write.csv(
  all_reps_b,
  here("Results", "tables", "all_reps_boin12.csv"),
  row.names = FALSE
)

write.csv(
  all_allocations_b,
  here("Results", "tables", "all_allocations_boin12.csv"),
  row.names = FALSE
)