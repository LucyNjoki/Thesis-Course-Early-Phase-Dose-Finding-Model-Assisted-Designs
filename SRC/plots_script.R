# Load Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(scales)
library(ggrepel)
library(purrr)
library(tibble)
library(readr)
library(stringr)
library(purrr)
library(tibble)
library(logistf)


# Settings ----
designs <- c("TEPI", "STEIN", "BOIN12", "RWR")

design_cols <- c(
  TEPI   = "#56B4E9",
  BOIN12 = "#E69F00",
  STEIN  = "#009E73",
  RWR    = "#CC79A7"
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

theme_thesis <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", colour = "grey80"),
    plot.margin = margin(8, 8, 8, 8)
  )

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

# Import Simulated Data ----
all_reps <- readRDS(here("Results","Simulated Data", "all_reps_final.rds"))

all_allocations <- readRDS(here("Results","Simulated Data", "all_allocations_final.rds"))

truth_tables <- read_csv(here("Results", "truth_tables_all_tepi_stein_rwr.csv"))

oc_summary <- all_reps %>%
  group_by(design, scenario, dgp) %>%
  summarise(
    mean_bias_idx = mean(bias_idx, na.rm = TRUE),
    mean_bias_log = mean(bias_log, na.rm = TRUE),
    p_below = mean(below_true, na.rm = TRUE),
    p_exact = mean(exact_true, na.rm = TRUE),
    p_within1 = mean(within_one, na.rm = TRUE),
    mean_selected_idx = mean(selected_idx, na.rm = TRUE),
    mean_selected_log = mean(selected_log, na.rm = TRUE),
    mean_n_tox = mean(n_tox, na.rm = TRUE),
    mean_n_eff = mean(n_eff, na.rm = TRUE),
    mean_n_neutral = mean(n_neutral, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  oc_summary,
  here("Results", "tables", "oc_summary.csv"),
  row.names = FALSE
)


# Bias and OBD recovery plot ----
bias_df <- oc_summary %>%
  select(design, scenario, dgp, mean_bias_idx, p_below, p_exact, p_within1) %>%
  pivot_longer(
    cols = c(mean_bias_idx, p_below, p_exact, p_within1),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = c("A", "B", "C", "D")),
    metric = factor(
      metric,
      levels = c("mean_bias_idx", "p_below", "p_exact", "p_within1"),
      labels = c(
        "Mean bias",
        "Below the\ntrue OBD",
        "Exact OBD\nselection",
        "Within one dose\nlevel of the OBD"
      )
    ),
    dgp = factor(dgp, levels = dgp_labels_all),
    design = factor(design, levels = c("TEPI", "BOIN12", "STEIN", "RWR"))
  ) %>% 
  mutate(
    scenario = str_wrap(scenario, width = 8),
    dgp = str_wrap(dgp, width = 14),
    metric = str_wrap(metric, width = 12)
  )

p_bias <- ggplot(bias_df, aes(x = design, y = value, fill = design)) +
  geom_col(width = 0.55, colour = "white") +
  facet_grid(metric ~ scenario + dgp, scales = "free_y") +
  scale_fill_manual(values = design_cols) +
  labs(
    title = "Bias and Accuracy of OBD Selection Across Designs",
    x = NULL,
    y = NULL,
    fill = "Design"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "top",
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(8, 8, 8, 8)
  )

print(p_bias)

ggsave(
  filename = here("Results", "figures", "bias_and_obd_recovery.png"),
  plot = p_bias,
  width = 18,
  height = 10,
  dpi = 300
)


# Distribution of selected dose relative to true OBD----
selection_df <- all_reps %>%
  mutate(delta = selected_idx - true_obd_idx) %>%
  count(design, scenario, dgp, delta) %>%
  group_by(design, scenario, dgp) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    dgp = factor(dgp, levels = dgp_labels_all),
    scenario = factor(scenario, levels = c("A", "B", "C", "D"), labels = c("A", "B", "C", "D")),
    design = factor(design, levels = c("TEPI", "BOIN12", "STEIN", "RWR"))
  ) %>% 
  mutate(
    scenario = str_wrap(scenario, width = 8),
    dgp = str_wrap(dgp, width = 14)
  )

p_selection <- ggplot(selection_df, aes(x = delta, y = prop, fill = design)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    colour = "white"
  ) +
  facet_grid(dgp ~ scenario) +
  scale_fill_manual(values = design_cols) +
  scale_x_continuous(
    breaks = seq(
      min(selection_df$delta, na.rm = TRUE),
      max(selection_df$delta, na.rm = TRUE),
      by = 1
    )
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Distribution of Selected Dose Relative to the True OBD",
    x = "Selected dose index - true OBD index",
    y = "Proportion of simulations",
    fill = "Design"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 10),
    plot.margin = margin(8, 8, 8, 8)
  )

print(p_selection)

ggsave(
  filename = here("Results", "figures", "selected_dose_relative_to_obd.png"),
  plot = p_selection,
  width = 11,
  height = 7,
  dpi = 300
)



# Probability of selecting within the therapeutic range----

therapeutic_range <- map_dfr(scenario_names, function(sc) {
  tibble(
    scenario = sc,
    mined_log = lords_scenarios[[sc]]$mined_log,
    obd_log   = lords_scenarios[[sc]]$obd_log,
    mtd_log   = lords_scenarios[[sc]]$mtd_log
  )
})

all_reps_therapeutic <- all_reps %>%
  left_join(therapeutic_range, by = "scenario") %>%
  mutate(
    within_therapeutic_range =
      selected_log >= mined_log &
      selected_log <= mtd_log
  )

p_therapeutic <- all_reps_therapeutic %>%
  group_by(design, scenario, dgp) %>%
  summarise(
    n_trials = n(),
    n_within_therapeutic_range = sum(within_therapeutic_range, na.rm = TRUE),
    p_within_therapeutic_range = mean(within_therapeutic_range, na.rm = TRUE),
    p_within_therapeutic_range_percent = 100 * p_within_therapeutic_range,
    .groups = "drop"
  )

write.csv(
  p_therapeutic,
  here("Results", "tables", "probability_within_therapeutic_range.csv"),
  row.names = FALSE
)

within_range_plot_df <- p_therapeutic %>%
  mutate(
    design = factor(design, levels = c("TEPI", "BOIN12", "STEIN", "RWR")),
    scenario = factor(scenario, levels = c("A", "B", "C", "D"), labels = scenario_labels_short),
    dgp = factor(dgp, levels = dgp_labels_all)
  )

p_within_range <- ggplot(
  within_range_plot_df,
  aes(x = scenario, y = p_within_therapeutic_range_percent, fill = design)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, colour = "white") +
  geom_text(
    aes(label = paste0(round(p_within_therapeutic_range_percent, 1), "%")),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 2.8
  ) +
  facet_wrap(~ dgp, ncol = 1) +
  scale_fill_manual(values = design_cols) +
  scale_y_continuous(
    limits = c(0, 105),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "Probability of Selecting Within the Therapeutic Range",
    subtitle = "Percentage of simulated trials selecting a final dose between MinED and MTD",
    x = "Scenario",
    y = "Probability within therapeutic range",
    fill = "Design"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_within_range)

ggsave(
  filename = here("Results", "figures", "probability_within_therapeutic_range.png"),
  plot = p_within_range,
  width = 12,
  height = 9,
  dpi = 300
)



# Mode selected dose table ----
mode_summary <- all_reps %>%
  group_by(design, scenario, dgp) %>%
  summarise(
    mode_selected_log = as.numeric(names(which.max(table(selected_log)))),
    true_obd_log = unique(true_obd_log),
    mode_minus_true_obd = mode_selected_log - true_obd_log,
    n_trials = n(),
    .groups = "drop"
  )

print(mode_summary)

mode_table <- mode_summary %>%
  mutate(
    scenario_dgp = paste(scenario, dgp, sep = " | "),
    mode_selected_log = round(mode_selected_log, 2)
  ) %>%
  select(design, scenario_dgp, mode_selected_log) %>%
  pivot_wider(
    names_from = scenario_dgp,
    values_from = mode_selected_log
  )

print(mode_table)

write.csv(
  mode_table,
  here("Results", "tables", "mode_selected_table.csv"),
  row.names = FALSE
)

# True OBD per scenario ----
true_obd_table <- all_reps %>%
  group_by(scenario) %>%
  summarise(
    true_obd_idx = first(true_obd_idx),
    true_obd_log = first(true_obd_log),
    .groups = "drop"
  ) %>%
  mutate(
    scenario_label = recode(
      scenario,
      A = "A: Narrow therapeutic window",
      B = "B: Wide therapeutic window",
      C = "C: Safe / low toxicity",
      D = "D: Unsafe / high toxicity"
    )
  ) %>%
  select(scenario, scenario_label, true_obd_idx, true_obd_log) %>%
  arrange(scenario)

print(true_obd_table)

write.csv(
  true_obd_table,
  here("Results", "tables", "true_obd_table.csv"),
  row.names = FALSE
)


saveRDS(all_reps, file = here("Results","Simulated Data", "all_reps_final.rds"))


# Dose allocation probability ----

# Dose allocation probability ----

dose_allocation_probability <- all_allocations %>%
  filter(design %in% designs) %>%
  group_by(design, scenario, dgp, dose_idx, dose_log, dose_mg) %>%
  summarise(
    n_allocated = n(),
    .groups = "drop"
  ) %>%
  group_by(design, scenario, dgp) %>%
  mutate(
    total_allocated = sum(n_allocated),
    p_allocation = n_allocated / total_allocated,
    p_allocation_percent = 100 * p_allocation
  ) %>%
  ungroup()

write.csv(
  dose_allocation_probability,
  here("Results", "tables", "dose_allocation_probability.csv"),
  row.names = FALSE
)


# Scenario-specific MTD table ----
mtd_table <- tibble(
  scenario = names(lords_scenarios),
  mtd_log = map_dbl(lords_scenarios, "mtd_log")
)
library(tidyr)
# For each scenario, find the first unsafe dose index
(unsafe_ranges <- expand_grid(
  scenario = names(lords_scenarios),
  dose_idx = sort(unique(dose_allocation_probability$dose_idx))
) %>%
    left_join(
      dose_allocation_probability %>%
        distinct(dose_idx, dose_log),
      by = "dose_idx"
    ) %>%
    left_join(mtd_table, by = "scenario") %>%
    mutate(
      unsafe = dose_log > mtd_log
    ) %>%
    group_by(scenario) %>%
    summarise(
      first_unsafe_idx = if (any(unsafe)) min(dose_idx[unsafe]) else NA_integer_,
      .groups = "drop"
    ) %>%
    mutate(
      xmin = first_unsafe_idx - 0.5,
      xmax = Inf
    ))
allocation_plot_df <- dose_allocation_probability %>%
  mutate(
    design = factor(design, levels = c("TEPI", "BOIN12", "STEIN", "RWR")),
    scenario = factor(scenario, levels = c("A", "B", "C", "D"), labels = scenario_labels_long),
    dgp = factor(dgp, levels = dgp_labels_all)
  ) %>% 
  mutate(dgp = str_wrap(as.character(dgp), width = 8),
         scenario = str_wrap(as.character(scenario), width = 8))

unsafe_ranges_plot <- unsafe_ranges %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("A", "B", "C", "D"),
      labels = scenario_labels_long
    ),
    scenario = str_wrap(as.character(scenario), width = 8)
  ) %>%
  tidyr::crossing(
    dgp = factor(c(
      "Base case: logit/logit",
      "Misspecification 1: probit/probit",
      "Misspecification 2: quadratic efficacy + logistic toxicity"
    ), levels = c(
      "Base case: logit/logit",
      "Misspecification 1: probit/probit",
      "Misspecification 2: quadratic efficacy + logistic toxicity"
    ))
  ) %>% 
  mutate(dgp = str_wrap(as.character(dgp), width = 8))

(p_allocation <- ggplot(
  allocation_plot_df,
  aes(
    x = factor(dose_idx),
    y = p_allocation_percent,
    fill = design
  )
) +
    geom_rect(
      data = unsafe_ranges_plot %>% filter(!is.na(xmin)),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "#FDE0DD",
      alpha = 0.35
    ) +
    geom_col(
      position = position_dodge(width = 0.75),
      width = 0.65,
      colour = "white"
    ) +
    facet_grid(dgp ~ scenario) +
    scale_fill_manual(values = design_cols) +
    scale_y_continuous(
      labels = label_percent(scale = 1),
      breaks = seq(0, 100, by = 20),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = "Probability of Dose Allocation During the Trial",
      subtitle = "Percentage of enrolled patients allocated to each dose across simulated trials",
      x = "Dose level",
      y = "Allocation probability (%)",
      fill = "Design",
      caption = "Shaded regions indicate dose levels above the scenario-specific MTD and are therefore considered unsafe."
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 12),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9.5),
      legend.key.size = unit(0.35, "cm"),
      legend.spacing.x = unit(0.15, "cm"),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      plot.caption = element_text(size = 9, hjust = 0)
    ))

ggsave(
  filename = here("Results", "figures", "dose_allocation_probability.png"),
  plot = p_allocation,
  width = 10,
  height = 8,
  dpi = 300
)

# write.csv(
#   dose_allocation_probability,
#   here("Results", "tables", "dose_allocation_probability.csv"),
#   row.names = FALSE
# )
# 
# allocation_plot_df <- dose_allocation_probability %>%
#   mutate(
#     design = factor(design, levels = c("TEPI", "BOIN12", "STEIN", "RWR")),
#     scenario = factor(scenario, levels = c("A", "B", "C", "D"), labels = scenario_labels_long),
#     dgp = factor(dgp, levels = dgp_labels_all)
#   ) %>% 
#   mutate(
#     scenario = str_wrap(scenario, width = 8),
#     dgp = str_wrap(dgp, width = 10)
#   )
# 
# p_allocation <- ggplot(
#   allocation_plot_df,
#   aes(
#     x = factor(dose_idx),
#     y = p_allocation_percent,
#     fill = design
#   )
# ) +
#   geom_col(
#     position = position_dodge(width = 0.75),
#     width = 0.65,
#     colour = "white"
#   ) +
#   facet_grid(dgp ~ scenario) +
#   scale_fill_manual(values = design_cols) +
#   scale_y_continuous(
#     labels = label_percent(scale = 1),
#     breaks = seq(0, 100, by = 20),
#     expand = expansion(mult = c(0, 0.02))
#   ) +
#   labs(
#     title = "Probability of Dose Allocation During the Trial",
#     subtitle = "Percentage of enrolled patients allocated to each dose across simulated trials",
#     x = "Dose level",
#     y = "Allocation probability (%)",
#     fill = "Design"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(face = "bold", size = 15),
#     plot.subtitle = element_text(size = 12),
#     legend.position = "top",
#     legend.title = element_text(face = "bold", size = 10),
#     legend.text = element_text(size = 9.5),
#     legend.key.size = unit(0.35, "cm"),
#     legend.spacing.x = unit(0.15, "cm"),
#     strip.text = element_text(face = "bold", size = 10),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
#     axis.title = element_text(face = "bold"),
#     axis.text.x = element_text(size = 10),
#     axis.text.y = element_text(size = 10)
#   )
# 
# print(p_allocation)

ggsave(
  filename = here("Results", "figures", "dose_allocation_probability.png"),
  plot = p_allocation,
  width = 12,
  height = 8,
  dpi = 300
)


# Expected outcomes from patient-level data ----

expected_outcomes_patient_level <- all_allocations %>%
  filter(design %in% designs) %>%
  group_by(design, scenario, dgp, rep) %>%
  summarise(
    n_patients = n(),
    
    # Toxicity outcome
    n_tox = sum(tox == 1, na.rm = TRUE),
    
    # Success = efficacy without toxicity
    n_success = sum(tox == 0 & eff == 1, na.rm = TRUE),
    
    # Neutral = no toxicity and no efficacy
    n_neutral = sum(tox == 0 & eff == 0, na.rm = TRUE),
    
    toxicity_rate = n_tox / n_patients,
    success_rate = n_success / n_patients,
    neutral_rate = n_neutral / n_patients,
    
    .groups = "drop"
  ) %>%
  group_by(design, scenario, dgp) %>%
  summarise(
    n_trials = n(),
    mean_n_patients = mean(n_patients, na.rm = TRUE),
    
    expected_toxicities = mean(n_tox, na.rm = TRUE),
    sd_toxicities = sd(n_tox, na.rm = TRUE),
    expected_toxicity_rate_percent = 100 * mean(toxicity_rate, na.rm = TRUE),
    
    expected_successes = mean(n_success, na.rm = TRUE),
    sd_successes = sd(n_success, na.rm = TRUE),
    expected_success_rate_percent = 100 * mean(success_rate, na.rm = TRUE),
    
    expected_neutral = mean(n_neutral, na.rm = TRUE),
    sd_neutral = sd(n_neutral, na.rm = TRUE),
    expected_neutral_rate_percent = 100 * mean(neutral_rate, na.rm = TRUE),
    
    .groups = "drop"
  )

# expected_outcomes_patient_level <- all_allocations %>%
#   filter(design %in% designs) %>%
#   group_by(design, scenario, dgp, rep) %>%
#   summarise(
#     n_patients = n(),
#     n_tox = sum(tox, na.rm = TRUE),
#     n_eff = sum(eff, na.rm = TRUE),
#     n_neutral = sum(outcome == "neutral", na.rm = TRUE),
#     toxicity_rate = n_tox / n_patients,
#     efficacy_rate = n_eff / n_patients,
#     neutral_rate = n_neutral / n_patients,
#     .groups = "drop"
#   ) %>%
#   group_by(design, scenario, dgp) %>%
#   summarise(
#     n_trials = n(),
#     mean_n_patients = mean(n_patients, na.rm = TRUE),
#     
#     expected_toxicities = mean(n_tox, na.rm = TRUE),
#     sd_toxicities = sd(n_tox, na.rm = TRUE),
#     expected_toxicity_rate_percent = 100 * mean(toxicity_rate, na.rm = TRUE),
#     
#     expected_efficacies = mean(n_eff, na.rm = TRUE),
#     sd_efficacies = sd(n_eff, na.rm = TRUE),
#     expected_efficacy_rate_percent = 100 * mean(efficacy_rate, na.rm = TRUE),
#     
#     expected_neutral = mean(n_neutral, na.rm = TRUE),
#     sd_neutral = sd(n_neutral, na.rm = TRUE),
#     expected_neutral_rate_percent = 100 * mean(neutral_rate, na.rm = TRUE),
#     
#     .groups = "drop"
#   )

write.csv(
  expected_outcomes_patient_level,
  here("Results", "tables", "expected_outcomes_patient_level.csv"),
  row.names = FALSE
)

# expected_plot_df <- expected_outcomes_patient_level %>%
#   mutate(
#     design = factor(design, levels = c("TEPI", "BOIN12", "STEIN", "RWR")),
#     scenario = factor(scenario, levels = c("A", "B", "C", "D"), labels = scenario_labels_short),
#     dgp = factor(dgp, levels = dgp_labels_all)
#   ) %>% 
#   mutate(dgp = str_wrap(as.character(dgp), width = 8))
# 
# # Toxicities
# p_expected_toxicities <- ggplot(
#   expected_plot_df,
#   aes(x = scenario, y = expected_toxicities, fill = design)
# ) +
#   geom_col(position = position_dodge(width = 0.75), width = 0.65, colour = "white") +
#   geom_text(
#     aes(label = round(expected_toxicities, 1)),
#     position = position_dodge(width = 0.75),
#     vjust = -0.35,
#     size = 3
#   ) +
#   facet_wrap(~ dgp, ncol = 1) +
#   scale_fill_manual(values = design_cols) +
#   labs(
#     title = "Expected Number of Toxicities",
#     subtitle = "Average number of toxicities per simulated trial",
#     x = "Scenario",
#     y = "Expected toxicities",
#     fill = "Design"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(face = "bold", size = 15),
#     plot.subtitle = element_text(size = 11),
#     legend.position = "top",
#     legend.title = element_text(face = "bold"),
#     strip.text = element_text(face = "bold"),
#     axis.title = element_text(face = "bold"),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank()
#   )
# 
# print(p_expected_toxicities)
# 
# ggsave(
#   filename = here("Results", "figures", "expected_toxicities_patient_level.png"),
#   plot = p_expected_toxicities,
#   width = 11,
#   height = 8,
#   dpi = 300
# )
# 
# # Efficacies
# p_expected_efficacies <- ggplot(
#   expected_plot_df,
#   aes(x = scenario, y = expected_successes, fill = design)
# ) +
#   geom_col(position = position_dodge(width = 0.75), width = 0.65, colour = "white") +
#   geom_text(
#     aes(label = round(expected_successes, 1)),
#     position = position_dodge(width = 0.75),
#     vjust = -0.35,
#     size = 3
#   ) +
#   facet_wrap(~ dgp, ncol = 1) +
#   scale_fill_manual(values = design_cols) +
#   labs(
#     title = "Expected Number of Efficacies",
#     subtitle = "Average number of efficacy responses per simulated trial",
#     x = "Scenario",
#     y = "Expected efficacies",
#     fill = "Design"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(face = "bold", size = 15),
#     plot.subtitle = element_text(size = 11),
#     legend.position = "top",
#     legend.title = element_text(face = "bold"),
#     strip.text = element_text(face = "bold"),
#     axis.title = element_text(face = "bold"),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank()
#   )
# 
# print(p_expected_efficacies)
# 
# ggsave(
#   filename = here("Results", "figures", "expected_efficacies_patient_level.png"),
#   plot = p_expected_efficacies,
#   width = 11,
#   height = 8,
#   dpi = 300
# )
# 
# # Neutral outcomes
# p_expected_neutral <- ggplot(
#   expected_plot_df,
#   aes(x = scenario, y = expected_neutral, fill = design)
# ) +
#   geom_col(position = position_dodge(width = 0.75), width = 0.65, colour = "white") +
#   geom_text(
#     aes(label = round(expected_neutral, 1)),
#     position = position_dodge(width = 0.75),
#     vjust = -0.35,
#     size = 3
#   ) +
#   facet_wrap(~ dgp, ncol = 1) +
#   scale_fill_manual(values = design_cols) +
#   labs(
#     title = "Expected Number of Neutral Outcomes",
#     subtitle = "Average number of neutral outcomes per simulated trial",
#     x = "Scenario",
#     y = "Expected neutral outcomes",
#     fill = "Design"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(face = "bold", size = 15),
#     plot.subtitle = element_text(size = 11),
#     legend.position = "top",
#     legend.title = element_text(face = "bold"),
#     strip.text = element_text(face = "bold"),
#     axis.title = element_text(face = "bold"),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank()
#   )
# 
# print(p_expected_neutral)
# 
# ggsave(
#   filename = here("Results", "figures", "expected_neutral_patient_level.png"),
#   plot = p_expected_neutral,
#   width = 11,
#   height = 8,
#   dpi = 300
# )

# Joint toxicity / efficacy comparison
expected_side_by_side_df <- expected_outcomes_patient_level %>%
  mutate(
    scenario = factor(scenario, levels = c("A", "B", "C", "D"), labels = scenario_labels_long),
    dgp = factor(dgp, levels = dgp_labels_all),
    design = factor(design, levels = c("TEPI", "BOIN12", "STEIN", "RWR"))
  ) %>%
  select(design, scenario, dgp, expected_toxicities, expected_successes, expected_neutral) %>%
  pivot_longer(
    cols = c(expected_toxicities, expected_successes, expected_neutral),
    names_to = "outcome_type",
    values_to = "expected_count"
  ) %>%
  mutate(
    outcome_type = factor(
      outcome_type,
      levels = c("expected_toxicities", "expected_neutral", "expected_successes"),
      labels = c("Toxicity", "Neutral", "Success")
    )
  ) %>% 
  mutate(dgp = str_wrap(as.character(dgp), width = 8),
         scenario = str_wrap(as.character(scenario), width = 8))

# print(p_expected_side_by_side)
p_expected_side_by_side <- ggplot(
  expected_side_by_side_df,
  aes(x = scenario, y = expected_count, fill = outcome_type)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, colour = "white") +
  geom_text(
    aes(label = round(expected_count, 1)),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 2.8
  ) +
  facet_grid(dgp ~ design) +
  scale_fill_manual(
    values = c(
      "Toxicity" = "#D73027",
      "Neutral"  = "#666666",
      "Success"  = "#2CA25F"
    )
  ) +
  labs(
    title = "Expected Outcomes by Design and Scenario",
    subtitle = "Average numbers per simulated 60-patient trial",
    x = "Scenario",
    y = "Expected count",
    fill = "Outcome"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_expected_side_by_side)
ggsave(
  filename = here("Results", "figures", "expected_toxicity_efficacy_side_by_side.png"),
  plot = p_expected_side_by_side,
  width = 13,
  height = 9,
  dpi = 300
)

# Prepare plotting data ----
expected_outcomes_long <- expected_outcomes_patient_level %>%
  select(design, scenario, dgp, expected_toxicities, expected_successes) %>%
  pivot_longer(
    cols = c(expected_toxicities, expected_successes),
    names_to = "outcome_type",
    values_to = "expected_count"
  ) %>%
  mutate(
    outcome_type = recode(
      outcome_type,
      expected_toxicities = "Toxicity",
      expected_successes = "Success"
    ),
    design = factor(design, levels = c("TEPI", "BOIN12", "STEIN", "RWR")),
    scenario = factor(
      scenario,
      levels = c("A", "B", "C", "D"),
      labels = scenario_labels_long
    ),
    dgp = factor(dgp, levels = dgp_labels_all)
  )

# Plot ----
p_expected_side_by_side <- ggplot(
  expected_side_by_side_df,
  aes(x = design, y = expected_count, fill = outcome_type)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    colour = "white"
  ) +
  geom_text(
    aes(label = round(expected_count, 1)),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 2.8
  ) +
  facet_grid(dgp ~ scenario) +
  scale_fill_manual(values = c("Toxicity" = "#D73027", "Success" = "#2CA25F")) +
  labs(
    title = "Expected Toxicity and Success Outcomes by Design and Scenario",
    subtitle = "Average numbers per simulated 60-patient trial",
    x = "Design",
    y = "Expected count",
    fill = "Outcome"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9)
  )

print(p_expected_side_by_side)

ggsave(
  filename = here("Results", "figures", "expected_toxicity_efficacy_side_by_side.png"),
  plot = p_expected_side_by_side,
  width = 13,
  height = 9,
  dpi = 300
)


# Model-assisted only: TEPI, STEIN, BOIN12 from all_allocations----
## 1. LORDs dose grid----

dose_log <- c(-1.20, -0.23, 0.92, 2.02, 3.00, 3.69, 4.38, 5.08, 5.77)
dose_mg  <- c(0.3, 0.8, 2.5, 7.5, 20, 40, 80, 160, 320)

## 2. LORDs scenarios ----

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

## 3. Settings ----

rep_keep <- 1
dgp_keep <- "base"
dgp_label <- dgp_settings[[dgp_keep]]$label

## 4. Helper functions ----

# Inverse link
inv_link <- function(eta, link = c("logit", "probit")) {
  link <- match.arg(link)
  if (link == "logit") {
    plogis(eta)
  } else {
    pnorm(eta)
  }
}

# Conditional efficacy under continuation-ratio linear form
pE_cond_cr_linear <- function(d, theta, link = "logit") {
  inv_link(theta["theta1"] + theta["theta2"] * d, link = link)
}

# Conditional efficacy under quadratic form
pE_cond_quadratic <- function(d, theta, dU = max(dose_log)) {
  dmax <- (qlogis(0.95) - theta["theta1"]) / theta["theta2"]
  
  p_left <- plogis(theta["theta1"] + theta["theta2"] * d)
  p_right <- ((0.8 - 0.95) * (d - dmax)^2) / ((dU - dmax)^2) + 0.95
  
  out <- ifelse(d <= dmax, p_left, p_right)
  pmin(pmax(out, 0), 1)
}

# Toxicity model
pT_model <- function(d, theta, tox_link = "logit") {
  inv_link(theta["theta3"] + theta["theta4"] * d, link = tox_link)
}

# Truth probabilities
truth_probs <- function(d, theta, dgp) {
  pE_cond <- switch(
    dgp$eff_model,
    cr_linear = pE_cond_cr_linear(d, theta, link = dgp$eff_link),
    quadratic = pE_cond_quadratic(d, theta, dU = max(dose_log))
  )
  
  pT <- pT_model(d, theta, tox_link = dgp$tox_link)
  
  tibble(
    dose_log = d,
    neutral  = (1 - pE_cond) * (1 - pT),
    success  = pE_cond * (1 - pT),
    toxicity = pT
  )
}

# Empirical success at each dose
empirical_summary <- function(trial, d_grid) {
  trial %>%
    mutate(success = tox == 0 & eff == 1) %>%
    group_by(dose_log) %>%
    summarise(
      success_hat = mean(success, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(tibble(dose_log = d_grid), by = "dose_log") %>%
    arrange(match(dose_log, d_grid))
}

# Fit toxicity and conditional efficacy models
fit_pair <- function(trial, link = c("logit", "probit"), method = c("glm", "firth")) {
  link <- match.arg(link)
  method <- match.arg(method)
  
  dat_eff <- trial %>% filter(tox == 0)
  
  if (nrow(trial) < 2 || nrow(dat_eff) < 2) {
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
    tox_fit <- try(logistf::logistf(tox ~ dose_log, data = trial), silent = TRUE)
    eff_fit <- try(logistf::logistf(eff ~ dose_log, data = dat_eff), silent = TRUE)
  }
  
  if (inherits(tox_fit, "try-error")) tox_fit <- NULL
  if (inherits(eff_fit, "try-error")) eff_fit <- NULL
  
  list(tox_fit = tox_fit, eff_fit = eff_fit)
}

# Safe prediction
predict_response_safe <- function(fit, newdat) {
  if (is.null(fit)) return(rep(NA_real_, nrow(newdat)))
  
  out <- try(
    predict(fit, newdata = newdat, type = "response"),
    silent = TRUE
  )
  
  if (inherits(out, "try-error")) return(rep(NA_real_, nrow(newdat)))
  as.numeric(out)
}

# Fitted success curve
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

## 5. Reconstruct curves from all_allocations ----

example_curves_assisted <- map_dfr(c("TEPI", "STEIN", "BOIN12"), function(des) {
  map_dfr(scenario_names, function(sc) {
    
    trial <- all_allocations %>%
      filter(
        design == des,
        scenario == sc,
        dgp == dgp_label,
        rep == rep_keep
      ) %>%
      select(dose_log, tox, eff)
    
    if (nrow(trial) == 0) return(NULL)
    
    sc_obj  <- lords_scenarios[[sc]]
    dgp_obj <- dgp_settings[[dgp_keep]]
    
    curve_true <- truth_probs(dose_log, sc_obj$theta, dgp_obj) %>%
      transmute(
        dose_log,
        success_hat = success,
        method = "Truth"
      )
    
    curve_emp <- empirical_summary(trial, dose_log) %>%
      mutate(method = "Empirical")
    
    fit_logit <- fit_pair(trial, link = "logit", method = "glm")
    curve_logit <- predict_success_curve(fit_logit, dose_log, "Logit GLM")
    
    fit_probit <- fit_pair(trial, link = "probit", method = "glm")
    curve_probit <- predict_success_curve(fit_probit, dose_log, "Probit GLM")
    
    fit_firth <- fit_pair(trial, link = "logit", method = "firth")
    curve_firth <- predict_success_curve(fit_firth, dose_log, "Firth logit")
    
    bind_rows(
      curve_true,
      curve_emp,
      curve_logit,
      curve_probit,
      curve_firth
    ) %>%
      mutate(
        design = des,
        scenario = sc,
        dgp = dgp_label
      )
  })
})

## 6. Format for plotting ----

example_curves_assisted <- example_curves_assisted %>%
  mutate(
    scenario = factor(scenario, levels = c("A", "B", "C", "D")),
    design = factor(design, levels = c("TEPI", "BOIN12", "STEIN"))
  )

## 7. Plot ----

p_example_curves_assisted <- ggplot() +
  geom_line(
    data = subset(example_curves_assisted, method %in% c("Truth", "Logit GLM", "Probit GLM", "Firth logit")),
    aes(x = dose_log, y = success_hat, colour = method, group = method),
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  geom_point(
    data = subset(example_curves_assisted, method %in% c("Truth", "Logit GLM", "Probit GLM", "Firth logit")),
    aes(x = dose_log, y = success_hat, colour = method, group = method),
    size = 1.1,
    na.rm = TRUE
  ) +
  geom_smooth(
    data = subset(example_curves_assisted, method == "Empirical" & !is.na(success_hat)),
    aes(x = dose_log, y = success_hat, colour = method, group = method),
    method = "loess",
    span = 0.9,
    se = FALSE,
    linewidth = 1,
    na.rm = TRUE
  ) +
  geom_point(
    data = subset(example_curves_assisted, method == "Empirical"),
    aes(x = dose_log, y = success_hat, colour = method, group = method),
    size = 1.1,
    na.rm = TRUE
  ) +
  facet_wrap(~scenario, ncol = 2) +
  scale_x_continuous(breaks = dose_log) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual(values = c(
    Truth = "#C77CFF",
    Empirical = "grey40",
    `Logit GLM` = "#1B9E77",
    `Probit GLM` = "#0072B2",
    `Firth logit` = "#D95F02"
  )) +
  labs(
    title = "Model-assisted designs: estimated success probability",
    subtitle = paste("TEPI, STEIN, BOIN12; reconstructed from all_allocations, replicate", rep_keep),
    x = "Dose (log scale)",
    y = "Estimated success probability"
  ) +
  theme_thesis +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 10)
  )

print(p_example_curves_assisted)

ggsave(
  filename = here("Results", "figures", "example_curves_assisted.png"),
  plot = p_example_curves_assisted,
  width = 10,
  height = 8,
  dpi = 300
)

# Mode selected dose summary ----
mode_summary <- all_reps %>%
  group_by(design, scenario, dgp) %>%
  summarise(
    mode_selected_log = as.numeric(names(which.max(table(selected_log)))),
    true_obd_log = first(true_obd_log),
    mode_minus_true_obd = mode_selected_log - true_obd_log,
    n_trials = n(),
    .groups = "drop"
  ) %>%
  mutate(
    mode_selected_log = round(mode_selected_log, 2),
    true_obd_log = round(true_obd_log, 2),
    mode_minus_true_obd = round(mode_minus_true_obd, 2)
  )

print(mode_summary)

write.csv(
  mode_summary,
  here("Results", "tables", "mode_summary.csv"),
  row.names = FALSE
)

# Wide mode table

mode_table <- mode_summary %>%
  mutate(
    scenario_dgp = paste(scenario, dgp, sep = " | ")
  ) %>%
  select(design, scenario_dgp, mode_selected_log) %>%
  pivot_wider(
    names_from = scenario_dgp,
    values_from = mode_selected_log
  )

print(mode_table)

write.csv(
  mode_table,
  here("Results", "tables", "mode_selected_table.csv"),
  row.names = FALSE
)

# True OBD table ----

true_obd_table <- all_reps %>%
  group_by(scenario) %>%
  summarise(
    true_obd_log = first(true_obd_log),
    .groups = "drop"
  ) %>%
  mutate(
    true_obd_log = round(true_obd_log, 2),
    scenario_label = recode(
      scenario,
      A = "A (Narrow)",
      B = "B (Wide)",
      C = "C (Safe)",
      D = "D (Unsafe)"
    )
  ) %>%
  select(scenario, scenario_label, true_obd_log)

print(true_obd_table)

write.csv(
  true_obd_table,
  here("Results", "tables", "true_obd_table.csv"),
  row.names = FALSE
)


# Base-case-only table

mode_table_base <- mode_summary %>%
  filter(dgp == "Base case: logit/logit") %>%
  select(design, scenario, mode_selected_log) %>%
  pivot_wider(
    names_from = scenario,
    values_from = mode_selected_log
  ) %>%
  left_join(
    true_obd_table %>%
      select(scenario, true_obd_log) %>%
      pivot_wider(
        names_from = scenario,
        values_from = true_obd_log
      ) %>%
      mutate(design = "True OBD"),
    by = "design"
  )

print(mode_table_base)

write.csv(
  mode_table_base,
  here("Results", "tables", "mode_selected_base_case_with_true_obd.csv"),
  row.names = FALSE
)

# Therapeutic range table ----

therapeutic_range <- tibble(
  scenario = names(lords_scenarios),
  mined_log = purrr::map_dbl(lords_scenarios, "mined_log"),
  mtd_log   = purrr::map_dbl(lords_scenarios, "mtd_log")
)

all_reps_therapeutic <- all_reps %>%
  left_join(therapeutic_range, by = "scenario") %>%
  mutate(
    within_therapeutic_range =
      selected_log >= mined_log &
      selected_log <= mtd_log
  )

p_therapeutic <- all_reps_therapeutic %>%
  group_by(design, scenario, dgp) %>%
  summarise(
    n_trials = n(),
    n_within_therapeutic_range = sum(within_therapeutic_range, na.rm = TRUE),
    p_within_therapeutic_range = mean(within_therapeutic_range, na.rm = TRUE),
    p_within_therapeutic_range_percent = 100 * p_within_therapeutic_range,
    .groups = "drop"
  )

print(p_therapeutic)

write.csv(
  p_therapeutic,
  here("Results", "tables", "probability_within_therapeutic_range.csv"),
  row.names = FALSE
)

# Heatmap ----

heatmap_table <- p_therapeutic %>%
  filter(design %in% c("TEPI", "BOIN12", "STEIN")) %>%
  mutate(
    design = recode(design,
                    "BOIN12" = "BOIN12"
    ),
    scenario_dgp = case_when(
      scenario == "A" & dgp == "Base case: logit/logit" ~ "A: Narrow\nBase case",
      scenario == "A" & dgp == "Misspecification 1: probit/probit" ~ "A: Narrow\nMisspec 1",
      scenario == "A" & dgp == "Misspecification 2: quadratic efficacy + logistic toxicity" ~ "A: Narrow\nMisspec 2",
      
      scenario == "B" & dgp == "Base case: logit/logit" ~ "B: Wide\nBase case",
      scenario == "B" & dgp == "Misspecification 1: probit/probit" ~ "B: Wide\nMisspec 1",
      scenario == "B" & dgp == "Misspecification 2: quadratic efficacy + logistic toxicity" ~ "B: Wide\nMisspec 2",
      
      scenario == "C" & dgp == "Base case: logit/logit" ~ "C: Safe\nBase case",
      scenario == "C" & dgp == "Misspecification 1: probit/probit" ~ "C: Safe\nMisspec 1",
      scenario == "C" & dgp == "Misspecification 2: quadratic efficacy + logistic toxicity" ~ "C: Safe\nMisspec 2",
      
      scenario == "D" & dgp == "Base case: logit/logit" ~ "D: Unsafe\nBase case",
      scenario == "D" & dgp == "Misspecification 1: probit/probit" ~ "D: Unsafe\nMisspec 1",
      scenario == "D" & dgp == "Misspecification 2: quadratic efficacy + logistic toxicity" ~ "D: Unsafe\nMisspec 2"
    )
  )

# Heatmap plot

p_heatmap <- ggplot(
  heatmap_table,
  aes(
    x = scenario_dgp,
    y = design,
    fill = p_within_therapeutic_range_percent
  )
) +
  geom_tile(colour = "white") +
  geom_text(
    aes(label = paste0(round(p_within_therapeutic_range_percent, 1), "%")),
    size = 3
  ) +
  scale_fill_gradient(
    low = "#F7F7F7",
    high = "#3B4CC0",
    limits = c(0, 100),
    name = "%"
  ) +
  labs(
    title = "Probability of Selecting Within the Therapeutic Range (MinED–MTD)",
    subtitle = "Percentage of simulated trials selecting a dose within the scenario-specific therapeutic range",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p_heatmap)

ggsave(
  filename = here("Results", "figures", "therapeutic_range_heatmap.png"),
  plot = p_heatmap,
  width = 9.5,
  height = 6,
  dpi = 300
)


# base-case-only mode table (slide-friendly format)

mode_table_slide <- mode_summary %>%
  filter(dgp == "Base case: logit/logit", design %in% c("TEPI", "BOIN12", "STEIN")) %>%
  mutate(
    design = factor(design, levels = c("TEPI", "BOIN12", "STEIN"))
  ) %>%
  select(design, scenario, mode_selected_log) %>%
  pivot_wider(
    names_from = scenario,
    values_from = mode_selected_log
  )

true_obd_row <- true_obd_table %>%
  select(scenario, true_obd_log) %>%
  pivot_wider(
    names_from = scenario,
    values_from = true_obd_log
  ) %>%
  mutate(design = "True OBD") %>%
  select(design, everything())

mode_table_slide <- bind_rows(mode_table_slide, true_obd_row)

print(mode_table_slide)

write.csv(
  mode_table_slide,
  here("Results", "tables", "mode_selected_slide_table.csv"),
  row.names = FALSE
)

oc_summary_with_rmse <- all_reps %>%
  group_by(design, scenario, dgp) %>%
  summarise(
    mean_bias_idx = mean(selected_idx - true_obd_idx, na.rm = TRUE),
    rmse_idx = sqrt(mean((selected_idx - true_obd_idx)^2, na.rm = TRUE)),
    mean_bias_log = mean(selected_log - true_obd_log, na.rm = TRUE),
    rmse_log = sqrt(mean((selected_log - true_obd_log)^2, na.rm = TRUE)),
    .groups = "drop"
  )

print(oc_summary_with_rmse)

write.csv(
  oc_summary_with_rmse,
  here("Results", "tables", "oc_summary_with_rmse.csv"),
  row.names = FALSE
)
