###############################################################################
# Analysis code for 2021 sea urchin heatwaves experiment at the Quadra Island # 
# Ecological Observatory - Hakai Institute                                    #
#                                                                             #
# Demipyramid growth (calcified mouth part)                                   #
#                                                                             #
# Author: Nathan B. Spindel                                                   #
###############################################################################

# Load CRAN libraries and utility functions  ----------------------------------
libs <- c(
  "tidyverse", 
  "janitor",
  "parallel",
  "rstan",
  "brms",
  "cmdstanr",
  "posterior",
  "tidybayes",
  "emmeans",
  "glmmTMB",
  "DHARMa",
  "ggrepel",
  "scales",
  "Hmisc",
  "readxl"
)

lapply(libs, require, character.only = TRUE)

source("Utility_Functions.R")

# Load/wrangle data -----------------------------------------------------------
data_demipyramid <- read.csv("./Data/Demipyramid_Growth.csv") %>%
  mutate(
    tank_id = as.factor(tank_id),
    pit_tag_id = as.factor(pit_tag),
    pco2_uatm = factor(pco2_uatm, levels = c("600", "1200"), labels = c("600 μatm", "1200 μatm"))
  ) %>%
  mutate(
    demi_growth_total_um = demi_final_um - demi_tag_um,
    demi_growth_per_day_um = demi_growth_total_um / total_days
  ) %>%
  filter(
    between(demi_growth_per_day_um,0.1,5)
  ) %>%
  mutate(
    temp_num = as.numeric(case_when(
      temp_deg_c == 10 ~ 10,
      temp_deg_c == 13 ~ 13,
      temp_deg_c == 16 ~ 16,
      temp_deg_c == 17 ~ 17,
      temp_deg_c == 18 ~ 18,
      temp_deg_c == 20 ~ 20,
      temp_deg_c == "18-14" ~ 16,
      temp_deg_c == "21-18" ~ 20
    ))
  )

path_initial <- "./Data/Biometrics_Initial.xlsx"
path_final   <- "./Data/Biometrics_Final.xlsx"

# Read and clean both datasets
bio_initial <- read_excel(path_initial) %>%
  dplyr::mutate(
    pit_tag_id = as.factor(pit_tag_id),
    pco2_uatm = as.factor(ifelse(pco2_h_l == "l", "600 μatm", "1200 μatm"))
    ) %>%
  janitor::clean_names() %>%
  rename_with(~ paste0("initial_", .x), 
              .cols = -c(pit_tag_id, temp_deg_c, pco2_h_l, pco2_uatm, tank_id)) %>%
  rename(initial_date = initial_date)  # if date column already named 'date'

bio_final <- read_excel(path_final) %>%
  dplyr::mutate(
    pit_tag_id = as.factor(pit_tag_id),
    pco2_uatm = as.factor(ifelse(pco2_h_l == "l", "600 μatm", "1200 μatm"))
  ) %>%
  janitor::clean_names() %>%
  rename_with(~ paste0("final_", .x), 
              .cols = -c(pit_tag_id, temp_deg_c, pco2_h_l, pco2_uatm, tank_id)) %>%
  rename(final_date = final_date)  # if date column already named 'date'

# Join the two datasets by unique identifiers and treatment keys
bio_joined <- full_join(
  bio_initial,
  bio_final,
  by = c("pit_tag_id", "temp_deg_c", "pco2_uatm", "tank_id")
) %>%
  select(
    initial_date,
    final_date,
    temp_deg_c,
    pco2_uatm,
    tank_id,
    pit_tag_id,
    initial_diam_mm,
    initial_ht_mm,
    initial_wet_mass_g,
    final_diam_mm,
    final_ht_mm,
    final_sex
  ) %>%
  mutate(
    tank_id = as.factor(tank_id)
  )

bio_demi_joined <- full_join(
  bio_joined,
  data_demipyramid,
  by = c("pit_tag_id", "temp_deg_c", "pco2_uatm", "tank_id")
)

skel_growth <- bio_demi_joined %>%
  mutate(
    temp_label = factor(case_when(
      temp_deg_c == "18-14" ~ "La Niña (18-14 °C)",
      temp_deg_c == "21-18" ~ "El Niño (21-18 °C)"
    ), 
    levels = c("La Niña (18-14 °C)", "El Niño (21-18 °C)"),
    labels = c("La Niña, 600 μatm", "El Niño, 600 μatm")),
    
    # Dates and elapsed time
    initial_date = as.Date(initial_date),
    final_date   = as.Date(final_date),
    delta_days   = as.numeric(difftime(final_date, initial_date, units = "days")),
    
    # --- Compute spheroid volumes (mL) ---
    initial_a_mm = 0.5 * initial_diam_mm,
    initial_c_mm = 0.5 * initial_ht_mm,
    initial_volume_ml = (4/3 * pi * initial_a_mm^2 * initial_c_mm) / 1000,
    
    final_a_mm = 0.5 * final_diam_mm,
    final_c_mm = 0.5 * final_ht_mm,
    final_volume_ml = (4/3 * pi * final_a_mm^2 * final_c_mm) / 1000,
    
    # --- Absolute deltas ---
    delta_diam_mm   = final_diam_mm   - initial_diam_mm,
    delta_ht_mm     = final_ht_mm     - initial_ht_mm,
    delta_volume_ml = final_volume_ml - initial_volume_ml,
    delta_demipyramid_um = demi_final_um - demi_tag_um,
    
    # --- Absolute growth rates (time-standardized) ---
    diam_growth_mm_day   = delta_diam_mm / delta_days,
    ht_growth_mm_day     = delta_ht_mm   / delta_days,
    volume_growth_ml_day = delta_volume_ml / delta_days,
    demipyramid_growth_um_day = delta_demipyramid_um / delta_days,
    
    diam_growth_mm_week   = diam_growth_mm_day * 7,
    ht_growth_mm_week     = ht_growth_mm_day * 7,
    volume_growth_ml_week = volume_growth_ml_day * 7,
    demipyramid_growth_um_week = demipyramid_growth_um_day * 7,
    
    diam_growth_mm_month   = diam_growth_mm_day * 30.44,
    ht_growth_mm_month     = ht_growth_mm_day * 30.44,
    volume_growth_ml_month = volume_growth_ml_day * 30.44,
    demipyramid_growth_um_month = demipyramid_growth_um_day * 30.44,
    
    # --- Percent changes (relative to initial) ---
    perc_change_diam   = (delta_diam_mm   / initial_diam_mm)   * 100,
    perc_change_ht     = (delta_ht_mm     / initial_ht_mm)     * 100,
    perc_change_volume = (delta_volume_ml / initial_volume_ml) * 100,
    perc_change_demipyramid = (delta_demipyramid_um / demi_tag_um) * 100,
    
    # --- Percent change rates (time-standardized) ---
    perc_change_diam_day   = perc_change_diam / delta_days,
    perc_change_ht_day     = perc_change_ht   / delta_days,
    perc_change_volume_day = perc_change_volume / delta_days,
    perc_change_demipyramid_day = perc_change_demipyramid / delta_days,
    
    perc_change_diam_week   = perc_change_diam_day * 7,
    perc_change_ht_week     = perc_change_ht_day * 7,
    perc_change_volume_week = perc_change_volume_day * 7,
    perc_change_demipyramid_week = perc_change_demipyramid_day * 7,
    
    perc_change_diam_month   = perc_change_diam_day * 30.44,
    perc_change_ht_month     = perc_change_ht_day * 30.44,
    perc_change_volume_month = perc_change_volume_day * 30.44,
    perc_change_demipyramid_month = perc_change_demipyramid_day * 30.44
  )

data_skel_growth_const_temp <- skel_growth %>%
  filter(temp_deg_c %nin% c("21-18", "18-14"))

data_skel_growth_dynamic_temp <- skel_growth %>%
  filter(temp_deg_c %in% c("21-18", "18-14"))

data_skel_growth_const_dynamic_temp <- skel_growth %>%
  filter(pco2_uatm == "600 μatm")

# Set general computing params ------------------------------------------------
cores        <- 26
chains       <- 4
iterations   <- 18000
warmup       <- 2000
control_list <- list(max_treedepth = 15, adapt_delta = 0.95)
file_refit   <- "always"
backend      <- "cmdstanr"
seed         <- 123
brms_fam <- Gamma(link = "log")

# unlink("./Model_Output/demipyramid_temp_oa_*", force = TRUE)


priors <- c( 
  prior(normal(0,1), class = "b"),
  prior(student_t(3, 0.3, 2.5), class = "Intercept"),
  prior(inv_gamma(4.6, 1.1), lb = 0, class = "lscale"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sdgp"),
  prior(gamma(0.01,0.01),lb = 0, class = "shape")
)
# Priors ----
priors_common <- c(
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 2, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01,0.01),lb = 0, class = "shape")
)

priors_gp <- c(
  prior(inv_gamma(4.6, 1.1), class = "lscale", coef = "gptemp_modelpco2_uatm1200"),
  prior(inv_gamma(4.6, 1.1), class = "lscale", coef = "gptemp_modelpco2_uatm600"),
  prior(R2D2(mean_R2 = 0.1, prec_R2 = 30, cons_D2 = 0.5), class = "sdgp")
)

priors_gp_noby <- c(
  prior(inv_gamma(4.6, 1.1), class = "lscale"),
  prior(R2D2(mean_R2 = 0.1, prec_R2 = 30, cons_D2 = 0.5), class = "sdgp")
)

priors_dynamic <- c( 
  prior(normal(0,0.1), class = "b"),
  prior(student_t(3, 0.2, 5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01,0.01),lb = 0, class = "shape")
)

# Candidate model formulae for constant temperatures ----
# response <- "delta_demipyramid_um"

# Full model
form_full <- brms::bf(
  delta_demipyramid_um ~
    log(initial_volume_ml) +
    gp(temp_num, k = 5, c = 5/4, by = pco2_uatm) +
    (1 | tank_id)
)

# Reduced models chosen to test specific hypotheses with principled nesting
form_gp_noby <- brms::bf(
  delta_demipyramid_um ~
    log(initial_volume_ml) +
    gp(temp_num, k = 5, c = 5/4) +
    (1 | tank_id)
)
# (Removes the pCO2-specific GP surfaces; retains nonparametric thermal response.)

form_lin_int <- brms::bf(
  delta_demipyramid_um ~
    log(initial_volume_ml) + temp_num * pco2_uatm +
    (1 | tank_id)
)
# (Parametric T × pCO2 interaction; no GP.)

form_lin_add <- brms::bf(
  delta_demipyramid_um ~
    log(initial_volume_ml) + temp_num + pco2_uatm +
    (1 | tank_id)
)
# (Additive parametric effects; no GP.)

form_size_only <- brms::bf(
  delta_demipyramid_um ~
    log(initial_volume_ml) + (1 | tank_id)
)
# (Allometry + random tank only.)

model_set <- list(
  full         = form_full,
  gp_noby      = form_gp_noby,
  lin_int      = form_lin_int,
  lin_add      = form_lin_add,
  size_only    = form_size_only
)

# Fit candidate constant temperature models ----
# Compact working df
d <- data_skel_growth_const_temp %>%
  select(delta_demipyramid_um, initial_volume_ml, temp_num, pco2_uatm, tank_id) %>%
  mutate(temp_num = as.numeric(temp_num)) %>%
  drop_na()

fit_all <- list()

for (nm in names(model_set)) {
  message("Refitting model: ", nm)
  
  fit_all[[nm]] <- brm(
    formula      = model_set[[nm]],
    data         = d,
    family       = brms_fam,
    prior        = priors_for(model_set[[nm]]),
    sample_prior = "yes",
    save_pars    = save_pars(all = TRUE),
    init         = 0,
    iter         = iterations,
    warmup       = warmup,
    chains       = chains,
    cores        = cores,
    control      = control_list,
    backend      = backend,
    seed         = seed,
    file_refit   = "always",  # force same number of draws for all
    file         = paste0("./Model_Output/demipyramid_temp_oa_", nm)
  )
}

brms_formula_full_dyn <- bf(delta_demipyramid_um ~ 
                              log(initial_volume_ml) + 
                              temp_num + 
                              (1|tank_id))

# fit model
fit_demipyramid_temp_dynamic_temp <- brm(
  formula = brms_formula_full_dyn, 
  data = data_skel_growth_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors_common,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  warmup = warmup,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/demipyramid_temp_dynamic_full"
)

brms_formula_full_const_dyn <- bf(delta_demipyramid_um ~
                                    log(initial_volume_ml) +
                                    temp_deg_c +
                                    (1|tank_id))

# fit model
fit_demipyramid_const_temp_dynamic_temp <- brm(
  formula = brms_formula_full_const_dyn,
  data = data_skel_growth_const_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors_common,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  warmup = warmup,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/demipyramid_temp_const_dynamic_full"
)

# Model selection based on LOO stacking weights ----
fit_all <- lapply(fit_all, add_criterion, criterion = "loo", reloo = FALSE)
loos    <- lapply(fit_all, \(x) x$criteria$loo)
sapply(loos, function(x) attr(x, "dims")[1])
# Should print identical numbers 
# full   gp_noby   lin_int   lin_add size_only 
# 64000     64000     64000     64000     64000

# Stacking and pseudo-BMA (Bayesian bootstrap) weights
# w_stack <- loo::model_weights(loos, method = "stacking")
w_stack <- brms::loo_model_weights(loos, method = "stacking")
w_stack
# Method: stacking
# 
#   weight
# full      1.000 # <- favored
# gp_noby   0.000 
# lin_int   0.000 
# lin_add   0.000 
# size_only 0.000

w_pbma  <- brms::loo_model_weights(loos, method = "pseudobma")
w_pbma
# Method: pseudo-BMA+ with Bayesian bootstrap
# 
#   weight
# full      0.664 # <- favored
# gp_noby   0.048 
# lin_int   0.103 
# lin_add   0.141 
# size_only 0.045 

# Pareto-k summaries
pk <- lapply(loos, \(x) x$diagnostics$pareto_k)
sapply(pk, \(v) table(cut(v,
                          breaks = c(-Inf, 0.5, 0.7, 1, Inf),
                          labels = c("good (≤0.5)", "ok (≤0.7)", "bad (≤1)", "very bad (>1)")
)))
#               full gp_noby lin_int lin_add size_only
# good (≤0.5)    311     311     311     311       311
# ok (≤0.7)        0       0       0       0         0
# bad (≤1)         0       0       0       0         0
# very bad (>1)    0       0       0       0         0

selected_model <- fit_all$full
selected_model$formula
# delta_demipyramid_um ~ log(initial_volume_ml) + gp(temp_num, k = 5, c = 5/4, by = pco2_uatm) + (1 | tank_id) 
bayes_R2(selected_model)
# Estimate  Est.Error      Q2.5     Q97.5
# R2 0.3074025 0.04344842 0.2216749 0.3902301

# Posterior predictive check for favored model
pp_check(selected_model, ndraws = 50)

# Save session info for provenance
sink("./Model_Output/sessionInfo_Demipyramid.txt"); print(sessionInfo()); sink()

# Figure 5. Visualize conditional effects on Demipyramid growth ----
fit_demi_const_temp <- readRDS("./Model_Output/demipyramid_temp_oa_full.rds")
selected_model <- fit_demi_const_temp
selected_model$formula
# delta_demipyramid_um ~ log(initial_volume_ml) + gp(temp_num, k = 5, c = 5/4, by = pco2_uatm) + (1 | tank_id) 

ce_95_pco2_comparison <- conditional_effects(selected_model, robust = F)

df_95_skel_constant_plot <- ce_95_pco2_comparison$`temp_num:pco2_uatm` %>%
  dplyr::mutate(pco2_uatm = factor(pco2_uatm, levels = c("600 μatm", "1200 μatm")))

df_95_skel_allom_constant_plot <- ce_95_pco2_comparison$initial_volume_ml

ce_80_pco2_comparison <- conditional_effects(selected_model, robust = F, prob = 0.80)

df_80_skel_constant_plot <- ce_80_pco2_comparison$`temp_num:pco2_uatm` %>%
  dplyr::mutate(pco2_uatm = factor(pco2_uatm, levels = c("600 μatm", "1200 μatm")))

df_80_skel_allom_constant_plot <- ce_80_pco2_comparison$initial_volume_ml

dynamic_mod <- readRDS("./Model_Output/demipyramid_temp_const_dynamic_full.rds")
dynamic_mod$formula
# delta_demipyramid_um ~ log(initial_volume_ml) + temp_deg_c + (1 | tank_id) 
ce_95_skel_dynamic_temp <- conditional_effects(dynamic_mod, 
                                               robust = F, prob = 0.95)

df_95_skel_dynamic_plot <- ce_95_skel_dynamic_temp$temp_deg_c %>%
  filter(temp_deg_c %in% c("21-18", "18-14")) %>%
  dplyr::mutate(
    temp_plot = ifelse(temp_deg_c == "18-14", 16, 20)
  ) %>%
  dplyr::mutate(
    temp_label = ifelse(temp_plot == 16, "La Niña, 600 μatm", "El Niño, 600 μatm")
  ) %>%
  dplyr::mutate(temp_label = factor(temp_label, levels = c("La Niña, 600 μatm", "El Niño, 600 μatm")))

ce_80_skel_dynamic_temp <- conditional_effects(dynamic_mod, 
                                               robust = F, prob = 0.80)

df_80_skel_dynamic_plot <- ce_80_skel_dynamic_temp$temp %>%
  filter(temp_deg_c %in% c("21-18", "18-14")) %>%
  dplyr::mutate(
    temp_plot = ifelse(temp_deg_c == "18-14", 16, 20)
  ) %>%
  dplyr::mutate(
    temp_label = ifelse(temp_plot == 16, "La Niña, 600 μatm", "El Niño, 600 μatm")
  ) %>%
  dplyr::mutate(temp_label = factor(temp_label, levels = c("La Niña, 600 μatm", "El Niño, 600 μatm")))

empirical_means_skel <- data_skel_growth_const_temp  %>%
  group_by(tank_id, temp_num, pco2_uatm) %>%
  dplyr::summarise(
    emp_mean_skel = mean(delta_demipyramid_um)
  ) %>%
  dplyr::filter(!is.na(emp_mean_skel)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    pco2_uatm = factor(pco2_uatm, levels = c("600 μatm", "1200 μatm"))) %>%
  dplyr::select(tank_id, temp_num, pco2_uatm, emp_mean_skel)

empirical_means_skel_dyn <- data_skel_growth_dynamic_temp %>%
  group_by(tank_id, temp_num, pco2_uatm) %>%
  dplyr::summarise(
    emp_mean_skel = mean(delta_demipyramid_um)
  ) %>%
  dplyr::filter(!is.na(emp_mean_skel)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    pco2_uatm = ifelse(temp_num == 16, "La Niña, 600 μatm", "El Niño, 600 μatm")
  ) %>%
  dplyr::select(tank_id, temp_num, pco2_uatm, emp_mean_skel)

df_empirical_means <- data.frame()
df_empirical_means <- rbind(empirical_means_skel, empirical_means_skel_dyn) #%>%
# dplyr::mutate(pco2_uatm = factor(pco2_uatm, levels = c("600 μatm", "1200 μatm"))) %>%
# filter(!is.na(emp_mean_skel))

p_skel_pco2_temp <-  ggplot() +
  geom_line(data = df_95_skel_constant_plot,
            aes(x = temp_num, y = estimate__,
                color = pco2_uatm, linetype = pco2_uatm),
            linewidth = 1.5, show.legend = F) +
  geom_ribbon(data = df_95_skel_constant_plot,
              aes(x = temp_num, y = estimate__,
                  ymin = lower__, ymax = upper__,
                  fill = pco2_uatm),
              alpha = 0.2) +
  geom_ribbon(data = df_80_skel_constant_plot,
              aes(x = temp_num, y = estimate__,
                  ymin = lower__, ymax = upper__,
                  fill = pco2_uatm),
              alpha = 0.4) +
  geom_linerange(data = df_95_skel_dynamic_plot ,
                 aes(x = temp_plot, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = temp_label),
                 linewidth = 1) +
  geom_linerange(data = df_80_skel_dynamic_plot ,
                 aes(x = temp_plot, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = temp_label),
                 linewidth = 3) +
  geom_jitter(data = df_empirical_means,
              aes(x = temp_num, y = emp_mean_skel,
                  shape = pco2_uatm, fill = pco2_uatm),
              size = 4, color = "black", width = 0.1) +
  scale_x_continuous(breaks = c(10,13,16,17,18,20)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_color_manual(values = c(
    "blue", # 600 uatm
    "red", # 1200 uatm
    "ivory4", # La Nina
    "gray20" # El Nino
  )) +
  scale_shape_manual(values = c(
    21, # 600 uatm
    24, # 1200 uatm
    22, # La Nina
    23 # El Nino
  )) +
  scale_fill_manual(values = c(
    "blue", # 600 uatm
    "red", # 1200 uatm
    "ivory4", # La Nina
    "gray20" # El Nino
  )) +
  scale_linetype_manual(values = c(
    1, # 600 uatm 
    3 # 1200 uatm
  )) +
  gg_options() +
  labs(
    x = "temperature (°C)",
    y = expression("demipyramid growth (μm)"),
    shape = NULL,
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.82,0.9)) 

p_skel_pco2_temp

ggsave("./Figures/Skeleton_Jaw.pdf", p_skel_pco2_temp,
       units = "in", width = 8, height = 8)
ggsave("./Figures/Skeleton_Jaw.png", p_skel_pco2_temp,
       units = "in", width = 8, height = 8)
ggsave("./Figures/Skeleton_Jaw.jpg", p_skel_pco2_temp,
       units = "in", width = 8, height = 8)

# Post hoc on full model for constant temps -----
fit_demi_const_temp_oa <- readRDS("./Model_Output/demipyramid_temp_oa_full.rds")
selected_model <- fit_demi_const_temp_oa
summary(selected_model)
# Family: gamma 
# Links: mu = log 
# Formula: delta_demipyramid_um ~ log(initial_volume_ml) + gp(temp_num, k = 5, c = 5/4, by = pco2_uatm) + (1 | tank_id) 
# Data: d (Number of observations: 311) 
# Draws: 4 chains, each with iter = 18000; warmup = 2000; thin = 1;
# total post-warmup draws = 64000
# 
# Gaussian Process Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sdgp(gptemp_numpco2_uatm1200μatm)       0.28      0.12     0.05     0.50 1.00    19804    25799
# sdgp(gptemp_numpco2_uatm600μatm)        0.21      0.12     0.01     0.45 1.00    19715    33059
# lscale(gptemp_numpco2_uatm1200μatm)     0.27      0.12     0.12     0.56 1.00    57380    39506
# lscale(gptemp_numpco2_uatm600μatm)      0.39      0.24     0.12     1.06 1.00    41529    42927
# 
# Multilevel Hyperparameters:
#   ~tank_id (Number of levels: 24) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.07      0.04     0.00     0.17 1.00    20335    28941
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                6.17      0.32     5.54     6.77 1.00    40218    47613
# loginitial_volume_ml    -0.42      0.08    -0.57    -0.27 1.00    87569    47475
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# shape     4.72      0.37     4.02     5.48 1.00   100358    45615
# 
# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
ems <- emmeans(selected_model, ~ pco2_uatm | temp_num,
               at = list(temp_num = c(10,13,16,17,18,20)))

cont_pco2 <- contrast(ems, method = "revpairwise", by = "pco2_uatm", 
                      type = "response")

emmeans::hpd.summary(cont_pco2, 
                     prob = 0.95, 
                     point.est = mean)
# pco2_uatm = 1200 μatm:
#   contrast                ratio lower.HPD upper.HPD
# temp_num13 / temp_num10 1.127     0.908      1.36
# temp_num16 / temp_num10 1.155     0.890      1.45
# temp_num16 / temp_num13 1.026     0.865      1.20
# temp_num17 / temp_num10 1.071     0.845      1.32
# temp_num17 / temp_num13 0.953     0.780      1.13
# temp_num17 / temp_num16 0.929     0.862      1.00
# temp_num18 / temp_num10 0.974     0.771      1.18
# temp_num18 / temp_num13 0.869     0.685      1.05
# temp_num18 / temp_num16 0.848     0.724      1.00
# temp_num18 / temp_num17 0.911     0.836      1.00
# temp_num20 / temp_num10 0.860     0.653      1.06
# temp_num20 / temp_num13 0.768     0.583      1.01
# temp_num20 / temp_num16 0.752     0.566      1.00
# temp_num20 / temp_num17 0.807     0.639      1.01
# temp_num20 / temp_num18 0.883     0.749      1.02
# 
# pco2_uatm = 600 μatm:
#   contrast                ratio lower.HPD upper.HPD
# temp_num13 / temp_num10 1.121     0.948      1.37
# temp_num16 / temp_num10 1.127     0.908      1.41
# temp_num16 / temp_num13 1.006     0.873      1.14
# temp_num17 / temp_num10 1.102     0.891      1.36
# temp_num17 / temp_num13 0.985     0.827      1.14
# temp_num17 / temp_num16 0.979     0.922      1.03
# temp_num18 / temp_num10 1.072     0.867      1.31
# temp_num18 / temp_num13 0.960     0.776      1.13
# temp_num18 / temp_num16 0.954     0.841      1.05
# temp_num18 / temp_num17 0.974     0.910      1.02
# temp_num20 / temp_num10 1.017     0.799      1.25
# temp_num20 / temp_num13 0.913     0.693      1.10
# temp_num20 / temp_num16 0.908     0.710      1.07
# temp_num20 / temp_num17 0.927     0.763      1.07
# temp_num20 / temp_num18 0.950     0.825      1.05
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95 

cont_posterior_pco2 <- gather_emmeans_draws(cont_pco2) 

cont_temp <- contrast(ems, method = "revpairwise", by = "temp_num", 
                      type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_temp, 
                     prob = 0.95, 
                     point.est = mean)
# temp_num = 10:
#   contrast             ratio lower.HPD upper.HPD
# 600 μatm / 1200 μatm  1.50      1.13      1.88
# 
# temp_num = 13:
#   contrast             ratio lower.HPD upper.HPD
# 600 μatm / 1200 μatm  1.49      1.17      1.83
# 
# temp_num = 16:
#   contrast             ratio lower.HPD upper.HPD
# 600 μatm / 1200 μatm  1.46      1.19      1.75
# 
# temp_num = 17:
#   contrast             ratio lower.HPD upper.HPD
# 600 μatm / 1200 μatm  1.54      1.29      1.79
# 
# temp_num = 18:
#   contrast             ratio lower.HPD upper.HPD
# 600 μatm / 1200 μatm  1.65      1.37      1.93
# 
# temp_num = 20:
#   contrast             ratio lower.HPD upper.HPD
# 600 μatm / 1200 μatm  1.79      1.32      2.28
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95

cont_posterior_temp <- gather_emmeans_draws(cont_temp) %>%
  mutate(
    temp_num = as.factor(temp_num)
  )

ems_overall <- emmeans(selected_model, ~ pco2_uatm | temp_num)

cont_overall <- contrast(ems_overall, method = "revpairwise", 
                         type = "response", ratios = TRUE)

cont_posterior_overall <- gather_emmeans_draws(cont_overall) %>%
  mutate(
    temp_num = as.factor("mean")
  ) 

cont_pairwise_and_overall_posterior_pco2 <- data.frame()
cont_pairwise_and_overall_posterior_pco2 <- rbind(cont_posterior_pco2, 
                                                  cont_posterior_overall)


summary_emm_const_pco2 <- cont_pairwise_and_overall_posterior_pco2 %>%
  group_by(contrast, pco2_uatm) %>%
  mean_hdi(.value, .width = c(0.95)) %>%
  mutate(
    .value = round(.value, digits = 2),
    .lower = round(.lower, digits = 2),
    .upper = round(.upper, digits = 2)
  ) %>%
  mutate(
    HDI = paste(.value, " (", .lower, ", ", .upper, ")", sep = "")
  )

cont_pairwise_and_overall_posterior_temp <- data.frame()
cont_pairwise_and_overall_posterior_temp <- rbind(cont_posterior_temp, 
                                                  cont_posterior_overall)

summary_emm_const_temp <- cont_pairwise_and_overall_posterior_temp %>%
  group_by(contrast, temp_num) %>%
  mean_hdi(.value, .width = c(0.95)) %>%
  mutate(
    .value = round(.value, digits = 2),
    .lower = round(.lower, digits = 2),
    .upper = round(.upper, digits = 2)
  ) %>%
  mutate(
    HDI = paste(.value, " (", .lower, ", ", .upper, ")", sep = "")
  )

plot_emm_fit <- ggplot(cont_pairwise_and_overall_posterior_temp, 
                       aes(x = .value, y = contrast)) +
  geom_vline(xintercept = 1, linetype = 2, color = "black", linewidth = 1, 
             alpha = 0.7) +
  stat_pointinterval(point_interval = "mean_hdi", .width = c(0.95, 0.80),
                     interval_size_range = c(1,2), point_size = 4) +
  geom_label_repel(data = summary_emm_const_temp,
                   aes(x = .value, label = HDI), 
                   fill = "white", nudge_x = -0.4, nudge_y = 1) +
  facet_grid(rows = vars(temp_num), switch = "y") +
  labs(
    x = "ratio",
    y = "temperature (°C)",
    title = "Jaw Growth",
    subtitle = expression("pCO"[2]*" effect (1200 µatm vs 600 µatm)")
  ) +
  gg_options() +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

plot_emm_fit

ggsave("./Figures/Supplemental/PWC_Skel_Jaw_pco2.jpg", plot_emm_fit,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Skel_Jaw_pco2.pdf", plot_emm_fit,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Skel_Jaw_pco2.png", plot_emm_fit,
       width = 6, height = 8, units = "in")

# Post hoc constant vs dynamic temp comparisons --------------------------------------------------------

fit_demi_const_dynamic_temp <- readRDS("./Model_Output/demipyramid_temp_const_dynamic_full.rds")
selected_model <- fit_demi_const_dynamic_temp
summary(selected_model)

ems <- emmeans(selected_model, ~ temp_deg_c)

cont <- contrast(ems, method = "revpairwise", by = NULL, 
                 type = "response", ratios = TRUE)

emmeans::hpd.summary(cont, 
                     prob = 0.95, 
                     point.est = mean)
# contrast          ratio lower.HPD upper.HPD
# 13 / 10           1.254     0.895     1.628
# 16 / 10           1.205     0.869     1.569
# 16 / 13           0.972     0.706     1.258
# 17 / 10           1.265     0.905     1.644
# 17 / 13           1.020     0.736     1.318
# 17 / 16           1.060     0.771     1.373
# 18 / 10           0.966     0.694     1.254
# 18 / 13           0.779     0.562     1.015
# 18 / 16           0.810     0.588     1.052
# 18 / 17           0.772     0.552     0.999
# (18-14) / 10      1.011     0.758     1.274
# (18-14) / 13      0.815     0.616     1.028
# (18-14) / 16      0.847     0.639     1.061
# (18-14) / 17      0.808     0.607     1.013
# (18-14) / 18      1.058     0.793     1.330
# 20 / 10           1.135     0.801     1.496
# 20 / 13           0.916     0.646     1.198
# 20 / 16           0.952     0.674     1.244
# 20 / 17           0.907     0.646     1.194
# 20 / 18           1.188     0.837     1.559
# 20 / (18-14)      1.130     0.837     1.444
# (21-18) / 10      0.691     0.517     0.873
# (21-18) / 13      0.558     0.418     0.702
# (21-18) / 16      0.579     0.436     0.727
# (21-18) / 17      0.552     0.416     0.695
# (21-18) / 18      0.723     0.547     0.912
# (21-18) / (18-14) 0.688     0.541     0.832
# (21-18) / 20      0.617     0.458     0.792
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95 

cont_posterior <- gather_emmeans_draws(cont) 

cont_posterior_filtered <- cont_posterior %>%
  filter(
    contrast %in% c(
      "(18-14) - 16",
      "(21-18) - 20"
    )
  ) %>%
  mutate(
    contrast_formatted = ifelse(contrast == "(18-14) - 16",
                                "La Niña vs 16°C", "El Niño vs 20°C")
  )

summary_emm_mhw <- cont_posterior_filtered %>%
  group_by(contrast_formatted) %>%
  mean_hdi(.value, .width = c(0.95)) %>%
  mutate(
    .value = round(.value, digits = 2),
    .lower = round(.lower, digits = 2),
    .upper = round(.upper, digits = 2)
  ) %>%
  mutate(
    HDI = paste(.value, " (", .lower, ", ", .upper, ")", sep = "")
  )

plot_emm_fit_tpc_mhw <- ggplot(cont_posterior_filtered, 
                               aes(x = .value, y = contrast_formatted)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black", size = 1, 
             alpha = 0.7) +
  stat_pointinterval(point_interval = "mean_hdi", .width = c(0.95, 0.80),
                     interval_size_range = c(1,2), point_size = 4) +
  geom_label_repel(data = summary_emm_mhw,
                   aes(x = .value, label = HDI), 
                   fill = "white", nudge_y = -0.4, nudge_x = 0.2) +
  facet_grid(rows = vars(contrast_formatted), switch = "y") +
  labs(
    x = NULL,
    y = NULL,
    title = "Jaw Growth",
    subtitle = "constant mean vs dynamic"
  ) +
  gg_options() +
  scale_x_continuous(breaks = pretty_breaks(n = 6), limits = c(-1.0,0.4)) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = c("#DEDEDE"))
  )

plot_emm_fit_tpc_mhw

ggsave("./Figures/Supplemental/PWC_Skel_Jaw_const_dynamic_temp.jpg", plot_emm_fit_tpc_mhw,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Skel_Jaw_const_dynamic_temp.pdf", plot_emm_fit_tpc_mhw,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Skel_Jaw_const_dynamic_temp.png", plot_emm_fit_tpc_mhw,
       width = 6, height = 8, units = "in")