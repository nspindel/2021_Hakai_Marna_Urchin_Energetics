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
  "ggrepel",
  "scales",
  "Hmisc"
)

lapply(libs, require, character.only = TRUE)

source("Utility_Functions.R")

# Load/wrangle data -----------------------------------------------------------
data_demipyramid <- read.csv("./Data/Demipyramid_Growth.csv") %>%
  mutate(
    tank_id = as.factor(tank),
    pit_tag_id = as.factor(pit_tag),
    pco2_uatm = factor(pco2, levels = c("600", "1200"), labels = c("600 μatm", "1200 μatm")),
    start_diam_mm = as.numeric(start_diam_mm),
    start_ht_mm = as.numeric(start_ht_mm),
    end_diam_mm = as.numeric(end_diam_mm),
    end_ht_mm = as.numeric(end_ht_mm)
  ) %>%
  mutate(
    demi_growth_total_um = demi_final_um - demi_tag_um,
    demi_growth_per_day_um = demi_growth_total_um / total_days,
    start.semi.major.axis.mm = 0.5 * start_diam_mm,
    start.semi.minor.axis.mm = 0.5 * start_ht_mm,
    start.spheroid.volume.ml = (4/3 * pi * start.semi.major.axis.mm ^2 * start.semi.minor.axis.mm)/1000, # # Equation for volume of oblate spheroid = 4/3pia^2c in mL
    end.semi.major.axis.mm = 0.5 * end_diam_mm,
    end.semi.minor.axis.mm = 0.5 * end_ht_mm,
    end.spheroid.volume.ml = (4/3 * pi * end.semi.major.axis.mm ^2 * end.semi.minor.axis.mm)/1000, # # Equation for volume of oblate spheroid = 4/3pia^2c in mL
  ) %>% 
  mutate(
    start_vol_ml = start.spheroid.volume.ml,
    end_vol_ml = end.spheroid.volume.ml
  ) %>%
  mutate(
    delta_diam_mm = end_diam_mm-start_diam_mm,
    delta_ht_mm = end_ht_mm-start_ht_mm,
    delta_vol_ml = end_vol_ml-start_vol_ml
  ) %>%
  mutate(
    diam_change_per_day_mm = delta_diam_mm/total_days,
    ht_change_per_day_mm = delta_ht_mm/total_days,
    vol_change_per_day_ml = delta_vol_ml/total_days,
  ) %>%
  filter(
    between(demi_growth_per_day_um,0.1,5)
  ) %>%
  filter(abs(delta_diam_mm)<10) %>%
  mutate(
    temp_model = as.numeric(case_when(
      temp == 10 ~ 10,
      temp == 13 ~ 13,
      temp == 16 ~ 16,
      temp == 17 ~ 17,
      temp == 18 ~ 18,
      temp == 20 ~ 20,
      temp == "18-14" ~ 16,
      temp == "21-18" ~ 20
    ))
  )

data_demipyramid_const_temp_oa <- data_demipyramid %>%
  filter(
    temp != "21-18" 
    & temp != "18-14"
    & !is.na(demi_growth_per_day_um)
  ) %>%
  select(
    tank_id,
    pit_tag_id,
    start_diam_mm,
    temp,
    temp_model,
    pco2_uatm,
    demi_growth_per_day_um
  )

data_demipyramid_const_dynamic_temp <- data_demipyramid %>%
  filter(
    pco2_uatm == "600 μatm"
    & !is.na(demi_growth_per_day_um)
  ) %>% 
  mutate(
    temp_label = factor(case_when(
      temp == "18-14" ~ "La Niña (18-14 °C)",
      temp == "21-18" ~ "El Niño (21-18 °C)"
    ), 
    levels = c("La Niña (18-14 °C)", "El Niño (21-18 °C)"),
    labels = c("La Niña, 600 μatm", "El Niño, 600 μatm"))
  ) %>%
  select(
    tank_id,
    pit_tag_id,
    start_diam_mm,
    temp,
    temp_label,
    pco2_uatm,
    demi_growth_per_day_um
  )

data_demipyramid_dynamic_temp <- data_demipyramid_const_dynamic_temp %>%
  filter(temp %in% c("21-18", "18-14"))

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

data_demipyramid <- read.csv("./Data/Demipyramid_Growth.csv") %>%
  mutate(
    tank_id = as.factor(tank_id),
    pit_tag_id = as.factor(pit_tag),
    pco2_uatm = factor(pco2_uatm, levels = c("600", "1200"), labels = c("600 μatm", "1200 μatm"))
  ) %>%
  select(
    tank_id,
    temp_deg_c,
    pco2_uatm,
    pit_tag_id,
    demi_tag_um,
    demi_final_um
  )

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

unlink("./Model_Output/demipyramid_temp_oa_*", force = TRUE)


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
# Full model
form_full <- brms::bf(
  demi_growth_per_day_um ~
    log(start_diam_mm) +
    gp(temp_model, k = 5, c = 5/4, by = pco2_uatm) +
    (1 | tank_id)
)

# Reduced models chosen to test specific hypotheses with principled nesting
form_gp_noby <- brms::bf(
  demi_growth_per_day_um ~
    log(start_diam_mm) +
    gp(temp_model, k = 5, c = 5/4) +
    (1 | tank_id)
)
# (Removes the pCO2-specific GP surfaces; retains nonparametric thermal response.)

form_lin_int <- brms::bf(
  demi_growth_per_day_um ~
    log(start_diam_mm) + temp_model * pco2_uatm +
    (1 | tank_id)
)
# (Parametric T × pCO2 interaction; no GP.)

form_lin_add <- brms::bf(
  demi_growth_per_day_um ~
    log(start_diam_mm) + temp_model + pco2_uatm +
    (1 | tank_id)
)
# (Additive parametric effects; no GP.)

form_size_only <- brms::bf(
  demi_growth_per_day_um ~
    log(start_diam_mm) + (1 | tank_id)
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
d <- data_skel_const_temp_oa %>%
  select(demi_growth_per_day_um, start_diam_mm, temp_model, pco2_uatm, tank_id) %>%
  mutate(temp_model = as.numeric(temp_model)) %>%
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

brms_formula_full_dyn <- bf(demi_growth_per_day_um ~ 
                              log(start_diam_mm) + 
                              temp_model + 
                              (1|tank_id))

# fit model
fit_demipyramid_temp_dynamic_temp <- brm(
  formula = brms_formula_full_dyn, 
  data = data_fce_dynamic_temp,
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

brms_formula_full_const_dyn <- bf(fce ~
                                    log(init_diam_mm) +
                                    temp_deg_c +
                                    (1|tank_id))

# fit model
fit_fce_const_temp_dynamic_temp <- brm(
  formula = brms_formula_full_const_dyn,
  data = data_fce_const_dynamic_temp,
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
  file = "./Model_Output/fce_temp_const_dynamic_full"
)


# Model skeletal growth -----
brms_fam <- Gamma(link = "log")

brms_formula_full <- bf(demi_growth_per_day_um ~ 
                          log(start_diam_mm) * temp_model * pco2_uatm +
                          gp(temp_model, k = 5, c = 5/4, by = pco2_uatm) + 
                          (1|tank_id))

# default_priors <- default_prior(brms_formula_full, data = data_skel_const_temp_oa, family = brms_fam)
priors <- c( 
  prior(normal(0,1), class = "b"),
  prior(student_t(3, 0.3, 2.5), class = "Intercept"),
  prior(inv_gamma(4.6, 1.1), lb = 0, class = "lscale"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sdgp"),
  prior(gamma(0.01,0.01),lb = 0, class = "shape")
)

fit_skel_const_temp_oa <- brm(
  formula = brms_formula_full, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_oa_full"
)

brms_formula_2 <- bf(demi_growth_per_day_um ~ 
                       log(start_diam_mm) + 
                       temp_model * pco2_uatm +
                       gp(temp_model, k = 5, c = 5/4, by = pco2_uatm) + 
                       (1|tank_id))

# fit model
fit_skel_const_temp_oa_2 <- brm(
  formula = brms_formula_2, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_oa_2"
)

brms_formula_3 <- bf(demi_growth_per_day_um ~ 
                       log(start_diam_mm) + 
                       temp_model + 
                       pco2_uatm +
                       gp(temp_model, k = 5, c = 5/4, by = pco2_uatm) + 
                       (1|tank_id))

# fit model
fit_skel_const_temp_oa_3 <- brm(
  formula = brms_formula_3, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_oa_3"
)

brms_formula_4 <- bf(demi_growth_per_day_um ~ 
                       log(start_diam_mm) * temp_model + 
                       pco2_uatm +
                       gp(temp_model, k = 5, c = 5/4, by = pco2_uatm) + 
                       (1|tank_id))

# fit model
fit_skel_const_temp_oa_4 <- brm(
  formula = brms_formula_4, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_oa_4"
)

brms_formula_5 <- bf(demi_growth_per_day_um ~ 
                       log(start_diam_mm) *  pco2_uatm +
                       temp_model + 
                       gp(temp_model, k = 5, c = 5/4, by = pco2_uatm) + 
                       (1|tank_id))

# fit model
fit_skel_const_temp_oa_5 <- brm(
  formula = brms_formula_5, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_oa_5"
)

brms_formula_6 <- bf(demi_growth_per_day_um ~ 
                       log(start_diam_mm) * temp_model +
                        gp(temp_model, k = 5, c = 5/4) + 
                       (1|tank_id))

# fit model
fit_skel_const_temp_oa_6 <- brm(
  formula = brms_formula_6, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_oa_6"
)

brms_formula_7 <- bf(demi_growth_per_day_um ~ 
                       log(start_diam_mm) * pco2_uatm + 
                       (1|tank_id))
priors2 <- c( 
  prior(normal(0,1), class = "b"),
  prior(student_t(3, 0.3, 2.5), class = "Intercept"),
  # prior(inv_gamma(4.6, 1.1), lb = 0, class = "lscale"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  # prior(student_t(3, 0, 2.5), lb = 0, class = "sdgp"),
  prior(gamma(0.01,0.01),lb = 0, class = "shape")
)


# fit model
fit_skel_const_temp_oa_7 <- brm(
  formula = brms_formula_7, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors2,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_oa_7"
)

brms_formula_8 <- bf(demi_growth_per_day_um ~ 
                       log(start_diam_mm) + 
                         (1|tank_id))

# fit model
fit_skel_const_temp_oa_8 <- brm(
  formula = brms_formula_8, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors2,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_oa_8"
)

brms_formula_9 <- bf(demi_growth_per_day_um ~ 
                       (1|tank_id))

priors3 <- c( 
  # prior(normal(0,1), class = "b"),
  prior(student_t(3, 0.3, 2.5), class = "Intercept"),
  # prior(inv_gamma(4.6, 1.1), lb = 0, class = "lscale"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  # prior(student_t(3, 0, 2.5), lb = 0, class = "sdgp"),
  prior(gamma(0.01,0.01),lb = 0, class = "shape")
)


# fit model
fit_skel_const_temp_oa_9 <- brm(
  formula = brms_formula_9, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors3,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_oa_9"
)

# Model selection via LOO stacking weight -------------------------------------
# Compute model weights via stacking of predictive distributions
# 
# Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian 
# model evaluation using leave-one-out cross-validation and WAIC. 
# Statistics and Computing. 27(5), 1413–1432. 
# doi:10.1007/s11222-016-9696-4 (journal version, preprint arXiv:1507.04544).
#
# Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2019). 
# Pareto smoothed importance sampling. preprint arXiv:1507.02646
#
# Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018) Using 
# stacking to average Bayesian predictive distributions. Bayesian 
# Analysis, advance publication, doi:10.1214/17-BA1091. (online).

loo_model_weights(
  fit_skel_const_temp_oa,
  fit_skel_const_temp_oa_2,
  fit_skel_const_temp_oa_3,
  fit_skel_const_temp_oa_4,
  fit_skel_const_temp_oa_5,
  fit_skel_const_temp_oa_6,
  fit_skel_const_temp_oa_7,
  fit_skel_const_temp_oa_8,
  fit_skel_const_temp_oa_9
)

# 2024-12-13
# Method: stacking
# 
#                           weight
# fit_skel_const_temp_oa   0.069 
# fit_skel_const_temp_oa_2 0.016 
# fit_skel_const_temp_oa_3 0.008 
# fit_skel_const_temp_oa_4 0.000 
# fit_skel_const_temp_oa_5 0.506 <- favored
# fit_skel_const_temp_oa_6 0.000 
# fit_skel_const_temp_oa_7 0.401 
# fit_skel_const_temp_oa_8 0.000 
# fit_skel_const_temp_oa_9 0.000

selected_model <- fit_skel_const_temp_oa_5

# Output summary stats for final model
summary(selected_model)

# Calculate Bayesian R-squared for final model
bayes_R2(selected_model)
#     Estimate  Est.Error      Q2.5     Q97.5
# R2 0.3327557 0.04240869 0.2498381 0.4138415

# Model diagnostics
# Pareto-k-diagnostic
plot(brms::loo(selected_model))

# Predictive performance loo
brms::loo(selected_model)

# Posterior predictive checks -----------------------------------------------------------------
ndraws = 200

ppc_cpue_dens_overlay <- pp_check(selected_model, ndraws = ndraws) + 
  xlim(0,20) +
  gg_options();ppc_cpue_dens_overlay

ppc_cpue_loo_ribbon <- pp_check(selected_model, ndraws = ndraws, 
                                type = "loo_ribbon") + 
  gg_options();ppc_cpue_loo_ribbon

ppc_cpue_loo_pit_qq <- pp_check(selected_model, ndraws = ndraws, 
                                type = "loo_pit_qq") + 
  gg_options();ppc_cpue_loo_pit_qq

ppc_multiplot <- cowplot::plot_grid(
  ppc_cpue_dens_overlay +
    theme(legend.position = c(0.8,0.85)),
  ppc_cpue_loo_ribbon +
    theme(legend.position = c(0.8,0.85)),
  ppc_cpue_loo_pit_qq,
  ncol = 3,
  align = "hv",
  axis = "tblr",
  labels = "AUTO"
);ppc_multiplot

ggsave("./PP_Checks/Skeletal_Growth_Const_Temp_OA_ppc.pdf", ppc_multiplot,
       width = 16, height = 6, units = "in")

# Sensitivity to priors -------------------------------------------------------
informative_priors <- c( 
  prior(normal(0,0.1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(inv_gamma(4.6, 1.1), lb = 0, class = "lscale"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sdgp"),
  prior(gamma(3, 1),lb = 0, class = "shape")
)

fit_informative <- brm(
  formula = brms_formula_full,
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = informative_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "C:/Users/nates/Documents/Marna_Urchin_Heatwaves_2021/Marna_Urchin_Heatwaves_2021/Model_Output/skel_temp_oa_full_PRIOR_INFORMATIVE"
)

vague_priors <- c( 
  prior(normal(0,1), class = "b"),
  prior(student_t(3, 1.1, 2.5), class = "Intercept"),
  prior(inv_gamma(4.6, 1.1), lb = 0, class = "lscale"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sdgp"),
  prior(gamma(3, 1),lb = 0, class = "shape")
)

fit_vague <- brm(
  formula = brms_formula_full,
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = vague_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "C:/Users/nates/Documents/Marna_Urchin_Heatwaves_2021/Marna_Urchin_Heatwaves_2021/Model_Output/skel_temp_oa_full_PRIOR_VAGUE"
)

diffuse_priors <- c( 
  prior(normal(0,10), class = "b"),
  prior(student_t(3, 1.1, 2.5), class = "Intercept"),
  prior(inv_gamma(4.6, 1.1), lb = 0, class = "lscale"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sdgp"),
  prior(gamma(2, 2),lb = 0, class = "shape")
)

brm(
  formula = brms_formula_full,
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = diffuse_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "C:/Users/nates/Documents/Marna_Urchin_Heatwaves_2021/Marna_Urchin_Heatwaves_2021/Model_Output/skel_temp_oa_full_PRIOR_DIFFUSE"
)

pp_check(fit_informative)
pp_check(fit_vague)
pp_check(fit_diffuse)

# Hypothesis tests ----
fit_skel_const_temp_oa_5 <- readRDS("./Model_Output/skeleton_temp_oa_5.rds")
selected_model <- fit_skel_const_temp_oa_5
selected_model$formula

# 1) body size
brms_formula_10 <- bf(demi_growth_per_day_um ~ 
                       pco2_uatm +
                       temp_model + 
                       gp(temp_model, k = 5, c = 5/4, by = pco2_uatm) + 
                       (1|tank_id))

# fit model
fit_skel_const_temp_oa_10 <- brm(
  formula = brms_formula_10, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", 
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_oa_10"
)

bayes_factor(selected_model, fit_skel_const_temp_oa_10)
# Estimated Bayes factor in favor of selected_model over fit_skel_const_temp_oa_10: 6605844.2944
# BF >> 100, decisive evidence in support of hypothesis that body size shapes outcomes

# 2) temperature
brms_formula_11 <- bf(demi_growth_per_day_um ~ 
                        log(start_diam_mm) * pco2_uatm +
                        # temp_model + 
                        # gp(temp_model, k = 5, c = 5/4, by = pco2_uatm) + 
                        (1|tank_id))

# fit model
fit_skel_const_temp_oa_11 <- brm(
  formula = brms_formula_11, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors2,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", 
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_oa_11"
)

bayes_factor(selected_model, fit_skel_const_temp_oa_11)
# Estimated Bayes factor in favor of selected_model over fit_skel_const_temp_oa_11: 0.00043
# BF < 1, no evidence in support of hypothesis that temperature shapes outcomes

# 3) pco2
fit_skel_const_temp_oa_6$formula

bayes_factor(selected_model, fit_skel_const_temp_oa_6)
# Estimated Bayes factor in favor of selected_model over fit_skel_const_temp_oa_6: 92490.44291
# BF >> 100, decisive evidence in support of hypothesis that pco2 shapes outcomes

# Predictions -----------------------------------------------------------------
fit_skel_const_temp_oa_5 <- readRDS("./Model_Output/skeleton_temp_oa_5.rds")
selected_model <- fit_skel_const_temp_oa_5
selected_model$formula

ce_95_pco2_comparison <- conditional_effects(selected_model, robust = F)
int_conditions <- list(temp_model = c(10,13,16,17,18,20))
foo <- conditional_effects(selected_model, robust = F,
                                             int_conditions = int_conditions,
                                             effects = "temp_model:pco2_uatm")
df_foo <- foo$`temp_model:pco2_uatm`


df_95_skel_constant_plot <- ce_95_pco2_comparison$`temp_model:pco2_uatm` 

ce_80_pco2_comparison <- conditional_effects(selected_model, robust = F, prob = 0.80)

df_80_skel_constant_plot <- ce_80_pco2_comparison$`temp_model:pco2_uatm` 

dynamic_mod <- skeleton_temp_const_dynamic_2 <- readRDS("./Model_Output/skeleton_temp_const_dynamic_2.rds")
dynamic_mod$formula

ce_95_skel_dynamic_temp <- conditional_effects(dynamic_mod, 
                                               robust = F, prob = 0.95)

df_95_skel_dynamic_plot <- ce_95_skel_dynamic_temp$temp %>%
  dplyr::mutate(
    temp_plot = ifelse(temp == "18-14", 16, 20)
  ) %>%
  dplyr::mutate(
    temp_label = ifelse(temp_plot == 16, "La Niña, 600 μatm", "El Niño, 600 μatm")
  ) %>%
  dplyr::mutate(temp_label = factor(temp_label, levels = c("La Niña, 600 μatm", "El Niño, 600 μatm")))


ce_80_skel_dynamic_temp <- conditional_effects(dynamic_mod, 
                                               robust = F, prob = 0.80)

df_80_skel_dynamic_plot <- ce_80_skel_dynamic_temp$temp %>%
  dplyr::mutate(
    temp_plot = ifelse(temp == "18-14", 16, 20)
  ) %>%
  dplyr::mutate(
    temp_label = ifelse(temp_plot == 16, "La Niña, 600 μatm", "El Niño, 600 μatm")
  ) %>%
  dplyr::mutate(temp_label = factor(temp_label, levels = c("La Niña, 600 μatm", "El Niño, 600 μatm")))

empirical_means_skel <- data_skel_const_temp_oa  %>%
  group_by(tank_id, temp_model, pco2_uatm) %>%
  dplyr::summarise(
    emp_mean_skel = mean(demi_growth_per_day_um)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(temp_plot = temp_model) %>%
  dplyr::select(tank_id, temp_plot, pco2_uatm, emp_mean_skel)

empirical_means_skel_dyn <- data_skel_dynamic_temp %>%
  group_by(tank_id, temp, temp_label, pco2_uatm) %>%
  dplyr::summarise(
    emp_mean_skel = mean(demi_growth_per_day_um)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    temp_plot = ifelse(temp_label == "La Niña, 600 μatm", 16, 20),
    pco2_uatm = temp_label
  ) %>%
  dplyr::select(tank_id, temp_plot, pco2_uatm, emp_mean_skel)

df_empirical_means <- data.frame()
df_empirical_means <- rbind(empirical_means_skel, empirical_means_skel_dyn) 

p_skel_pco2_temp <-  ggplot() +
  geom_line(data = df_95_skel_constant_plot,
            aes(x = temp_model, y = estimate__,
                color = pco2_uatm, linetype = pco2_uatm),
            linewidth = 1.5, show.legend = F) +
  geom_ribbon(data = df_95_skel_constant_plot,
              aes(x = temp_model, y = estimate__,
                  ymin = lower__, ymax = upper__,
                  fill = pco2_uatm),
              alpha = 0.2) +
  geom_ribbon(data = df_80_skel_constant_plot,
              aes(x = temp_model, y = estimate__,
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
              aes(x = temp_plot, y = emp_mean_skel,
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
    y = expression("demipyramid growth (μm day"^-1*")"),
    shape = NULL,
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.8,0.9)); p_skel_pco2_temp

ggsave("./Figures/Skeleton.pdf", p_skel_pco2_temp,
       units = "in", width = 8, height = 8)
ggsave("./Figures/Skeleton.png", p_skel_pco2_temp,
       units = "in", width = 8, height = 8)
ggsave("./Figures/Skeleton.jpg", p_skel_pco2_temp,
       units = "in", width = 8, height = 8)

# Pairwise comparisons --------------------------------------------------------
selected_model <- fit_skel_const_temp_oa_5 <- readRDS("./Model_Output/skeleton_temp_oa_5.rds")
selected_model$formula


ems <- emmeans(selected_model, ~ temp_model + pco2_uatm | log(start_diam_mm),
               at = list(temp_model = c(10,13,16,17,18,20)))

cont_pco2 <- contrast(ems, method = "revpairwise", by = "pco2_uatm", 
                      type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_pco2, 
                     prob = 0.95, 
                     point.est = mean)

cont_posterior_pco2 <- gather_emmeans_draws(cont_pco2) 

cont_temp <- contrast(ems, method = "revpairwise", by = "temp_model", 
                      type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_temp, 
                     prob = 0.95, 
                     point.est = mean)

cont_posterior_temp <- gather_emmeans_draws(cont_temp) %>%
  mutate(
    temp_model = as.factor(temp_model)
  )

ems_overall <- emmeans(selected_model, ~ temp_model + pco2_uatm | log(start_diam_mm))

cont_overall <- contrast(ems_overall, method = "revpairwise", 
                         type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_overall, 
                     prob = 0.95, 
                     point.est = mean)


cont_posterior_overall <- gather_emmeans_draws(cont_overall) %>%
  mutate(
    temp_model = as.factor("mean")
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
  group_by(contrast, temp_model) %>%
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
  geom_vline(xintercept = 0, linetype = 2, color = "black", linewidth = 1, 
             alpha = 0.7) +
  stat_pointinterval(point_interval = "mean_hdi", .width = c(0.95, 0.80),
                     interval_size_range = c(1,2), point_size = 4) +
  geom_label_repel(data = summary_emm_const_temp,
                   aes(x = .value, label = HDI), 
                   fill = "white", nudge_x = -0.4, nudge_y = 1) +
  facet_grid(rows = vars(temp_model), switch = "y") +
  labs(
    x = NULL,
    y = "temperature (°C)",
    title = "Skeletal Growth",
    subtitle = expression("pCO"[2]*" effect (1200 µatm vs 600 µatm)")
  ) +
  gg_options() +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  );plot_emm_fit

ggsave("./Figures/Supplemental/pwc_skeleton_pco2.pdf", plot_emm_fit,
       width = 6, height = 8, units = "in")

# CONSTANT & DYNAMIC TEMP Prior predictive check ------------------------------------
priors_const_dynamic <- c(
  prior(normal(0,1), class = "b"),
  prior(student_t(3, 0.4, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01, 0.01), lb = 0, class = "shape")
)

brms_formula_full <- bf(demi_growth_per_day_um ~ 
                          log(start_diam_mm) * temp + 
                          (1|tank_id))

# default_priors <- get_prior(brms_formula_full, data_skel_const_temp_dynamic_temp, family = brms_fam)

# fit model
brm(
  formula = brms_formula_full, 
  data = data_skel_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors_const_dynamic,
  sample_prior = "only",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/skeleton_temp_const_dynamic_PRIORS"
)

# default_priors <- get_prior(brms_formula_full, data_agg_feed_const_temp_dynamic_temp, family = brms_fam)

# Model constant vs dynamic temperature ------------------------------------------------------------
brms_fam <- Gamma(link = "log")

brms_formula_full <- bf(demi_growth_per_day_um ~ 
                          log(start_diam_mm) * temp +
                          (1|tank_id))

# default_priors <- default_prior(brms_formula_full, data = data_skel_dynamic_temp, family = brms_fam)

priors <- c( 
  prior(normal(0,10), class = "b"),
  prior(student_t(3, 0.3, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01,0.01),lb = 0, class = "shape")
)


# fit model
fit_skel_const_temp_dynamic_temp <- brm(
  formula = brms_formula_full, 
  data = data_skel_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_const_dynamic_full"
)

brms_formula_2 <- bf(demi_growth_per_day_um ~ 
                       log(start_diam_mm) + 
                       temp + 
                       (1|tank_id))

# fit model
fit_skel_const_temp_dynamic_temp_2 <- brm(
  formula = brms_formula_2, 
  data = data_skel_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_const_dynamic_2"
)

brms_formula_3 <- bf(demi_growth_per_day_um ~ 
                       temp + 
                       (1|tank_id))

# fit model
fit_skel_const_temp_dynamic_temp_3 <- brm(
  formula = brms_formula_3, 
  data = data_skel_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_const_dynamic_3"
)

brms_formula_4 <- bf(demi_growth_per_day_um ~ (1|tank_id))

priors_const_dynamic_null <- c(
  prior(student_t(3, 0.3, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01, 0.01), lb = 0, class = "shape")
)

# fit model
fit_skel_const_temp_dynamic_temp_4 <- brm(
  formula = brms_formula_4, 
  data = data_skel_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors_const_dynamic_null,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = "on_change", #file_refit,
  backend = backend,
  seed = 123,
  file = "./Model_Output/skeleton_temp_const_dynamic_4"
)

# Model selection via LOO stacking weight -------------------------------------
# Compute model weights via stacking of predictive distributions
# 
# Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian 
# model evaluation using leave-one-out cross-validation and WAIC. 
# Statistics and Computing. 27(5), 1413–1432. 
# doi:10.1007/s11222-016-9696-4 (journal version, preprint arXiv:1507.04544).
#
# Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2019). 
# Pareto smoothed importance sampling. preprint arXiv:1507.02646
#
# Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018) Using 
# stacking to average Bayesian predictive distributions. Bayesian 
# Analysis, advance publication, doi:10.1214/17-BA1091. (online).

loo_model_weights(
  fit_skel_const_temp_dynamic_temp,
  fit_skel_const_temp_dynamic_temp_2,
  fit_skel_const_temp_dynamic_temp_3,
  fit_skel_const_temp_dynamic_temp_4
  )

# 2024-12-16
# Method: stacking
# 
#                                     weight
# fit_skel_const_temp_dynamic_temp   0.332 
# fit_skel_const_temp_dynamic_temp_2 0.620 <- favored
# fit_skel_const_temp_dynamic_temp_3 0.048 
# fit_skel_const_temp_dynamic_temp_4 0.000
selected_model <- fit_skel_const_temp_dynamic_temp_2

# Output summary stats for final model
summary(selected_model)

# Calculate Bayesian R-squared for final model
bayes_R2(selected_model)
#     Estimate  Est.Error      Q2.5    Q97.5
# R2 0.2621525 0.07616975 0.1093392 0.4051911

# Model diagnostics
# Pareto-k-diagnostic
plot(brms::loo(selected_model))

# Predictive performance loo
brms::loo(selected_model)

# Posterior predictive checks -----------------------------------------------------------------
ndraws = 200

ppc_cpue_dens_overlay <- pp_check(selected_model, ndraws = ndraws) + 
  xlim(0,20) +
  gg_options();ppc_cpue_dens_overlay

ppc_cpue_loo_ribbon <- pp_check(selected_model, ndraws = ndraws, 
                                type = "loo_ribbon") + 
  gg_options();ppc_cpue_loo_ribbon

ppc_cpue_loo_pit_qq <- pp_check(selected_model, ndraws = ndraws, 
                                type = "loo_pit_qq") + 
  gg_options();ppc_cpue_loo_pit_qq

ppc_multiplot <- cowplot::plot_grid(
  ppc_cpue_dens_overlay +
    theme(legend.position = c(0.8,0.85)),
  ppc_cpue_loo_ribbon +
    theme(legend.position = c(0.8,0.85)),
  ppc_cpue_loo_pit_qq,
  ncol = 3,
  align = "hv",
  axis = "tblr",
  labels = "AUTO"
);ppc_multiplot

ggsave("./PP_Checks/Skeletal_Growth_Const_Temp_Dynamic_Temp_ppc.pdf", ppc_multiplot,
       width = 16, height = 6, units = "in")

# Sensitivity to priors -------------------------------------------------------
informative_priors <- c(
  prior(normal(0,0.1), class = "b"),
  prior(student_t(3, 0.4, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01, 0.01), lb = 0, class = "shape")
)

fit_informative <- brm(
  formula = brms_formula_full,
  data = data_skel_const_temp_dynamic_temp,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = informative_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/skel_const_temp_dynamic_temp_full_PRIOR_INFORMATIVE"
)

vague_priors <- c(
  prior(normal(0,1), class = "b"),
  prior(student_t(3, 0.4, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01, 0.01), lb = 0, class = "shape")
)

fit_vague <- brm(
  formula = brms_formula_full,
  data = data_skel_const_temp_dynamic_temp,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = vague_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/skel_const_temp_dynamic_temp_full_PRIOR_VAGUE"
)

diffuse_priors <- c(
  prior(normal(0,10), class = "b"),
  prior(student_t(3, 0.4, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01, 0.01), lb = 0, class = "shape")
)

brm(
  formula = brms_formula_full,
  data = data_skel_const_temp_dynamic_temp,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = diffuse_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/skel_const_temp_dynamic_temp_full_PRIOR_DIFFUSE"
)

# Pairwise comparisons --------------------------------------------------------
selected_model <- skeleton_temp_const_dynamic_2 <- readRDS("./Model_Output/skeleton_temp_const_dynamic_2.rds")
selected_model$formula

ems <- emmeans(selected_model, ~ log(start_diam_mm) + temp,
               at = list(temp = c("16","18-14", "20", "21-18")))

emmip(ems, ~ temp, CIs = T)
# cont_temp <- contrast(ems, method = "revpairwise", by = "temp_deg_c", 
#                       type = "response", ratios = TRUE)
cont_temp <- contrast(ems, method = "revpairwise",
                      type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_temp, 
                     prob = 0.95, 
                     point.est = mean)

cont_posterior_temp <- gather_emmeans_draws(cont_temp)

summary_emm_const_temp <- cont_posterior_temp %>%
  group_by(contrast) %>%
  mean_hdi(.value, .width = c(0.95)) %>%
  mutate(
    .value = round(.value, digits = 2),
    .lower = round(.lower, digits = 2),
    .upper = round(.upper, digits = 2)
  ) %>%
  mutate(
    HDI = paste(.value, " (", .lower, ", ", .upper, ")", sep = "")
  )

plot_emm_fit <- ggplot(cont_posterior_temp, 
                       aes(x = .value, y = contrast)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black", linewidth = 1, 
             alpha = 0.7) +
  stat_pointinterval(point_interval = "mean_hdi", .width = c(0.95, 0.80),
                     interval_size_range = c(1,2), point_size = 4) +
  geom_label_repel(data = summary_emm_const_temp,
                   aes(x = .value, label = HDI), 
                   fill = "white", nudge_x = -0.4, nudge_y = 1) +
  facet_grid(rows = vars(contrast), switch = "y") +
  labs(
    x = NULL,
    y = "temperature (°C)",
    title = "Skeletal Growth Rate"
  ) +
  gg_options() +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  );plot_emm_fit

ggsave("./Figures/Supplemental/pwc_consumption_dynamic_vs_constant.pdf", plot_emm_fit,
       width = 6, height = 18, units = "in")
# Test growth analysis ----
library(tidyverse)
library(glmmTMB)

# Read and clean data
bio_growth <- read_csv("./Data/Data_Final/Biometrics_Growth_PercentChange.csv") %>%
  mutate(
    temp_deg_c = factor(temp_deg_c),
    pco2_h_l   = factor(pco2_h_l),
    tank_id    = factor(tank_id)
  ) %>%
  filter(!temp_deg_c %in% c("21-18", "18-14")) %>%
  drop_na(perc_change_diam, initial_diam_mm)   # ensure no NAs in model variables

# Fit the model on the filtered dataset
m1 <- glmmTMB(
  perc_change_diam ~ log(initial_diam_mm) + temp_deg_c * pco2_h_l + (1 | tank_id),
  data = bio_growth,
  family = t_family(link = "identity")
)

m1_disp <- glmmTMB(
  perc_change_diam ~ log(initial_diam_mm) + 
    temp_deg_c * pco2_h_l + 
    (1 | tank_id),
  # Model variance (dispersion) as a function of treatment and size
  dispformula = ~ log(initial_diam_mm) + temp_deg_c * pco2_h_l,
  data = bio_growth,
  family = t_family(link = "identity")
)

AIC(m1, m1_disp)
# df      AIC
# m1      16 1817.221 <- simpler, dispersion model-free formula favored
# m1_disp 28 1833.923
summary(m1)
# Family: t  ( identity )
# Formula:          
#   perc_change_diam ~ log(initial_diam_mm) + temp_deg_c * pco2_h_l +  
#   (1 | tank_id)
# Data: bio_growth
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 1817.2    1878.3    -892.6    1785.2       321 
# 
# Random effects:
#   
#   Conditional model:
#   Groups  Name        Variance Std.Dev.
# tank_id (Intercept) 0.2644   0.5142  
# Number of obs: 337, groups:  tank_id, 24
# 
# Dispersion estimate for t family (sigma^2): 6.86 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             72.6518     6.6766  10.882  < 2e-16 ***
#   log(initial_diam_mm)   -17.1836     1.7153 -10.018  < 2e-16 ***
#   temp_deg_c13            -3.5578     1.0163  -3.501 0.000464 ***
#   temp_deg_c16            -6.9086     0.9820  -7.035 1.99e-12 ***
#   temp_deg_c17            -5.3061     1.0012  -5.300 1.16e-07 ***
#   temp_deg_c18            -5.1484     0.9915  -5.193 2.07e-07 ***
#   temp_deg_c20            -3.1419     0.9776  -3.214 0.001309 ** 
#   pco2_h_ll               -4.8279     1.0107  -4.777 1.78e-06 ***
#   temp_deg_c13:pco2_h_ll   4.0386     1.4250   2.834 0.004595 ** 
#   temp_deg_c16:pco2_h_ll   8.1339     1.3705   5.935 2.94e-09 ***
#   temp_deg_c17:pco2_h_ll   7.0019     1.4225   4.922 8.56e-07 ***
#   temp_deg_c18:pco2_h_ll   3.7932     1.3986   2.712 0.006685 ** 
#   temp_deg_c20:pco2_h_ll   4.2279     1.4166   2.984 0.002841 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# 2. MODEL DIAGNOSTICS
res <- simulateResiduals(m1)
plot(res)                      # Residual diagnostics plot
testDispersion(res)            # Check for over/under dispersion
testUniformity(res)            # Check uniformity of residuals
testZeroInflation(res)         # Check if zero-inflation is needed

# 3. ANOVA TABLE (TYPE II OR III)
# Type III ANOVA (for interaction models)
anova_m1 <- car::Anova(m1, type = 3)
print(anova_m1)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: perc_change_diam
# Chisq Df Pr(>Chisq)    
# (Intercept)          118.409  1  < 2.2e-16 ***
#   log(initial_diam_mm) 100.356  1  < 2.2e-16 ***
#   temp_deg_c            57.482  5  4.023e-11 ***
#   pco2_h_l              22.816  1  1.783e-06 ***
#   temp_deg_c:pco2_h_l   42.135  5  5.532e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# verify alignment
nrow(bio_growth)             # should match
nobs(m1)                     # should match
identical(nrow(bio_growth), nobs(m1))  # should return TRUE

# Extract residuals and fitted values
bio_growth <- bio_growth %>%
  mutate(
    # model-predicted fitted values for visualization
    fitted_perc_change_diam = predict(m1, newdata = bio_growth, type = "response"),
    
    # residuals: computed from fitted data only (no newdata argument)
    resid_perc_change_diam  = residuals(m1, type = "response"),
    
    # standardized (z-scored) residuals: the size-corrected growth index
    size_corrected_growth_index = scale(resid_perc_change_diam)[, 1]
  )

# Inspect residual distribution
hist(bio_growth$resid_perc_change_diam, breaks = 20,
     main = "Residual Growth (size-corrected)", xlab = "Residual % change in diameter")

# 4. POST-HOC COMPARISONS
emm_treat <- emmeans(m1, ~ temp_deg_c * pco2_h_l)
pairs(emm_treat, simple = "each", combine = TRUE, adjust = "tukey")

# Summarize marginal means
summary(emm_treat)
# temp_deg_c pco2_h_l  emmean    SE  df asymp.LCL asymp.UCL
# 10         h         5.3865 0.719 Inf     3.977     6.796
# 13         h         1.8288 0.723 Inf     0.412     3.246
# 16         h        -1.5221 0.660 Inf    -2.816    -0.228
# 17         h         0.0804 0.694 Inf    -1.281     1.441
# 18         h         0.2381 0.679 Inf    -1.093     1.570
# 20         h         2.2446 0.667 Inf     0.937     3.552
# 10         l         0.5586 0.724 Inf    -0.860     1.978
# 13         l         1.0394 0.694 Inf    -0.320     2.399
# 16         l         1.7839 0.658 Inf     0.495     3.073
# 17         l         2.2545 0.719 Inf     0.845     3.664
# 18         l        -0.7966 0.697 Inf    -2.162     0.569
# 20         l         1.6446 0.726 Inf     0.222     3.067
# 
# Confidence level used: 0.95 

# Plot interaction of treatment × pCO₂
emmip(m1, pco2_h_l ~ temp_deg_c,
      CIs = TRUE, type = "response") +
  theme_bw() +
  labs(
    x = "Temperature Treatment (°C)",
    y = "Predicted Δ Volume (mL)",
    color = "pCO₂ Treatment"
  )

