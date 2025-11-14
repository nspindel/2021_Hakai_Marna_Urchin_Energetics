###############################################################################
# Analysis code for 2021 sea urchin heatwaves experiment at the Quadra Island # 
# Ecological Observatory - Hakai Institute                                    #
#                                                                             #
# Metabolic rate                                                              #
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
# Load & wrangle data: ---------------------------------------------------------
data_metabolism <- read.csv("./Data/Respirometry.csv") %>%
  mutate(
    trial_id = case_when(
      trial_day %in% c(19:26) ~ 1,
      trial_day %in% c(41:43) ~ 2,
      trial_day %in% c(66:71) ~ 3
    )
  )

data_metabolism_const_temp_oa <- data_metabolism %>%
  filter(
    temp_deg_c != "21-18" 
    & temp_deg_c != "18-14"
    & !is.na(afdm_g)
  ) %>%
  mutate(
    temp_deg_c = as.numeric(temp_deg_c),
    pco2_uatm = as.factor(pco2_uatm)
  )

data_metabolism_const_temp_dynamic_temp <- data_metabolism %>%
  filter(
    pco2_uatm == "600"
    & !is.na(afdm_g)
  ) %>% 
  mutate(
    temp_deg_c = as.factor(temp_deg_c),
    pco2_uatm = as.factor(pco2_uatm),
    temp_model = case_when(
      temp_deg_c == 10 ~ 10,
      temp_deg_c == 13 ~ 13,
      temp_deg_c == 16 ~ 16,
      temp_deg_c == 17 ~ 17,
      temp_deg_c == 18 ~ 18,
      temp_deg_c == 20 ~ 20,
      temp_deg_c == "18-14" ~ 16,
      temp_deg_c == "21-18" ~ 20
    ),
    temp_label = factor(case_when(
      temp_model == 16 ~ "La Niña (18-14 °C)",
      temp_model == 20 ~ "El Niño (21-18 °C)"
    ), levels = c("La Niña (18-14 °C)", "El Niño (21-18 °C)"))
  ) %>%
  select(
    afdm_g,
    temp_model,
    temp_deg_c,
    temp_label,
    pco2_uatm,
    tank_id,
    absolute_rate_vo2_mg_o2_per_hour
  )

data_wet_mass_rmr_dynamic_temp <- data_metabolism %>%
  filter(
    temp_deg_c %in% c("21-18", "18-14")
    & !is.na(wet_mass_g)
  ) %>% 
  mutate(
    temp_deg_c = as.factor(temp_deg_c),
    pco2_uatm = as.factor(pco2_uatm),
    temp_model = case_when(
      temp_deg_c == 10 ~ 10,
      temp_deg_c == 13 ~ 13,
      temp_deg_c == 16 ~ 16,
      temp_deg_c == 17 ~ 17,
      temp_deg_c == 18 ~ 18,
      temp_deg_c == 20 ~ 20,
      temp_deg_c == "18-14" ~ 16,
      temp_deg_c == "21-18" ~ 20
    ),
    temp_label = factor(case_when(
      temp_model == 16 ~ "La Niña (18-14 °C)",
      temp_model == 20 ~ "El Niño (21-18 °C)"
    ), levels = c("La Niña (18-14 °C)", "El Niño (21-18 °C)"))
  ) %>%
  select(
    trial_id,
    trial_day,
    wet_mass_g,
    temp_model,
    temp_deg_c,
    temp_label,
    pco2_uatm,
    tank_id,
    absolute_rate_vo2_mg_o2_per_hour
  )

data_afdm_rmr_dynamic_temp <- data_metabolism %>%
  filter(
    temp_deg_c %in% c("21-18", "18-14")
    & !is.na(afdm_g)
  ) %>% 
  mutate(
    temp_deg_c = as.factor(temp_deg_c),
    pco2_uatm = as.factor(pco2_uatm),
    temp_model = case_when(
      temp_deg_c == 10 ~ 10,
      temp_deg_c == 13 ~ 13,
      temp_deg_c == 16 ~ 16,
      temp_deg_c == 17 ~ 17,
      temp_deg_c == 18 ~ 18,
      temp_deg_c == 20 ~ 20,
      temp_deg_c == "18-14" ~ 16,
      temp_deg_c == "21-18" ~ 20
    ),
    temp_label = factor(case_when(
      temp_model == 16 ~ "La Niña (18-14 °C)",
      temp_model == 20 ~ "El Niño (21-18 °C)"
    ), levels = c("La Niña (18-14 °C)", "El Niño (21-18 °C)"))
  ) %>%
  select(
    trial_id,
    trial_day,
    afdm_g,
    temp_model,
    temp_deg_c,
    temp_label,
    pco2_uatm,
    tank_id,
    absolute_rate_vo2_mg_o2_per_hour
  )

# Set general computing params ------------------------------------------------
cores        <- 26
chains       <- 4
iterations   <- 18000
warmup       <- 2000
control_list <- list(max_treedepth = 15, adapt_delta = 0.95)
file_refit   <- "on_change"
backend      <- "cmdstanr"
seed         <- 123
brms_fam     <- Gamma(link = "log")

unlink("./Model_Output/metabolism_temp_oa_*", force = TRUE)

# Priors ----
priors_common <- c(
  prior(normal(0, 10), class = "b"),
  prior(student_t(3, 1.1, 5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01, 0.01), lb = 0, class = "shape")
)

priors_gp <- c(
  prior(inv_gamma(4.6, 1.1), class = "lscale", coef = "gptemp_deg_cpco2_uatm1200"),
  prior(inv_gamma(4.6, 1.1), class = "lscale", coef = "gptemp_deg_cpco2_uatm600"),
  prior(R2D2(mean_R2 = 0.1, prec_R2 = 30, cons_D2 = 0.5), class = "sdgp")
)

priors_gp_noby <- c(
  prior(inv_gamma(4.6, 1.1), class = "lscale"),
  prior(R2D2(mean_R2 = 0.1, prec_R2 = 30, cons_D2 = 0.5), class = "sdgp")
)

priors_dynamic <- c( 
  prior(normal(0,10), class = "b"),
  prior(student_t(3, 0.2, 5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01, 0.01),lb = 0, class = "shape")
)

# Candidate model formulae for constant temperatures ----
# Full model
form_full <- brms::bf(
  absolute_rate_vo2_mg_o2_per_hour ~
    log(afdm_g) +
    gp(temp_deg_c, k = 5, c = 5/4, by = pco2_uatm) +
    (1 | tank_id)
)

# Reduced models chosen to test specific hypotheses with principled nesting
form_gp_noby <- brms::bf(
  absolute_rate_vo2_mg_o2_per_hour ~
    log(afdm_g) +
    gp(temp_deg_c, k = 5, c = 5/4) +
    (1 | tank_id)
)
# (Removes the pCO2-specific GP surfaces; retains nonparametric thermal response.)

form_lin_int <- brms::bf(
  absolute_rate_vo2_mg_o2_per_hour ~
    log(afdm_g) + temp_deg_c * pco2_uatm +
    (1 | tank_id)
)
# (Parametric T × pCO2 interaction; no GP.)

form_lin_add <- brms::bf(
  absolute_rate_vo2_mg_o2_per_hour ~
    log(afdm_g) + temp_deg_c + pco2_uatm +
    (1 | tank_id)
)
# (Additive parametric effects; no GP.)

form_size_only <- brms::bf(
  absolute_rate_vo2_mg_o2_per_hour ~
    log(afdm_g) + (1 | tank_id)
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
d <- data_metabolism_const_temp_oa %>%
  select(absolute_rate_vo2_mg_o2_per_hour, afdm_g, temp_deg_c, pco2_uatm, tank_id) %>%
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
    file         = paste0("./Model_Output/metabolism_temp_oa_", nm)
  )
}

brms_formula_dynamic <- brms::bf(absolute_rate_vo2_mg_o2_per_hour ~ 
                                   log(afdm_g) +
                                   temp_deg_c +
                                   (1|tank_id))
# fit model
fit_metabolism_dynamic_temp_afdm_rmr_full <- brm(
  formula = brms_formula_dynamic, 
  data = data_afdm_rmr_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors_dynamic,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  warmup = warmup,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  seed = seed,
  file = "./Model_Output/metabolism_dynamic_temp_afdm_rmr_full"
)

# Model selection based on LOO stacking weights ----
fit_all <- lapply(fit_all, add_criterion, criterion = "loo", reloo = FALSE)
loos    <- lapply(fit_all, \(x) x$criteria$loo)
sapply(loos, function(x) attr(x, "dims")[1])
# Should print identical numbers (e.g., 64000 for all)
# full   gp_noby   lin_int   lin_add size_only 
# 64000     64000     64000     64000     64000

# Stacking and pseudo-BMA (Bayesian bootstrap) weights
w_stack <- brms::loo_model_weights(loos, method = "stacking")
w_stack
# Method: stacking
# 
#   weight
# full      0.000 
# gp_noby   0.002 
# lin_int   0.000 
# lin_add   0.915 # <- favored
# size_only 0.082

w_pbma  <- brms::loo_model_weights(loos, method = "pseudobma")
w_pbma
# Method: pseudo-BMA+ with Bayesian bootstrap
# 
#   weight
# full      0.040 
# gp_noby   0.107 
# lin_int   0.219 
# lin_add   0.500 # <- favored
# size_only 0.135

# Pareto-k summaries
pk <- lapply(loos, \(x) x$diagnostics$pareto_k)
sapply(pk, \(v) table(cut(v,
                          breaks = c(-Inf, 0.5, 0.7, 1, Inf),
                          labels = c("good (≤0.5)", "ok (≤0.7)", "bad (≤1)", "very bad (>1)")
)))
# full gp_noby lin_int lin_add size_only null
# good (≤0.5)    281     281     281     281       281  282
# ok (≤0.7)        0       0       0       0         0    0
# bad (≤1)         0       0       0       0         0    0
# very bad (>1)    1       1       1       1         1    0

bayes_R2(fit_all$lin_add)
# Estimate  Est.Error      Q2.5     Q97.5
# R2 0.8212704 0.01653791 0.7846591 0.8492003

# Posterior predictive check for favored model
pp_check(fit_all$lin_add, ndraws = 50)

# Save session info for provenance
sink("./Model_Output/sessionInfo_RMR.txt"); print(sessionInfo()); sink()

# Visualize conditional effects for constant ----
# selected_model <- fit_all[[winner]]
selected_model <- fit_all[["full"]]

ce_95_const <- conditional_effects(selected_model, robust = F)
ce_80_const <- conditional_effects(selected_model, robust = F, prob = 0.8)

data <- d %>% filter(afdm_g > 0.1) # exclude erroneous measurement

df_ce_temp_pco2_plot_95 <- ce_95_const$`temp_deg_c:pco2_uatm` %>%
  dplyr::mutate(
    plot_label = factor(pco2_uatm, labels = c("600 μatm", "1200 μatm"))
  )

df_ce_temp_pco2_plot_80 <- ce_80_const$`temp_deg_c:pco2_uatm` %>%
  dplyr::mutate(
    plot_label = factor(pco2_uatm, labels = c("600 μatm", "1200 μatm"))
  )

fit_metabolism_dynamic <- readRDS("./Model_Output/metabolism_dynamic_temp_afdm_rmr_full.rds")

ce_95_dyn <- conditional_effects(fit_metabolism_dynamic, robust = F)
ce_80_dyn <- conditional_effects(fit_metabolism_dynamic, prob = 0.8, robust = F)

df_dynamic_95 <- ce_95_dyn$temp_deg_c %>%
  mutate(
    temp_plot = as.numeric(ifelse(temp_deg_c == "21-18", 20, 16)),
    plot_label = ifelse(temp_deg_c == "18-14", "La Niña, 600 μatm", "El Niño, 600 μatm")
  ) %>%
  dplyr::mutate(plot_label = factor(plot_label, levels = c("La Niña, 600 μatm", "El Niño, 600 μatm"))
  )

df_dynamic_80 <- ce_80_dyn$temp_deg_c %>%
  mutate(
    temp_plot = as.numeric(ifelse(temp_deg_c == "21-18", 20, 16)),
    plot_label = ifelse(temp_deg_c == "18-14", "La Niña, 600 μatm", "El Niño, 600 μatm")
  ) %>%
  dplyr::mutate(plot_label = factor(plot_label, levels = c("La Niña, 600 μatm", "El Niño, 600 μatm"))
  )

empirical_means_rmr <- d %>%
  group_by(tank_id, temp_deg_c, pco2_uatm) %>%
  dplyr::summarise(
    emp_mean_rmr = mean(absolute_rate_vo2_mg_o2_per_hour)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    plot_label = factor(pco2_uatm, labels = c("600 μatm", "1200 μatm"))
  ) %>% 
  dplyr::select(tank_id, temp_deg_c, pco2_uatm, plot_label, emp_mean_rmr)

empirical_means_rmr_dyn <- data_afdm_rmr_dynamic_temp %>%
  dplyr::filter(temp_deg_c %in% c("18-14", "21-18")) %>%
  dplyr::group_by(tank_id, temp_deg_c, pco2_uatm) %>%
  dplyr::summarise(
    emp_mean_rmr = mean(absolute_rate_vo2_mg_o2_per_hour)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    temp_deg_c = as.numeric(dplyr::case_when(
      temp_deg_c == "18-14" ~ 16,
      temp_deg_c == "21-18" ~ 20
    )),
    plot_label = factor(temp_deg_c, labels = c("La Niña, 600 μatm", "El Niño, 600 μatm"))
  ) %>% 
  dplyr::select(tank_id, temp_deg_c, pco2_uatm, plot_label, emp_mean_rmr)

df_empirical_means <- rbind(empirical_means_rmr, empirical_means_rmr_dyn) %>%
  dplyr::mutate(temp_deg_c = as.numeric(temp_deg_c))

pd <- position_dodge(width = 0.2)

p_temp_pco2_rmr <- ggplot() +
  geom_line(data = ce_95_const$`temp_deg_c:pco2_uatm`,
            aes(x = temp_deg_c, y = estimate__,
                color = pco2_uatm, linetype = pco2_uatm),
            linewidth = 1.5, position = pd, show.legend = F) +
  geom_ribbon(data = df_ce_temp_pco2_plot_95,
              aes(x = temp_deg_c, y = estimate__,
                  ymin = lower__, ymax = upper__,
                  fill = plot_label), 
              alpha = 0.2) +
  geom_ribbon(data = df_ce_temp_pco2_plot_80,
              aes(x = temp_deg_c, y = estimate__,
                  ymin = lower__, ymax = upper__,
                  fill = plot_label), 
              alpha = 0.4) +
  geom_point(data = df_empirical_means,
             aes(x = temp_deg_c, y = emp_mean_rmr,
                 shape = plot_label, fill = plot_label),
             size = 4, color = "black", position = pd) +
  geom_linerange(data = df_dynamic_95,
                 aes(x = temp_plot, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = plot_label
                 ),
                 linewidth = 1, show.legend = F) +
  geom_linerange(data = df_dynamic_80,
                 aes(x = temp_plot, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = plot_label
                 ),
                 linewidth = 3, show.legend = F) +
  
  scale_color_manual(values = c(
    "blue", # 600 uatm
    "red", # 1200 uatm
    "ivory4", # La Nina 18-14
    "black" # El Nino 21-18
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
    "black" # El Nino
  )) +
  scale_linetype_manual(values = c(
    1, # 600 uatm 
    3 # 1200 uatm
  )) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_x_continuous(breaks = c(10,13,16,17,18,20)) +
  labs(
    x = "temperature (°C)",
    y = expression("respiration rate (mgO"["2"]*"h"^-1*"g"^-1*" ash-free dry mass)"),
    shape = NULL,
    color = NULL,
    title = NULL,
    fill = NULL
  ) +
  gg_options() +
  theme(
    legend.position = c(0.2,0.6)
  )

p_temp_pco2_rmr



# Post hoc on favored model (Linear additive (stacking-weighted winner)) ----
fit_lin_add_metabolism <- readRDS("./Model_Output/metabolism_temp_oa_lin_add.rds")
selected_model <- fit_lin_add_metabolism
selected_model$formula
# absolute_rate_vo2_mg_o2_per_hour ~ log(afdm_g) + temp_deg_c + pco2_uatm + (1 | tank_id) 
stopifnot(inherits(selected_model, "brmsfit"))
stopifnot(exists("d"))

d <- data_metabolism_const_temp_oa %>%
  select(absolute_rate_vo2_mg_o2_per_hour, afdm_g, temp_deg_c, pco2_uatm, tank_id) %>%
  drop_na()

d <- d %>%
  mutate(
    pco2_uatm = factor(pco2_uatm),          # enforce factor coding
    logafdm_g = if ("logafdm_g" %in% names(d)) d$logafdm_g else log(afdm_g)
  )

levels(d$pco2_uatm)                         # confirm treatment levels & baseline

# Population-level = marginalize random effects (re_formula = NA)
emm_pop_link  <- emmeans(selected_model, ~ log(afdm_g) + temp_deg_c + pco2_uatm,
                         re_formula = NA, 
                         type = "link") # on link scale (log mu)
                         # at = list(temp_deg_c = c(10,13,16,17,18,20)))   

emmeans::hpd.summary(emm_pop_link, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")
# afdm_g temp_deg_c pco2_uatm response lower.HPD upper.HPD
# 4.71       15.7 600           2.49      1.83      3.21
# 4.71       15.7 1200          3.94      2.84      5.08
# 
# Point estimate displayed: median 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95

# Differences on response scale (additive in mu)
pairs_pop_diff <- contrast(emm_pop_link, "revpairwise")

emmeans::hpd.summary(pairs_pop_diff, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")

# contrast                                                                                                                       ratio
# afdm_g4.71379432624113 temp_deg_c15.677304964539 pco2_uatm1200 / afdm_g4.71379432624113 temp_deg_c15.677304964539 pco2_uatm600  1.61
# lower.HPD upper.HPD
# 1.02      2.28
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95

# Do the effects of pCO₂ differ across temperature? Possible nonlinear attenuation at 20 °C. Use GP (by pCO₂) model to explore pairwise contrasts by temperature ---- 
fit_full_metabolism <- readRDS("./Model_Output/metabolism_temp_oa_full.rds")

ems <- emmeans(fit_full_metabolism, ~ pco2_uatm | temp_deg_c,
               at = list(temp_deg_c = c(10,13,16,17,18,20)))

cont_pco2 <- contrast(ems, method = "revpairwise", by = "pco2_uatm", ratios = TRUE)

emmeans::hpd.summary(cont_pco2, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")

# pco2_uatm = 600:
#   contrast                    ratio lower.HPD upper.HPD
# temp_deg_c13 / temp_deg_c10 0.812     0.488      1.17
# temp_deg_c16 / temp_deg_c10 1.053     0.518      1.69
# temp_deg_c16 / temp_deg_c13 1.298     0.847      1.80
# temp_deg_c17 / temp_deg_c10 1.350     0.657      2.13
# temp_deg_c17 / temp_deg_c13 1.686     0.935      2.51
# temp_deg_c17 / temp_deg_c16 1.290     0.997      1.52
# temp_deg_c18 / temp_deg_c10 1.808     0.816      2.88
# temp_deg_c18 / temp_deg_c13 2.291     0.952      3.64
# temp_deg_c18 / temp_deg_c16 1.751     0.985      2.45
# temp_deg_c18 / temp_deg_c17 1.341     0.994      1.63
# temp_deg_c20 / temp_deg_c10 3.019     0.892      5.39
# temp_deg_c20 / temp_deg_c13 3.892     0.911      7.10
# temp_deg_c20 / temp_deg_c16 3.003     0.913      5.43
# temp_deg_c20 / temp_deg_c17 2.262     0.939      3.75
# temp_deg_c20 / temp_deg_c18 1.644     0.957      2.37
# 
# pco2_uatm = 1200:
#   contrast                    ratio lower.HPD upper.HPD
# temp_deg_c13 / temp_deg_c10 0.953     0.637      1.24
# temp_deg_c16 / temp_deg_c10 1.102     0.629      1.62
# temp_deg_c16 / temp_deg_c13 1.156     0.853      1.54
# temp_deg_c17 / temp_deg_c10 1.222     0.720      1.84
# temp_deg_c17 / temp_deg_c13 1.291     0.863      1.86
# temp_deg_c17 / temp_deg_c16 1.109     0.970      1.29
# temp_deg_c18 / temp_deg_c10 1.363     0.790      2.12
# temp_deg_c18 / temp_deg_c13 1.453     0.849      2.27
# temp_deg_c18 / temp_deg_c16 1.244     0.924      1.70
# temp_deg_c18 / temp_deg_c17 1.113     0.959      1.33
# temp_deg_c20 / temp_deg_c10 1.572     0.831      2.64
# temp_deg_c20 / temp_deg_c13 1.698     0.821      3.03
# temp_deg_c20 / temp_deg_c16 1.456     0.812      2.47
# temp_deg_c20 / temp_deg_c17 1.290     0.837      1.97
# temp_deg_c20 / temp_deg_c18 1.144     0.881      1.53
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95 

cont_posterior_pco2 <- gather_emmeans_draws(cont_pco2) 

cont_temp <- contrast(ems, method = "revpairwise", by = "temp_deg_c", ratios = TRUE)

emmeans::hpd.summary(cont_temp, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")

# temp_deg_c = 10:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600 1.678     0.782      2.71
# 
# temp_deg_c = 13:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600 1.991     0.874      3.17
# 
# temp_deg_c = 16:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600 1.762     0.919      2.63
# 
# temp_deg_c = 17:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600 1.506     0.883      2.16
# 
# temp_deg_c = 18:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600 1.255     0.728      1.80
# 
# temp_deg_c = 20:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600 0.915     0.403      1.54
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95

cont_posterior_temp <- gather_emmeans_draws(cont_temp) %>%
  mutate(
    temp_deg_c = as.factor(temp_deg_c)
  )

ems_overall <- emmeans(fit_all$full, ~ pco2_uatm | temp_deg_c)

cont_overall <- contrast(ems_overall, method = "revpairwise", 
                         type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_overall, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")

# temp_deg_c = 15.6666666666667:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600  1.62      1.18      2.09
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95

cont_posterior_overall <- gather_emmeans_draws(cont_overall) %>%
  mutate(
    temp_deg_c = as.factor("mean")
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

# write.xlsx(as.data.frame(summary_emm_const_pco2), "./Pairwise_Comparisons_Emmeans_pco2.xlsx", 
#            row.names = F)

cont_pairwise_and_overall_posterior_temp <- data.frame()
cont_pairwise_and_overall_posterior_temp <- rbind(cont_posterior_temp, 
                                                  cont_posterior_overall)

summary_emm_const_temp <- cont_pairwise_and_overall_posterior_temp %>%
  group_by(contrast, temp_deg_c) %>%
  mean_hdi(.value, .width = c(0.95)) %>%
  mutate(
    .value = round(.value, digits = 2),
    .lower = round(.lower, digits = 2),
    .upper = round(.upper, digits = 2)
  ) %>%
  mutate(
    HDI = paste(.value, " (", .lower, ", ", .upper, ")", sep = "")
  )

# write.xlsx(as.data.frame(summary_emm_const_temp), "./Pairwise_Comparisons_Emmeans_temp.xlsx", 
#            row.names = F)

plot_emm_fit <- ggplot(cont_pairwise_and_overall_posterior_temp, 
                       aes(x = .value, y = contrast)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black", linewidth = 1, 
             alpha = 0.7) +
  stat_pointinterval(point_interval = "mean_hdi", .width = c(0.95, 0.80),
                     interval_size_range = c(1,2), point_size = 4) +
  geom_label_repel(data = summary_emm_const_temp,
                   aes(x = .value, label = HDI), 
                   fill = "white", nudge_x = -0.4, nudge_y = 1) +
  facet_grid(rows = vars(temp_deg_c), switch = "y") +
  labs(
    x = NULL,
    y = "temperature (°C)",
    title = NULL,
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

ggsave("./Figures/Supplemental/PWC_Metabolism.jpg", plot_emm_fit,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Metabolism.pdf", plot_emm_fit,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Metabolism.png", plot_emm_fit,
       width = 6, height = 8, units = "in")

# Model RMR in El Nino and La Nina ----
brms_formula_full_dyn_wet <- bf(absolute_rate_vo2_mg_o2_per_hour ~ 
                          log(wet_mass_g) +
                          gp(trial_day, k = 5, c = 5/4, by = temp_label) +
                          (1|tank_id))

# default_priors <- default_prior(object = brms_formula_full, 
#                                 data = data_wet_mass_rmr_dynamic_temp, 
#                                 family = brms_fam)

priors_dyn_wet <- c( 
  prior(normal(0,10), class = "b"),
  prior(student_t(3, 0.6, 5), class = "Intercept"),
  prior(inv_gamma(1.882105, 0.123406), class = "lscale", coef = "gptrial_daytemp_labelElNiño21M18°C"),
  prior(inv_gamma(1.896485, 0.126358), class = "lscale", coef = "gptrial_daytemp_labelLaNiña18M14°C"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(R2D2(mean_R2 = 0.1, prec_R2 = 30, cons_D2 = 0.5), class = "sdgp"),
  prior(gamma(0.01, 0.01),lb = 0, class = "shape")
)

# fit model
fit_metabolism_rmr_dyn_wet_mass <- brm(
  formula = brms_formula_full_dyn_wet, 
  data = data_wet_mass_rmr_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors_dyn_wet,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  warmup = warmup,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  seed = seed,
  file = "./Model_Output/metabolism_dyn_wet_mass_full"
)

# Visualize conditional effects in El Nino and La Nina ----
# fit_metabolism_rmr_dyn_wet_mass <- readRDS("./Model_Output/metabolism_dyn_wet_mass_full.rds")
ce_95_dyn_wet_mass <- conditional_effects(fit_metabolism_rmr_dyn_wet_mass, robust = F)
ce_80_dyn_wet_mass <- conditional_effects(fit_metabolism_rmr_dyn_wet_mass, prob = 0.8, robust = F)

pd <- position_dodge(width = 2)

empirical_means_rmr_dyn_wet_mass <- data_wet_mass_rmr_dynamic_temp %>%
  dplyr::mutate(
    plot_day = case_when(
      trial_id == 1 ~ 22,
      trial_id == 2 ~ 42,
      trial_id == 3 ~ 68
    )
  ) %>%
  dplyr::group_by(plot_day, tank_id, temp_label) %>% 
  dplyr::summarise(
    emp_mean_rmr = mean(absolute_rate_vo2_mg_o2_per_hour),
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(plot_day, tank_id, temp_label, emp_mean_rmr) 

p_dyn_wet_mass_rmr <- ggplot() +
  geom_line(data = ce_95_dyn_wet_mass$`trial_day:temp_label`,
            aes(x = trial_day, y = estimate__,
                color = temp_label, linetype = temp_label),
            linewidth = 1.5, show.legend = F) +
  geom_ribbon(data = ce_95_dyn_wet_mass$`trial_day:temp_label`,
              aes(x = trial_day, y = estimate__,
                  ymin = lower__, ymax = upper__,
                  fill = temp_label), 
              alpha = 0.2) +
  geom_ribbon(data = ce_80_dyn_wet_mass$`trial_day:temp_label`,
              aes(x = trial_day, y = estimate__,
                  ymin = lower__, ymax = upper__,
                  fill = temp_label), 
              alpha = 0.4) +
  geom_point(data = empirical_means_rmr_dyn_wet_mass,
             aes(x = plot_day, y = emp_mean_rmr,
                 shape = temp_label, fill = temp_label),
             size = 4, color = "black", position = pd) +
  # geom_point(data = data_wet_mass_rmr_dynamic_temp,
  #            aes(x = trial_day, y = absolute_rate_vo2_mg_o2_per_hour,
  #                shape = temp_label, fill = temp_label),
  #            size = 1) +
  scale_color_manual(values = c(
    "ivory4", # La Nina 18-14
    "black" # El Nino 21-18
  )) +
  scale_shape_manual(values = c(
    22, # La Nina
    23 # El Nino
  )) +
  scale_fill_manual(values = c(
    "ivory4", # La Nina
    "black" # El Nino
  )) +
  scale_linetype_manual(values = c(
    1, # La Nina
    3 # El Nino
  )) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6),
    labels = scales::number_format(accuracy = 1)
  ) +
  labs(
    x = "days of exposure",
    y = expression("respiration rate (mgO"["2"]*"h"^-1*"g"^-1*"wet mass)"),
    y = NULL,
    shape = NULL,
    color = NULL,
    title = NULL,
    fill = NULL
  ) +
  gg_options() +
  theme(
    legend.position = c(0.6,0.6)
  )

p_dyn_wet_mass_rmr

# Figure 1: a) (constant) RMR x temp x pco2 ;b) (El Nino & La Nina) RMR x temp x time ----
p_respiration_multi <- cowplot::plot_grid(
  egg::tag_facet(p_temp_pco2_rmr, tag_pool = c("a")), 
  egg::tag_facet(p_dyn_wet_mass_rmr, tag_pool = c("b"))
)

p_respiration_multi

ggsave("./Figures/Respiration_Multi.png", p_respiration_multi,
       units = "in", width = 12, height = 6)
ggsave("./Figures/Respiration_Multi.pdf", p_respiration_multi,
       units = "in", width = 12, height = 6)
ggsave("./Figures/Respiration_Multi.jpg", p_respiration_multi,
       units = "in", width = 12, height = 6)

# Model constant and El Nino/La Nina jointly at pco2 = 600 uatm ----
brms_formula_full_const_dyn <- bf(absolute_rate_vo2_mg_o2_per_hour ~ 
                          log(afdm_g) +
                          temp_deg_c +
                          (1|tank_id))

priors_const_dyn <- c( 
  prior(normal(0,10), class = "b"),
  prior(student_t(3, 0.5, 5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01, 0.01),lb = 0, class = "shape")
)

# fit model
fit_metabolism_dynamic_const_temp_afdm_rmr_full <- brm(
  formula = brms_formula_full_const_dyn, 
  data = data_metabolism_const_temp_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors_const_dyn,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  warmup = warmup,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  seed = seed,
  file = "./Model_Output/metabolism_dynamic_const_temp_afdm_rmr_full"
)

# Post hoc pairwise comparisons (constant & El Nino/La Nina) by temperature ----
fit_metabolism_dynamic_const_temp_afdm_rmr_full <- readRDS("./Model_Output/metabolism_dynamic_const_temp_afdm_rmr_full.rds")
ems <- emmeans(fit_metabolism_dynamic_const_temp_afdm_rmr_full, ~ temp_deg_c + log(afdm_g))

cont_temp <- contrast(ems, "pairwise", ratios = TRUE, type = "response")

emmeans::hpd.summary(cont_temp, 
                     prob = 0.95)
# contrast                                                        ratio lower.HPD upper.HPD
# 10 afdm_g4.82194915254237 / 13 afdm_g4.82194915254237           0.966    0.4698     1.634
# 10 afdm_g4.82194915254237 / 16 afdm_g4.82194915254237           0.945    0.4635     1.601
# 10 afdm_g4.82194915254237 / 17 afdm_g4.82194915254237           0.776    0.3775     1.309
# 10 afdm_g4.82194915254237 / 18 afdm_g4.82194915254237           0.375    0.1825     0.633
# 10 afdm_g4.82194915254237 / (18-14 afdm_g4.82194915254237)      1.049    0.5854     1.671
# 10 afdm_g4.82194915254237 / 20 afdm_g4.82194915254237           0.153    0.0744     0.260
# 10 afdm_g4.82194915254237 / (21-18 afdm_g4.82194915254237)      1.464    0.8009     2.324
# 13 afdm_g4.82194915254237 / 16 afdm_g4.82194915254237           0.978    0.4791     1.655
# 13 afdm_g4.82194915254237 / 17 afdm_g4.82194915254237           0.803    0.3814     1.358
# 13 afdm_g4.82194915254237 / 18 afdm_g4.82194915254237           0.388    0.1879     0.653
# 13 afdm_g4.82194915254237 / (18-14 afdm_g4.82194915254237)      1.086    0.5816     1.704
# 13 afdm_g4.82194915254237 / 20 afdm_g4.82194915254237           0.159    0.0759     0.269
# 13 afdm_g4.82194915254237 / (21-18 afdm_g4.82194915254237)      1.516    0.8374     2.414
# 16 afdm_g4.82194915254237 / 17 afdm_g4.82194915254237           0.823    0.4041     1.394
# 16 afdm_g4.82194915254237 / 18 afdm_g4.82194915254237           0.397    0.1878     0.666
# 16 afdm_g4.82194915254237 / (18-14 afdm_g4.82194915254237)      1.111    0.6130     1.775
# 16 afdm_g4.82194915254237 / 20 afdm_g4.82194915254237           0.162    0.0766     0.274
# 16 afdm_g4.82194915254237 / (21-18 afdm_g4.82194915254237)      1.551    0.8304     2.436
# 17 afdm_g4.82194915254237 / 18 afdm_g4.82194915254237           0.483    0.2330     0.817
# 17 afdm_g4.82194915254237 / (18-14 afdm_g4.82194915254237)      1.355    0.7323     2.143
# 17 afdm_g4.82194915254237 / 20 afdm_g4.82194915254237           0.198    0.0905     0.331
# 17 afdm_g4.82194915254237 / (21-18 afdm_g4.82194915254237)      1.890    1.0280     2.991
# 18 afdm_g4.82194915254237 / (18-14 afdm_g4.82194915254237)      2.803    1.5145     4.459
# 18 afdm_g4.82194915254237 / 20 afdm_g4.82194915254237           0.409    0.1992     0.691
# 18 afdm_g4.82194915254237 / (21-18 afdm_g4.82194915254237)      3.916    2.1555     6.260
# (18-14 afdm_g4.82194915254237) / 20 afdm_g4.82194915254237      0.146    0.0762     0.231
# (18-14 afdm_g4.82194915254237) / (21-18 afdm_g4.82194915254237) 1.396    0.8565     2.033
# 20 afdm_g4.82194915254237 / (21-18 afdm_g4.82194915254237)      9.571    5.0757    15.378 # <- Contrast of interest
# 
# Point estimate displayed: median 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95 