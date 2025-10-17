###############################################################################
# Analysis code for 2021 sea urchin heatwaves experiment at the Quadra Island # 
# Ecological Observatory - Hakai Institute                                    #
#                                                                             #
# Metabolic rate                                                              #
#                                                                             #
# Author: Nathan B. Spindel                                                   #
###############################################################################
# Load CRAN libraries and utility functions  ----------------------------------
# source("Packages.R")

# source("Functions.R")
libs <- c(
  "tidyverse", 
  "janitor",
  "parallel",
  "rstan",
  "brms",
  "cmdstanr",
  "posterior"
  )
lapply(libs, require, character.only = TRUE)

options(scipen = 100000)

# Helper: add GP priors only if GP term present
priors_for <- function(formula) {
  # Collapse possible multiline bf() into one string
  f_chr <- paste(deparse(formula), collapse = " ")
  
  # Case 1: full GP with by = pco2_uatm (two GP surfaces)
  if (grepl("gp\\s*\\([^,]+,\\s*[^,]*by\\s*=\\s*pco2_uatm", f_chr)) {
    message("→ Using GP priors (by-pCO2 specific surfaces)")
    c(priors_common, priors_gp)
    
    # Case 2: single GP term (no by)
  } else if (grepl("gp\\(", f_chr)) {
    message("→ Using GP priors (single surface, no by-factor)")
    c(priors_common, priors_gp_noby)
    
    # Case 3: no GP term
  } else {
    message("→ Using common (non-GP) priors")
    priors_common
  }
}

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
file_refit   <- "always"
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

# Candidate model formulae ----
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
# Fit candidate models ----
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
# Add LOO criteria and compute weights ----
fit_all <- lapply(fit_all, add_criterion, criterion = "loo", reloo = FALSE)
loos    <- lapply(fit_all, \(x) x$criteria$loo)
sapply(loos, function(x) attr(x, "dims")[1])
# Should print identical numbers (e.g., 64000 for all)

# Stacking and pseudo-BMA (Bayesian bootstrap) weights
# w_stack <- loo::model_weights(loos, method = "stacking")
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
# full      0.041 
# gp_noby   0.119 
# lin_int   0.215 
# lin_add   0.479 
# size_only 0.147

# Bayes factors for targeted hypotheses ----
# Bridge sampling is most stable when models are reasonably close in complexity.
# Compare (a) full vs gp_noby: does pCO2-specific thermal surface improve fit?
# Compare (b) lin_int vs lin_add: is T × pCO2 interaction warranted?
bf_full_vs_gp_noby <- bayes_factor(fit_all$full, fit_all$gp_noby)
bf_linint_vs_linadd <- bayes_factor(fit_all$lin_int, fit_all$lin_add)

interpret_bf_manual <- function(bf) {
  if (bf < 1) "Evidence for H0"
  else if (bf < 3) "Anecdotal"
  else if (bf < 10) "Moderate"
  else if (bf < 30) "Strong"
  else if (bf < 100) "Very strong"
  else "Extreme"
}

bf_tbl <- tibble(
  comparison = c("full vs gp_noby (by-pCO2 GP?)",
                 "lin_int vs lin_add (T×pCO2?)"),
  BF10 = c(as.numeric(bf_full_vs_gp_noby$bf),
           as.numeric(bf_linint_vs_linadd$bf)),
  interpretation = sapply(BF10, interpret_bf_manual)
) %>%
  print()
# # A tibble: 2 × 3
# comparison                      BF10 interpretation 
# <chr>                          <dbl> <chr>          
# 1 full vs gp_noby (by-pCO2 GP?) 0.129  Evidence for H0
# 2 lin_int vs lin_add (T×pCO2?)  0.0113 Evidence for H0

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

bayes_R2(fit_all$gp_noby)
# Estimate  Est.Error      Q2.5     Q97.5
# R2 0.8194751 0.01652719 0.7832296 0.8478939

# Posterior predictive check for full model
pp_check(fit_all$full, ndraws = 50)

# Save session info for provenance
sink("./Model_Output/sessionInfo_RMR.txt"); print(sessionInfo()); sink()

# Visualize conditional effects ----
# selected_model <- fit_all[[winner]]
selected_model <- fit_all[["full"]]

ce_95 <- conditional_effects(selected_model, robust = F)
ce_80 <- conditional_effects(selected_model, robust = F, prob = 0.8)

data <- d %>% filter(afdm_g > 0.1) # exclude erroneous measurement

df_ce_temp_pco2_plot_95 <- ce_95$`temp_deg_c:pco2_uatm` %>%
  dplyr::mutate(
    plot_label = factor(pco2_uatm, labels = c("600 μatm", "1200 μatm"))
  )

df_ce_temp_pco2_plot_80 <- ce_80$`temp_deg_c:pco2_uatm` %>%
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
  geom_line(data = ce_95$`temp_deg_c:pco2_uatm`,
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
    y = expression("respiration rate (mgO"["2"]*"h"^-1*"kg"^-1*")"),
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


