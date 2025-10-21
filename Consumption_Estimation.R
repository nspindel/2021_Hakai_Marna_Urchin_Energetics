###############################################################################
# Analysis code for 2021 sea urchin heatwaves experiment at the Quadra Island # 
# Ecological Observatory - Hakai Institute                                    #
# Bayesian Thermal Performance Curve (TPC) estimation using `brms`            #
# Estimation of consumption rate.                                             #
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
  "scales"
)

lapply(libs, require, character.only = TRUE)

source("Utility_Functions.R")

# Load/wrangle data -----------------------------------------------------------
  data_agg_feeding <- read.csv("./Data/Aggregated_Feeding.csv") %>%
  dplyr::mutate(
    trial_id = as.factor(trial_id),
    tank_id = as.factor(tank_id),
    pco2_uatm = as.factor(pco2_uatm),
    days_feeding = hours_feeding / 24,
    per_cap_feed_per_day = (initial_food_dry_mass_g - leftover_food_dry_mass_g) / number_individuals / days_feeding
  ) %>%
  dplyr::mutate(
    per_cap_feed_per_day = ifelse(per_cap_feed_per_day < 0, 0.0001, per_cap_feed_per_day)
  )

data_agg_feeding_const_temp_oa <- data_agg_feeding %>%
  filter(
    temp_deg_c != "21-18" 
    & temp_deg_c != "18-14"
    & !is.na(per_cap_feed_per_day)
  ) %>%
  mutate(
    temp_deg_c = as.numeric(temp_deg_c)
  ) %>%
  select(
    trial_id,
    tank_id,
    temp_deg_c,
    pco2_uatm,
    per_cap_feed_per_day
  )

data_agg_feed_const_temp_dynamic_temp <- data_agg_feeding %>%
  filter(
    pco2_uatm == "600"
    & !is.na(per_cap_feed_per_day)
  ) %>% 
  mutate(
    temp_deg_c = factor(temp_deg_c, 
                           levels = c("10","13","16","18-14","17","18","20","21-18"),
                           labels = c("10","13","16","18-14 La Niña","17","18","20","21-18 El Niño"))  
  ) %>%
  select(
    trial_id,
    tank_id,
    temp_deg_c,
    pco2_uatm,
    per_cap_feed_per_day
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

unlink("./Model_Output/consumption_temp_oa_*", force = TRUE)


# Priors ----
priors_common <- c(
  prior(normal(0, 10), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd")
  # prior(gamma(1, 1), lb = 0, class = "shape")
)

priors_common_no_b <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd")
  # prior(gamma(1, 1), lb = 0, class = "shape")
)

priors_gp <- c(
  prior(inv_gamma(4.6, 1.1), class = "lscale", coef = "gptemp_deg_cpco2_uatm1200"),
  prior(inv_gamma(4.6, 1.1), class = "lscale", coef = "gptemp_deg_cpco2_uatm600"),
  prior(student_t(3, 0, 1), class = "sdgp")
)

priors_gp_noby <- c(
  prior(inv_gamma(4.6, 1.1), class = "lscale"),
  prior(student_t(3, 0, 1), class = "sdgp")
)

# Add priors for the shape sub-model (distributional Gamma)
priors_shape_distr <- c(
  prior(normal(2, 1), class = "Intercept", dpar = "shape"),
  prior(normal(0, 0.5), class = "b", dpar = "shape")
)

# Candidate model formulae for constant temperatures ----
# Full GP (distributional Gamma)
form_full_distr <- brms::bf(
  per_cap_feed_per_day ~ gp(temp_deg_c, k = 10, c = 5/4, by = pco2_uatm) + (1 | tank_id),
  shape ~ 1 + pco2_uatm + temp_deg_c
)

# Reduced GP (single surface, no by)
form_gp_noby <- brms::bf(
  per_cap_feed_per_day ~ gp(temp_deg_c, k = 10, c = 5/4) + (1 | tank_id)
)

# Linear interaction
form_lin_int <- brms::bf(
  per_cap_feed_per_day ~ temp_deg_c * pco2_uatm + (1 | tank_id)
)

# Linear additive
form_lin_add <- brms::bf(
  per_cap_feed_per_day ~ temp_deg_c + pco2_uatm + (1 | tank_id)
)

# Null model
form_null <- brms::bf(
  per_cap_feed_per_day ~ (1 | tank_id)
)

model_set <- list(
  full_distr   = form_full_distr,
  gp_noby      = form_gp_noby,
  lin_int      = form_lin_int,
  lin_add      = form_lin_add,
  null_mod     = form_null
)

# Fit candidate constant temperature models ----
# Compact working df
d <- data_agg_feeding_const_temp_oa %>%
  select(per_cap_feed_per_day, temp_deg_c, pco2_uatm, tank_id) %>%
  drop_na()

# tweak helper function for priors
priors_for <- function(formula) {
  f_chr <- paste(deparse(formula), collapse = " ")
  
  has_gp    <- grepl("gp\\(", f_chr)
  has_by    <- grepl("by\\s*=\\s*pco2_uatm", f_chr)
  is_distr  <- grepl("shape\\s*~", f_chr)
  
  # Detect fixed effects in main (mean) structure
  terms_obj <- try(suppressWarnings(terms(formula)), silent = TRUE)
  has_param <- !inherits(terms_obj, "try-error") &&
    length(attr(terms_obj, "term.labels")) > 0
  
  # Logic tree
  if (has_gp && has_by && is_distr) {
    message("→ Using distributional GP priors (by-pCO2 surfaces + shape sub-model)")
    # exclude class="shape" from priors_common_no_b before combining
    priors_main <- priors_common_no_b[!grepl("class = \"shape\"", deparse(priors_common_no_b))]
    c(priors_main, priors_gp, priors_shape_distr)
    
  } else if (has_gp && has_by) {
    message("→ Using GP priors (by-pCO2 surfaces)")
    c(priors_common_no_b, priors_gp)
    
  } else if (has_gp && !has_by) {
    message("→ Using GP priors (single surface, no by-factor)")
    if (has_param) {
      c(priors_common, priors_gp_noby)
    } else {
      c(priors_common_no_b, priors_gp_noby)
    }
    
  } else if (has_param) {
    message("→ Using common priors with parametric effects")
    priors_common
    
  } else {
    message("→ Using intercept-only priors (null model)")
    priors_common_no_b
  }
}

lapply(model_set, function(m) {
  message("\nModel:\n", deparse(m))
  invisible(priors_for(m))
})

fit_all <- list()

for (nm in names(model_set)) {
  message("\nRefitting model: ", nm)
  
  fit_all[[nm]] <- brm(
    formula      = model_set[[nm]],
    data         = d,
    family       = Gamma(link = "log"),
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
    file_refit   = "always",
    file         = paste0("./Model_Output/consumption_temp_oa_", nm)
  )
}

# Model selection based on LOO stacking weights ----
fit_all <- lapply(fit_all, add_criterion, criterion = "loo", reloo = FALSE)
loos    <- lapply(fit_all, \(x) x$criteria$loo)
sapply(loos, function(x) attr(x, "dims")[1])
# Should print identical numbers (e.g., 64000 for all)
# full  gp_noby  lin_int  lin_add null_mod 
# 64000    64000    64000    64000    64000 

# Stacking and pseudo-BMA (Bayesian bootstrap) weights
# w_stack <- loo::model_weights(loos, method = "stacking")
w_stack <- brms::loo_model_weights(loos, method = "stacking")
w_stack
# Method: stacking
# 
#   weight
#   full_distr 1.000 # favored
#   gp_noby    0.000 
#   lin_int    0.000 
#   lin_add    0.000 
#   null_mod   0.000

w_pbma  <- brms::loo_model_weights(loos, method = "pseudobma")
w_pbma
# Method: pseudo-BMA+ with Bayesian bootstrap
# 
#   weight
#   full_distr 1.000 # favored
#   gp_noby    0.000 
#   lin_int    0.000 
#   lin_add    0.000 
#   null_mod   0.000 

# Pareto-k summaries
pk <- lapply(loos, \(x) x$diagnostics$pareto_k)
sapply(pk, \(v) table(cut(v,
                          breaks = c(-Inf, 0.5, 0.7, 1, Inf),
                          labels = c("good (≤0.5)", "ok (≤0.7)", "bad (≤1)", "very bad (>1)")
)))
# full_distr gp_noby lin_int lin_add null_mod
# good (≤0.5)          191     189     189     189      189
# ok (≤0.7)              0       2       2       0        0
# bad (≤1)               1       1       1       3        3
# very bad (>1)          0       0       0       0        0

selected_model <- fit_all$full

bayes_R2(selected_model)
# Estimate  Est.Error      Q2.5     Q97.5
# R2 0.4520927 0.07946571 0.2861106 0.5872066

# Posterior predictive check for favored model
pp_check(selected_model, ndraws = 50)
# Save session info for provenance
sink("./Model_Output/sessionInfo_Consumption.txt"); print(sessionInfo()); sink()




brms_formula_dynamic <- brms::bf(per_cap_feed_per_day ~ 
                                   log(afdm_g) +
                                   temp_deg_c +
                                   (1|tank_id))
# fit model
fit_consumption_dynamic_temp_afdm_rmr_full <- brm(
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
  file = "./Model_Output/consumption_dynamic_temp_afdm_rmr_full"
)

brms_formula_full <- bf(per_cap_feed_per_day ~ 
                          temp_deg_c * pco2_uatm +
                          gp(temp_deg_c, k = 5, c = 5/4, by = pco2_uatm) + 
                          (1|tank_id))

# default_priors <- default_prior(brms_formula_full, data_agg_feeding_const_temp_oa, brms_fam)
fit_feed_const_temp_oa <- brm(
  formula = brms_formula_full, 
  data = data_agg_feeding_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  seed = seed,
  backend = backend,
  file = "./Model_Output/consumption_temp_oa_full"
)

brms_formula_2 <- bf(per_cap_feed_per_day ~ 
                       temp_deg_c + 
                       pco2_uatm +
                       gp(temp_deg_c, k = 5, c = 5/4, by = pco2_uatm) + 
                       (1|tank_id))

# fit model
fit_feed_const_temp_oa_2 <- brm(
  formula = brms_formula_2, 
  data = data_agg_feeding_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  seed = seed,
  backend = backend,
  file = "./Model_Output/consumption_temp_oa_2"
)

brms_formula_3 <- bf(per_cap_feed_per_day ~ 
                       temp_deg_c + 
                       gp(temp_deg_c, k = 5, c = 5/4) + 
                       (1|tank_id))

# fit model
fit_feed_const_temp_oa_3 <- brm(
  formula = brms_formula_3, 
  data = data_agg_feeding_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  seed = seed,
  backend = backend,
  file = "./Model_Output/consumption_temp_oa_3"
)

brms_formula_4 <- bf(per_cap_feed_per_day ~ 
                       pco2_uatm + 
                       gp(temp_deg_c, k = 5, c = 5/4, by = pco2_uatm) + 
                       (1|tank_id))

# fit model
fit_feed_const_temp_oa_4 <- brm(
  formula = brms_formula_4, 
  data = data_agg_feeding_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  seed = seed,
  backend = backend,
  file = "./Model_Output/consumption_temp_oa_4"
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

brms::loo_model_weights(
  fit_feed_const_temp_oa,
  fit_feed_const_temp_oa_2,
  fit_feed_const_temp_oa_3,
  fit_feed_const_temp_oa_4
)

# 2024-12-11
# Method: stacking
# 
#                          weight
# fit_feed_const_temp_oa   1.000 <- favored
# fit_feed_const_temp_oa_2 0.000 
# fit_feed_const_temp_oa_3 0.000 
# fit_feed_const_temp_oa_4 0.000 

selected_model <- fit_feed_const_temp_oa

# Output summary stats for final model
summary(selected_model)

# Calculate Bayesian R-squared for final model
bayes_R2(selected_model)
#     Estimate  Est.Error      Q2.5     Q97.5
# R2 0.5133133 0.06469998 0.3695859 0.6202064


# Posterior predictive checks -----------------------------------------------------------------
# Pareto-k-diagnostic
plot(brms::loo(selected_model))

# Predictive performance loo
brms::loo(selected_model)

# Posterior predictive checks
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

ggsave("./PP_Checks/Consumption_Agg_Const_Temp_OA_ppc.pdf", ppc_multiplot,
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
  data = data_agg_feeding_const_temp_oa,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = informative_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/consumption_temp_oa_full_PRIOR_INFORMATIVE"
)

vague_priors <- c(
  prior(normal(0,1), class = "b"),
  prior(student_t(3, 0.6, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sds"),
  prior(gamma(3, 0.5), lb = 0, class = "shape"),
  prior(normal(0,1), class = "b", dpar = "hu"),
  prior(logistic(0, 1), class = "Intercept", dpar = "hu"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sds", dpar = "hu")
)

fit_vague <- brm(
  formula = brms_formula_full,
  data = data_agg_feeding_const_temp_oa,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = vague_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/consumption_temp_oa_full_PRIOR_VAGUE"
)

diffuse_priors <- c(
  prior(normal(0,10), class = "b"),
  prior(student_t(3, 0.6, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sds"),
  prior(gamma(3, 0.5), lb = 0, class = "shape"),
  prior(normal(0,1), class = "b", dpar = "hu"),
  prior(logistic(0, 1), class = "Intercept", dpar = "hu"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sds", dpar = "hu")
)

brm(
  formula = brms_formula_full,
  data = data_agg_feeding_const_temp_oa,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = diffuse_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/consumption_temp_oa_full_PRIOR_DIFFUSE"
)

pp_check(fit_informative)
pp_check(fit_vague)
pp_check(fit_diffuse)

# Hypothesis tests ----
fit_feed_const_temp_oa <- readRDS("./Model_Output/consumption_temp_oa_full.rds")
selected_model <- fit_feed_const_temp_oa
selected_model$formula

# temperature
brms_formula_5 <- bf(per_cap_feed_per_day ~ 
                          pco2_uatm + temp_deg_c +
                       
                          (1|tank_id))
# default_prior(brms_formula_5, data_agg_feeding_const_temp_oa, brms_fam)
priors2 <- c( 
  prior(normal(0,1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  # prior(inv_gamma(4.6, 1.1), lb = 0, class = "lscale"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  # prior(student_t(3, 0, 2.5), lb = 0, class = "sdgp"),
  prior(gamma(3, 1),lb = 0, class = "shape")
)

fit_feed_const_temp_oa_5 <- brm(
  formula = brms_formula_5, 
  data = data_agg_feeding_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors2,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/consumption_temp_oa_5"
)

brms_formula_6 <- bf(per_cap_feed_per_day ~ 
                       pco2_uatm +
                       
                       (1|tank_id))

fit_feed_const_temp_oa_6 <- brm(
  formula = brms_formula_6, 
  data = data_agg_feeding_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors2,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/consumption_temp_oa_6"
)

bayes_factor(fit_feed_const_temp_oa_5, fit_feed_const_temp_oa_6)
# Estimated Bayes factor in favor of fit_feed_const_temp_oa_5 over fit_feed_const_temp_oa_6: 21.89841
# BF > 20, strong evidence for temperature effect

brms_formula_7 <- bf(per_cap_feed_per_day ~ 
                       temp_deg_c +
                       
                       (1|tank_id))

fit_feed_const_temp_oa_7 <- brm(
  formula = brms_formula_7, 
  data = data_agg_feeding_const_temp_oa,
  family = brms_fam,
  cores = cores,
  prior = priors2,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/consumption_temp_oa_7"
)
bayes_factor(fit_feed_const_temp_oa_5, fit_feed_const_temp_oa_7)
# Estimated Bayes factor in favor of fit_feed_const_temp_oa_5 over fit_feed_const_temp_oa_7: 283.06023
# BF > 20, strong evidence for pco2 effect


# Predictions -----------------------------------------------------------------
ce_95_pco2_comparison <- brms::conditional_effects(fit_feed_const_temp_oa_2, robust = F)

df_ce_temp_pco2_plot_95 <- ce_95_pco2_comparison$`temp_deg_c:pco2_uatm` %>%
  dplyr::mutate(
    pco2_plot = factor(pco2_uatm, labels = c("600 μatm", "1200 μatm"))
  )

ce_80_pco2_comparison <- brms::conditional_effects(fit_feed_const_temp_oa_2, robust = F, prob = 0.8)

df_ce_temp_pco2_plot_80 <- ce_80_pco2_comparison$`temp_deg_c:pco2_uatm` %>%
  dplyr::mutate(
    pco2_plot = factor(pco2_uatm, labels = c("600 μatm", "1200 μatm"))
  )

empirical_means_feed <- data_agg_feeding_const_temp_oa %>%
  dplyr::group_by(tank_id, temp_deg_c, pco2_uatm) %>%
  dplyr::summarise(
    emp_mean_feed = mean(per_cap_feed_per_day)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    plot_label = factor(pco2_uatm, labels = c("600 μatm", "1200 μatm"))
  ) %>% 
  dplyr::select(tank_id, temp_deg_c, plot_label, emp_mean_feed)

empirical_means_feed_dynamic <- data_agg_feed_const_temp_dynamic_temp %>%
  dplyr::filter(temp_deg_c %in% c("18-14 La Niña", "21-18 El Niño")) %>%
  dplyr::group_by(tank_id, temp_deg_c) %>%
  dplyr::summarise(
    emp_mean_feed = mean(per_cap_feed_per_day)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    temp_deg_c = as.numeric(dplyr::case_when(
      temp_deg_c == "18-14 La Niña" ~ 16,
      temp_deg_c == "21-18 El Niño" ~ 20
    )),
    plot_label = factor(temp_deg_c, labels = c("18-14 La Niña, 600 μatm", "21-18 El Niño, 600 μatm"))
  ) %>% 
  dplyr::select(tank_id, temp_deg_c, plot_label, emp_mean_feed)

df_empirical_means <- rbind(empirical_means_feed, empirical_means_feed_dynamic)

pd <- position_jitterdodge(dodge.width = 0.1)

p_agg_feed <- ggplot() +
  # geom_jitter(data = data_agg_feeding_const_temp_oa,
  #             aes(x = temp_deg_c, y = per_cap_feed_per_day,
  #                 color = pco2_uatm, shape = pco2_uatm),
  #             alpha = 0.4) +
  geom_line(data = df_ce_temp_pco2_plot_95,
            aes(x = temp_deg_c, y = estimate__,
                color = pco2_plot, linetype = pco2_plot),
            linewidth = 1.5, show.legend = F) +
  geom_ribbon(data = df_ce_temp_pco2_plot_95,
              aes(x = temp_deg_c, y = estimate__,
                  ymin = lower__, ymax = upper__, 
                  fill = pco2_plot),
              alpha = 0.2) +
  geom_ribbon(data = df_ce_temp_pco2_plot_80,
              aes(x = temp_deg_c, y = estimate__,
                  ymin = lower__, ymax = upper__, 
                  fill = pco2_plot),
              alpha = 0.4) +
  geom_linerange(data = dynamic_modeled_means_95,
                 aes(x = temp_plot, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = temp_plot_label), 
                 linewidth = 1) +
  geom_linerange(data = dynamic_modeled_means_80,
                 aes(x = temp_plot, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = temp_plot_label), 
                 linewidth = 3) +
  geom_point(data = df_empirical_means,
             aes(x = temp_deg_c, y = emp_mean_feed,
                 shape = plot_label, fill = plot_label),
             size = 4, color = "black", position = pd) +
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
    y = expression("consumption (g dry mass urchin"^-1*"day"^-1*")"),
    shape = NULL,
    color = NULL,
    fill = NULL
  ) +
  theme(
    legend.position = c(0.8,0.15)
  ); p_agg_feed

ggsave("./Figures/Per_Cap_Consumption.pdf", p_agg_feed,
       units = "in", width = 8, height = 8)
ggsave("./Figures/Per_Cap_Consumption.png", p_agg_feed,
       units = "in", width = 8, height = 8)
ggsave("./Figures/Per_Cap_Consumption.jpg", p_agg_feed,
       units = "in", width = 8, height = 8)

# Pairwise comparisons --------------------------------------------------------
fit_feed_const_temp_oa <- readRDS("./Model_Output/consumption_temp_oa_full.rds")
selected_model <- fit_feed_const_temp_oa
selected_model$formula

ems <- emmeans(selected_model, ~ pco2_uatm | temp_deg_c,
               at = list(temp_deg_c = c(10,13,16,17,18,20)))

cont_pco2 <- contrast(ems, method = "revpairwise", by = "pco2_uatm", 
                      type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_pco2, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")

cont_temp <- contrast(ems, method = "revpairwise", by = "temp_deg_c", 
                      type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_temp, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")

cont_posterior_temp <- gather_emmeans_draws(cont_temp) %>%
  mutate(
    temp_deg_c = as.factor(temp_deg_c)
  )

ems_overall <- emmeans(fit_feed_const_temp_oa_2, ~ pco2_uatm | temp_deg_c)

cont_overall <- contrast(ems_overall, method = "revpairwise", 
                         type = "response", ratios = TRUE)

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
    title = "Per Capita Consumption",
    subtitle = expression("pCO"[2]*" effect (1200 µatm vs 600 µatm)")
  ) +
  gg_options() +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  );plot_emm_fit

ggsave("./Figures/Supplemental/pwc_consumption_pco2.pdf", plot_emm_fit,
       width = 6, height = 8, units = "in")

# CONSTANT & DYNAMIC TEMP Prior predictive check ------------------------------------
brms_formula_full <- bf(per_cap_feed_per_day ~ 
                          temp_deg_c + 
                          (1|tank_id))

# default_priors <- get_prior(brms_formula_full, data_agg_feed_const_temp_dynamic_temp, family = brms_fam)

priors_const_dynamic <- c( 
  prior(normal(0,1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(3, 1),lb = 0, class = "shape")
)

# fit model
prior_check <- brm(
  formula = brms_formula_full, 
  data = data_agg_feed_const_temp_dynamic_temp,
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
  file = "./HPC/HPC_Model_Output/consumption_temp_const_dynamic_PRIORS"
)

# (1) [full model] --------------------------------------------
brms_formula_full <- bf(per_cap_feed_per_day ~ 
                          temp_deg_c + 
                          (1|tank_id))

priors_const_dynamic <- c(
  prior(normal(0,10), class = "b"),
  prior(student_t(3, 0.9, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(gamma(0.01, 0.01), class = "shape")
)

# fit model
fit_feed_const_temp_dynamic_temp <- brm(
  formula = brms_formula_full, 
  data = data_agg_feed_const_temp_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors_const_dynamic,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/consumption_temp_const_dynamic_full"
)

# (2) [no temp effect] --------------------------------------------------------
brms_formula_2 <- bf(per_cap_feed_per_day ~ (1|tank_id))

priors_const_dynamic <- c(
  prior(student_t(3, 0.9, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01, 0.01), lb = 0, class = "shape")
)

# fit model
fit_feed_const_temp_dynamic_temp_2 <- brm(
  formula = brms_formula_2, 
  data = data_agg_feed_const_temp_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors_const_dynamic,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/consumption_temp_const_dynamic_2"
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

brms::loo_model_weights(
  fit_feed_const_temp_dynamic_temp,
  fit_feed_const_temp_dynamic_temp_2
  )

# 2024-10-21
# Method: stacking
# 
#                                     weight
# fit_feed_const_temp_dynamic_temp   1.000 <- favored
# fit_feed_const_temp_dynamic_temp_2 0.000 

selected_model <- fit_feed_const_temp_dynamic_temp

# Output summary stats for final model
summary(selected_model)

# Calculate Bayesian R-squared for final model
bayes_R2(selected_model)
#     Estimate Est.Error           Q2.5     Q97.5
# R2 0.4156406 0.1693348 0.000008674612 0.5084325
# Model diagnostics
# Pareto-k-diagnostic
plot(brms::loo(selected_model))

# Predictive performance loo
brms::loo(selected_model)

# Posterior predictive checks -----------------------------------------------------------------
ndraws = 200

ppc_cpue_dens_overlay <- pp_check(selected_model, ndraws = ndraws) + 
  # xlim(0,20) +
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

ggsave("./PP_Checks/Consumption_Agg_Const_Temp_Dynamic_Temp_ppc.pdf", ppc_multiplot,
       width = 16, height = 6, units = "in")

# Sensitivity to priors -------------------------------------------------------
informative_priors <- c(
  prior(normal(0,0.1), class = "b"),
  prior(student_t(3, 0.9, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(gamma(0.01, 0.01), class = "shape")
)

fit_informative <- brm(
  formula = brms_formula_full,
  data = data_agg_feed_const_temp_dynamic_temp,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = informative_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/consumption_temp_const_dynamic_full_PRIOR_INFORMATIVE"
)

vague_priors <- c(
  prior(normal(0,1), class = "b"),
  prior(student_t(3, 0.9, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(gamma(0.01, 0.01), class = "shape")
)

fit_vague <- brm(
  formula = brms_formula_full,
  data = data_agg_feed_const_temp_dynamic_temp,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = vague_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/consumption_temp_const_dynamic_full_PRIOR_VAGUE"
)

vague_priors <- c(
  prior(normal(0,10), class = "b"),
  prior(student_t(3, 0.9, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(gamma(0.01, 0.01), class = "shape")
)

brm(
  formula = brms_formula_full,
  data = data_agg_feed_const_temp_dynamic_temp,
  family = brms_fam,
  cores = parallel::detectCores(),
  prior = vague_priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  iter = iterations,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/consumption_temp_const_dynamic_full_PRIOR_DIFFUSE"
)

pp_check(fit_informative)
pp_check(fit_vague)
pp_check(fit_diffuse)

# Predictions -----------------------------------------------------------------
ce_95_const_dynamic_comparison <- brms::conditional_effects(selected_model, robust = F)

dynamic_modeled_means_95 <- ce_95_const_dynamic_comparison$temp_deg_c %>%
  dplyr::mutate(
    temp_plot = as.numeric(dplyr::case_when(
      temp_deg_c == "10" ~ 10,
      temp_deg_c == "13" ~ 13,
      temp_deg_c == "16" ~ 16,
      temp_deg_c == "18-14 La Niña" ~ 16,
      temp_deg_c == "17" ~ 17,
      temp_deg_c == "18" ~ 18,
      temp_deg_c == "20" ~ 20,
      temp_deg_c == "21-18 El Niño" ~ 20
    ))
  ) %>% 
  dplyr::filter(temp_deg_c %in% c("18-14 La Niña", "21-18 El Niño")) %>%
  dplyr::mutate(
    temp_plot_label = factor(temp_deg_c, labels = c("18-14 La Niña, 600 μatm", "21-18 El Niño, 600 μatm"))
  )


ce_80_const_dynamic_comparison <- conditional_effects(selected_model, robust = F, prob = 0.8)

dynamic_modeled_means_80 <- ce_80_const_dynamic_comparison$temp_deg_c %>%
  dplyr::mutate(
    temp_plot = as.numeric(dplyr::case_when(
      temp_deg_c == "10" ~ 10,
      temp_deg_c == "13" ~ 13,
      temp_deg_c == "16" ~ 16,
      temp_deg_c == "18-14 La Niña" ~ 16,
      temp_deg_c == "17" ~ 17,
      temp_deg_c == "18" ~ 18,
      temp_deg_c == "20" ~ 20,
      temp_deg_c == "21-18 El Niño" ~ 20
    ))
  ) %>% 
  dplyr::filter(temp_deg_c %in% c("18-14 La Niña", "21-18 El Niño")) %>%
  dplyr::mutate(
    temp_plot_label = factor(temp_deg_c, labels = c("18-14 La Niña, 600 μatm", "21-18 El Niño, 600 μatm"))
  )

empirical_means_feed_const_dyn <- data_agg_feed_const_temp_dynamic_temp %>%
  group_by(tank_id, temp_deg_c) %>%
  dplyr::summarise(
    emp_mean_feed_const_dyn = mean(per_cap_feed_per_day)
  )

p_agg_feed_const_dynamic_point_estimates <- ggplot() +
  geom_jitter(data = data_agg_feed_const_temp_dynamic_temp,
              aes(x = temp_deg_c, y = per_cap_feed_per_day,
                  color = temp_deg_c),
              alpha = 0.4, width = 0.25
  ) +
  geom_linerange(data = ce_95_const_dynamic_comparison$temp_deg_c,
                 aes(x = temp_deg_c, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = temp_deg_c),
                 alpha = 0.6, linewidth = 1) +
  geom_linerange(data = ce_80_const_dynamic_comparison$temp_deg_c,
                 aes(x = temp_deg_c, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = temp_deg_c),
                 alpha = 0.9, linewidth = 2.5) +
  geom_jitter(data = empirical_means_feed_const_dyn,
              aes(x = temp_deg_c, y = emp_mean_feed_const_dyn,
                  fill = temp_deg_c),
              size = 4, shape = 21, color = "black", width = 0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_color_viridis(discrete = T, option = "plasma")+
  scale_fill_viridis(discrete = T, option = "plasma")+
  gg_options() +
  labs(
    x = "temperature (°C)",
    y = expression("consumption (g dry mass urchin"^-1*"day"^-1*")")
  ) +
  theme(
    legend.position = "none"
  ); p_agg_feed_const_dynamic_point_estimates

# Pairwise comparisons --------------------------------------------------------
fit_feed_const_temp_dynamic_temp$formula
# per_cap_feed_per_day ~ temp_deg_c + (1 | tank_id) 

ems <- emmeans(fit_feed_const_temp_dynamic_temp, ~ temp_deg_c,
               at = list(temp_deg_c = c("16","18-14 La Niña", "20", "21-18 El Niño")))
# ems <- emmeans(fit_feed_const_temp_dynamic_temp, ~ temp_deg_c)

emmip(ems, ~ temp_deg_c, CIs = T)
# cont_temp <- contrast(ems, method = "revpairwise", by = "temp_deg_c", 
#                       type = "response", ratios = TRUE)
cont_temp <- contrast(ems, method = "revpairwise", 
                      type = "response", ratios = TRUE)

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
    title = "Per Capita Consumption"
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
