###############################################################################
# Analysis code for 2021 sea urchin heatwaves experiment at the Quadra Island # 
# Ecological Observatory - Hakai Institute                                    #
#                                                                             #
# Gonad mass                                                                  #
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
data_gonad <- read.csv("./Data/Derived_Gonad_Production.csv") %>%
  filter(!is.na(single_gonad_dry_mass_g)) %>%
  mutate(
    semi_major_axis_mm = try(0.5 * final_diam_mm),
    semi_minor_axis_mm = try(0.5 * final_ht_mm),
    spheroid_volume_ml = 
      # Equation for volume of oblate spheroid = 4/3pia^2c in mL
      try((4/3 * pi * semi_major_axis_mm ^2 * semi_minor_axis_mm)/1000),
    pco2_uatm = factor(pco2_uatm, levels = c("600", "1200"), labels = c("600 μatm", "1200 μatm"))
  ) %>% 
  mutate(
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
    final_total_gonad_dry_mass_g = single_gonad_dry_mass_g*5,
    estimated_dry_gonad_production_g = final_total_gonad_dry_mass_g - estimated_initial_dry_gonad
  )

data_gonad_const_temp_oa <- data_gonad %>%
  filter(
    temp_deg_c != "21-18" 
    & temp_deg_c != "18-14"
    & !is.na(estimated_dry_gonad_production_g)
  ) %>%
  mutate(
    temp_deg_c = as.numeric(temp_deg_c)
  ) %>%
  select(
    spheroid_volume_ml,
    temp_deg_c,
    pco2_uatm,
    tank_id,
    estimated_dry_gonad_production_g
  )

data_gonad_const_temp_dynamic_temp <- data_gonad %>%
  filter(
    pco2_uatm == "600 μatm"
    & !is.na(estimated_dry_gonad_production_g)
  ) %>% 
  mutate(
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
    spheroid_volume_ml,
    single_gonad_dry_mass_g,
    temp_model,
    temp_deg_c,
    temp_label,
    pco2_uatm,
    tank_id,
    estimated_dry_gonad_production_g
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
brms_fam     <- gaussian()

# unlink("./Model_Output/gonad_temp_oa_*", force = TRUE)

# Estimate initial gonad mass based on wild source population ------------------
dat_gonad_init_wild <- read.csv("./Data/Initial_Wild_Morphometrics.csv")

ggplot(data = dat_gonad_init_wild,
       aes(x = spheroid_volume_ml, y = total_dry_gonad_mass_g)) +
  geom_point() +
  geom_smooth()

brms_formula_init_wild <- bf(total_dry_gonad_mass_g ~ log(spheroid_volume_ml))

# default_priors <- brms::default_prior(brms_formula_init_wild, dat_gonad_init_wild, family = brms_fam)

priors <- c( 
  prior(normal(0,1), class = "b"),
  prior(student_t(3, 0.1, 2.5), class = "Intercept"),
  prior(gamma(0.01, 0.01),lb = 0, class = "shape")
)

# fit model
fit_gonad_init_wild <- brm(
  formula = brms_formula_init_wild, 
  data = dat_gonad_init_wild,
  family = brms_fam,
  cores = cores,
  prior = priors,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  warmup = warmup,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/gonad_init_wild"
)

# predict initial gonad mass for dynamic treatment subjects
dat_init_dynamic <- data_gonad_const_temp_dynamic_temp %>%
  dplyr::filter(temp_deg_c %in% c("18-14", "21-18")) %>%
  dplyr::select(tank_id, temp_deg_c, temp_model, spheroid_volume_ml, single_gonad_dry_mass_g)

int_conditions <- list(spheroid_volume_ml = c(dat_init_dynamic$spheroid_volume_ml))

ce_init_gonad_95 <- brms::conditional_effects(
  fit_gonad_init_wild,
  int_conditions = int_conditions,
  robust = F
)

df_init_gonad_95 <- ce_init_gonad_95$spheroid_volume_ml %>%
  dplyr::mutate(
    init_est_dry_gonad_mass_g = estimate__,
    init_est_dry_gonad_mass_g_lci = lower__,
    init_est_dry_gonad_mass_g_uci = upper__
  ) %>%
  dplyr::select(spheroid_volume_ml, init_est_dry_gonad_mass_g, init_est_dry_gonad_mass_g_lci, init_est_dry_gonad_mass_g_uci)

# Priors ----
priors_common <- c(
  prior(normal(0, 10), class = "b"),
  prior(student_t(3, 1.1, 5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sigma")
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
  prior(student_t(3, 0, 2.5), lb = 0, class = "sigma")
)

# Candidate model formulae for constant temperatures ----
# Full model
form_full <- brms::bf(
  estimated_dry_gonad_production_g ~
    log(spheroid_volume_ml) +
    gp(temp_deg_c, k = 5, c = 5/4, by = pco2_uatm) +
    (1 | tank_id)
)

# Reduced models chosen to test specific hypotheses with principled nesting
form_gp_noby <- brms::bf(
  estimated_dry_gonad_production_g ~
    log(spheroid_volume_ml) +
    gp(temp_deg_c, k = 5, c = 5/4) +
    (1 | tank_id)
)
# (Removes the pCO2-specific GP surfaces; retains nonparametric thermal response.)

form_lin_int <- brms::bf(
  estimated_dry_gonad_production_g ~
    log(spheroid_volume_ml) + temp_deg_c * pco2_uatm +
    (1 | tank_id)
)
# (Parametric T × pCO2 interaction; no GP.)

form_lin_add <- brms::bf(
  estimated_dry_gonad_production_g ~
    log(spheroid_volume_ml) + temp_deg_c + pco2_uatm +
    (1 | tank_id)
)
# (Additive parametric effects; no GP.)

form_size_only <- brms::bf(
  estimated_dry_gonad_production_g ~
    log(spheroid_volume_ml) + (1 | tank_id)
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
d <- data_gonad_const_temp_oa %>%
  select(estimated_dry_gonad_production_g, spheroid_volume_ml, temp_deg_c, pco2_uatm, tank_id) %>%
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
    file         = paste0("./Model_Output/gonad_temp_oa_", nm)
  )
}

brms_formula_full_const_dyn <- bf(estimated_dry_gonad_production_g ~ 
                                    log(spheroid_volume_ml) + 
                                    temp_deg_c + 
                                    (1|tank_id))

priors_const_dynamic <- c( 
  prior(normal(0,10), class = "b"),
  prior(student_t(3, 1.2, 5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sigma")
)

# fit model
fit_gonad_const_temp_dynamic_temp <- brm(
  formula = brms_formula_full_const_dyn, 
  data = data_gonad_const_temp_dynamic_temp,
  family = brms_fam,
  cores = cores,
  prior = priors_const_dynamic,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE),
  init = 0,
  iter = iterations,
  warmup = warmup,
  control = control_list,
  file_refit = file_refit,
  backend = backend,
  file = "./Model_Output/gonad_temp_const_dynamic_full"
)

# Model selection based on LOO stacking weights ----
fit_all <- lapply(fit_all, add_criterion, criterion = "loo", reloo = FALSE)
loos    <- lapply(fit_all, \(x) x$criteria$loo)
sapply(loos, function(x) attr(x, "dims")[1])
# Should print identical numbers (e.g., 64000 for all)
# full   gp_noby   lin_int   lin_add size_only 
# 64000     64000     64000     64000     64000

# Stacking and pseudo-BMA (Bayesian bootstrap) weights
# w_stack <- loo::model_weights(loos, method = "stacking")
w_stack <- brms::loo_model_weights(loos, method = "stacking")
w_stack
# Method: stacking
# 
#   weight
# full      0.000 
# gp_noby   1.000 # <- favored
# lin_int   0.000 
# lin_add   0.000 
# size_only 0.000

w_pbma  <- brms::loo_model_weights(loos, method = "pseudobma")
w_pbma
# Method: pseudo-BMA+ with Bayesian bootstrap
# 
#   weight
# full      0.059 
# gp_noby   0.794 # <- favored
# lin_int   0.036 
# lin_add   0.043 
# size_only 0.068

# Pareto-k summaries
pk <- lapply(loos, \(x) x$diagnostics$pareto_k)
sapply(pk, \(v) table(cut(v,
                          breaks = c(-Inf, 0.5, 0.7, 1, Inf),
                          labels = c("good (≤0.5)", "ok (≤0.7)", "bad (≤1)", "very bad (>1)")
)))
# full gp_noby lin_int lin_add size_only
# good (≤0.5)    345     345     345     345       345
# ok (≤0.7)        0       0       0       0         0
# bad (≤1)         0       0       0       0         0
# very bad (>1)    0       0       0       0         0

selected_model <- fit_all$gp_noby
selected_model$formula
# estimated_dry_gonad_production_g ~ log(spheroid_volume_ml) + gp(temp_deg_c, k = 5, c = 5/4) + (1 | tank_id) 
bayes_R2(selected_model)
# Estimate  Est.Error      Q2.5     Q97.5
# R2 0.08604576 0.02931259 0.03517894 0.1489722

# Posterior predictive check for favored model
pp_check(selected_model, ndraws = 50)

# Save session info for provenance
sink("./Model_Output/sessionInfo_gonad.txt"); print(sessionInfo()); sink()

# Post hoc on favored model (GP with no pCO2 (stacking-weighted winner)) ----
selected_model <- readRDS("./Model_Output/gonad_temp_oa_gp_noby.rds")
selected_model$formula
# estimated_dry_gonad_production_g ~ log(spheroid_volume_ml) + gp(temp_deg_c, k = 5, c = 5/4) + (1 | tank_id) 

d <- data_gonad_const_temp_oa %>%
  select(estimated_dry_gonad_production_g, spheroid_volume_ml, temp_deg_c, pco2_uatm, tank_id) %>%
  drop_na()

stopifnot(inherits(selected_model, "brmsfit"))
stopifnot(exists("d"))

d <- d %>%
  mutate(
    pco2_uatm = factor(pco2_uatm),          # enforce factor coding
    logvolume = if ("logvolume" %in% names(d)) d$logvolume else log(spheroid_volume_ml)
  )

levels(d$pco2_uatm)                         # confirm treatment levels & baseline

# Population-level = marginalize random effects (re_formula = NA)
emm_pop_link  <- emmeans(selected_model, ~ log(spheroid_volume_ml) + temp_deg_c,
                         re_formula = NA, 
                         type = "link", # on link scale (log mu)
at = list(temp_deg_c = c(10,13,16,17,18,20)))

emmeans::hpd.summary(emm_pop_link, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")
# spheroid_volume_ml temp_deg_c emmean lower.HPD upper.HPD
# 33.8         10   1.19     0.937      1.44
# 33.8         13   1.61     1.389      1.83
# 33.8         16   1.80     1.600      1.99
# 33.8         17   1.70     1.539      1.87
# 33.8         18   1.55     1.389      1.71
# 33.8         20   1.26     1.019      1.51
# 
# Point estimate displayed: mean 
# HPD interval probability: 0.95

pairs_pop_diff <- contrast(emm_pop_link, "revpairwise")

emmeans::hpd.summary(pairs_pop_diff, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")
# contrast                                                                                          estimate lower.HPD upper.HPD
# spheroid_volume_ml33.8182941277219 temp_deg_c13 - spheroid_volume_ml33.8182941277219 temp_deg_c10   0.4163    0.1302   0.69887
# spheroid_volume_ml33.8182941277219 temp_deg_c16 - spheroid_volume_ml33.8182941277219 temp_deg_c10   0.6036    0.2580   0.94180
# spheroid_volume_ml33.8182941277219 temp_deg_c16 - spheroid_volume_ml33.8182941277219 temp_deg_c13   0.1873   -0.0370   0.40423
# spheroid_volume_ml33.8182941277219 temp_deg_c17 - spheroid_volume_ml33.8182941277219 temp_deg_c10   0.5076    0.1918   0.82826
# spheroid_volume_ml33.8182941277219 temp_deg_c17 - spheroid_volume_ml33.8182941277219 temp_deg_c13   0.0913   -0.1604   0.34422
# spheroid_volume_ml33.8182941277219 temp_deg_c17 - spheroid_volume_ml33.8182941277219 temp_deg_c16  -0.0960   -0.1837  -0.00708
# spheroid_volume_ml33.8182941277219 temp_deg_c18 - spheroid_volume_ml33.8182941277219 temp_deg_c10   0.3585    0.0622   0.65862
# spheroid_volume_ml33.8182941277219 temp_deg_c18 - spheroid_volume_ml33.8182941277219 temp_deg_c13  -0.0578   -0.3365   0.22940
# spheroid_volume_ml33.8182941277219 temp_deg_c18 - spheroid_volume_ml33.8182941277219 temp_deg_c16  -0.2451   -0.4310  -0.05727
# spheroid_volume_ml33.8182941277219 temp_deg_c18 - spheroid_volume_ml33.8182941277219 temp_deg_c17  -0.1491   -0.2530  -0.04653
# spheroid_volume_ml33.8182941277219 temp_deg_c20 - spheroid_volume_ml33.8182941277219 temp_deg_c10   0.0701   -0.2392   0.38129
# spheroid_volume_ml33.8182941277219 temp_deg_c20 - spheroid_volume_ml33.8182941277219 temp_deg_c13  -0.3462   -0.6839   0.00516
# spheroid_volume_ml33.8182941277219 temp_deg_c20 - spheroid_volume_ml33.8182941277219 temp_deg_c16  -0.5335   -0.8679  -0.18392
# spheroid_volume_ml33.8182941277219 temp_deg_c20 - spheroid_volume_ml33.8182941277219 temp_deg_c17  -0.4374   -0.7214  -0.15571
# spheroid_volume_ml33.8182941277219 temp_deg_c20 - spheroid_volume_ml33.8182941277219 temp_deg_c18  -0.2884   -0.4945  -0.09366
# 
# Point estimate displayed: mean 
# HPD interval probability: 0.95 

# Do the effects of pCO₂ differ across temperature? Explore using full GP by pCO2 ----
fit_gonad_const_temp_oa_full <- readRDS("./Model_Output/gonad_temp_oa_full.rds")
selected_model <- fit_gonad_const_temp_oa_full
selected_model$formula
# estimated_dry_gonad_production_g ~ log(spheroid_volume_ml) + gp(temp_deg_c, k = 5, c = 5/4, by = pco2_uatm) + (1 | tank_id) 

ems <- emmeans(selected_model, ~ pco2_uatm * temp_deg_c,
               at = list(temp_deg_c = c(10,13,16,17,18,20)))

cont_pco2 <- contrast(ems, method = "revpairwise", by = "pco2_uatm", ratios = TRUE)

emmeans::hpd.summary(cont_pco2, 
                     prob = 0.95, 
                     point.est = mean)

# pco2_uatm = 600 μatm:
#   contrast                    estimate lower.HPD upper.HPD
# temp_deg_c13 - temp_deg_c10   0.2489   -0.0632    0.6090
# temp_deg_c16 - temp_deg_c10   0.3901   -0.0799    0.8951
# temp_deg_c16 - temp_deg_c13   0.1413   -0.1134    0.4355
# temp_deg_c17 - temp_deg_c10   0.3244   -0.1083    0.7848
# temp_deg_c17 - temp_deg_c13   0.0755   -0.2149    0.3956
# temp_deg_c17 - temp_deg_c16  -0.0657   -0.1870    0.0377
# temp_deg_c18 - temp_deg_c10   0.2188   -0.1637    0.6253
# temp_deg_c18 - temp_deg_c13  -0.0300   -0.3786    0.3128
# temp_deg_c18 - temp_deg_c16  -0.1713   -0.4379    0.0580
# temp_deg_c18 - temp_deg_c17  -0.1056   -0.2548    0.0267
# temp_deg_c20 - temp_deg_c10   0.0144   -0.3615    0.3959
# temp_deg_c20 - temp_deg_c13  -0.2345   -0.7162    0.1999
# temp_deg_c20 - temp_deg_c16  -0.3758   -0.8869    0.0800
# temp_deg_c20 - temp_deg_c17  -0.3100   -0.7261    0.0591
# temp_deg_c20 - temp_deg_c18  -0.2045   -0.4821    0.0433
# 
# pco2_uatm = 1200 μatm:
#   contrast                    estimate lower.HPD upper.HPD
# temp_deg_c13 - temp_deg_c10   0.1476   -0.0893    0.4748
# temp_deg_c16 - temp_deg_c10   0.2271   -0.1250    0.6940
# temp_deg_c16 - temp_deg_c13   0.0795   -0.1218    0.3389
# temp_deg_c17 - temp_deg_c10   0.2029   -0.1375    0.6467
# temp_deg_c17 - temp_deg_c13   0.0553   -0.1966    0.3369
# temp_deg_c17 - temp_deg_c16  -0.0242   -0.1264    0.0515
# temp_deg_c18 - temp_deg_c10   0.1599   -0.1635    0.5674
# temp_deg_c18 - temp_deg_c13   0.0124   -0.2965    0.3102
# temp_deg_c18 - temp_deg_c16  -0.0671   -0.2812    0.0991
# temp_deg_c18 - temp_deg_c17  -0.0430   -0.1644    0.0459
# temp_deg_c20 - temp_deg_c10   0.0668   -0.2644    0.4279
# temp_deg_c20 - temp_deg_c13  -0.0808   -0.4899    0.2730
# temp_deg_c20 - temp_deg_c16  -0.1603   -0.5730    0.1523
# temp_deg_c20 - temp_deg_c17  -0.1361   -0.4757    0.1111
# temp_deg_c20 - temp_deg_c18  -0.0931   -0.3233    0.0746
# 
# Point estimate displayed: mean 
# HPD interval probability: 0.95

cont_posterior_pco2 <- gather_emmeans_draws(cont_pco2) 

cont_temp <- contrast(ems, method = "revpairwise", by = "temp_deg_c", 
                      type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_temp, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")

# temp_deg_c = 10:
#   contrast              estimate lower.HPD upper.HPD
# 1200 μatm - 600 μatm  0.102302    -0.332     0.539
# 
# temp_deg_c = 13:
#   contrast              estimate lower.HPD upper.HPD
# 1200 μatm - 600 μatm  0.000988    -0.364     0.357
# 
# temp_deg_c = 16:
#   contrast              estimate lower.HPD upper.HPD
# 1200 μatm - 600 μatm -0.060759    -0.454     0.334
# 
# temp_deg_c = 17:
#   contrast              estimate lower.HPD upper.HPD
# 1200 μatm - 600 μatm -0.019196    -0.353     0.316
# 
# temp_deg_c = 18:
#   contrast              estimate lower.HPD upper.HPD
# 1200 μatm - 600 μatm  0.043418    -0.254     0.347
# 
# temp_deg_c = 20:
#   contrast              estimate lower.HPD upper.HPD
# 1200 μatm - 600 μatm  0.154732    -0.273     0.566
# 
# Point estimate displayed: mean 
# HPD interval probability: 0.95

cont_posterior_temp <- gather_emmeans_draws(cont_temp) %>%
  mutate(
    temp_deg_c = as.factor(temp_deg_c)
  )

ems_overall <- emmeans(selected_model, ~ pco2_uatm + temp_deg_c + log(spheroid_volume_ml))

cont_overall <- contrast(ems_overall, method = "revpairwise", 
                         type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_overall, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")

# contrast                                                                                                                                         estimate lower.HPD upper.HPD
# 1200 μatm temp_deg_c15.5739130434783 spheroid_volume_ml33.8182941277219 - 600 μatm temp_deg_c15.5739130434783 spheroid_volume_ml33.8182941277219  -0.0688    -0.475      0.34
# 
# Point estimate displayed: mean 
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
    title = "Gonad Production",
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

ggsave("./Figures/Supplemental/PWC_Gonad.jpg", plot_emm_fit,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Gonad.pdf", plot_emm_fit,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Gonad.png", plot_emm_fit,
       width = 6, height = 8, units = "in")


# Figure 3: Visualize conditional effects on gonad production ----
full_model <- readRDS("./Model_Output/gonad_temp_oa_full.rds")
full_model$formula
# estimated_dry_gonad_production_g ~ log(spheroid_volume_ml) + gp(temp_deg_c, k = 5, c = 5/4, by = pco2_uatm) + (1 | tank_id) 
selected_model <- full_model

ce_95_gonad_constant <- brms::conditional_effects(selected_model, 
                                                  robust = F, prob = 0.95)

ce_80_gonad_constant <- brms::conditional_effects(selected_model, 
                                                  robust = F, prob = 0.80)

fit_gonad_const_temp_dynamic_temp <- readRDS("./Model_Output/gonad_temp_const_dynamic_full.rds")

# int_conditions <- list(temp_model = c(16, 20))
ce_95_gonad_const_temp_dynamic_temp <- brms::conditional_effects(
  fit_gonad_const_temp_dynamic_temp, 
  # int_conditions = int_conditions,
  robust = F, prob = 0.95)

df_dynamic_plot_95 <- ce_95_gonad_const_temp_dynamic_temp$temp_deg_c %>%
  dplyr::filter(temp_deg_c %in% c("18-14", "21-18")) %>%
  dplyr::mutate(plot_label = ifelse(temp_deg_c == "18-14", "La Niña, 600 μatm", "El Niño, 600 μatm")) %>%
  dplyr::mutate(temp_plot = ifelse(temp_deg_c == "18-14", 16, 20)) %>%
  dplyr::mutate(plot_label = factor(plot_label, levels = c("La Niña, 600 μatm", "El Niño, 600 μatm"))
  )

ce_80_gonad_const_temp_dynamic_temp <- brms::conditional_effects(
  fit_gonad_const_temp_dynamic_temp, 
  # int_conditions = int_conditions,
  robust = F, prob = 0.80)

df_dynamic_plot_80 <- ce_80_gonad_const_temp_dynamic_temp$temp_deg_c %>%
  dplyr::filter(temp_deg_c %in% c("18-14", "21-18")) %>%
  dplyr::mutate(plot_label = ifelse(temp_deg_c == "18-14", "La Niña, 600 μatm", "El Niño, 600 μatm")) %>%
  dplyr::mutate(temp_plot = ifelse(temp_deg_c == "18-14", 16, 20)) %>%
  dplyr::mutate(plot_label = factor(plot_label, levels = c("La Niña, 600 μatm", "El Niño, 600 μatm"))
  )

empirical_means_gonad <- data_gonad_const_temp_oa %>%
  dplyr::filter(temp_deg_c %nin% c("18-14", "21-18")) %>%
  dplyr::group_by(tank_id, temp_deg_c, pco2_uatm) %>%
  dplyr::summarise(
    emp_mean_gonad = mean(estimated_dry_gonad_production_g)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(tank_id, temp_deg_c, pco2_uatm, emp_mean_gonad)

empirical_means_gonad_dyn <- data_gonad_const_temp_dynamic_temp %>%
  dplyr::filter(temp_deg_c %in% c("18-14", "21-18")) %>%
  dplyr::group_by(tank_id, temp_deg_c) %>%
  dplyr::summarise(
    emp_mean_gonad = mean(estimated_dry_gonad_production_g)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    temp_deg_c = as.numeric(dplyr::case_when(
      temp_deg_c == "18-14" ~ 16,
      temp_deg_c == "21-18" ~ 20
    )),
    pco2_uatm = factor(temp_deg_c, labels = c("La Niña, 600 μatm", "El Niño, 600 μatm"))
  ) %>% 
  dplyr::select(tank_id, temp_deg_c, pco2_uatm, emp_mean_gonad)

df_empirical_means <- rbind(empirical_means_gonad, empirical_means_gonad_dyn) %>%
  dplyr::mutate(temp_deg_c = as.numeric(temp_deg_c))

p_gonad_pco2_temp <- ggplot() +
  geom_line(data = ce_95_gonad_constant$`temp_deg_c:pco2_uatm`,
            aes(x = temp_deg_c, y = estimate__,
                color = pco2_uatm, linetype = pco2_uatm),
            linewidth = 1.5, show.legend = F) +
  geom_ribbon(data = ce_95_gonad_constant$`temp_deg_c:pco2_uatm`,
              aes(x = temp_deg_c, y = estimate__,
                  ymin = lower__, ymax = upper__,
                  fill = pco2_uatm),
              alpha = 0.2) +
  geom_ribbon(data = ce_80_gonad_constant$`temp_deg_c:pco2_uatm`,
              aes(x = temp_deg_c, y = estimate__,
                  ymin = lower__, ymax = upper__,
                  fill = pco2_uatm),
              alpha = 0.4) +
  geom_linerange(data = df_dynamic_plot_95 ,
                 aes(x = temp_plot, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = plot_label),
                 linewidth = 1) +
  geom_linerange(data = df_dynamic_plot_80 ,
                 aes(x = temp_plot, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = plot_label),
                 linewidth = 3) +
  geom_jitter(data = df_empirical_means,
              aes(x = temp_deg_c, y = emp_mean_gonad,
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
    y = "gonad mass produced (g dry)",
    shape = NULL,
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.2,0.8)); p_gonad_pco2_temp

ggsave("./Figures/Gonad.jpg", p_gonad_pco2_temp,
       units = "in", width = 8, height = 8)
ggsave("./Figures/Gonad.pdf", p_gonad_pco2_temp,
       units = "in", width = 8, height = 8)
ggsave("./Figures/Gonad.png", p_gonad_pco2_temp,
       units = "in", width = 8, height = 8)
