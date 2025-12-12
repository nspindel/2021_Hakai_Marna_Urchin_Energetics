###############################################################################
# Analysis code for 2021 sea urchin heatwaves experiment at the Quadra Island # 
# Ecological Observatory - Hakai Institute                                    #
# Bayesian Thermal Performance Curve (TPC) estimation using `brms`            #
# Estimation of food conversion efficiency.                                   #
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
  "readxl",
  "boot"
)
lapply(libs, require, character.only = TRUE)

source("Utility_Functions.R")

# Load & wrangle data ---------------------------------------------------------
dat_raw <- read_excel("./Data/Assimilation_Efficiency.xlsx", sheet = 1)

dat_processed <- dat_raw %>%
  mutate(
    tank_id = as.factor(tank_id),
    temp    = as.numeric(temp),
    pco2    = ifelse(pco2 == "l", "600 μatm", "1200 μatm"),
    pco2    = factor(pco2, levels = c("600 μatm","1200 μatm")),
    temp_pco2_interaction = paste0(temp, "_", pco2),
    bin_id  = as.factor(bin_id),
    pit_id  = as.factor(pit_id),
    
    semi_major_axis_mm = 0.5 * diam_mm,
    semi_minor_axis_mm = 0.5 * ht_mm,
    
    # Equation for volume of oblate spheroid = 4/3*pi*a^2*c in mL
    spheroid_volume_ml = (4/3 * pi * semi_major_axis_mm^2 * semi_minor_axis_mm) / 1000,
    
    leftover_dry_food_g   = food_wb_plus_food_g - food_wb_mass_g,
    dry_food_consumed_g   = food_mass_g - leftover_dry_food_g,
    dry_feces_egested_g   = feces_wb_plus_feces_g - feces_wb_mass_g,
    hours_feeding         = as.numeric(hours_feeding),
    minutes_feeding       = as.numeric(minutes_feeding),
    
    dry_food_consumed_per_day_g  = dry_food_consumed_g / days_feeding,
    dry_feces_egested_per_day_g  = dry_feces_egested_g / days_feeding,
    
    assimilation_efficiency_percent =
      ((dry_food_consumed_per_day_g - dry_feces_egested_per_day_g) / dry_food_consumed_per_day_g) * 100,
    assimilation_efficiency_proportion =
      ((dry_food_consumed_per_day_g - dry_feces_egested_per_day_g) / dry_food_consumed_per_day_g)
  ) %>%
  # robust guardrail for Beta regression (even if later data edits introduce 0/1)
  mutate(
    assimilation_efficiency_proportion =
      pmin(pmax(assimilation_efficiency_proportion, 1e-6), 1 - 1e-6)
  )

# Quick plots ---------------------------------------------------------------
indiv_feed_plot <- ggplot(dat_processed, aes(x = temp, y = dry_food_consumed_per_day_g,
                                             fill = pco2, color = pco2)) +
  geom_jitter(alpha = 0.2) +
  geom_smooth() +
  scale_color_manual(values = c("blue","red")) +
  scale_fill_manual(values  = c("blue","red")) +
  gg_options() +
  labs(
    x = "Temperature (°C)",
    y = expression("Consumption (g dry mass urchin"^-1*"day"^-1*")"),
    title = "Sea Urchin Condos"
  ) +
  guides(
    color = guide_legend(title = expression("pCO"[2]),
                         override.aes = list(fill = c("blue","red"), color = c("blue","red"))),
    fill  = guide_legend(title = expression("pCO"[2]),
                         override.aes = list(fill = c("blue","red"), color = c("blue","red")))
  ) +
  theme(legend.position = "top")
indiv_feed_plot

indiv_assim_eff_plot <- ggplot(dat_processed, aes(x = temp, y = assimilation_efficiency_percent,
                                                  fill = pco2, color = pco2)) +
  geom_jitter(alpha = 0.2) +
  geom_smooth() +
  scale_color_manual(values = c("blue","red")) +
  scale_fill_manual(values  = c("blue","red")) +
  gg_options() +
  labs(
    x = "Temperature (°C)",
    y = "Assimilation Efficiency (%)",
    title = "Sea Urchin Condos"
  ) +
  guides(
    color = guide_legend(title = expression("pCO"[2]),
                         override.aes = list(fill = c("blue","red"), color = c("blue","red"))),
    fill  = guide_legend(title = expression("pCO"[2]),
                         override.aes = list(fill = c("blue","red"), color = c("blue","red")))
  ) +
  theme(legend.position = "top")
indiv_assim_eff_plot

hist(dat_processed$assimilation_efficiency_proportion)

# Set general computing params ------------------------------------------------
cores        <- 26
chains       <- 4
iterations   <- 18000
warmup       <- 2000
control_list <- list(max_treedepth = 15, adapt_delta = 0.95)
backend      <- "cmdstanr"
seed         <- 123

# Beta regression with logit mean and log precision (phi)
brms_fam <- Beta(link = "logit", link_phi = "log")

# Compact working df ----------------------------------------------------------
d <- dat_processed %>%
  dplyr::select(assimilation_efficiency_proportion, spheroid_volume_ml, temp, pco2, tank_id) %>%
  drop_na()

# Candidate model formulae (ALL include phi submodel) --------------------------
# NOTE: distributional submodels must be specified inside bf(..., phi ~ ...)
#       Do NOT define a separate bf(phi ~ ...) and add it with `+`.

form_full_by <- bf(
  assimilation_efficiency_proportion ~
    log(spheroid_volume_ml) +
    gp(temp, k = 6, c = 5/4, by = pco2) +
    (1 | tank_id),
  phi ~ 1 + temp + pco2 + (1 | tank_id)
)

form_gp_noby <- bf(
  assimilation_efficiency_proportion ~
    log(spheroid_volume_ml) +
    gp(temp, k = 6, c = 5/4) +
    (1 | tank_id),
  phi ~ 1 + temp + pco2 + (1 | tank_id)
)

form_lin_int <- bf(
  assimilation_efficiency_proportion ~
    log(spheroid_volume_ml) +
    temp * pco2 +
    (1 | tank_id),
  phi ~ 1 + temp + pco2 + (1 | tank_id)
)

form_lin_add <- bf(
  assimilation_efficiency_proportion ~
    log(spheroid_volume_ml) +
    temp + pco2 +
    (1 | tank_id),
  phi ~ 1 + temp + pco2 + (1 | tank_id)
)

form_null <- bf(
  assimilation_efficiency_proportion ~
    log(spheroid_volume_ml) +
    (1 | tank_id),
  phi ~ 1 + temp + pco2 + (1 | tank_id)
)

model_set <- list(
  full_by   = form_full_by,
  gp_noby   = form_gp_noby,
  lin_int   = form_lin_int,
  lin_add   = form_lin_add,
  null_mod  = form_null
)

# Priors: "default + replace only what I care about" --------------------------
# Strategy:
#   1) pull default_prior() for the specific candidate model
#   2) convert to brms prior object
#   3) append overrides LAST (last matching prior wins)

priors_override <- c(
  # Mean model (mu) priors 
  prior(normal(0, 10), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "mu"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd", dpar = "mu"),
  
  # Phi (precision) submodel priors (dpar = "phi")
  prior(normal(log(70), 0.75), class = "Intercept", dpar = "phi"),
  prior(normal(0, 0.5),        class = "b",         dpar = "phi"),
  prior(student_t(3, 0, 1), lb = 0, class = "sd",   dpar = "phi")
)

priors_for <- function(formula, data, family) {
  
  dp <- default_prior(formula, data = data, family = family)
  
  # Identify which default priors to REPLACE.
  # override:
  #   - mu: Intercept, b, sd
  #   - phi: Intercept, b, sd
  dp_trim <- dp %>%
    dplyr::filter(
      !(
        # mean model (mu)
        (dpar == ""    & class %in% c("b", "Intercept", "sd")) |
          # phi submodel
          (dpar == "phi" & class %in% c("b", "Intercept", "sd"))
      )
    )
  
  # Override set 
  overrides <- priors_override
  
  c(dp_trim, overrides)
}


# Sanity check: show GP defaults being pulled (optional)
lapply(model_set, function(fm) {
  dp <- default_prior(fm, data = d, family = brms_fam)
  dp %>% dplyr::filter(class %in% c("sdgp", "lscale")) %>% print()
})

# Fit candidate models ---------------------------------------------------------
fit_all <- list()

for (nm in names(model_set)) {
  message("\nFitting model: ", nm)
  
  fit_all[[nm]] <- brm(
    formula      = model_set[[nm]],
    data         = d,
    family       = brms_fam,
    prior        = priors_for(model_set[[nm]], data = d, family = brms_fam),
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
    file         = paste0("./Model_Output/assimilation_efficiency_temp_oa_", nm)
  )
}

# Model selection based on LOO stacking weights --------------------------------
fit_all <- lapply(fit_all, add_criterion, criterion = "loo", reloo = FALSE)
loos    <- lapply(fit_all, \(x) x$criteria$loo)

w_stack <- brms::loo_model_weights(loos, method = "stacking")
w_stack
# Method: stacking
# 
#   weight
# full_by  0.937 # <- favored
# gp_noby  0.043 
# lin_int  0.000 
# lin_add  0.019 
# null_mod 0.000 

w_pbma  <- brms::loo_model_weights(loos, method = "pseudobma")
w_pbma
# Method: pseudo-BMA+ with Bayesian bootstrap
# 
#   weight
# full_by  0.357 # <- favored
# gp_noby  0.123 
# lin_int  0.200 
# lin_add  0.190 
# null_mod 0.130 

# Pareto-k summaries
pk <- lapply(loos, \(x) x$diagnostics$pareto_k)
sapply(pk, \(v) table(cut(v,
                          breaks = c(-Inf, 0.5, 0.7, 1, Inf),
                          labels = c("good (≤0.5)", "ok (≤0.7)", "bad (≤1)", "very bad (>1)")
)))
# full_by gp_noby lin_int lin_add null_mod
# good (≤0.5)        73      72      74      73       73
# ok (≤0.7)           5       6       4       5        5
# bad (≤1)            2       2       1       2        2
# very bad (>1)       0       0       1       0        0

# Choose preferred model (example: stacking winner)
selected_name  <- names(which.max(w_stack))
selected_name
# [1] "full_by"
selected_model <- fit_all[[selected_name]]
selected_model
# Family: beta 
# Links: mu = logit; phi = log 
# Formula: assimilation_efficiency_proportion ~ log(spheroid_volume_ml) + gp(temp, k = 6, c = 5/4, by = pco2) + (1 | tank_id) 
# phi ~ 1 + temp + pco2 + (1 | tank_id)
# Data: d (Number of observations: 80) 
# Draws: 4 chains, each with iter = 18000; warmup = 2000; thin = 1;
# total post-warmup draws = 64000
# 
# Gaussian Process Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat
# sdgp(gptemppco2600μatm)        1.40      0.86     0.31     3.58 1.00
# sdgp(gptemppco21200μatm)       0.89      0.75     0.05     2.89 1.00
# lscale(gptemppco2600μatm)      0.40      0.12     0.22     0.69 1.00
# lscale(gptemppco21200μatm)     0.47      0.19     0.23     0.95 1.00
# Bulk_ESS Tail_ESS
# sdgp(gptemppco2600μatm)       28216    22058
# sdgp(gptemppco21200μatm)      20945    24281
# lscale(gptemppco2600μatm)     52292    39538
# lscale(gptemppco21200μatm)    47626    36991
# 
# Multilevel Hyperparameters:
#   ~tank_id (Number of levels: 10) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
# sd(Intercept)         0.25      0.22     0.01     0.82 1.00    17046
# sd(phi_Intercept)     1.06      0.38     0.42     1.92 1.00    21405
# Tail_ESS
# sd(Intercept)        28903
# sd(phi_Intercept)    22305
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
# Intercept                 3.98      0.97     1.97     5.84 1.00    37195
# phi_Intercept            -0.08      1.72    -3.34     3.55 1.00    35651
# logspheroid_volume_ml    -0.45      0.21    -0.87    -0.05 1.00    72341
# phi_temp                  0.26      0.11     0.03     0.47 1.00    35890
# phi_pco21200μatm          0.30      0.42    -0.56     1.11 1.00    54117
# Tail_ESS
# Intercept                28492
# phi_Intercept            38840
# logspheroid_volume_ml    49851
# phi_temp                 38522
# phi_pco21200μatm         49636
# 
# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

bayes_R2(selected_model)
# Estimate  Est.Error      Q2.5     Q97.5
# R2 0.2953068 0.09556297 0.1172975 0.4804096

# Posterior predictive checks (use density overlays to diagnose peak/tails)
pp_check(selected_model, type = "dens_overlay", ndraws = 200) +
  coord_cartesian(xlim = c(0, 1))
pp_check(selected_model, type = "ecdf_overlay", ndraws = 200) +
  coord_cartesian(xlim = c(0, 1))

