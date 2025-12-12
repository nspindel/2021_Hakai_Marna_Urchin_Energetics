###############################################################################
# Analysis code for 2021 sea urchin heatwaves experiment at the Quadra Island # 
# Ecological Observatory - Hakai Institute                                    #
# Bayesian Thermal Performance Curve (TPC) estimation using `brms`            #
# Estimation of consumption rate.                                             #
#                                                                             #
# Author: Nathan B. Spindel                                                   #
##############################################################################
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
file_refit   <- "on_change"
backend      <- "cmdstanr"
seed         <- 123
brms_fam     <- Gamma(link = "log")

# unlink("./Model_Output/consumption_temp_oa_*", force = TRUE)


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
  per_cap_feed_per_day ~ gp(temp_deg_c, k = 6, c = 5/4, by = pco2_uatm) + (1 | tank_id)
  # shape ~ 1 + pco2_uatm + temp_deg_c
)

# Reduced GP (single surface, no by)
form_gp_noby <- brms::bf(
  per_cap_feed_per_day ~ gp(temp_deg_c, k = 6, c = 5/4) + (1 | tank_id)
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
    file_refit   = "on_change",
    file         = paste0("./Model_Output/consumption_temp_oa_", nm)
  )
}

brms_formula_full_const_dyn <- bf(per_cap_feed_per_day ~ 
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
  formula = brms_formula_full_const_dyn, 
  data = data_agg_feed_const_temp_dynamic_temp,
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
  file = "./Model_Output/consumption_temp_const_dynamic_full"
)

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
# good (≤0.5)          191     188     189     189      189
# ok (≤0.7)              0       1       2       0        0
# bad (≤1)               1       3       1       3        3
# very bad (>1)          0       0       0       0        0

selected_model <- fit_all$full
selected_model
# Family: gamma 
# Links: mu = log 
# Formula: per_cap_feed_per_day ~ gp(temp_deg_c, k = 6, c = 5/4, by = pco2_uatm) + (1 | tank_id) 
# Data: d (Number of observations: 192) 
# Draws: 4 chains, each with iter = 18000; warmup = 2000; thin = 1;
# total post-warmup draws = 64000
# 
# Gaussian Process Hyperparameters:
#   Estimate Est.Error l-95% CI u-95% CI
# sdgp(gptemp_deg_cpco2_uatm600)        0.71      0.40     0.23     1.68
# sdgp(gptemp_deg_cpco2_uatm1200)       0.62      0.68     0.11     1.99
# lscale(gptemp_deg_cpco2_uatm600)      0.31      0.15     0.12     0.66
# lscale(gptemp_deg_cpco2_uatm1200)     0.40      0.31     0.12     1.27
# Rhat Bulk_ESS Tail_ESS
# sdgp(gptemp_deg_cpco2_uatm600)    1.00     1887     4914
# sdgp(gptemp_deg_cpco2_uatm1200)   1.01      551      164
# lscale(gptemp_deg_cpco2_uatm600)  1.00    30758    32046
# lscale(gptemp_deg_cpco2_uatm1200) 1.01      551      163
# 
# Multilevel Hyperparameters:
#   ~tank_id (Number of levels: 24) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.09      0.06     0.00     0.24 1.00    21305    24149
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.01      0.35    -0.90     0.55 1.01      637      161
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# shape     2.78      0.28     2.27     3.36 1.00     2451     2982
# 
# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
bayes_R2(selected_model)
# Estimate  Est.Error      Q2.5     Q97.5
# R2 0.4544242 0.07923617 0.2889653 0.5879063

# Posterior predictive check for favored model
pp_check(selected_model, ndraws = 50)
# Save session info for provenance
sink("./Model_Output/sessionInfo_Consumption.txt"); print(sessionInfo()); sink()

# Figure 2: Visualize conditional effects on consumption ----
ce_95_pco2_comparison <- brms::conditional_effects(selected_model, robust = F)

df_ce_temp_pco2_plot_95 <- ce_95_pco2_comparison$`temp_deg_c:pco2_uatm` %>%
  dplyr::mutate(
    pco2_plot = factor(pco2_uatm, labels = c("600 μatm", "1200 μatm"))
  )

ce_80_pco2_comparison <- brms::conditional_effects(selected_model, robust = F, prob = 0.8)

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

fit_feed_const_temp_dynamic_temp <- readRDS("./Model_Output/consumption_temp_const_dynamic_full.rds")

ce_95_const_dynamic_comparison <- brms::conditional_effects(fit_feed_const_temp_dynamic_temp, robust = F)

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


ce_80_const_dynamic_comparison <- conditional_effects(fit_feed_const_temp_dynamic_temp, robust = F, prob = 0.8)

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
  ) 

p_agg_feed

ggsave("./Figures/Per_Cap_Consumption.pdf", p_agg_feed,
       units = "in", width = 8, height = 8)
ggsave("./Figures/Per_Cap_Consumption.png", p_agg_feed,
       units = "in", width = 8, height = 8)
ggsave("./Figures/Per_Cap_Consumption.jpg", p_agg_feed,
       units = "in", width = 8, height = 8)

# Pairwise comparisons --------------------------------------------------------
selected_model <- readRDS("./Model_Output/consumption_temp_oa_full_distr.rds")
selected_model$formula
# per_cap_feed_per_day ~ gp(temp_deg_c, k = 10, c = 5/4, by = pco2_uatm) + (1 | tank_id) 

ems <- emmeans(selected_model, ~ pco2_uatm | temp_deg_c,
               at = list(temp_deg_c = c(10,13,16,17,18,20)))

cont_pco2 <- contrast(ems, method = "revpairwise", by = "pco2_uatm", 
                      type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_pco2, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")
# pco2_uatm = 600:
#   contrast                    ratio lower.HPD upper.HPD
# temp_deg_c13 / temp_deg_c10 1.810     1.063      2.67
# temp_deg_c16 / temp_deg_c10 2.318     1.339      3.41
# temp_deg_c16 / temp_deg_c13 1.295     0.925      1.67
# temp_deg_c17 / temp_deg_c10 2.275     1.380      3.22
# temp_deg_c17 / temp_deg_c13 1.284     0.849      1.74
# temp_deg_c17 / temp_deg_c16 0.989     0.849      1.12
# temp_deg_c18 / temp_deg_c10 2.256     1.395      3.18
# temp_deg_c18 / temp_deg_c13 1.279     0.800      1.79
# temp_deg_c18 / temp_deg_c16 0.987     0.735      1.22
# temp_deg_c18 / temp_deg_c17 0.995     0.872      1.12
# temp_deg_c20 / temp_deg_c10 2.409     1.335      3.60
# temp_deg_c20 / temp_deg_c13 1.357     0.833      1.93
# temp_deg_c20 / temp_deg_c16 1.056     0.681      1.48
# temp_deg_c20 / temp_deg_c17 1.068     0.715      1.46
# temp_deg_c20 / temp_deg_c18 1.071     0.786      1.42
# 
# pco2_uatm = 1200:
#   contrast                    ratio lower.HPD upper.HPD
# temp_deg_c13 / temp_deg_c10 1.479     0.967      2.12
# temp_deg_c16 / temp_deg_c10 1.753     0.990      2.51
# temp_deg_c16 / temp_deg_c13 1.189     0.896      1.49
# temp_deg_c17 / temp_deg_c10 1.720     0.991      2.40
# temp_deg_c17 / temp_deg_c13 1.174     0.828      1.54
# temp_deg_c17 / temp_deg_c16 0.986     0.872      1.09
# temp_deg_c18 / temp_deg_c10 1.653     0.987      2.29
# temp_deg_c18 / temp_deg_c13 1.134     0.767      1.52
# temp_deg_c18 / temp_deg_c16 0.953     0.762      1.15
# temp_deg_c18 / temp_deg_c17 0.964     0.860      1.06
# temp_deg_c20 / temp_deg_c10 1.486     0.929      2.15
# temp_deg_c20 / temp_deg_c13 1.021     0.652      1.42
# temp_deg_c20 / temp_deg_c16 0.864     0.563      1.17
# temp_deg_c20 / temp_deg_c17 0.875     0.606      1.15
# temp_deg_c20 / temp_deg_c18 0.904     0.679      1.13
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95 
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95
cont_temp <- contrast(ems, method = "revpairwise", by = "temp_deg_c", 
                      type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_temp, 
                     prob = 0.95, 
                     point.est = mean,
                     type = "response")
# temp_deg_c = 10:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600  2.14     1.160      3.26
# 
# temp_deg_c = 13:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600  1.74     1.099      2.43
# 
# temp_deg_c = 16:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600  1.59     1.107      2.09
# 
# temp_deg_c = 17:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600  1.58     1.165      2.04
# 
# temp_deg_c = 18:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600  1.53     1.113      1.99
# 
# temp_deg_c = 20:
#   contrast                     ratio lower.HPD upper.HPD
# pco2_uatm1200 / pco2_uatm600  1.31     0.831      1.86
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95 
cont_posterior_temp <- gather_emmeans_draws(cont_temp) %>%
  mutate(
    temp_deg_c = as.factor(temp_deg_c)
  )
cont_posterior_pco2 <- gather_emmeans_draws(cont_pco2) 

ems_overall <- emmeans(selected_model, ~ pco2_uatm | temp_deg_c)

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
  )

plot_emm_fit

ggsave("./Figures/Supplemental/PWC_Consumption.jpg", plot_emm_fit,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Consumption.pdf", plot_emm_fit,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Consumption.png", plot_emm_fit,
       width = 6, height = 8, units = "in")

fit_consumption_dynamic_temp <- readRDS("./Model_Output/consumption_temp_const_dynamic_full.rds")
selected_model <- fit_consumption_dynamic_temp
selected_model$formula
# fce ~ log(init_diam_mm) + temp_deg_c + (1 | tank_id) 

ems <- emmeans(selected_model, ~ temp_deg_c)

cont <- contrast(ems, method = "revpairwise", by = NULL, 
                 type = "response", ratios = TRUE)

emmeans::hpd.summary(cont, 
                     prob = 0.95, 
                     point.est = mean)
# contrast                          ratio lower.HPD upper.HPD
# 13 / 10                           2.298     1.152      3.57
# 16 / 10                           2.430     1.211      3.79
# 16 / 13                           1.099     0.553      1.71
# (18-14 La Niña) / 10              3.587     2.038      5.29
# (18-14 La Niña) / 13              1.622     0.902      2.40
# (18-14 La Niña) / 16              1.535     0.868      2.29
# 17 / 10                           3.053     1.532      4.73
# 17 / 13                           1.380     0.699      2.15
# 17 / 16                           1.306     0.674      2.05
# 17 / (18-14 La Niña)              0.867     0.487      1.29
# 18 / 10                           2.443     1.243      3.80
# 18 / 13                           1.105     0.561      1.72
# 18 / 16                           1.046     0.524      1.62
# 18 / (18-14 La Niña)              0.694     0.393      1.03
# 18 / 17                           0.831     0.417      1.29
# 20 / 10                           2.904     1.475      4.53
# 20 / 13                           1.313     0.651      2.03
# 20 / 16                           1.243     0.634      1.96
# 20 / (18-14 La Niña)              0.825     0.462      1.22
# 20 / 17                           0.988     0.513      1.55
# 20 / 18                           1.234     0.607      1.91
# (21-18 El Niño) / 10              3.315     1.869      4.91
# (21-18 El Niño) / 13              1.498     0.844      2.21
# (21-18 El Niño) / 16              1.419     0.799      2.12
# (21-18 El Niño) / (18-14 La Niña) 0.941     0.609      1.32
# (21-18 El Niño) / 17              1.127     0.650      1.68
# (21-18 El Niño) / 18              1.408     0.778      2.08
# (21-18 El Niño) / 20              1.186     0.666      1.76
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95
cont_posterior <- gather_emmeans_draws(cont) 

cont_posterior_filtered <- cont_posterior %>%
  filter(
    contrast %in% c(
      "(18-14 La Niña) - 16",
      "(21-18 El Niño) - 20"
    )
  ) %>%
  mutate(
    contrast_formatted = ifelse(contrast == "(18-14 La Niña) - 16",
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
    title = "Per Capita Consumption",
    subtitle = "constant mean vs dynamic"
  ) +
  gg_options() +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = c("#DEDEDE"))
  );plot_emm_fit_tpc_mhw

ggsave("./Figures/Supplemental/PWC_Consumption_const_dynamic_temp.jpg", plot_emm_fit_tpc_mhw,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Consumption_const_dynamic_temp.pdf", plot_emm_fit_tpc_mhw,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_Consumption_const_dynamic_temp.png", plot_emm_fit_tpc_mhw,
       width = 6, height = 8, units = "in")