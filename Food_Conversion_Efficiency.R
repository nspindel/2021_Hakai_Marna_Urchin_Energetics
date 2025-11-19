###############################################################################
# Analysis code for 2021 sea urchin heatwaves experiment at the Quadra Island # 
# Ecological Observatory - Hakai Institute                                    #
#                                                                             #
# Food conversion efficiency (FCE)                                            #
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
data_fce <- read.csv("./Data/FCE_Data.csv") %>%
mutate(
  tank_id = as.factor(tank_id),
  pit_tag_id = as.factor(pit_tag_id),
  pco2_uatm = factor(pco2_h_l, levels = c("600 µatm", "1200 µatm")),
  feed_est_total_g = daily_consumption_g_dry_est*total_experiment_days,
  delta_gonad_mass_g = total_wet_gonad_mass_g - init_nad_est
) %>% 
  mutate(    
    fce = (delta_gonad_mass_g/feed_est_total_g),#*100, # express as %
    temp_model = as.numeric(case_when(
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

data_fce_const_temp_oa <- data_fce %>%
  filter(
    temp_deg_c != "21-18" 
    & temp_deg_c != "18-14"
    & fce > 0
  ) %>%
  select(
    tank_id,
    pit_tag_id,
    init_diam_mm,
    init_ht_mm,
    temp_deg_c,
    temp_model,
    pco2_uatm,
    fce
  )

data_fce_const_dynamic_temp <- data_fce %>%
  filter(
    pco2_uatm == "600 µatm"
    & fce > 0
  ) 

data_fce_dynamic_temp <- data_fce_const_dynamic_temp %>%
  filter(
    temp_deg_c %in% c("21-18", "18-14")
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
brms_fam     <- Beta(link = "logit", link_phi = "log")

# unlink("./Model_Output/fce_temp_oa_*", force = TRUE)

# Priors ----
priors_common <- c(
  prior(normal(0, 0.1), class = "b"),
  prior(student_t(3, 2, 2.5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01,0.01),lb = 0, class = "phi")
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
  prior(normal(0,0.1), class = "b"),
  prior(student_t(3, 0.2, 5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), lb = 0, class = "sd"),
  prior(gamma(0.01,0.01),lb = 0, class = "phi")
)

# Candidate model formulae for constant temperatures ----
# Full model
form_full <- brms::bf(
  fce ~
    log(init_diam_mm) +
    gp(temp_deg_c, k = 5, c = 5/4, by = pco2_uatm) +
    (1 | tank_id)
)

# Reduced models chosen to test specific hypotheses with principled nesting
form_gp_noby <- brms::bf(
  fce ~
    log(init_diam_mm) +
    gp(temp_deg_c, k = 5, c = 5/4) +
    (1 | tank_id)
)
# (Removes the pCO2-specific GP surfaces; retains nonparametric thermal response.)

form_lin_int <- brms::bf(
  fce ~
    log(init_diam_mm) + temp_deg_c * pco2_uatm +
    (1 | tank_id)
)
# (Parametric T × pCO2 interaction; no GP.)

form_lin_add <- brms::bf(
  fce ~
    log(init_diam_mm) + temp_deg_c + pco2_uatm +
    (1 | tank_id)
)
# (Additive parametric effects; no GP.)

form_size_only <- brms::bf(
  fce ~
    log(init_diam_mm) + (1 | tank_id)
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
d <- data_fce_const_temp_oa %>%
  select(fce, init_diam_mm, temp_deg_c, pco2_uatm, tank_id) %>%
  mutate(temp_deg_c = as.numeric(temp_deg_c)) %>%
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
    file         = paste0("./Model_Output/fce_temp_oa_", nm)
  )
}

brms_formula_full_dyn <- bf(fce ~ 
                                    log(init_diam_mm) + 
                                    temp_deg_c + 
                                    (1|tank_id))

# fit model
fit_fce_temp_dynamic_temp <- brm(
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
  file = "./Model_Output/fce_temp_dynamic_full"
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
# full      0.896 # <- favored
# gp_noby   0.020 
# lin_int   0.036 
# lin_add   0.035 
# size_only 0.013

# Pareto-k summaries
pk <- lapply(loos, \(x) x$diagnostics$pareto_k)
sapply(pk, \(v) table(cut(v,
                          breaks = c(-Inf, 0.5, 0.7, 1, Inf),
                          labels = c("good (≤0.5)", "ok (≤0.7)", "bad (≤1)", "very bad (>1)")
)))
#               full gp_noby lin_int lin_add size_only
# good (≤0.5)    346     346     346     346       346
# ok (≤0.7)        0       0       0       0         0
# bad (≤1)         0       0       0       0         0
# very bad (>1)    0       0       0       0         0

selected_model <- fit_all$full
selected_model$formula
# fce ~ log(init_diam_mm) + gp(temp_deg_c, k = 5, c = 5/4, by = pco2_uatm) + (1 | tank_id) 
bayes_R2(selected_model)
# Estimate  Est.Error      Q2.5     Q97.5
# R2 0.4128564 0.03632055 0.338815 0.4803027

# Posterior predictive check for favored model
pp_check(selected_model, ndraws = 50)

# Save session info for provenance
sink("./Model_Output/sessionInfo_FCE.txt"); print(sessionInfo()); sink()


# Figure 4: Visualize conditional effects on FCE ----
fit_fce_const_temp_oa <- readRDS("./Model_Output/fce_temp_oa_full.rds")
selected_model <- fit_fce_const_temp_oa
selected_model$formula
# fce ~ log(init_diam_mm) + gp(temp_deg_c, k = 5, c = 5/4, by = pco2_uatm) + (1 | tank_id) 

ce_95_fce_constant <- conditional_effects(selected_model, robust = F, prob = 0.95)
ce_95_fce_constant_plot <- ce_95_fce_constant$`temp_deg_c:pco2_uatm` %>%
  mutate(
    pco2_uatm = factor(pco2_uatm, levels = c("600 µatm", "1200 µatm"))
    )

ce_80_fce_constant <- conditional_effects(selected_model, robust = F, prob = 0.80)
ce_80_fce_constant_plot <- ce_80_fce_constant$`temp_deg_c:pco2_uatm` %>%
  mutate(
    pco2_uatm = factor(pco2_uatm, levels = c("600 µatm", "1200 µatm"))
  )

fit_fce_dynamic_temp <- readRDS("./Model_Output/fce_temp_dynamic_full.rds")

ce_95_fce_dynamic_temp <- conditional_effects(fit_fce_dynamic_temp, 
                                                         robust = F, prob = 0.95)
ce_80_fce_dynamic_temp <- conditional_effects(fit_fce_dynamic_temp, 
                                                           robust = F, prob = 0.80)

plot_temps_dynamic_95 <- ce_95_fce_dynamic_temp$temp_deg_c %>%
  filter(temp_deg_c %in% c("18-14", "21-18")) %>%
  mutate(
    temp_label = ifelse(temp_deg_c == "18-14", "La Niña, 600 μatm", "El Niño, 600 μatm")) %>%
  mutate(temp_label = factor(temp_label, levels = c("La Niña, 600 μatm", "El Niño, 600 μatm"))) %>%
  mutate(
    temp_plot = case_when(
      temp_label == "La Niña, 600 μatm" ~ 16,
      temp_label == "El Niño, 600 μatm" ~ 20,
    )
  )

plot_temps_dynamic_80 <- ce_80_fce_dynamic_temp$temp_deg_c %>%
  filter(temp_deg_c %in% c("18-14", "21-18")) %>%
  mutate(
    temp_label = ifelse(temp_deg_c == "18-14", "La Niña, 600 μatm", "El Niño, 600 μatm")) %>%
  mutate(temp_label = factor(temp_label, levels = c("La Niña, 600 μatm", "El Niño, 600 μatm"))) %>%
  mutate(
    temp_plot = case_when(
      temp_label == "La Niña, 600 μatm" ~ 16,
      temp_label == "El Niño, 600 μatm" ~ 20,
    )
  )

empirical_means_fce <- data_fce_const_temp_oa %>%
  group_by(tank_id, temp_deg_c, pco2_uatm) %>%
  dplyr::summarise(
    emp_mean_fce = mean(fce)
  ) %>%
  dplyr::mutate(temp_plot = as.numeric(temp_deg_c))

empirical_means_fce_dyn <- data_fce_dynamic_temp %>%
  filter(temp_deg_c %in% c("18-14", "21-18") & fce > 0) %>%
  dplyr::mutate(
    pco2_uatm = factor(temp_deg_c, labels = c("La Niña, 600 μatm", "El Niño, 600 μatm"))
  ) %>%
  dplyr::group_by(tank_id, temp_deg_c, pco2_uatm) %>%
  dplyr::summarise(
    emp_mean_fce = mean(fce)
  ) %>%
  dplyr::mutate(
    temp_plot = ifelse(temp_deg_c == "18-14", 16, 20)
  )

df_empirical_means <- data.frame()
df_empirical_means <- rbind(empirical_means_fce, empirical_means_fce_dyn) #%>%
  # dplyr::mutate(
  #   plot_label = factor(dplyr::case_when(
  #     temp_deg_c %in% c("10","13","16","17","18","20") ~ pco2_uatm,
  #     temp_deg_c == "18-14" ~ "La Niña, 600 μatm",
  #     temp_deg_c == "21-18" ~ "El Niño, 600 μatm"
  #   ), levels = c("600 μatm", "1200 μatm", "La Niña, 600 μatm", "El Niño, 600 μatm"))
  # )

p_fce_pco2_temp <- ggplot() +
  geom_line(data = ce_95_fce_constant_plot,
            aes(x = temp_deg_c, y = estimate__,
                color = pco2_uatm, linetype = pco2_uatm),
            linewidth = 1.5, show.legend = F) +
  geom_ribbon(data = ce_95_fce_constant_plot,
              aes(x = temp_deg_c, y = estimate__,
                  ymin = lower__, ymax = upper__,
                  fill = pco2_uatm),
              alpha = 0.2) +
  geom_ribbon(data = ce_80_fce_constant_plot,
              aes(x = temp_deg_c, y = estimate__,
                  ymin = lower__, ymax = upper__,
                  fill = pco2_uatm),
              alpha = 0.4) +
  geom_linerange(data = plot_temps_dynamic_95 ,
                 aes(x = temp_plot, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = temp_label),
                 linewidth = 1) +
  geom_linerange(data = plot_temps_dynamic_80 ,
                 aes(x = temp_plot, y = estimate__,
                     ymin = lower__, ymax = upper__,
                     color = temp_label),
                 linewidth = 3) +
  geom_jitter(data = df_empirical_means,
              aes(x = temp_plot, y = emp_mean_fce,
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
    y = "food conversion efficiency",
    shape = NULL,
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.8,0.8))

p_fce_pco2_temp

ggsave("./Figures/FCE.jpg", p_fce_pco2_temp,
       units = "in", width = 8, height = 8)
ggsave("./Figures/FCE.pdf", p_fce_pco2_temp,
       units = "in", width = 8, height = 8)
ggsave("./Figures/FCE.png", p_fce_pco2_temp,
       units = "in", width = 8, height = 8)

# Post hoc on full model for constant temps -----
fit_fce_const_temp_oa <- readRDS("./Model_Output/fce_temp_oa_full.rds")
selected_model <- fit_fce_const_temp_oa
selected_model$formula
# fce ~ log(init_diam_mm) + gp(temp_deg_c, k = 5, c = 5/4, by = pco2_uatm) + (1 | tank_id) 
ems <- emmeans(selected_model, ~ pco2_uatm | temp_deg_c,
               at = list(temp_deg_c = c(10,13,16,17,18,20)))

cont_pco2 <- contrast(ems, method = "revpairwise", by = "pco2_uatm", 
                      type = "response")

emmeans::hpd.summary(cont_pco2, 
                     prob = 0.95, 
                     point.est = mean)
# pco2_uatm = 600 µatm:
#   contrast                    odds.ratio lower.HPD upper.HPD
# temp_deg_c13 / temp_deg_c10      0.915     0.741     1.096
# temp_deg_c16 / temp_deg_c10      0.614     0.491     0.745
# temp_deg_c16 / temp_deg_c13      0.672     0.572     0.779
# temp_deg_c17 / temp_deg_c10      0.509     0.409     0.610
# temp_deg_c17 / temp_deg_c13      0.558     0.461     0.661
# temp_deg_c17 / temp_deg_c16      0.829     0.777     0.881
# temp_deg_c18 / temp_deg_c10      0.418     0.336     0.503
# temp_deg_c18 / temp_deg_c13      0.459     0.370     0.555
# temp_deg_c18 / temp_deg_c16      0.683     0.593     0.771
# temp_deg_c18 / temp_deg_c17      0.823     0.764     0.883
# temp_deg_c20 / temp_deg_c10      0.299     0.222     0.380
# temp_deg_c20 / temp_deg_c13      0.328     0.244     0.416
# temp_deg_c20 / temp_deg_c16      0.489     0.364     0.617
# temp_deg_c20 / temp_deg_c17      0.589     0.454     0.720
# temp_deg_c20 / temp_deg_c18      0.714     0.587     0.836
# 
# pco2_uatm = 1200 µatm:
#   contrast                    odds.ratio lower.HPD upper.HPD
# temp_deg_c13 / temp_deg_c10      1.015     0.849     1.196
# temp_deg_c16 / temp_deg_c10      0.915     0.726     1.113
# temp_deg_c16 / temp_deg_c13      0.902     0.774     1.025
# temp_deg_c17 / temp_deg_c10      0.865     0.684     1.046
# temp_deg_c17 / temp_deg_c13      0.854     0.709     1.009
# temp_deg_c17 / temp_deg_c16      0.946     0.889     1.004
# temp_deg_c18 / temp_deg_c10      0.821     0.653     1.013
# temp_deg_c18 / temp_deg_c13      0.812     0.653     1.005
# temp_deg_c18 / temp_deg_c16      0.899     0.786     1.010
# temp_deg_c18 / temp_deg_c17      0.949     0.883     1.009
# temp_deg_c20 / temp_deg_c10      0.776     0.597     1.003
# temp_deg_c20 / temp_deg_c13      0.768     0.581     1.005
# temp_deg_c20 / temp_deg_c16      0.851     0.655     1.035
# temp_deg_c20 / temp_deg_c17      0.898     0.728     1.051
# temp_deg_c20 / temp_deg_c18      0.944     0.816     1.057
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log odds ratio scale 
# HPD interval probability: 0.95 
cont_posterior_pco2 <- gather_emmeans_draws(cont_pco2) 

cont_temp <- contrast(ems, method = "revpairwise", by = "temp_deg_c", 
                      type = "response", ratios = TRUE)

emmeans::hpd.summary(cont_temp, 
                     prob = 0.95, 
                     point.est = mean)

# temp_deg_c = 10:
#   contrast             odds.ratio lower.HPD upper.HPD
# 1200 µatm / 600 µatm      0.418     0.324     0.522
# 
# temp_deg_c = 13:
#   contrast             odds.ratio lower.HPD upper.HPD
# 1200 µatm / 600 µatm      0.463     0.365     0.565
# 
# temp_deg_c = 16:
#   contrast             odds.ratio lower.HPD upper.HPD
# 1200 µatm / 600 µatm      0.621     0.517     0.733
# 
# temp_deg_c = 17:
#   contrast             odds.ratio lower.HPD upper.HPD
# 1200 µatm / 600 µatm      0.708     0.593     0.822
# 
# temp_deg_c = 18:
#   contrast             odds.ratio lower.HPD upper.HPD
# 1200 µatm / 600 µatm      0.817     0.675     0.963
# 
# temp_deg_c = 20:
#   contrast             odds.ratio lower.HPD upper.HPD
# 1200 µatm / 600 µatm      1.089     0.795     1.393
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log odds ratio scale 
# HPD interval probability: 0.95 

cont_posterior_temp <- gather_emmeans_draws(cont_temp) %>%
  mutate(
    temp_deg_c = as.factor(temp_deg_c)
  )

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
    title = "Food Conversion Efficiency",
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

ggsave("./Figures/Supplemental/PWC_FCE_pco2.jpg", plot_emm_fit,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_FCE_pco2.pdf", plot_emm_fit,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_FCE_pco2.png", plot_emm_fit,
       width = 6, height = 8, units = "in")

# Post hoc constant vs dynamic temp comparisons --------------------------------------------------------
fit_fce_dynamic_temp <- readRDS("./Model_Output/fce_temp_const_dynamic_full.rds")
selected_model <- fit_fce_dynamic_temp
selected_model$formula
# fce ~ log(init_diam_mm) + temp_deg_c + (1 | tank_id) 

ems <- emmeans(selected_model, ~ temp_deg_c)

cont <- contrast(ems, method = "revpairwise", by = NULL, 
                 type = "response", ratios = TRUE)

emmeans::hpd.summary(cont, 
                     prob = 0.95, 
                     point.est = mean)
# contrast          odds.ratio lower.HPD upper.HPD
# 13 / 10                1.055     0.852      1.26
# 16 / 10                1.023     0.833      1.22
# 16 / 13                0.979     0.724      1.25
# 17 / 10                1.027     0.837      1.22
# 17 / 13                0.983     0.729      1.25
# 17 / 16                1.014     0.758      1.30
# 18 / 10                1.001     0.812      1.19
# 18 / 13                0.958     0.707      1.22
# 18 / 16                0.988     0.728      1.26
# 18 / 17                0.983     0.729      1.25
# (18-14) / 10           0.958     0.786      1.15
# (18-14) / 13           0.917     0.677      1.17
# (18-14) / 16           0.945     0.696      1.20
# (18-14) / 17           0.941     0.693      1.20
# (18-14) / 18           0.966     0.715      1.23
# 20 / 10                0.975     0.787      1.16
# 20 / 13                0.933     0.686      1.19
# 20 / 16                0.962     0.714      1.23
# 20 / 17                0.958     0.711      1.22
# 20 / 18                0.983     0.727      1.25
# 20 / (18-14)           1.027     0.761      1.30
# (21-18) / 10           0.935     0.757      1.11
# (21-18) / 13           0.895     0.657      1.15
# (21-18) / 16           0.923     0.683      1.18
# (21-18) / 17           0.919     0.677      1.17
# (21-18) / 18           0.943     0.702      1.20
# (21-18) / (18-14)      0.985     0.734      1.25
# (21-18) / 20           0.968     0.725      1.23
# 
# Point estimate displayed: mean 
# Results are back-transformed from the log odds ratio scale 
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
    title = "Food Conversion Efficiency",
    subtitle = "constant mean vs dynamic"
  ) +
  gg_options() +
  scale_x_continuous(breaks = pretty_breaks(n = 6), limits = c(-1.0,0.4)) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = c("#DEDEDE"))
  );plot_emm_fit_tpc_mhw

ggsave("./Figures/Supplemental/PWC_FCE_const_dynamic_temp.jpg", plot_emm_fit_tpc_mhw,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_FCE_const_dynamic_temp.pdf", plot_emm_fit_tpc_mhw,
       width = 6, height = 8, units = "in")
ggsave("./Figures/Supplemental/PWC_FCE_const_dynamic_temp.png", plot_emm_fit_tpc_mhw,
       width = 6, height = 8, units = "in")