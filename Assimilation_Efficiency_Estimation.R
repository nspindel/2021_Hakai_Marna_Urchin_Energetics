###############################################################################
# Analysis code for 2021 sea urchin heatwaves experiment at the Quadra Island # 
# Ecological Observatory - Hakai Institute                                    #
# Bayesian Thermal Performance Curve (TPC) estimation using `brms`            #
# Estimation of food conversion efficiency.                                   #
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
    temp = as.numeric(temp),
    pco2 = ifelse(pco2 == "l", "600 μatm", "1200 μatm"),
    pco2 = factor(pco2, levels = c("600 μatm","1200 μatm")),
    temp_pco2_interaction = paste0(temp, "_", pco2),
    bin_id = as.factor(bin_id),
    pit_id = as.factor(pit_id),
    semi_major_axis_mm = try(0.5 * diam_mm),
    semi_minor_axis_mm = try(0.5 * ht_mm),
    spheroid_volume_ml = 
      # Equation for volume of oblate spheroid = 4/3pia^2c in mL
      try((4/3 * pi * semi_major_axis_mm ^2 * semi_minor_axis_mm)/1000), 
    leftover_dry_food_g = food_wb_plus_food_g - food_wb_mass_g,
    dry_food_consumed_g = food_mass_g - leftover_dry_food_g,
    dry_feces_egested_g = feces_wb_plus_feces_g - feces_wb_mass_g,
    hours_feeding = as.numeric(hours_feeding),
    minutes_feeding = as.numeric(minutes_feeding),
    dry_food_consumed_per_day_g = dry_food_consumed_g/days_feeding,
    dry_food_consumed_per_hour_g = dry_food_consumed_g/hours_feeding,
    dry_food_consumed_per_min_g = dry_food_consumed_g/minutes_feeding,
    dry_feces_egested_per_day_g = dry_feces_egested_g/days_feeding,
    dry_feces_egested_per_hour_g = dry_feces_egested_g/hours_feeding,
    dry_feces_egested_per_min_g = dry_feces_egested_g/minutes_feeding,
    assimilation_efficiency_percent = 
      ((dry_food_consumed_per_day_g - dry_feces_egested_per_day_g)/dry_food_consumed_per_day_g)*100
)

feed_means <- dat_processed %>%
  data.frame()%>%
  group_by(temp, pco2)%>%
  do(boot_stat(.$dry_food_consumed_per_day_g))
feed_means

indiv_feed_plot <- ggplot(data = dat_processed, aes(x = temp, y = dry_food_consumed_per_day_g,  
                                                    fill = pco2,color = pco2))+
  geom_jitter(alpha = 0.2)+
  # geom_pointrange(data = feed_means, aes(x = temp, y = y, ymin = ymin, ymax = ymax),
  #                 position = position_jitter(width = 0.12))+
  geom_smooth()+
  scale_color_manual(values = c("blue","red"))+
  scale_fill_manual(values = c("blue","red"))+
  gg_options()+
  labs(
    x= "Temperature (°C)",
    y= expression("Consumption (g dry mass urchin"^-1*"day"^-1*")"),
    title = "Sea Urchin Condos"
  )+
  guides(
    color = guide_legend(
      title = expression("pCO"[2]),
      override.aes = list(fill = c("blue","red"), color = c("blue","red"))
    ),
    fill = guide_legend(
      title = expression("pCO"[2]),
      override.aes = list(fill = c("blue","red"), color = c("blue","red"))
    )
  )+
  theme(
    legend.position = "top"
  );indiv_feed_plot

assim_eff_means <- dat_processed %>%
  data.frame()%>%
  group_by(temp, pco2)%>%
  do(boot_stat(.$assimilation_efficiency_percent))
assim_eff_means

indiv_assim_eff_plot <- ggplot(data = dat_processed, aes(x = temp, y = assimilation_efficiency_percent,  
                                                    fill = pco2,color = pco2))+
  geom_jitter(alpha = 0.2)+
  # geom_pointrange(data = assim_eff_means, aes(x = temp, y = y, ymin = ymin, ymax = ymax),
  #                 position = position_jitter(width = 0.12))+
  geom_smooth()+
  scale_color_manual(values = c("blue","red"))+
  scale_fill_manual(values = c("blue","red"))+
  gg_options()+
  labs(
    x= "Temperature (°C)",
    y= "Assimilation Efficiency (%)",
    title = "Sea Urchin Condos"
  )+
  guides(
    color = guide_legend(
      title = expression("pCO"[2]),
      override.aes = list(fill = c("blue","red"), color = c("blue","red"))
    ),
    fill = guide_legend(
      title = expression("pCO"[2]),
      override.aes = list(fill = c("blue","red"), color = c("blue","red"))
    )
  )+
  theme(
    legend.position = "top"
  );indiv_assim_eff_plot


# 
# # brms_form <- bf(dry_food_consumed_per_day_g ~ log(spheroid_volume_ml) + 
#                   # t2(temp, by = pco2) + (1|tank_id))
# brms_form <- bf(dry_food_consumed_per_day_g~t2(temp,by = temp_pco2_interaction))                  
# 
# 
# brms_fam <- Gamma(link = "log")
# 
# # Get all parameters and parameters classes to define priors on
# (prior <- get_prior(formula = brms_form, 
#                     data = dat_processed, 
#                     family = brms_fam
# ))
# 
# # Generic weakly informative prior 
# # (https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations)
# prior$prior[2:11] <- "normal(0, 5)"
# 
# # Fit the model
# fit_condo <- brm(
#   formula = brms_form,
#   data = dat_processed,
#   family = brms_fam,
#   cores = getOption("mc.cores", 12),
#   prior = prior,
#   sample_prior = "yes",
#   iter = 2000,
#   control = list(
#     max_treedepth = 15,
#     adapt_delta = 0.9
#   ),
#   backend = "cmdstanr"
# )
# 
# # R-squared for Bayesian regression models
# bayes_R2(fit_condo)
# # LOO cross validation
# loo(fit_condo)
# # Posterior predictive check
# pp_check(fit_condo, ndraws = 100)
# # Save as RDS
# saveRDS(fit_gonad, file = "./Model_Fits_RDS/brms_gonad.rds")
