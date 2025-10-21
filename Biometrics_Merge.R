###############################################################################
# Merge code for biometrics files                                             #
#                                                                             #
# Author: Nathan B. Spindel                                                   #
###############################################################################

# Load libraries
libs <- c("tidyverse", "readxl", "janitor")
lapply(libs, require, character.only = TRUE)

# Define file paths (adjust if needed)
path_initial <- "./Data/Biometrics_Initial.xlsx"
path_final   <- "./Data/Biometrics_Final.xlsx"

# Read and clean both datasets
bio_initial <- read_excel(path_initial) %>%
  dplyr::mutate(pit_tag_id = as.factor(pit_tag_id)) %>%
  janitor::clean_names() %>%
  rename_with(~ paste0("initial_", .x), 
              .cols = -c(pit_tag_id, temp_deg_c, pco2_h_l, tank_id)) %>%
  rename(initial_date = initial_date)  # if date column already named 'date'

bio_final <- read_excel(path_final) %>%
  dplyr::mutate(pit_tag_id = as.factor(pit_tag_id)) %>%
  janitor::clean_names() %>%
  rename_with(~ paste0("final_", .x), 
              .cols = -c(pit_tag_id, temp_deg_c, pco2_h_l, tank_id)) %>%
  rename(final_date = final_date)  # if date column already named 'date'

# Join the two datasets by unique identifiers and treatment keys
bio_joined <- full_join(
  bio_initial,
  bio_final,
  by = c("pit_tag_id", "temp_deg_c", "pco2_h_l", "tank_id")
) %>%
  mutate()



# ============================================================
# 1. GROWTH AND PERCENT-CHANGE CALCULATIONS
# ============================================================
bio_growth <- bio_joined %>%
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
    
    # --- Absolute growth rates (time-standardized) ---
    diam_growth_mm_day   = delta_diam_mm / delta_days,
    ht_growth_mm_day     = delta_ht_mm   / delta_days,
    volume_growth_ml_day = delta_volume_ml / delta_days,
    
    diam_growth_mm_week   = diam_growth_mm_day * 7,
    ht_growth_mm_week     = ht_growth_mm_day * 7,
    volume_growth_ml_week = volume_growth_ml_day * 7,
    
    diam_growth_mm_month   = diam_growth_mm_day * 30.44,
    ht_growth_mm_month     = ht_growth_mm_day * 30.44,
    volume_growth_ml_month = volume_growth_ml_day * 30.44,
    
    # --- Percent changes (relative to initial) ---
    perc_change_diam   = (delta_diam_mm   / initial_diam_mm)   * 100,
    perc_change_ht     = (delta_ht_mm     / initial_ht_mm)     * 100,
    perc_change_volume = (delta_volume_ml / initial_volume_ml) * 100,
    
    # --- Percent change rates (time-standardized) ---
    perc_change_diam_day   = perc_change_diam / delta_days,
    perc_change_ht_day     = perc_change_ht   / delta_days,
    perc_change_volume_day = perc_change_volume / delta_days,
    
    perc_change_diam_week   = perc_change_diam_day * 7,
    perc_change_ht_week     = perc_change_ht_day * 7,
    perc_change_volume_week = perc_change_volume_day * 7,
    
    perc_change_diam_month   = perc_change_diam_day * 30.44,
    perc_change_ht_month     = perc_change_ht_day * 30.44,
    perc_change_volume_month = perc_change_volume_day * 30.44
  )

# ============================================================
# 2. SUMMARY (MEAN ± SD ± n) BY treatment × pCO₂ × tank
# ============================================================
bio_summary <- bio_growth %>%
  group_by(temp_deg_c, pco2_h_l, tank_id) %>%
  summarise(
    n = n(),
    across(
      .cols = c(
        delta_diam_mm, delta_ht_mm, delta_volume_ml,
        diam_growth_mm_day, ht_growth_mm_day, volume_growth_ml_day,
        diam_growth_mm_week, ht_growth_mm_week, volume_growth_ml_week,
        diam_growth_mm_month, ht_growth_mm_month, volume_growth_ml_month,
        perc_change_diam, perc_change_ht, perc_change_volume,
        perc_change_diam_day, perc_change_ht_day, perc_change_volume_day,
        perc_change_diam_week, perc_change_ht_week, perc_change_volume_week,
        perc_change_diam_month, perc_change_ht_month, perc_change_volume_month
      ),
      .fns = list(mean = ~mean(.x, na.rm = TRUE),
                  sd   = ~sd(.x,   na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# ============================================================
# 3. ANOVA-READY LONG FORMAT SUMMARY
# ============================================================
bio_summary_long <- bio_summary %>%
  pivot_longer(
    cols = -c(temp_deg_c, pco2_h_l, tank_id, n),
    names_to  = c("metric", ".value"),
    names_pattern = "(.*)_(mean|sd)$"
  )

# The resulting columns:
#   temp_deg_c | pco2_h_l | tank_id | n | metric | mean | sd
# Each row = one growth/percent metric for one tank replicate.

# ============================================================
# 4. EXPORT OUTPUTS
# ============================================================
write.csv(bio_growth,        "./Data/Data_Final/Biometrics_Growth_PercentChange.csv", row.names = F)

# ============================================================
#  Statistical Analysis Module — Mixed ANOVA with glmmTMB
#  Author: N. B. Spindel
#  Date:   2025-10-13
# ============================================================

# ---- Load packages ----
libs <- c("tidyverse", "glmmTMB", "performance", "emmeans", "DHARMa", "car", "Hmisc")
lapply(libs, require, character.only = TRUE)

# ---- Read processed dataset ----
bio_growth <- read_csv("./Data/Data_Final/Biometrics_Growth_PercentChange.csv")

# ---- Ensure correct data structure ----
bio_growth <- bio_growth %>%
  mutate(
    temp_deg_c = factor(temp_deg_c),
    pco2_h_l   = factor(pco2_h_l),
    tank_id    = factor(tank_id)
  ) %>% filter(temp_deg_c %nin% c("21-18", "18-14"))

# ============================================================
# 1. MODEL FORMULATION
# ============================================================
# Example response variables (choose one or iterate over several):
#   delta_diam_mm, delta_ht_mm, delta_volume_ml,
#   diam_growth_mm_day, volume_growth_ml_day,
#   perc_change_volume, etc.

# Example: model absolute change in spheroid volume (Δvolume, mL)
m1 <- glmmTMB(
  perc_change_diam ~ log(initial_diam_mm) + 
    temp_deg_c * pco2_h_l + 
    (1 | tank_id),
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
performance::r2_nakagawa(m1, m1_disp)

# ---- Model summary ----
summary(m1)

# ============================================================
# 2. MODEL DIAGNOSTICS
# ============================================================
# --- Residual diagnostics using DHARMa ---
res <- simulateResiduals(m1)
plot(res)                      # Residual diagnostics plot
testDispersion(res)            # Check for over/under dispersion
testUniformity(res)            # Check uniformity of residuals
testZeroInflation(res)         # Check if zero-inflation is needed

# ============================================================
# 3. ANOVA TABLE (TYPE II OR III)
# ============================================================
# Type III ANOVA (for interaction models)
anova_m1 <- car::Anova(m1, type = 3)
print(anova_m1)
# ---- Extract residuals and fitted values ----
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


# ============================================================
# 4. POST-HOC COMPARISONS
# ============================================================
# Estimated marginal means (EMMs) and pairwise contrasts
emm_treat <- emmeans(m1, ~ temp_deg_c * pco2_h_l)
pairs(emm_treat, simple = "each", combine = TRUE, adjust = "tukey")

# Summarize marginal means
summary(emm_treat)

# Plot interaction of treatment × pCO₂
emmip(m1, pco2_h_l ~ temp_deg_c,
      CIs = TRUE, type = "response") +
  theme_bw() +
  labs(
    x = "Temperature Treatment (°C)",
    y = "Predicted Δ Volume (mL)",
    color = "pCO₂ Treatment"
  )

# ============================================================
# 5. LOOP OVER MULTIPLE GROWTH METRICS (OPTIONAL)
# ============================================================
# Define the variables of interest
growth_vars <- c(
  "delta_diam_mm", "delta_ht_mm", "delta_volume_ml",
  "diam_growth_mm_day", "ht_growth_mm_day", "volume_growth_ml_day",
  "perc_change_diam", "perc_change_ht", "perc_change_volume"
)

# Function to fit glmmTMB model and extract ANOVA summary
fit_growth_model <- function(var) {
  f <- as.formula(paste(var, "~ temp_deg_c * pco2_h_l + (1 | tank_id)"))
  model <- glmmTMB(f, data = bio_growth, family = t_family(link = "identity"))
  list(
    variable = var,
    anova = car::Anova(model, type = 3),
    summary = summary(model)
  )
}

# Apply to all variables
model_results <- lapply(growth_vars, fit_growth_model)

# Extract Type III ANOVA tables into tidy summary
anova_summary <- map_dfr(model_results, function(x) {
  as.data.frame(x$anova) %>%
    mutate(
      term = rownames(.),
      response = x$variable
    ) %>%
    select(response, term, everything())
})

# Write ANOVA results
write_csv(anova_summary, "Biometrics_Growth_ANOVA_Results.csv")

# ============================================================
# 6. MODEL COMPARISON OR VISUALIZATION (OPTIONAL)
# ============================================================
# Compare models with or without interaction
m0 <- glmmTMB(delta_volume_ml ~ temp_deg_c + pco2_h_l + (1 | tank_id),
              data = bio_growth, family = t_family(link = "identity"))
AIC(m0, m1)

# Evaluate model fit (R² conditional and marginal)
performance::r2_nakagawa(m1)
