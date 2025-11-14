###############################################################################
# Mapping script for environmental context                                    #
#                                                                             #
###############################################################################

# Load required libraries --------
library(marmap)   # NOAA ETOPO1 and GEBCO data downloader
library(terra)    # raster handling
library(raster)   # for compatibility
library(sf)
library(Hmisc)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(lubridate)
library(janitor)
library(here)
library(forcats)

source("./Utility_Functions.R")

# Generate map of S. purpuratus range overlaid with relevant temperature and pCO2 data --------
# 1. Define geographic extent
# Covers: Southern Alaska (~65°N) to Ensenada (~28°N)
lon_min <- -155
lon_max <- -115
lat_min <- 28
lat_max <- 65

# 2. Download ETOPO1 bathymetry (1 arc-minute resolution)
bathy <- getNOAA.bathy(lon1 = lon_min, 
                       lon2 = lon_max, 
                       lat1 = lat_min, 
                       lat2 = lat_max, 
                       resolution = 1)

# 3. Convert to raster
bathy_raster <- marmap::as.raster(bathy)

# 4. Convert to terra::SpatRaster
bathy_terra <- terra::rast(bathy_raster)

# 5. (Optional) Save as GeoTIFF for reuse
terra::writeRaster(bathy_terra, filename = "./Data/composite_bath_alaska_to_mexico.tif", overwrite = TRUE)

# 6. Reload to confirm
bathy_loaded <- terra::rast("./Data/composite_bath_alaska_to_mexico.tif")

# 7. Plot quicklook
plot(bathy_loaded, main = "NOAA ETOPO1 Bathymetry: Alaska to Mexico")

# Load bathymetry raster from GeoTIFF 
bathy_raster <- rast("./Data/composite_bath_alaska_to_mexico.tif")

roi_extent <- ext(lon_min, lon_max, lat_min, lat_max)

bathy_clipped <- crop(bathy_raster, roi_extent)

# Generate mask for depth between 0 and 100 m (i.e., values between -100 and 0)
depth_mask <- classify(bathy_clipped, 
                       rcl = matrix(c(-Inf, -100.01, NA, -100, 0, 1, 0, Inf, NA), 
                                    ncol = 3, byrow = TRUE))

# Convert raster to polygons
depth_sf <- as.polygons(depth_mask, dissolve = TRUE, na.rm = TRUE)
depth_sf <- st_as_sf(depth_sf)
depth_sf <- st_make_valid(depth_sf)  # ensure geometry validity

# Load land layer for base map
land <- ne_countries(scale = "medium", returnclass = "sf")

# Define inset bounding boxes (WGS84)
inset_boxes <- list(
  bamfield  = st_bbox(c(xmin = -145, xmax = -135,   
                        ymin = 45.0, ymax = 50.0), 
                      crs = 4326),
  newport = st_bbox(c(xmin = -145, xmax = -135, 
                      ymin = 39, ymax = 44), 
                    crs = 4326)
) %>% purrr::map(st_as_sfc)

# Define site coordinates
site_data <- data.frame(
  site = c("Bamfield", "Newport", "Quadra"),
  lat = c(48.5, 44.0, 50.116),
  lon = c(-125.25, -124, -125.222)
)

# Convert to sf
site_sf <- st_as_sf(site_data, coords = c("lon", "lat"), crs = 4326)

# Global mapping from variable suffix to readable labels
canonical_var_names <- c(
  "AlkalinityTot" = "Alkalinity (µmol/kg)",
  "OmegaAragSat"  = "Ω Aragonite",
  "pH"            = "pH (total scale)",
  "Salinity"      = "Salinity (PSU)",
  "TCO2"          = "TCO₂ (µmol/kg)",
  "WaterTemp"     = "Temperature (°F)",
  "CO2"           = "pCO₂ (µatm)"
)

# Get filepaths for time series
all_files <- load_site_files("Newport_OR")

# Apply loader function
ocean_data_newport <- map_dfr(all_files, load_ocean_file) 

ocean_pco2_newport <- ocean_data_newport %>% 
  filter(variable == "pCO₂ (µatm)")
ocean_temp_newport <- ocean_data_newport %>% 
  filter(variable == "Temperature (°F)", depth_ft %in% c("-3 ft", "-23 ft")) %>%
  mutate(temp_deg_c = f_to_c(value))

merged_df <- full_join(ocean_temp_newport, 
                       ocean_pco2_newport, 
                       by = c("datetime", "depth_ft"))

scale_factor_newport <- max(merged_df$value.y, na.rm = TRUE) /
  max(merged_df$temp_deg_c, na.rm = TRUE)

# Compute y-axis (pCO2) limits that accommodate the full temperature range
temp_breaks <- c(0,2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22)
pco2_breaks <- c(0,200,400,600,800,1000,1200, 1400, 1600, 1800)


# Clean and parse datetime
df <- merged_df %>%
  mutate(datetime = ymd_hms(datetime),
         date = as.Date(datetime),
         month = as.Date(floor_date(datetime, "month")),
         depth_ft = as.numeric(str_remove(depth_ft, " ft"))) %>%
  dplyr::rename(pCO2 = value.y)

# Define time span for full water year (April to March)
start_date <- as.Date(paste0(year(min(df$date)), "-04-01"))
end_date <- start_date + years(1) - days(1)

# Create full daily and monthly date sequences
all_dates <- tibble(date = seq.Date(start_date, end_date, by = "day"))
all_months <- tibble(month = seq.Date(start_date, end_date, by = "month"))

# Summarize daily
daily_summary <- df %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(mean_pCO2 = mean(pCO2, na.rm = TRUE), .groups = "drop") %>%
  right_join(all_dates, by = "date") %>%
  dplyr::mutate(source = "Daily", time = date)

# Summarize monthly
monthly_summary <- df %>%
  dplyr::group_by(month) %>%
  dplyr::summarise(mean_pCO2 = mean(pCO2, na.rm = TRUE), .groups = "drop") %>%
  right_join(all_months, by = "month") %>%
  dplyr::mutate(source = "Monthly", time = month)

# Combine
plot_data <- bind_rows(daily_summary, monthly_summary)

# Define month order
month_order <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep",
                 "Oct", "Nov", "Dec", "Jan", "Feb", "Mar")

# Plot
ggplot(plot_data, aes(x = time, y = mean_pCO2)) +
  geom_point(data = filter(plot_data, source == "Daily"),
             shape = 1, size = 2) +  # no jitter
  geom_point(data = filter(plot_data, source == "Monthly"),
             shape = 16, size = 3, color = "black") +
  scale_x_date(
    breaks = seq.Date(from = start_date, to = end_date, by = "1 month"),
    labels = function(x) {
      m <- month(x)
      factor(month.abb[m], levels = month_order)
    },
    expand = c(0.01, 0.01),
    limits = c(start_date, end_date)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=6), limits = c(0,1500)) +
  geom_hline(yintercept = 600, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 1200, linetype = "dashed", color = "red") +
  labs(title = NULL,
       x = "Month", y = expression(pCO[2]~(μatm))) +
  # theme_bw(base_size = 9) +
  gg_options() +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Get filepaths for time series
all_files_bamfield <- load_site_files(region = "Bamfield")

# Apply the loader function
ocean_data_bamfield <- map_dfr(all_files_bamfield, load_ocean_file) 
# write.csv(ocean_data_bamfield, "./Data/bamfield_oceanography.csv", row.names = F)

ocean_pco2_bamfield <- ocean_data_bamfield %>% filter(variable == "pCO₂ (µatm)")

# Clean and parse datetime
bamfield_df <- ocean_pco2_bamfield %>%
  mutate(datetime = ymd_hms(datetime),
         date = as.Date(datetime),
         month = as.Date(floor_date(datetime, "month")),
         depth_ft = as.numeric(str_remove(depth_ft, " ft"))) %>%
  dplyr::rename(pCO2 = value)

# Summarize daily
daily_summary_bam <- bamfield_df %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(mean_pCO2 = mean(pCO2, na.rm = TRUE), .groups = "drop") %>%
  right_join(all_dates, by = "date") %>%
  dplyr::mutate(source = "Daily", time = date)

# Summarize monthly
monthly_summary_bam <- bamfield_df %>%
  dplyr::group_by(month) %>%
  dplyr::summarise(mean_pCO2 = mean(pCO2, na.rm = TRUE), .groups = "drop") %>%
  right_join(all_months, by = "month") %>%
  dplyr::mutate(source = "Monthly", time = month)

# Combine
plot_data_bam <- bind_rows(daily_summary_bam, monthly_summary_bam)

# Plot
ggplot(plot_data_bam, aes(x = time, y = mean_pCO2)) +
  geom_point(data = filter(plot_data_bam, source == "Daily"),
             shape = 1, size = 2) +  # no jitter
  geom_point(data = filter(plot_data_bam, source == "Monthly"),
             shape = 16, size = 3, color = "black") +
  scale_x_date(
    breaks = seq.Date(from = start_date, to = end_date, by = "1 month"),
    labels = function(x) {
      m <- month(x)
      factor(month.abb[m], levels = month_order)
    },
    expand = c(0.01, 0.01),
    limits = c(start_date, end_date)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=6), limits = c(0,1500)) +
  
  geom_hline(yintercept = 600, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 1200, linetype = "dashed", color = "red") +
  labs(title = NULL,
       x = "Month", y = expression(pCO[2]~(μatm))) +
  # theme_bw(base_size = 9) +
  gg_options() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



"purple"
"darkviolet"
"darkmagenta"
"mediumpurple4"
# Step 1: List and read all files, skip metadata
file_paths_quadra <- list.files(
  path = "./Data/Sensor_Data/QuadraBoL/",
  pattern = "^189920.*\\.csv$",
  full.names = TRUE
)

# Step 2: Define expected column aliases to handle variation
standard_colnames <- c("datetime", "temp", "pco2")

all_data_quadra <- file_paths_quadra %>%
  map_dfr(~ {
    df <- read_csv(.x, skip = 6, show_col_types = FALSE,
                   col_types = cols(.default = col_character()))  # Force all as character for now
    
    # Try flexible renaming
    names(df) <- str_replace_all(names(df), fixed(" "), "_")
    
    # Standardize known variable names
    names(df)[str_detect(names(df), regex("date.*gmt", ignore_case = TRUE))] <- "datetime"
    names(df)[str_detect(names(df), regex("TSG.*T.*deg", ignore_case = TRUE))] <- "temp"
    names(df)[str_detect(names(df), regex("SWpCO2.*TSG_T", ignore_case = TRUE))] <- "pco2"
    
    # If any required column is missing, warn and skip
    if (!all(standard_colnames %in% names(df))) {
      warning(paste("Skipping file:", .x, "- missing columns"))
      return(tibble())
    }
    
    df %>%
      dplyr::select(all_of(standard_colnames))
  })

# Step 3: Convert types and filter
df_quadra <- all_data_quadra %>%
  mutate(
    date = as.Date(datetime),
    day = day(date),
    month = month(date, label = T, abbr = T),
    year = year(date),
    temp = as.numeric(temp),
    pco2 = as.numeric(pco2)
  ) %>%
  drop_na(datetime, temp, pco2)

# Summarize daily
daily_summary_quadra <- df_quadra %>%
  dplyr::group_by(day, month, year) %>%
  dplyr::summarise(mean_pCO2 = mean(pco2, na.rm = TRUE), .groups = "drop") 

# Summarize monthly
monthly_summary_quadra <- df_quadra %>%
  dplyr::group_by(month, year) %>%
  dplyr::summarise(mean_pCO2 = mean(pco2, na.rm = TRUE), .groups = "drop")

# Combine
plot_data_quadra <- bind_rows(daily_summary_quadra, monthly_summary_quadra)
# Step 1: Ensure month is ordered April–March
month_levels <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep",
                  "Oct", "Nov", "Dec", "Jan", "Feb", "Mar")

# Assume your dataframe has these columns: day (int), month (factor or string), year (numeric), value (numeric)

# Step 1: Clean input and compute actual dates
df <- df_quadra %>%
  mutate(
    month = factor(month, levels = month.abb, ordered = TRUE),
    date = make_date(year, match(month, month.abb), day),
    
    # Determine "season year" starting in April
    season_year = if_else(month(date) >= 4, year, year - 1),
    
    # Reference date is April 1 of season year
    ref_date = as.Date(paste0(season_year, "-04-01")),
    
    # Compute day-of-season as days since April 1
    doy_season = as.numeric(difftime(date, ref_date, units = "days")) + 1
  )

# Step 2: Daily averages
daily_df <- df %>%
  dplyr::group_by(season_year, doy_season) %>%
  dplyr::summarize(value = mean(pco2, na.rm = TRUE), .groups = "drop")

# Step 3: Monthly averages (use 15th of each month as anchor)
monthly_df <- df %>%
  dplyr::mutate(month_day = 15) %>%
  dplyr::mutate(monthly_date = make_date(year, match(month, month.abb), month_day),
         season_year = if_else(month(monthly_date) >= 4, year, year - 1),
         ref_date = as.Date(paste0(season_year, "-04-01")),
         doy_season = as.numeric(difftime(monthly_date, ref_date, units = "days")) + 1) %>%
  dplyr::group_by(season_year, doy_season) %>%
  dplyr::summarize(value = mean(pco2, na.rm = TRUE), .groups = "drop")

# Step 4: Build month labels and breaks (April 1 to March 1)
month_starts <- seq.Date(as.Date("2000-04-01"), as.Date("2001-03-01"), by = "1 month")
month_labels <- format(month_starts, "%b")
month_breaks <- as.numeric(difftime(month_starts, as.Date("2000-04-01"), units = "days")) + 1

# Step 5: Plot
ggplot() +
  geom_point(data = daily_df, aes(x = doy_season, y = value), shape = 16, color = "black", size = 1, alpha = 0.5) +
  geom_point(data = monthly_df, aes(x = doy_season, y = value), shape = 21, color = "black", size = 2, fill = "black", alpha = 0.8) +
  geom_hline(yintercept = 600, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 1200, linetype = "dashed", color = "red") +
  
  scale_x_continuous(
    name = NULL,
    breaks = month_breaks,
    labels = month_labels,
    expand = c(0.01, 0.01),
    limits = c(1, 366)  # optional: March 31
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=6), limits = c(0,1500)) +
  labs(y = expression("pCO"[2]*" (µatm)")) +
  gg_options() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

df_quadra_filtered <- df_quadra # %>% 
  # filter(datetime > "2022-01-01")

pco2_range_quadra <- range(df_quadra_filtered$pco2, na.rm = TRUE)
temp_range_quadra <- range(df_quadra_filtered$temp, na.rm = TRUE)
scale_factor_quadra <- diff(pco2_range_quadra) / diff(temp_range_quadra)

# Step 4: Plot
p_quadra_pco2_temp <- ggplot(df_quadra_filtered, aes(x = datetime)) +
  geom_line(aes(y = pco2), color = "black") +
  geom_line(aes(y = temp * scale_factor_quadra), color = "red") +
  scale_y_continuous(
    name = expression("pCO"[2]*" (µatm)"),
    breaks = pco2_breaks,
    sec.axis = sec_axis(~ . / scale_factor_quadra, name = "Temperature (°C)", breaks = temp_breaks)
  ) +
  labs(x = NULL, title = NULL) +
  theme_bw(base_size = 9) +
  theme(
    axis.title.y.left = element_text(color = "black"),
    axis.title.y.right = element_text(color = "red"),
    axis.text.y.left = element_text(color = "black"),
    axis.text.y.right = element_text(color = "red"),
    axis.line.y.left = element_line(color = "black"),
    axis.line.y.right = element_line(color = "red"),
    axis.line.x = element_line(color = "black"),
    # Transparent background settings
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

p_quadra_pco2_temp

# 7. Plot base map + urchin range polygon
p_base_map <- ggplot() +
  geom_sf(data = depth_sf, fill = "mediumpurple", alpha = 0.9, color = NA) +
  geom_sf(data = land, fill = "grey85", color = "white", linewidth = 0.3) +
  
  geom_sf(data = site_sf, size = 3, color = "red") +
  geom_sf_text(data = site_sf, 
               aes(label = site), 
               nudge_x = 3, size = 4, fontface = "bold") +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  labs(
    title = NULL,
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 13)

p_base_map


# Combine all components into final figure
final_map <- cowplot::ggdraw() +
  cowplot::draw_plot(p_base_map) +
  cowplot::draw_plot(p_quadra_pco2_temp, 
                     x = 0.56, y = 0.7, 
                     width = 0.42, height = 0.22) +
  geom_segment(
    aes(x = 0.75, xend = 0.75, 
        y = 0.7, yend = 0.61),  # coordinates in plot space
    arrow = arrow(length = unit(0.2, "cm")),
    lineend = "round",
    color = "black",
    linewidth = 1.2
  ) +
  
  cowplot::draw_plot(p_bamfield_pco2_temp, 
                     x = 0.14, y = 0.4, 
                     width = 0.42, height = 0.22) +
  geom_segment(
    aes(x = 0.57, xend = 0.73, 
        y = 0.53, yend = 0.55),  # coordinates in plot space
    arrow = arrow(length = unit(0.2, "cm")),
    lineend = "round",
    color = "black",
    linewidth = 1.2
  ) +
  cowplot::draw_plot(p_newport_pco2_temp, 
                     x = 0.19, y = 0.15, 
                     width = 0.42, height = 0.22) +
  geom_segment(
    aes(x = 0.62, xend = 0.77, 
        y = 0.3, yend = 0.45),  # coordinates in plot space
    arrow = arrow(length = unit(0.2, "cm")),
    lineend = "round",
    color = "black",
    linewidth = 1.2
  )
  
final_map

ggsave("./Figures/Final/Map.pdf", final_map,
       units = "in", width = 8, height = 9)
ggsave("./Figures/Final/Map.png", final_map,
       units = "in", width = 8, height = 9)
ggsave("./Figures/Final/Map.jpg", final_map,
       units = "in", width = 8, height = 9)

ggsave("./Figures/Final/Base_Map.pdf", p_base_map,
       units = "in", width = 9, height = 10)
ggsave("./Figures/Final/Base_Map.png", p_base_map,
       units = "in", width = 9, height = 10)
ggsave("./Figures/Final/Base_Map.jpg", p_base_map,
       units = "in", width = 9, height = 10)


ggsave("./Figures/Final/Quadra_TS.pdf", p_quadra_pco2_temp,
       units = "in", width = 4, height = 2)
ggsave("./Figures/Final/Quadra_TS.png", p_quadra_pco2_temp,
       units = "in", width = 4, height = 2)
ggsave("./Figures/Final/Quadra_TS.jpg", p_quadra_pco2_temp,
       units = "in", width = 4, height = 2)

ggsave("./Figures/Final/Bamfield_TS.pdf", p_bamfield_pco2_temp,
       units = "in", width = 4, height = 2)
ggsave("./Figures/Final/Bamfield_TS.png", p_bamfield_pco2_temp,
       units = "in", width = 4, height = 2)
ggsave("./Figures/Final/Bamfield_TS.jpg", p_bamfield_pco2_temp,
       units = "in", width = 4, height = 2)

ggsave("./Figures/Final/Newport_TS.pdf", p_newport_pco2_temp,
       units = "in", width = 4, height = 2)
ggsave("./Figures/Final/Newport_TS.png", p_newport_pco2_temp,
       units = "in", width = 4, height = 2)
ggsave("./Figures/Final/Newport_TS.jpg", p_newport_pco2_temp,
       units = "in", width = 4, height = 2)

p_ts_multi <- cowplot::plot_grid(
  p_quadra_pco2_temp,
  p_bamfield_pco2_temp,
  p_newport_pco2_temp,
  align = "tblr",
  axis = "hv",
  ncol = 1
)
p_ts_multi

ggsave("./Figures/Final/TS_Multi.pdf", p_ts_multi,
       units = "in", width = 4, height = 6)
ggsave("./Figures/Final/TS_Multi.png", p_ts_multi,
       units = "in", width = 4, height = 6)
ggsave("./Figures/Final/TS_Multi.jpg", p_ts_multi,
       units = "in", width = 4, height = 6)

final_map <- cowplot::ggdraw() +
  cowplot::draw_plot(p_base_map) +
  cowplot::draw_plot(p_ts_multi, 
                     x = 0.16, y = 0.05, 
                     width = 0.4, height = 0.7)
final_map

ggsave("./Figures/Final/Map.pdf", final_map,
       units = "in", width = 7.2, height = 9)
ggsave("./Figures/Final/Map.png", final_map,
       units = "in", width = 7.2, height = 9, dpi = 600)
ggsave("./Figures/Final/Map.jpg", final_map,
       units = "in", width = 7.2, height = 9)

# Get metabolism plots -------------------
selected_model <- fit_metabolism_const_temp_oa <- readRDS("./Model_Output/metabolism_temp_oa_full.rds")
selected_model$formula

ce_95 <- conditional_effects(selected_model, robust = F)
ce_80 <- conditional_effects(selected_model, robust = F, prob = 0.8)

data <- data_metabolism_const_temp_oa %>% filter(afdm_g > 0.1) # exclude erroneous measurement

df_ce_temp_pco2_plot_95 <- ce_95$`temp_deg_c:pco2_uatm` %>%
  dplyr::mutate(
    plot_label = factor(pco2_uatm, labels = c("600 μatm", "1200 μatm"))
  )

df_ce_temp_pco2_plot_80 <- ce_80$`temp_deg_c:pco2_uatm` %>%
  dplyr::mutate(
    plot_label = factor(pco2_uatm, labels = c("600 μatm", "1200 μatm"))
  )

pd <- position_dodge(width = 0.2)

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

empirical_means_rmr <- data %>%
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

fit_metabolism_rmr_dyn_wet_mass <- readRDS("./Model_Output/metabolism_dyn_wet_mass_full.rds")
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
    # y = expression("respiration rate (mgO"["2"]*"h"^-1*")"),
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

p_respiration_multi <- cowplot::plot_grid(
  egg::tag_facet(p_temp_pco2_rmr, tag_pool = c("a")), 
  egg::tag_facet(p_dyn_wet_mass_rmr, tag_pool = c("b"))
)

p_respiration_multi

ggsave("./Figures/Final/Respiration_Multi.png", p_respiration_multi,
       units = "in", width = 12, height = 6)
ggsave("./Figures/Final/Respiration_Multi.pdf", p_respiration_multi,
       units = "in", width = 12, height = 6)
ggsave("./Figures/Final/Respiration_Multi.jpg", p_respiration_multi,
       units = "in", width = 12, height = 6)

# DEVEL ------
library(tidyverse)
library(lubridate)

# ---- Step 1: Load and parse data ----
# Update to your folder path
data_dir <- "./Data/pCO2_ph_ts_NCEI_0173932/1.1/data/0-data/Coastal_Shelf/"
files <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)

read_site_data <- function(file) {
  site <- tools::file_path_sans_ext(basename(file))
  read_tsv(file, comment = "#", na = "NaN", show_col_types = FALSE) %>%
    dplyr::select(datetime_utc, pCO2_sw) %>%
    mutate(
      site = site,
      datetime = ymd_hms(datetime_utc),
      date = as.Date(datetime)
    )
}

all_data <- map_dfr(files, read_site_data)

# ---- Step 2: Compute daily averages ----
daily_df <- all_data %>%
  dplyr::group_by(site, date) %>%
  dplyr::summarize(pCO2_sw = mean(pCO2_sw, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    year = year(date),
    month = month(date),
    day = day(date),
    season_year = if_else(month >= 4, year, year - 1),
    season_ref = as.Date(paste0(season_year, "-04-01")),
    doy_season = as.numeric(date - season_ref) + 1
  )

# ---- Step 3: Compute monthly averages ----
monthly_df <- all_data %>%
  dplyr::mutate(
    year = year(date),
    month = month(date),
    season_year = if_else(month >= 4, year, year - 1),
    mid_month = as.Date(paste0(year, "-", month, "-15")),
    season_ref = as.Date(paste0(season_year, "-04-01")),
    doy_season = as.numeric(mid_month - season_ref) + 1
  ) %>%
  dplyr::group_by(site, season_year, doy_season) %>%
  dplyr::summarize(pCO2_sw = mean(pCO2_sw, na.rm = TRUE), .groups = "drop")

# ---- Step 4: Rank sites by mean pCO2_sw ----
site_order <- all_data %>%
  dplyr::group_by(site) %>%
  dplyr::summarize(mean_pCO2 = mean(pCO2_sw, na.rm = TRUE)) %>%
  arrange(desc(mean_pCO2)) %>%
  pull(site)

daily_df$site <- factor(daily_df$site, levels = site_order)
monthly_df$site <- factor(monthly_df$site, levels = site_order)

# ---- Step 5: Define month labels for x-axis ----
month_starts <- seq.Date(as.Date("2000-04-01"), as.Date("2001-03-01"), by = "1 month")
month_labels <- format(month_starts, "%b")
month_breaks <- as.numeric(month_starts - as.Date("2000-04-01")) + 1

# ---- Step 6: Plot ----
ggplot() +
  geom_point(data = daily_df, aes(x = doy_season, y = pCO2_sw), shape = 1, size = 0.6, alpha = 0.5) +
  geom_point(data = monthly_df, aes(x = doy_season, y = pCO2_sw), shape = 16, size = 1.5, color = "blue") +
  scale_x_continuous(
    name = NULL,
    breaks = month_breaks,
    labels = month_labels,
    limits = c(1, 366),
    expand = c(0.01, 0.01)
  ) +
  labs(y = expression("Seawater " * pCO[2] * " (µatm)")) +
  facet_wrap(~ site, nrow = 1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=6)) +
  gg_options() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )
