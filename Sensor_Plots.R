###############################################################################
# Code to visualize data scraped from the Hakai Sensor Network associated     # 
# with manuscript "Thermal suppression of gametogenesis explains historical   #
# collapses in larval recruitment"                                            #
#                                                                             #
# Author: Nathan B. Spindel                                                   #
###############################################################################

# Load CRAN libraries ----
library(readxl)
library(tidyverse)
library(reshape2)
library(scales)
library(viridis)
library(zoo)

gg_options <-  function(){
  theme_bw(base_size = 16)+theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.background =  element_blank(),
    legend.background =  element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank())}

# Visualize Marna Lab and Limpet Temp Time Series Data ------------------------------
dat_raw_quadra_limpet_seafet <- 
  read_excel("./Data/Sensor_Data/Quadra_Island_Limpet_SeaFET.xlsx",
             skip = 3) 

dat_processed_quadra_limpet_seafet <- dat_raw_quadra_limpet_seafet %>%
  filter(
    !is.na(TSG_T_Avg)
  ) %>%
  mutate(date = as.Date(date))

temp_quadra_ts_plot_day <- ggplot(data = dat_processed_quadra_limpet_seafet, 
                                  aes(x = date, y = TSG_T_Avg))+
  scale_x_date(date_breaks = "5 day", date_labels = "%D")+
  geom_point(size = 1)+
  geom_line()+
  geom_smooth()+
  scale_y_continuous( 
    breaks = pretty_breaks(n = 6),
    limits = c(12,18))+
  labs(
    y = "temperature (daily mean deg C)",
    x = NULL,
    title = "Quadra Island Limpet SeaFET (depth ~ 20 m)"
  )+
  gg_options()+
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.75)
  )+
  guides(
    color = guide_legend(override.aes = list(size=5))
  );temp_quadra_ts_plot_day

ggsave("./Figures/Quadra_Limpet_SeaFET_Temperature_Time_Series.pdf", temp_quadra_ts_plot_day,
       width = 10, height = 6, units = "in")
ggsave("./Figures/Quadra_Limpet_SeaFET_Temperature_Time_Series.png", temp_quadra_ts_plot_day,
       width = 10, height = 6, units = "in")
ggsave("./Figures/Quadra_Limpet_SeaFET_Temperature_Time_Series.jpg", temp_quadra_ts_plot_day,
       width = 10, height = 6, units = "in")

dat_raw_marna <- 
  read_excel("./Data/Sensor_Data/Marna_Lab_Mesocosm_Temps.xlsx",
             skip = 3) 

# reshape data frame to long format
columns <- names(dat_raw_marna %>% 
                   select(!c("Walchem_Temperature_Tank7":"Walchem_Temperature_Tank9_Max")))

dat_raw_marna_melt <- reshape2::melt(dat_raw_marna, id.vars=columns,
                                     variable.name="sensor_category", value.name="temp")

# get metadata for each sensor
source("./Data/Sensor_Data/Sensor_Metadata.R")

dat_processed_marna_melt <- dat_raw_marna_melt %>%
  filter(!is.na(temp)) %>%
  mutate(
    date = as.Date(date),
    constant_dynamic = as.factor(ifelse(sensor_category %in% dynamic_tank_sensor_categories, 
                                        "dynamic", "constant"))
  ) %>%
  arrange(date) %>%
  mutate(trial_day = as.integer(factor(date, levels = sort(unique(date))))) %>%
  filter(grepl("Avg", sensor_category)) %>%
  mutate(
    treatment = factor(case_when(
      sensor_category %in% la_nina_tank_categories ~ "La Niña (18-14 °C)",
      sensor_category %in% el_nino_tank_categories ~ "El Niño (21-18 °C)",
      sensor_category %in% const_10_tanks ~ "10 °C",
      sensor_category %in% const_13_tanks ~ "13 °C (ambient)",
      sensor_category %in% const_16_tanks ~ "16 °C",
      sensor_category %in% const_17_tanks ~ "17 °C",
      sensor_category %in% const_18_tanks ~ "18 °C",
      sensor_category %in% const_20_tanks ~ "20 °C"
    ), levels = c("10 °C","13 °C (ambient)","16 °C","La Niña (18-14 °C)","17 °C","18 °C","20 °C","El Niño (21-18 °C)"))) %>%
  dplyr::filter(!is.na(treatment)) %>%
  dplyr::mutate(
    treatment = factor(treatment, labels = c("10 °C","13 °C (ambient)","16 °C","18-14 °C (La Niña)","17 °C","18 °C","20 °C","21-18 °C (El Niño)"))
  ) %>%
  dplyr::group_by(date, trial_day, treatment) %>%
  dplyr::summarise(
    mean_temp_c = mean(temp),
    sd_temp_c = sd(temp),
    n = n()
  ) %>%
  mutate(
    se_temp_c = sd_temp_c/sqrt(n)
  )

dat_processed_marna_melt_reps <- dat_raw_marna_melt %>%
  dplyr::filter(!is.na(temp)) %>%
  dplyr::mutate(
    date = as.Date(date),
    constant_dynamic = as.factor(ifelse(sensor_category %in% dynamic_tank_sensor_categories, 
                                        "dynamic", "constant"))
  ) %>%
  dplyr::filter(grepl("Avg", sensor_category)) %>%
  dplyr::mutate(
    treatment = factor(dplyr::case_when(
      sensor_category %in% la_nina_tank_categories ~ "La Niña (18-14 °C)",
      sensor_category %in% el_nino_tank_categories ~ "El Niño (21-18 °C)",
      sensor_category %in% const_10_tanks ~ "10 °C",
      sensor_category %in% const_13_tanks ~ "13 °C (ambient)",
      sensor_category %in% const_16_tanks ~ "16 °C",
      sensor_category %in% const_17_tanks ~ "17 °C",
      sensor_category %in% const_18_tanks ~ "18 °C",
      sensor_category %in% const_20_tanks ~ "20 °C"
    ), levels = c("10 °C","13 °C (ambient)","16 °C","La Niña (18-14 °C)","17 °C","18 °C","20 °C","El Niño (21-18 °C)"))) %>%
  dplyr::filter(!is.na(treatment)) %>%
  dplyr::mutate(
    treatment = factor(treatment, labels = c("10 °C","13 °C (ambient)","16 °C","18-14 °C (La Niña)","17 °C","18 °C","20 °C","21-18 °C (El Niño)"))
  )

# initial_temps <- dat_processed_marna_melt %>%
#   filter(date == "2021-09-21")
# 
# initial_temp_mean <- initial_temps %>%
#   group_by(date) %>%
#   dplyr::summarise(
#     mean_temp_c = mean(temp),
#     sd_temp_c = sd(temp),
#     n = n()
#   ) %>%
#   mutate(
#     se_temp_c = sd_temp_c/sqrt(n)
#   )

linetypes <- c(
  1, # "10 °C"
  1, # "13 °C (ambient)"
  1, # "16 °C"
  4, # "La Niña (18-14 °C)" 
  1, # "17 °C"
  1, # "18 °C"
  1, # "20 °C"
  4 # "El Niño (21-18 °C)"
)

temp_marna_ts_plot_day <- ggplot()+
  scale_x_date(date_breaks = "1 week", date_labels = "%D")+
  geom_rect(aes(xmin = as.Date(c("2021-09-21")), xmax = as.Date(c("2021-09-28")),
                ymin = -Inf, ymax = Inf), fill = "gray")+
  geom_point(data = dat_processed_marna_melt_reps, 
             aes(x = date, y = temp, 
                 color = treatment),
             size = 1, shape = 1)+
  geom_line(data = dat_processed_marna_melt, 
            aes(x = date, y = mean_temp_c,
                color = treatment, linetype = treatment), 
            linewidth = 1.5) +
  scale_linetype_manual(values = linetypes) +
  scale_color_viridis(discrete = T, option = "plasma")+
  scale_y_continuous( 
    breaks = pretty_breaks(n = 6),
    limits = c(9,22))+
  labs(
    y = "temperature (daily mean °C)",
    x = NULL,
    title = NULL,
    color = NULL,
    linetype = NULL
  )+
  gg_options()+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.75)
  )+
  guides(
    color = guide_legend(override.aes = list(size=5)),
    linetype = guide_legend(override.aes = list(linewidth = 0.5))
  );temp_marna_ts_plot_day

ggsave("./Figures/Marna_Sensor_Time_Series.pdf", temp_marna_ts_plot_day,
       width = 10, height = 6, units = "in")
ggsave("./Figures/Marna_Sensor_Time_Series.png", temp_marna_ts_plot_day,
       width = 10, height = 6, units = "in")
ggsave("./Figures/Marna_Sensor_Time_Series.jpg", temp_marna_ts_plot_day,
       width = 10, height = 6, units = "in")

# Visualize Marna Lab carbonate chemistry --------------------------------------
dat_raw <- read_excel("./Data/Sensor_Data/Urchins Fall 2021 - Adult & Larvae PCO2 - Data.xlsx",
                      skip = 1)

dat_summary <- dat_raw %>%
  dplyr::rename(
    date = "Date of Collection",
    treatment = "Treatment (PCO2 - Temp",
    TA_umol_per_kg = "TA (umol/kg)",
    TCO2_umol_per_kg = "TCO2 (umol/kg)",
    pCO2_uatm = "pCO2 (uatm)",
    pH_total = "pH (total)",
    pH_NBS = "pH (NBS) - do conversion for any samples dates after 2022-10-26",
    omega_arag = "Omega_arag",
    omega_ca = "Omega_ca"
  ) %>%
  select(
    date, treatment, TA_umol_per_kg, TCO2_umol_per_kg, pCO2_uatm, pH_total, pH_NBS, omega_arag, omega_ca
  ) %>%
  mutate(
    date = as.Date(date),
    treatment = as.factor(treatment)
  ) %>%  
  group_by(treatment) %>%
  dplyr::summarise(
    TA_umol_per_kg_avg = mean(TA_umol_per_kg),
    TCO2_umol_per_kg_avg = mean(TCO2_umol_per_kg),
    pCO2_uatm_avg = mean(pCO2_uatm),
    pH_total_avg = mean(pH_total),
    pH_NBS_avg = mean(pH_NBS),
    omega_arag_avg = mean(omega_arag),
    omega_ca_avg = mean(omega_ca)
  )

dat_carbo_ts <- dat_raw %>%
  dplyr::rename(
    tank_id = "Tank ID",
    date = "Date of Collection",
    treatment = "Treatment (PCO2 - Temp",
    TA_umol_per_kg = "TA (umol/kg)",
    TCO2_umol_per_kg = "TCO2 (umol/kg)",
    pCO2_uatm = "pCO2 (uatm)",
    pH_total = "pH (total)",
    pH_NBS = "pH (NBS) - do conversion for any samples dates after 2022-10-26",
    omega_arag = "Omega_arag",
    omega_ca = "Omega_ca"
  ) %>%
  dplyr::mutate(
    trt = as.factor(case_when(
      treatment == "Ambient 10" ~ "600 μatm, 10 °C",
      treatment == "Ambient 13" ~ "600 μatm, 13 °C",
      treatment == "Ambient 16" ~ "600 μatm, 16 °C",
      treatment == "Ambient 17" ~ "600 μatm, 17 °C",
      treatment == "Ambient 18" ~ "600 μatm, 18 °C",
      treatment == "Ambient 18-14 Heatwave" ~ "600 μatm, 18-14 °C (La Niña)",
      treatment == "Ambient 20" ~ "600 μatm, 20 °C",
      treatment == "Ambient 21-18 Heatwave" ~ "600 μatm, 21-18 °C (El Niño)",
      treatment == "High 10" ~ "1200 μatm, 10 °C",
      treatment == "High 13" ~ "1200 μatm, 13 °C",
      treatment == "High 16" ~ "1200 μatm, 16 °C",
      treatment == "High 17" ~ "1200 μatm, 17 °C",
      treatment == "High 18" ~ "1200 μatm, 18 °C",
      treatment == "High 20" ~ "1200 μatm, 20 °C"
    ))
  ) %>%
  dplyr::select(
    tank_id, date, trt, TA_umol_per_kg, TCO2_umol_per_kg, pCO2_uatm, pH_total, pH_NBS, omega_arag, omega_ca
  ) %>%
  dplyr::mutate(
    date = as.Date(date),
    trt = factor(trt, levels = c(
      "600 μatm, 10 °C",               
      "600 μatm, 13 °C",               
      "600 μatm, 16 °C",
      "600 μatm, 17 °C",
      "600 μatm, 18 °C",
      "600 μatm, 20 °C",
      "600 μatm, 18-14 °C (La Niña)",
      "600 μatm, 21-18 °C (El Niño)",
      "1200 μatm, 10 °C",
      "1200 μatm, 13 °C",
      "1200 μatm, 16 °C",
      "1200 μatm, 17 °C",
      "1200 μatm, 18 °C",
      "1200 μatm, 20 °C"
    ))
  ) 

# reshape data frame to long format
columns_carbo <- names(dat_carbo_ts %>% 
                   select(!c(TA_umol_per_kg, TCO2_umol_per_kg, pCO2_uatm, pH_total, pH_NBS, omega_arag, omega_ca)))

dat_carbo_melt <- reshape2::melt(dat_carbo_ts, id.vars = columns_carbo,
                                     variable.name = "carbo_param", value.name="value") %>%
  dplyr::filter(carbo_param == "pCO2_uatm")

dat_carbo_melt_summary <- dat_carbo_melt %>%
  dplyr::group_by(date, trt, carbo_param) %>%
  dplyr::summarise(
    mean_val = mean(value)
  )

levels(dat_carbo_melt_summary$trt)

linetypes_carbo <- c(
  1, # "600 μatm, 10 °C"              
  1, # "600 μatm, 13 °C"              
  1, # "600 μatm, 16 °C"             
  1, # "600 μatm, 17 °C"              
  1, # "600 μatm, 18 °C"
  1, # "600 μatm, 20 °C" 
  4, # "600 μatm, 18-14 °C (La Niña)"         
  4, # "600 μatm, 21-18 °C (El Niño)" 
  3, # "1200 μatm, 10 °C"             
  3, # "1200 μatm, 13 °C"
  3, # "1200 μatm, 16 °C"             
  3, # "1200 μatm, 17 °C"             
  3, # "1200 μatm, 18 °C"             
  3 # "1200 μatm, 20 °C"
  )

pco2_marna_ts_plot_day <- ggplot()+
  geom_rect(aes(xmin = as.Date(c("2021-10-04")), xmax = as.Date(c("2021-10-18")),
                ymin = -Inf, ymax = Inf), fill = "gray")+
  
  geom_point(data = dat_carbo_melt, 
             aes(x = date, y = value,
                 color = trt), 
              size = 1, shape = 1)+
  geom_line(data = dat_carbo_melt_summary,
            aes(x = date, y = mean_val,
                color = trt, linetype = trt),
            linewidth = 1.5) +
  scale_x_date(date_breaks = "1 week", date_labels = "%D") +
  scale_linetype_manual(values = linetypes_carbo) +
  scale_color_viridis(discrete = T, option = "plasma") +
  scale_fill_viridis(discrete = T, option = "plasma") +
  scale_y_continuous(breaks = scales::breaks_pretty(n=6)) +
  labs(
    y = expression("pCO"[2]*" (μatm)"),
    x = NULL,
    title = NULL,
    fill = NULL,
    color = NULL,
    linetype = NULL
  )+
  gg_options()+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.75)
  )+
  guides(
    color = guide_legend(override.aes = list(size=5)),
    linetype = guide_legend(override.aes = list(linewidth = 0.5))
  );pco2_marna_ts_plot_day

ggsave("./Figures/Marna_Carbo_Chem_Time_Series.pdf", pco2_marna_ts_plot_day,
       width = 10, height = 6, units = "in")
ggsave("./Figures/Marna_Carbo_Chem_Time_Series.png", pco2_marna_ts_plot_day,
       width = 10, height = 6, units = "in")
ggsave("./Figures/Marna_Carbo_Chem_Time_Series.jpg", pco2_marna_ts_plot_day,
       width = 10, height = 6, units = "in")

# Visualize pCO2 bouy time series ---------------------------------------------
data_directory <- "C:/Users/nates/Documents/Marna_Urchin_Heatwaves_2021/Marna_Urchin_Heatwaves_2021/Data/Remote_Sensing/pCO2_Bouy_Time_Series/Papa/"
# Sutton, A.J., C.L. Sabine, J.T. Mathis and A. Kozyr. 2014. A high-frequency 
# atmospheric and seawater pCO2 data set from 14 open ocean sites using a moored 
# autonomous system. ORNL/CDIAC-158, NDP-092A. http://cdiac.ornl.gov/oceans/Moorings/ndp092.html 
# Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, 
# US Department of Energy, Oak Ridge, Tennessee. doi: 10.3334/CDIAC/OTG.TSM_NDP092

### Get urchin and abalone data
files.list <- list.files(data_directory)

library(plyr)
library(dplyr)
library(zoo)
library(lubridate)
test <- read.csv(paste0(data_directory, "Papa_145W_50N_Jun2007_Jun2008.csv"), header = TRUE, skip = 4) %>%
  select(Date, pCO2.SW..sat...uatm.)

data_pco2_bouy_papa <- ldply(files.list, function(x) {
  dat = read.csv(paste0(data_directory,x), header = TRUE, skip = 4)
  dat = dat %>%
    select(Date, pCO2.SW..sat...uatm.)
  return(dat)
}
) %>%
  rename(
    pco2_uatm = pCO2.SW..sat...uatm.
  ) %>%
  mutate(
    Date = parse_date_time(Date,"%m-%d-%Y")
  ) %>%
  filter(
    pco2_uatm > 0
  )

summary_07_08 <- data_pco2_bouy_papa %>%
  filter(Date >= "2007-01-01" & Date <= "2008-01-01") %>%
  dplyr::summarise(
    n = n(),
    pco2_sd = sd(pco2_uatm),
    pco2_mean = mean(pco2_uatm)
  ) %>%
  mutate(
    pco2_se = pco2_sd/sqrt(n)
  )
  
summary_21_22 <- data_pco2_bouy_papa %>%
  filter(Date >= "2021-01-01" & Date <= "2022-01-01") %>%
  dplyr::summarise(
    n = n(),
    pco2_sd = sd(pco2_uatm),
    pco2_mean = mean(pco2_uatm)
  ) %>%
  mutate(
    pco2_se = pco2_sd/sqrt(n)
  )

# 9.12115% increase in average pco2 between 2007-8 and 2021-22

p_papa_mooring_ts <- ggplot() +
  geom_point(data = data_pco2_bouy_papa,
             aes(x = as.Date(Date), y = pco2_uatm)) +
  geom_smooth(data = data_pco2_bouy_papa,
              aes(x = as.Date(Date), y = pco2_uatm)) +
  scale_x_date(date_breaks = 'year', date_labels = "%Y") +
  geom_vline(xintercept = as.Date("2007-06-01"), linetype = 2, color = "red") +
  geom_hline(yintercept = 368.7311, linetype = 2, color = "red") +
  geom_vline(xintercept = as.Date("2021-09-01"), linetype = 2, color = "red") +
  geom_hline(yintercept = 402.3636, linetype = 2, color = "red") +
  labs(
    x = NULL,
    y = expression("pCO"[2]*" (µatm)"),
    title = "CO2 Time-series and Moorings Project, Pacific Ocean",
    subtitle = "Papa mooring - offshore of western Canada"
  ) +
  gg_options();p_papa_mooring_ts

ggsave("./Figures/pCO2_Mooring_TS.pdf", p_papa_mooring_ts,
       units = "in", width = 10, height = 6)


# Load necessary libraries
library(tidyverse)

# Set the directory containing the .txt files
data_dir <- "Data/pCO2_ph_ts_NCEI_0173932/1.1/data/0-data/AOI_Urchins"  # Modify this path as needed

# List all .txt files in the directory
file_list <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)

# Function to read a single file, skipping first 96 lines and adding source info
read_oceanographic_file <- function(filepath) {
  read_delim(
    file = filepath,
    delim = "\t",                 # Adjust if files use a different delimiter
    skip = 96,
    show_col_types = FALSE
  ) %>%
    mutate(source_file = basename(filepath))
}

# Read and combine all files
combined_data <- file_list %>%
  map_dfr(read_oceanographic_file)

# Pivot to long format
long_data <- combined_data %>%
  pivot_longer(
    cols = -c(source_file, datetime_utc),  # Retain these as metadata
    names_to = "variable",
    values_to = "value"
  )

# Preview the result
print(head(long_data))

ggplot(data = long_data, 
       aes(x = datetime_utc, y = value)) +
  geom_line(aes(color = variable)) +
  facet_grid(~source_file)




