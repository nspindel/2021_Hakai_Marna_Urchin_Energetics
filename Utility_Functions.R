###############################################################################
# Utility functions                                                           #
#                                                                             #
###############################################################################

# set number of digits for visualizations for the session
options(scipen = 100000)

# Function to set up consistent ggplot aesthetics:
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

# Flexible path management
load_site_files <- function(region = "Newport_OR") {
  data_path <- file.path(here("Data"), region)
  all_files <- list.files(path = data_path, pattern = "\\.csv$", full.names = TRUE)
  message("Loading ", length(all_files), " files from: ", data_path)
  return(all_files)
}

load_ocean_file <- function(file_path, site_prefix_regex = ".*-(H[12])_") {
  df <- read_csv(file_path, show_col_types = FALSE) %>%
    clean_names()
  
  # Rename first three columns by position
  names(df)[1:3] <- c("datetime", "depth_ft", "value")
  
  df <- df %>%
    mutate(
      datetime = ymd_hms(datetime, quiet = TRUE),
      value = as.numeric(value),
      variable = basename(file_path) %>%
        str_remove(site_prefix_regex) %>%
        str_remove("\\.csv$") %>%
        recode(!!!canonical_var_names)
    ) %>%
    drop_na(datetime, value)
  
  return(df)
}


f_to_c <- function(temp_f) {
  (temp_f - 32) * 5 / 9
}