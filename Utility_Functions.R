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
