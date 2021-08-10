######################################

# This script:
# - imports outputted model summaries from EMIS and TPP EHR backends
# - combined the model coefficients using inverse-variance-weighted averages
# - saves and plots outputs
# IMPORTANT: this script should only be run OFFLINE

######################################



# Preliminaries ----

## Import libraries
library('here')
library('tidyverse')
library('lubridate')
library('survival')


# import models ----
dummy_data <- TRUE

output_directory <- if(dummy_data){
 here("output", "combined")
} else{
  here("released_output", "output", "combined")
}

tidy_tpp <- read_csv(fs::path(output_directory, "tidy_tpp.csv"))
tidy_emis <- read_csv(fs::path(output_directory, "tidy_emis.csv"))

# Combine estimates ----

tidy_stack <- 
  bind_rows(
    tidy_tpp, 
    tidy_emis
  ) 

tidy_combined <- tidy_stack %>%
  group_by(term, variable, var_label, var_class, var_type, var_nlevels, contrasts, contrasts_type, reference_row, label) %>% 
  summarise(
    n_obs = sum(n_obs),
    n_event = sum(n_event),
    exposure = sum(exposure),
    estimate = weighted.mean(estimate, std.error^-2),
    std.error = sqrt(1/sum(std.error^-2)),
    statistic = estimate/std.error,
    p.value = 2 * pmin(pnorm(statistic), pnorm(-statistic)),
    conf.low = estimate + qnorm(0.025)*std.error,
    conf.high = estimate + qnorm(0.975)*std.error,
    
  ) %>% ungroup()

write_csv(tidy_combined, here("released_outputs", "combined", "meta_estimates.csv"))



