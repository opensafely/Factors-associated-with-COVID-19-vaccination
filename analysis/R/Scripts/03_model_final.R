######################################

# This script:
# - imports processed data
# - fit univariate and multivariable stratified cox model(s) using the coxph package
# - saves models

######################################


# Preliminaries ----

## Import libraries
library('here')
library('readr')
library('tidyr')
library('tidyverse')
library('lubridate')
library('survival')
library('gtsummary')
library('gt')
library('survminer')
library('glue')
library('fs')

## Create output directory
dir.create(here::here("output", "model"), showWarnings = FALSE, recursive=TRUE)

## Custom function
tidy_wald <- function(x, conf.int = TRUE, conf.level = .95, exponentiate = TRUE, ...) {
  
  # to use Wald CIs instead of profile CIs.
  ret <- broom::tidy(x, conf.int = FALSE, conf.level = conf.level, exponentiate = exponentiate)
  
  if(conf.int){
    ci <- confint.default(x, level = conf.level)
    if(exponentiate){ci = exp(ci)}
    ci <- as_tibble(ci, rownames = "term")
    names(ci) <- c("term", "conf.low", "conf.high")
    
    ret <- dplyr::left_join(ret, ci, by = "term")
  }
  ret
}

## Import processed data
data_tte <- read_rds(here::here("output", "data", "data_modelling.rds"))


# MODEL ----

## Stratified Cox PH model - adjusted; baseline demographics + comorbs, 
mod.strat.coxph.adj <- coxph(Surv(follow_up_time, covid_vax) ~
                               ageband + sex + ethnicity + imd + immunosuppression + ckd + 
                             chronic_respiratory_disease + diabetes + chronic_liver_disease + 
                             chronic_neuro_dis_inc_sig_learn_dis + chronic_heart_disease + asplenia + 
                             sev_mental_ill + morbid_obesity + strata(practice_id_latest_active_registration),
                             data = data_tte)


# Save outputs ----

## Save model
write_rds(mod.strat.coxph.adj, here::here("output", "model", "mod_strat_coxph_adj.rds"), compress="gz")

## Save a "tidy" copy of each model output. Create "dummy" emis/tpp outputs (identical) for use with combine script
#tidy_model <- broom.helpers::tidy_plus_plus(mod.strat.coxph.adj, tidy_fun = tidy_wald, exponentiate = FALSE)
tidy_model <- tidy_wald(mod.strat.coxph.adj, exponentiate = FALSE)

if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  for(backend in c("tpp", "emis")){
    write_csv(tidy_model, here("output", "model", glue("tidy_{backend}.csv")))
  }
} else {
  write_csv(tidy_model, here("output", "model", glue("tidy_{Sys.getenv('OPENSAFELY_BACKEND')}.csv")))
}

