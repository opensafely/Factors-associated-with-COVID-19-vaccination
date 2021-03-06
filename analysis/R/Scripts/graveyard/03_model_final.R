######################################

# This script:
# - imports processed data
# - fit univariate and multivariable stratified cox model(s) using the coxph package
# - saves models

######################################


# Preliminaries ----

## Import libraries
library('tidyverse')
library('lubridate')
library('survival')
library('gtsummary')
library('gt')
library('survminer')
#library('ehahelper')

## Create output directory
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  for(backend in c("tpp", "emis")){
    
    dir.create(here::here("output", {backend}, "models", "final"), showWarnings = FALSE, recursive=TRUE)
    
  }
} else {
  
  dir.create(here::here("output", "models", "final"), showWarnings = FALSE, recursive=TRUE)
  
  }

## Import processed data
data_tte <- read_rds(here::here("output", "data", "data_modelling.rds"))

# MODELS ----

## Stratified Cox PH model - adjusted; baseline demographics + comorbs, 
mod.strat.coxph.adj <- coxph(Surv(follow_up_time, covid_vax) ~
                               ageband + sex + ethnicity + imd + immunosuppression + ckd + 
                             chronic_respiratory_disease + diabetes + chronic_liver_disease + 
                             chronic_neuro_dis_inc_sig_learn_dis + chronic_heart_disease + asplenia + 
                             sev_mental_ill + morbid_obesity + strata(practice_id_latest_active_registration),
                             data = data_tte)

## Save model
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  for(backend in c("tpp", "emis")){
    write_rds(mod.strat.coxph.adj, here::here("output", {backend}, "models", "final", "mod_strat_coxph_adj.rds"), compress="gz")
  }
} else {
  write_rds(mod.strat.coxph.adj, here::here("output", "models", "final", "mod_strat_coxph_adj.rds"), compress="gz")
}


