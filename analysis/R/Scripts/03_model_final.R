######################################

# This script:
# - imports processed data
# - fits a stratified cox model using the coxph package
# - saves model summaries (tables and figures)

######################################


# Preliminaries ----

## Import libraries
library('tidyverse')
library('lubridate')
library('survival')
library('coxme')
library('gtsummary')
library('gt')
library('survminer')
#library('ehahelper')

## Create output directory
dir.create(here::here("output", "models", "final"), showWarnings = FALSE, recursive=TRUE)

## Import processed data
data_tte <- read_rds(here::here("output", "data", "data_modelling.rds"))

## Converts logical to integer so that model coefficients print nicely in gtsummary methods
data_cox <- data_tte %>%
  mutate(
    across(
      where(is.logical),
      ~.x*1L
    )
  )


# MODELS ----

## Stratified Cox PH model - adjusted; baseline demographics, comorbs, geographical, flu and shielding 
mod.strat.coxph.adj <- coxph(Surv(follow_up_time, covid_vax) ~
                               ageband + sex + ethnicity + morbid_obesity +
                               chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                               chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                               asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                               immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded +
                               shielded_since_feb_15 + strata(practice_id_latest_active_registration),
                             data = data_cox)

write_rds(mod.strat.coxph.adj, here::here("output", "models", "final", "mod_strat_coxph_adj.rds"), compress="gz")

