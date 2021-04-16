######################################

# This script:
# - imports processed data
# - fits a mixed effects cox model using the coxme package
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
#library('ehahelper')

## Create output directory
dir.create(here::here("output", "models"), showWarnings = FALSE, recursive=TRUE)

## Import processed data
data_tte <- read_rds(here::here("output", "data", "data_all.rds"))

## Converts logical to integer so that model coefficients print nicely in gtsummary methods
data_cox <- data_tte %>%
  mutate(
    across(
      where(is.logical),
      ~.x*1L
    )
  )


# MODELS ----

## Cox PH model - unadjusted
mod.coxph.unadj <- coxph(Surv(follow_up_time, covid_vax) ~ 1, data = data_cox)

write_rds(mod.coxph.unadj, here::here("output", "models", "mod_coxph_unadj.rds"), compress="gz")

## Cox PH model - adjusted; baseline demographics, comorbs, geographical, flu and shielding 
mod.coxph.adj <- coxph(Surv(follow_up_time, covid_vax) ~ 
                         ageband + sex + ethnicity + morbid_obesity +
                         chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                         chronic_kidney_disease_all_stages_1_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                         asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                         immunosuppression_medication + imd + stp + region + rural_urban + flu_vaccine + shielded +
                         shielded_since_feb_15,
                       data = data_cox)

write_rds(mod.coxph.adj, here::here("output", "models", "mod_coxph_adj.rds"), compress="gz")

# Cox model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & practice id
mod.coxph.adjb <- coxph(Surv(follow_up_time, covid_vax) ~
                         ageband + sex + ethnicity + morbid_obesity +
                         chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                         chronic_kidney_disease_all_stages_1_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                         asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                         immunosuppression_medication + imd + region + rural_urban + flu_vaccine + shielded +
                         shielded_since_feb_15 + strata(practice_id),
                       data = data_cox)

write_rds(mod.coxph.adjb, here::here("output", "models", "mod_coxph_adjb.rds"), compress="gz")

# Mixed effects Cox model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & pracice as random effect
mod.coxme.adj <- coxme(Surv(follow_up_time, covid_vax) ~
                          ageband + sex + ethnicity + morbid_obesity +
                          chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                          chronic_kidney_disease_all_stages_1_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                          asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                          immunosuppression_medication + imd + region + rural_urban + flu_vaccine + shielded +
                          shielded_since_feb_15 + (1 | practice_id),
                       sparse=c(50,.03),
                       data = data_cox)

write_rds(mod.coxme.adj, here::here("output", "models", "mod_coxme_adj.rds"), compress="gz")