######################################

# This script:
# - imports processed data
# - fits several cox models using the coxph package
# - fits a mixed effects cox model using the coxme package
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
dir.create(here::here("output", "models", "testing"), showWarnings = FALSE, recursive=TRUE)

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

# ## Exclude practices with less than 100 registered patients
# 
# ### Registered patients counts
# practice_counts <- data_cox %>% 
#   group_by(practice_id_latest_active_registration) %>%
#   summarise(`Number_of_registered_patients` = n())
# 
# ## Exclude 
# data_cox_stratification <- data_cox %>%
#   filter(practice_id_latest_active_registration %in% subset(practice_counts, Number_of_registered_patients >= 100)$practice_id_latest_active_registration)



# MODELS ----

## Cox PH model - unadjusted
mod.coxph.unadj <- coxph(Surv(follow_up_time, covid_vax) ~ 1, data = data_cox)

write_rds(mod.coxph.unadj, here::here("output", "models", "testing", "mod_coxph_unadj.rds"), compress="gz")

## Cox PH model - adjusted; baseline demographics, comorbs, geographical, flu and shielding
mod.coxph.adj <- coxph(Surv(follow_up_time, covid_vax) ~
                         ageband + sex + ethnicity + morbid_obesity +
                         chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                         chronic_kidney_disease_all_stages_1_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                         asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                         immunosuppression_medication + imd + stp + region + rural_urban + prior_covid + flu_vaccine + 
                         shielded + shielded_since_feb_15,
                       data = data_cox)

write_rds(mod.coxph.adj, here::here("output", "models", "testing", "mod_coxph_adj.rds"), compress="gz")

# Cox model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & practice id
# mod.coxph.adjb <- coxph(Surv(follow_up_time, covid_vax) ~
#                          ageband + sex + ethnicity + morbid_obesity +
#                          chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
#                          chronic_kidney_disease_all_stages_1_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
#                          asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
#                          immunosuppression_medication + imd + region + rural_urban + flu_vaccine + shielded +
#                          shielded_since_feb_15 + practice_id_latest_active_registration,
#                        data = data_cox)
# 
# write_rds(mod.coxph.adjb, here::here("output", "models", "testing", "mod_coxph_adjb.rds"), compress="gz")

# Stratified Cox model - unadjusted;
mod.strat.coxph.unadj <- coxph(Surv(follow_up_time, covid_vax) ~ strata(practice_id_latest_active_registration),
                        data = data_cox)

write_rds(mod.strat.coxph.unadj, here::here("output", "models", "testing", "mod_strat_coxph_unadj.rds"), compress="gz")

# Stratified Cox model - adjusted; baseline demographics, comorbs, geographical, flu, shielding
mod.strat.coxph.adj <- coxph(Surv(follow_up_time, covid_vax) ~
                               ageband + sex + ethnicity + morbid_obesity +
                               chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                               chronic_kidney_disease_all_stages_1_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                               asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                               immunosuppression_medication + imd + region + rural_urban + prior_covid + flu_vaccine + shielded +
                               shielded_since_feb_15 + strata(practice_id_latest_active_registration),
                             data = data_cox)

write_rds(mod.strat.coxph.adj, here::here("output", "models", "testing", "mod_strat_coxph_adj.rds"), compress="gz")

# Mixed effects Cox model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & pracice as random effect
# mod.coxme.adj <- coxme(Surv(follow_up_time, covid_vax) ~
#                          ageband + sex + morbid_obesity +
#                          chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
#                          chronic_kidney_disease_all_stages_1_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
#                          asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
#                          immunosuppression_medication + imd + region + rural_urban + flu_vaccine + shielded +
#                          shielded_since_feb_15 + (1 | practice_id_latest_active_registration),
#                        data = data_cox_stratification)
# 
# write_rds(mod.coxme.adj, here::here("output", "models", "testing", "mod_coxme_adj.rds"), compress="gz")
