######################################

# This script:
# - imports processed data
# - fits a number of stratified coxph and AFT models to several subsets of the data

######################################


# Preliminaries ----

## Import libraries
library('tidyverse')
library('lubridate')
library('survival')
library('gtsummary')
library('gt')
library("survminer")

## Create output directory
dir.create(here::here("output", "models", "cox_vs_aft"), showWarnings = FALSE, recursive=TRUE)

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

## Stratified Cox model
mod.stratcoxph.adj <- coxph(Surv(follow_up_time, covid_vax) ~
                              ageband + sex + ethnicity + morbid_obesity + chronic_heart_disease + 
                              diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                              chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + 
                              chronic_neuro_dis_inc_sig_learn_dis + asplenia + chronic_liver_disease + 
                              chronis_respiratory_disease + immunosuppression_diagnosis +
                              immunosuppression_medication + imd + region + rural_urban + prior_covid + 
                              flu_vaccine + shielded + shielded_since_feb_15 + 
                              strata(practice_id_latest_active_registration),
                            data = data_cox)

write_rds(mod.stratcoxph.adj, here::here("output", "models", "cox_vs_aft", "mod_stratcoxph_adj.rds"), compress="gz")

## AFT models
### Weibull
mod.aft.adj.weibulla <- survreg(Surv(follow_up_time, covid_vax) ~
                                 ageband + sex + ethnicity + morbid_obesity + chronic_heart_disease + 
                                 diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                                 chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + 
                                 chronic_neuro_dis_inc_sig_learn_dis + asplenia + chronic_liver_disease + 
                                 chronis_respiratory_disease + immunosuppression_diagnosis +
                                 immunosuppression_medication + imd + region + rural_urban + prior_covid + 
                                 flu_vaccine + shielded + shielded_since_feb_15,                               dist = "weibull",
                               data = data_cox)

write_rds(mod.aft.adj.weibulla, here::here("output", "models", "cox_vs_aft", "mod_aft_weibull_adja.rds"), compress="gz")


### Weibull
mod.aft.adj.weibull <- survreg(Surv(follow_up_time, covid_vax) ~
                            ageband + sex + ethnicity + morbid_obesity + chronic_heart_disease + 
                            diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                            chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + 
                            chronic_neuro_dis_inc_sig_learn_dis + asplenia + chronic_liver_disease + 
                            chronis_respiratory_disease + immunosuppression_diagnosis +
                            immunosuppression_medication + imd + region + rural_urban + prior_covid + 
                            flu_vaccine + shielded + shielded_since_feb_15 + 
                            practice_id_latest_active_registration,
                          dist = "weibull",
                          data = data_cox)

write_rds(mod.aft.adj.weibull, here::here("output", "models", "cox_vs_aft", "mod_aft_weibull_adj.rds"), compress="gz")

# ### Exponential
# mod.aft.adj.exponential <- survreg(Surv(follow_up_time, covid_vax) ~
#                                     ageband + sex + ethnicity + morbid_obesity + chronic_heart_disease + 
#                                     diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
#                                     chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + 
#                                     chronic_neuro_dis_inc_sig_learn_dis + asplenia + chronic_liver_disease + 
#                                     chronis_respiratory_disease + immunosuppression_diagnosis +
#                                     immunosuppression_medication + imd + region + rural_urban + prior_covid + 
#                                     flu_vaccine + shielded + shielded_since_feb_15 + 
#                                     practice_id_latest_active_registration,
#                                   dist = "exponential",
#                                   data = data_cox)
# 
# write_rds(mod.aft.adj.exponential, here::here("output", "models", "cox_vs_aft", "mod_aft_exponential_adj.rds"), compress="gz")
# 
# ### Loglogistic
# mod.aft.adj.loglogistic <- survreg(Surv(follow_up_time, covid_vax) ~
#                                     ageband + sex + ethnicity + morbid_obesity + chronic_heart_disease + 
#                                     diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
#                                     chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + 
#                                     chronic_neuro_dis_inc_sig_learn_dis + asplenia + chronic_liver_disease + 
#                                     chronis_respiratory_disease + immunosuppression_diagnosis +
#                                     immunosuppression_medication + imd + region + rural_urban + prior_covid + 
#                                     flu_vaccine + shielded + shielded_since_feb_15 + 
#                                     practice_id_latest_active_registration,
#                                   dist = "loglogistic",
#                                   data = data_cox)
# 
# write_rds(mod.aft.adj.loglogistic, here::here("output", "models", "cox_vs_aft", "mod_aft_loglogistic_adj.rds"), compress="gz")
# 
# ### Lognormal
# mod.aft.adj.lognormal <- survreg(Surv(follow_up_time, covid_vax) ~
#                                     ageband + sex + ethnicity + morbid_obesity + chronic_heart_disease + 
#                                     diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
#                                     chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + 
#                                     chronic_neuro_dis_inc_sig_learn_dis + asplenia + chronic_liver_disease + 
#                                     chronis_respiratory_disease + immunosuppression_diagnosis +
#                                     immunosuppression_medication + imd + region + rural_urban + prior_covid + 
#                                     flu_vaccine + shielded + shielded_since_feb_15 + 
#                                     practice_id_latest_active_registration,
#                                   dist = "lognormal",
#                                   data = data_cox)
# 
# write_rds(mod.aft.adj.lognormal, here::here("output", "models", "cox_vs_aft", "mod_aft_lognormal_adj.rds"), compress="gz")
# 
# ## AFT models with random effects
# ### Weibull
# mod.aft.re.adj.weibull <- survreg(Surv(follow_up_time, covid_vax) ~
#                                     ageband + sex + ethnicity + morbid_obesity + chronic_heart_disease + 
#                                     diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
#                                     chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + 
#                                     chronic_neuro_dis_inc_sig_learn_dis + asplenia + chronic_liver_disease + 
#                                     chronis_respiratory_disease + immunosuppression_diagnosis +
#                                     immunosuppression_medication + imd + region + rural_urban + prior_covid + 
#                                     flu_vaccine + shielded + shielded_since_feb_15 + 
#                                     frailty(practice_id_latest_active_registration),
#                                   dist = "weibull",
#                                   data = data_cox)
# 
# write_rds(mod.aft.re.adj.weibull, here::here("output", "models", "cox_vs_aft", "mod_aft_re_weibull_adj.rds"), compress="gz")
# 
# ### Exponential
# mod.aft.re.adj.exponential <- survreg(Surv(follow_up_time, covid_vax) ~
#                                         ageband + sex + ethnicity + morbid_obesity + chronic_heart_disease + 
#                                         diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
#                                         chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + 
#                                         chronic_neuro_dis_inc_sig_learn_dis + asplenia + chronic_liver_disease + 
#                                         chronis_respiratory_disease + immunosuppression_diagnosis +
#                                         immunosuppression_medication + imd + region + rural_urban + prior_covid + 
#                                         flu_vaccine + shielded + shielded_since_feb_15 + 
#                                         frailty(practice_id_latest_active_registration),
#                                       dist = "exponential",
#                                       data = data_cox)
# 
# write_rds(mod.aft.re.adj.exponential, here::here("output", "models", "cox_vs_aft", "mod_aft_re_exponential_adj.rds"), compress="gz")
# 
# ### Loglogistic
# mod.aft.re.adj.loglogistic <- survreg(Surv(follow_up_time, covid_vax) ~
#                                         ageband + sex + ethnicity + morbid_obesity + chronic_heart_disease + 
#                                         diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
#                                         chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + 
#                                         chronic_neuro_dis_inc_sig_learn_dis + asplenia + chronic_liver_disease + 
#                                         chronis_respiratory_disease + immunosuppression_diagnosis +
#                                         immunosuppression_medication + imd + region + rural_urban + prior_covid + 
#                                         flu_vaccine + shielded + shielded_since_feb_15 + 
#                                         frailty(practice_id_latest_active_registration),
#                                       dist = "loglogistic",
#                                       data = data_cox)
# 
# write_rds(mod.aft.re.adj.loglogistic, here::here("output", "models", "cox_vs_aft", "mod_aft_re_loglogistic_adj.rds"), compress="gz")
# 
# ### Lognormal
# mod.aft.re.adj.lognormal <- survreg(Surv(follow_up_time, covid_vax) ~
#                                       ageband + sex + ethnicity + morbid_obesity + chronic_heart_disease + 
#                                       diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
#                                       chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + 
#                                       chronic_neuro_dis_inc_sig_learn_dis + asplenia + chronic_liver_disease + 
#                                       chronis_respiratory_disease + immunosuppression_diagnosis +
#                                       immunosuppression_medication + imd + region + rural_urban + prior_covid + 
#                                       flu_vaccine + shielded + shielded_since_feb_15 + 
#                                       frailty(practice_id_latest_active_registration),
#                                     dist = "lognormal",
#                                     data = data_cox)
# 
# write_rds(mod.aft.re.adj.lognormal, here::here("output", "models", "cox_vs_aft", "mod_aft_re_lognormal_adj.rds"), compress="gz")
# 
# 
#   