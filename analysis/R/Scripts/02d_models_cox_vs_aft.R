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

## Stratified Cox model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & practice as strata
mod.stratcoxph.adj <- coxph(Surv(follow_up_time, covid_vax) ~
                              ageband + sex + ethnicity + morbid_obesity +
                              chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                              chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                              asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                              immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15 + 
                              rural_urban + region + strata(practice_id_latest_active_registration),
                            data = data_cox)

write_rds(mod.stratcoxph.adj, here::here("output", "models", "cox_vs_aft", "mod_stratcoxph_adj.rds"), compress="gz")

cox.zph(mod.stratcoxph.adj, transform = "identity")$table
cox.zph(mod.stratcoxph.adj, transform = "rank")$table
cox.zph(mod.stratcoxph.adj, transform = "log") $table
cox.zph(mod.stratcoxph.adj, transform = "km")$table

# ## AFT model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & practice as random effect
# mod.aft.re.adj <- survreg(Surv(follow_up_time, covid_vax) ~
#                             ageband + sex + ethnicity + morbid_obesity +
#                             chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
#                             chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
#                             asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
#                             immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15 + 
#                             rural_urban + region + frailty(practice_id_latest_active_registration),
#                           data = data_cox)
# 
# write_rds(mod.aft.re.adj, here::here("output", "models", "cox_vs_aft", "mod_aft_re_adj.rds"), compress="gz")
# 



  