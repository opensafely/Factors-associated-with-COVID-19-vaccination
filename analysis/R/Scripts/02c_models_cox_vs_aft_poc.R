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
dir.create(here::here("output", "models", "testing", "cox_vs_aft"), showWarnings = FALSE, recursive=TRUE)

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

## Run models on different subsets of data (i.e. 10/50/100 practices)
sample_size = c(10,50,100,250)
timings_all = data.frame(Method = c("Coxph",
                                    "Stratified coxph", 
                                    "AFT",
                                    "Stratified AFT",
                                    "AFT with RE"))

for (i in 1:length(sample_size)) {
  
  # Subset data
  data_sub <- data_cox %>%
    filter(practice_id_latest_active_registration %in% unique(data_cox$practice_id_latest_active_registration)[1:sample_size[i]])
  
  # Fit models and save model output
  ## Cox model - adjusted; baseline demographics, comorbs, geographical, flu, shielding
  mod.coxph.adj <- coxph(Surv(follow_up_time, covid_vax) ~
                           ageband + sex + ethnicity + morbid_obesity +
                           chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                           chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                           asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                           immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15 + 
                           rural_urban,
                         data = data_sub)
  
  write_rds(mod.coxph.adj, here::here("output", "models", "testing", "cox_vs_aft", 
                                      paste("mod_coxph_adj_", sample_size[i],".rds", sep = "")), compress="gz")
  
  ## Stratified Cox model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & practice as strata
  mod.stratcoxph.adj <- coxph(Surv(follow_up_time, covid_vax) ~
                                ageband + sex + ethnicity + morbid_obesity +
                                chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                                chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                                asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                                immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15 + 
                                strata(practice_id_latest_active_registration),
                              data = data_sub)
  
  write_rds(mod.stratcoxph.adj, here::here("output", "models", "testing", "cox_vs_aft", 
                                           paste("mod_stratcoxph_adj_", sample_size[i],".rds", sep = "")), compress="gz")
  
  # AFT model - adjusted; baseline demographics, comorbs, geographical, flu, shielding
  mod.aft.adj <- survreg(Surv(follow_up_time, covid_vax) ~
                           ageband + sex + ethnicity + morbid_obesity +
                           chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                           chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                           asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                           immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15,
                         dist = "lognormal",
                         data = data_sub)
  
  write_rds(mod.aft.adj, here::here("output", "models", "testing", "cox_vs_aft", 
                                    paste("mod_aft_adj_", sample_size[i],".rds", sep = "")), compress="gz")
  
  # Stratified AFT model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & practice as strata
  mod.strat.aft.adj <- survreg(Surv(follow_up_time, covid_vax) ~
                                 ageband + sex + ethnicity + morbid_obesity +
                                 chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                                 chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                                 asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                                 immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15 + 
                                 strata(practice_id_latest_active_registration),
                               data = data_sub)
  
  write_rds(mod.strat.aft.adj, here::here("output", "models", "testing", "cox_vs_aft", 
                                          paste("mod_strat_aft_adj_", sample_size[i],".rds", sep = "")), compress="gz")
  
  # AFT model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & practice as random effect
  mod.aft.re.adj <- survreg(Surv(follow_up_time, covid_vax) ~
                              ageband + sex + ethnicity + morbid_obesity +
                              chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                              chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                              asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                              immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15 + 
                              frailty(practice_id_latest_active_registration),
                            data = data_sub)
  
  write_rds(mod.aft.re.adj, here::here("output", "models", "testing", "cox_vs_aft", 
                                       paste("mod_aft_re_adj_", sample_size[i],".rds", sep = "")), compress="gz")
  
  
  # Tables
  ## Cox model - adjusted
  table_results_coxph.adj <- data.frame(summary(mod.coxph.adj)$coefficients) %>%
    rownames_to_column(var = "Variable") %>%
    mutate(LCI = round(exp.coef. - 1.96*se.coef., digits = 2),
           UCI = round(exp.coef. + 1.96*se.coef., digits = 2),
           `CoxPH HR (95% CI)` = paste(round(exp.coef., digits = 2),
                                       " (", LCI, " - ", UCI, ")", sep = "")) %>%
    select(Variable, `CoxPH HR (95% CI)`)
  
  ## Stratified Cox model
  table_results_mod.stratcoxph.adj <- data.frame(summary(mod.stratcoxph.adj)$coefficients) %>%
    rownames_to_column(var = "Variable") %>%
    mutate(LCI = round(exp.coef. - 1.96*se.coef., digits = 2),
           UCI = round(exp.coef. + 1.96*se.coef., digits = 2),
           `Strat CoxPH HR (95% CI)` = paste(round(exp.coef., digits = 2),
                                             " (", LCI, " - ", UCI, ")", sep = "")) %>%
    select(Variable, `Strat CoxPH HR (95% CI)`)
  
  # AFT model
  table_results_mod.aft.adj <- data.frame(summary(mod.aft.adj)$table) %>%
    rownames_to_column(var = "Variable") %>%
    filter(!(Variable %in% c("(Intercept)", "Log(scale)"))) %>%
    mutate(LCI = round(Value - 1.96*Std..Error, digits = 2),
           UCI = round(Value + 1.96*Std..Error, digits = 2),
           `AFT Time ratio (95% CI)` = paste(round(Value, digits = 2),
                                             " (", LCI, " - ", UCI, ")", sep = "")) %>%
    select(Variable, `AFT Time ratio (95% CI)`)
  
  # Stratified AFT model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & practice as strata
  table_results_mod.strat.aft.adj <- data.frame(summary(mod.strat.aft.adj)$table) %>%
    rownames_to_column(var = "Variable") %>%
    mutate(LCI = round(Value - 1.96*Std..Error, digits = 2),
           UCI = round(Value + 1.96*Std..Error, digits = 2),
           `Strat AFT Time ratio (95% CI)` = paste(round(Value, digits = 2),
                                                   " (", LCI, " - ", UCI, ")", sep = "")) %>%
    select(Variable, `Strat AFT Time ratio (95% CI)`)
  table_results_mod.strat.aft.adj <- table_results_mod.strat.aft.adj[2:62,]
  
  # AFT model with random effect
  table_results_mod.aft.re.adj <- data.frame(summary(mod.aft.re.adj)$table) %>%
    rownames_to_column(var = "Variable") %>%
    mutate(LCI = round(Value - 1.96*Std..Error, digits = 2),
           UCI = round(Value + 1.96*Std..Error, digits = 2),
           `AFT with RE Time ratio (95% CI)` = paste(round(Value, digits = 2),
                                                     " (", LCI, " - ", UCI, ")", sep = "")) %>%
    select(Variable, `AFT with RE Time ratio (95% CI)`)
  table_results_mod.aft.re.adj <- table_results_mod.aft.re.adj[2:62,]
  
  ## Combine tables
  table_results <- left_join(table_results_coxph.adj, table_results_mod.stratcoxph.adj, by = c("Variable")) %>%
    left_join(table_results_mod.aft.adj, by = c("Variable")) %>%
    left_join(table_results_mod.strat.aft.adj, by = c("Variable")) %>%
    left_join(table_results_mod.aft.re.adj, by = c("Variable"))
  
  write_csv(table_results, here::here("output", "models", "testing", "cox_vs_aft", 
                                      paste("table_results_", sample_size[i],".csv", sep = "")))
  
  
  # Timings
  ## Summary table
  timings <- data.frame(Method = c("Coxph",
                                   "Stratified coxph", 
                                   "AFT",
                                   "Stratified AFT",
                                   "AFT with RE"),
                        Time = NA)
  
  
  
  ## Cox model
  fit1 <- system.time(coxph(Surv(follow_up_time, covid_vax) ~
                              ageband + sex + ethnicity + morbid_obesity +
                              chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                              chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                              asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                              immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15                            data = data_sub))
  
  timings[1,2] <- fit1[3]
  
  ## Stratified Cox model
  fit2 <- system.time(coxph(Surv(follow_up_time, covid_vax) ~
                              ageband + sex + ethnicity + morbid_obesity +
                              chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                              chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                              asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                              immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15 + 
                              strata(practice_id_latest_active_registration),
                            data = data_sub))
  
  timings[2,2] <- fit2[3]
  
  # AFT model a
  fit3 <- system.time(survreg(Surv(follow_up_time, covid_vax) ~
                                ageband + sex + ethnicity + morbid_obesity +
                                chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                                chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                                asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                                immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15                              dist = "lognormal",
                              data = data_sub))
  
  timings[3,2] <- fit3[3]
  
  # AFT model b
  fit4 <- system.time(survreg(Surv(follow_up_time, covid_vax) ~
                                ageband + sex + ethnicity + morbid_obesity +
                                chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                                chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                                asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                                immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15 + 
                                strata(practice_id_latest_active_registration),
                              dist = "lognormal",
                              data = data_sub))
  timings[4,2] <- fit4[3]
  
  # AFT model c
  fit5 <- system.time(survreg(Surv(follow_up_time, covid_vax) ~
                                ageband + sex + ethnicity + morbid_obesity +
                                chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                                chronic_kidney_disease_all_stages_3_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                                asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                                immunosuppression_medication + imd + rural_urban + prior_covid + flu_vaccine + shielded + shielded_since_feb_15 + 
                                frailty(practice_id_latest_active_registration),
                              dist = "lognormal",
                              data = data_sub))
  timings[5,2] <- fit5[3]
  
  timings_all <- left_join(timings_all, timings, by = c("Method"))
  
  print(i)
  
}

## Plot timings
timings_plot <- timings_all %>%
  pivot_longer(!Method, names_to = "Number of practices", values_to = "Time") %>%
  mutate(`Number of practices` = ifelse(`Number of practices` == "Time.x", 10, `Number of practices`),
         `Number of practices` = ifelse(`Number of practices` == "Time.y", 50, `Number of practices`),
         `Number of practices` = ifelse(`Number of practices` == "Time.x.x", 100, `Number of practices`),
         `Number of practices` = ifelse(`Number of practices` == "Time.y.y", 250, `Number of practices`),
         `Number of practices` = as.numeric(`Number of practices`),
         Method = factor(Method, levels = c("Coxph",
                                            "Stratified coxph", 
                                            "AFT",
                                            "Stratified AFT",
                                            "AFT with RE")))

plot <- ggplot(timings_plot, aes(x = `Number of practices`, y = Time, colour = Method)) +
  geom_point() +
  geom_line()

ggsave(here::here("output", "models", "testing", "cox_vs_aft", "plot_models_fit_time.svg"),
       plot,
       units = "cm", width = 10, height = 8)


