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
dir.create(here::here("output", "models", "testing"), showWarnings = FALSE, recursive=TRUE)

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

## Exclude practices with less than 100 registered patients

### Registered patients counts
practice_counts <- data_cox %>% 
  group_by(practice_id) %>%
  summarise(`Number_of_registered_patients` = n())

## Exclude 
data_cox_stratification <- data_cox %>%
  filter(practice_id %in% subset(practice_counts, Number_of_registered_patients >= 100)$practice_id,
         practice_id %in% subset(practice_counts, Number_of_registered_patients >= 100)$practice_id[1:10])

print(length(unique(data_cox_stratification$practice_id)))
# MODELS ----

# Mixed effects Cox model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & pracice as random effect
mod.coxme.adj <- coxme(Surv(follow_up_time, covid_vax) ~
                         ageband + (1 | practice_id),
                       data = data_cox_stratification)

write_rds(mod.coxme.adj, here::here("output", "models", "testing", "mod_test_coxme_adj.rds"), compress="gz")
