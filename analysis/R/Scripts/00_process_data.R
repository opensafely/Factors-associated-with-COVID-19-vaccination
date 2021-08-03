######################################

# This script:
# - imports data extracted by the cohort extractor
# - combines ethnicity columns
# - calculates survival time
# - standardises some variables (eg convert to factor) and derives some new ones
# - saves processed one-row-per-patient dataset

######################################


# Preliminaries ----

## Import libraries
library('tidyverse')
library('lubridate')

## Custom functions
### tte
tte <- function(origin_date, event_date, censor_date, na.censor=FALSE){
  # returns time-to-event date or time to censor date, which is earlier
  
  if (na.censor)
    time <- event_date-origin_date
  else
    time <- pmin(event_date-origin_date, censor_date-origin_date, na.rm=TRUE)
  as.numeric(time)
}

### Factorise
fct_case_when <- function(...) {
  # uses dplyr::case_when but converts the output to a factor,
  # with factors ordered as they appear in the case_when's  ... argument
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]
  factor(dplyr::case_when(...), levels=levels)
}

## Output processed data to rds
dir.create(here::here("output", "data"), showWarnings = FALSE, recursive=TRUE)


# Process data ----

## Print variable names
read_csv(here::here("output", "input.csv"),
         n_max = 0,
         col_types = cols()) %>%
  names() %>%
  print()

## Read in data (don't rely on defaults)
data_extract0 <- read_csv(
  here::here("output", "input.csv"),
  col_types = cols_only(
    
    # Identifier
    patient_id = col_integer(),
    
    # Outcome
    covid_vax_1_date = col_date(format="%Y-%m-%d"),
    
    # Censoring
    death_date = col_date(format="%Y-%m-%d"),
    dereg_date = col_date(format="%Y-%m-%d"),
    
    # Demographic
    age = col_integer(),
    sex = col_character(),
    ethnicity = col_character(),
    ethnicity_other = col_date(format="%Y-%m-%d"),
    ethnicity_not_given = col_date(format="%Y-%m-%d"),
    ethnicity_not_stated = col_date(format="%Y-%m-%d"),
    ethnicity_no_record = col_date(format="%Y-%m-%d"),
    
    # Clinical measurements and comorbidities
    bmi = col_double(),
    bmi_stage_date = col_date(format="%Y-%m-%d"),
    sev_obesity = col_date(format="%Y-%m-%d"),
    chronic_heart_disease = col_logical(),
    diabetes = col_logical(),
    chronic_kidney_disease_diagnostic = col_date(format="%Y-%m-%d"),
    chronic_kidney_disease_all_stages = col_date(format="%Y-%m-%d"),
    chronic_kidney_disease_all_stages_3_5 = col_date(format="%Y-%m-%d"),
    sev_mental_ill = col_date(format="%Y-%m-%d"),
    learning_disability = col_date(format="%Y-%m-%d"),
    chronic_neuro_dis_inc_sig_learn_dis = col_date(format="%Y-%m-%d"),
    asplenia = col_logical(),
    chronic_liver_disease = col_logical(),
    chronic_respiratory_disease = col_date(format="%Y-%m-%d"),
    immunosuppression_diagnosis = col_date(format="%Y-%m-%d"),
    immunosuppression_medication = col_date(format="%Y-%m-%d"),
    
    # Geographical
    practice_id_at_start = col_double(),
    practice_id_at_end = col_double(),
    practice_id_at_death = col_double(),
    practice_id_at_dereg = col_double(),
    imd = col_character()
    # region = col_character(),
    # stp = col_character(),
    # rural_urban = col_character(),
    
    # Other
    # flu_vaccine = col_logical(),
    # shielded = col_logical(),
    # shielded_since_feb_15 = col_logical()
    # prior_covid_date = col_date(format="%Y-%m-%d")
  
  ),
  na = character() # more stable to convert to missing later
)

## Parse NAs
data_extract <- data_extract0 %>%
  mutate(across(
    .cols = where(is.character),
    .fns = ~na_if(.x, "")
  )) %>%
  mutate(across(
    .cols = c(where(is.numeric), -ends_with("_id")), #convert numeric+integer but not id variables
    .fns = ~na_if(.x, 0)
  )) %>%
  arrange(patient_id) %>%
  select(all_of((names(data_extract0))))



## Format columns (i.e, set factor levels)
data_processed <- data_extract %>%
  mutate(
    
    # Start date
    start_date = as.Date("2020-12-07", format = "%Y-%m-%d"),
    
    # End date
    end_date = as.Date("2021-03-17", format = "%Y-%m-%d"),
    
    # COVID vaccination
    covid_vax = as.integer(ifelse(is.na(covid_vax_1_date), 0, 1)),
    
    # Censoring
    censor_date = pmin(death_date, 
                       dereg_date, 
                       as.Date("2021-03-17", format = "%Y-%m-%d"), 
                       na.rm=TRUE),
    
    # Follow-up time
    follow_up_time = tte(start_date,
                         covid_vax_1_date,
                         censor_date),
    
    # Age
    ageband = cut(
      age,
      breaks = c(-Inf, 75, 80, 85, 90, 95, Inf),
      labels = c("70-74", "75-79", "80-84", "85-89", "90-94", "95+"),
      right = FALSE
    ),
    
    # Sex
    sex = fct_case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      #sex == "I" ~ "Inter-sex",
      #sex == "U" ~ "Unknown",
      TRUE ~ NA_character_
    ),
    
    # Ethnicity
    ethnicity =  ifelse(is.na(ethnicity) & !is.na(ethnicity_other), 17, 
                        ifelse(is.na(ethnicity) & !is.na(ethnicity_not_given), 18,
                               ifelse(is.na(ethnicity) & !is.na(ethnicity_not_stated), 19,
                                      ifelse(is.na(ethnicity) & !is.na(ethnicity_no_record), 20,
                                             ethnicity)))),
    
    ethnicity = ifelse(is.na(ethnicity), 20, ethnicity),
    
    ethnicity = fct_case_when(
      ethnicity == "1" ~ "White - British",
      ethnicity == "2" ~ "White - Irish",
      ethnicity == "3" ~ "White - Any other White background",
      ethnicity == "4" ~ "Mixed - White and Black Caribbean",
      ethnicity == "5" ~ "Mixed - White and Black African",
      ethnicity == "6" ~ "Mixed - White and Asian",
      ethnicity == "7" ~ "Mixed - Any other mixed background",
      ethnicity == "8" ~ "Asian or Asian British - Indian",
      ethnicity == "9" ~ "Asian or Asian British - Pakistani",
      ethnicity == "10" ~ "Asian or Asian British - Bangladeshi",
      ethnicity == "11" ~ "Asian or Asian British - Any other Asian background",
      ethnicity == "12" ~ "Black or Black British - Caribbean",
      ethnicity == "13" ~ "Black or Black British - African",
      ethnicity == "14" ~ "Black or Black British - Any other Black background",
      ethnicity == "15" ~ "Other ethnic groups - Chinese",
      ethnicity == "16" ~ "Other ethnic groups - Any other ethnic group",
      ethnicity == "17" ~ "Patients with any other ethnicity code",
      ethnicity == "18" ~ "Ethnicity not given - patient refused",
      ethnicity == "19" ~ "Ethnicity not stated",
      ethnicity == "20" ~ "Ethnicity not recorded",
      #TRUE ~ "Unknown",
      TRUE ~ NA_character_
    ),
    
    # Morbidly obese
    morbid_obesity_bmi = ifelse(bmi >= 40, 1, 0),
    morbid_obesity_date = ifelse(sev_obesity >= bmi_stage_date, 1, 0),
    morbid_obesity = ifelse(morbid_obesity_date == 1 & morbid_obesity_bmi == 0, 1, morbid_obesity_bmi),
    morbid_obesity = ifelse(is.na(morbid_obesity), 0, morbid_obesity),
    
    # ckd
    ckd = case_when(
      !is.na(chronic_kidney_disease_diagnostic) ~ TRUE,
      is.na(chronic_kidney_disease_all_stages) ~ FALSE,
      !is.na(chronic_kidney_disease_all_stages_3_5) & (chronic_kidney_disease_all_stages_3_5 >= chronic_kidney_disease_all_stages) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # Mental illness
    sev_mental_ill = !is.na(sev_mental_ill),
    
    # Learning disability
    learning_disability = !is.na(learning_disability),
    
    # CND inc LD
    chronic_neuro_dis_inc_sig_learn_dis = !is.na(chronic_neuro_dis_inc_sig_learn_dis),
    chronic_neuro_dis_inc_sig_learn_dis = chronic_neuro_dis_inc_sig_learn_dis | learning_disability,
    
    # CRD
    chronic_respiratory_disease = !is.na(chronic_respiratory_disease),
    
    # Immunosuppression
    immunosuppression_diagnosis = !is.na(immunosuppression_diagnosis),
    immunosuppression_medication = !is.na(immunosuppression_medication),
    immunosuppression = immunosuppression_diagnosis | immunosuppression_medication,
    
    # IMD
    imd = na_if(imd, "0"),
    imd = fct_case_when(
      imd == 1 ~ "1 most deprived",
      imd == 2 ~ "2",
      imd == 3 ~ "3",
      imd == 4 ~ "4",
      imd == 5 ~ "5 least deprived",
      #TRUE ~ "Unknown",
      TRUE ~ NA_character_
    ),
    
    # Practice id at death, dereg or end
    practice_id_latest_active_registration = ifelse(!is.na(death_date) & death_date < end_date, practice_id_at_death,
                                                    ifelse(!is.na(dereg_date) & dereg_date < end_date,
                                                           practice_id_at_dereg, practice_id_at_end)),
    
    # # Region
    # region = fct_case_when(
    #   region == "London" ~ "London",
    #   region == "East" ~ "East of England",
    #   region == "East Midlands" ~ "East Midlands",
    #   region == "North East" ~ "North East",
    #   region == "North West" ~ "North West",
    #   region == "South East" ~ "South East",
    #   region == "South West" ~ "South West",
    #   region == "West Midlands" ~ "West Midlands",
    #   region == "Yorkshire and The Humber" ~ "Yorkshire and the Humber",
    #   #TRUE ~ "Unknown",
    #   TRUE ~ NA_character_
    # ),
    # 
    # # stp
    # stp = as.factor(stp),
    
    # # Rural/urban
    # rural_urban = fct_case_when(
    #   rural_urban == 1 ~ "Urban - major conurbation",
    #   rural_urban == 2 ~ "Urban - minor conurbation",
    #   rural_urban == 3 ~ "Urban - city and town",
    #   rural_urban == 4 ~ "Urban - city and town in a sparse setting",
    #   rural_urban == 5 ~ "Rural - town and fringe",
    #   rural_urban == 6 ~ "Rural - town and fringe in a sparse setting",
    #   rural_urban == 7 ~ "Rural village and dispersed",
    #   rural_urban == 8 ~ "Rural village and dispersed in a sparse setting",
    #   #TRUE ~ "Unknown",
    #   TRUE ~ NA_character_
    # ),
    
    # # Prior covid
    # prior_covid = as.integer(ifelse(is.na(prior_covid_date), 0, 1))
    
  ) %>%
  filter(age >= 70,
         sex %in% c("Male", "Female"),
         !is.na(imd),
         !is.na(ethnicity)
         # !is.na(rural_urban),
  ) %>%
  select(patient_id, covid_vax, follow_up_time,
         ageband, sex, ethnicity, imd, immunosuppression, ckd, chronic_respiratory_disease,
         diabetes, chronic_liver_disease, chronic_neuro_dis_inc_sig_learn_dis, chronic_heart_disease,
         asplenia, sev_mental_ill, morbid_obesity, 
         practice_id_latest_active_registration, practice_id_at_death, practice_id_at_dereg, practice_id_at_end)

# Data for modelling
data_processed_modelling <- data_processed %>%
  mutate(practice_id_latest_active_registration = as.factor(practice_id_latest_active_registration)) %>%
  droplevels() %>%
  mutate(
    across(
      where(is.logical),
      ~.x*1L
    )
  )

# Save dataset as .rds files ----
write_rds(data_processed, here::here("output", "data", "data_all.rds"), compress="gz")
write_rds(data_processed_modelling, here::here("output", "data", "data_modelling.rds"), compress="gz")
