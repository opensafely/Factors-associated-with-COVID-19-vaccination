######################################

# This script:
# imports data extracted by the cohort extractor
# fills in unknown ethnicity from GP records with ethnicity from SUS (secondary care)
# tidies missing values
# re-orders date variables so no negative time differences (only actually does anything for dummy data)
# standardises some variables (eg convert to factor) and derives some new ones
# saves processed one-row-per-patient dataset
# saves one-row-per-patient dataset for vaccines and for hospital admissions

######################################


# Preliminaries ----

## Import libraries
library('tidyverse')
library('lubridate')

## Import custom user functions from lib
source(here::here("lib", "utility_functions.R"))

# Import globally defined repo variables from
gbl_vars <- jsonlite::fromJSON(
  txt="./analysis/global-variables.json"
)
gbl_vars$run_date =date(file.info(here::here("metadata","extract_all.log"))$ctime)
#list2env(gbl_vars, globalenv())

## Output processed data to rds
dir.create(here::here("output", "data"), showWarnings = FALSE, recursive=TRUE)


# Process ----

## Print variable names
read_csv(
  here::here("output", "input_all.csv"),
  n_max=0,
  col_types = cols()
) %>%
names() %>%
print()

## Format columns (don't rely on defaults)
data_extract0 <- read_csv(
  here::here("output", "input_all.csv"),
  col_types = cols_only(

    # identifiers
    patient_id = col_integer(),
    practice_id = col_integer(),

    # demographic / administrative
    stp = col_character(),
    region = col_character(),
    imd = col_character(),
    care_home_type = col_character(),
    care_home = col_logical(),

    registered_at_latest = col_logical(),
    has_follow_up_previous_year = col_logical(),

    age = col_integer(),
    sex = col_character(),
    ethnicity = col_character(),
    ethnicity_6_sus = col_character(),
    #ethnicity_16 = col_character(),

    # dates
    dereg_date = col_date(format="%Y-%m-%d"),

    prior_positive_test_date = col_date(format="%Y-%m-%d"),
    prior_primary_care_covid_case_date = col_date(format="%Y-%m-%d"),
    prior_covidadmitted_date = col_date(format="%Y-%m-%d"),

    admitted_0_date = col_date(format="%Y-%m-%d"),
    admitted_1_date = col_date(format="%Y-%m-%d"),
    admitted_2_date = col_date(format="%Y-%m-%d"),
    admitted_3_date = col_date(format="%Y-%m-%d"),
    admitted_4_date = col_date(format="%Y-%m-%d"),
    admitted_5_date = col_date(format="%Y-%m-%d"),

    discharged_0_date = col_date(format="%Y-%m-%d"),
    discharged_1_date = col_date(format="%Y-%m-%d"),
    discharged_2_date = col_date(format="%Y-%m-%d"),
    discharged_3_date = col_date(format="%Y-%m-%d"),
    discharged_4_date = col_date(format="%Y-%m-%d"),
    discharged_5_date = col_date(format="%Y-%m-%d"),

    primary_care_probable_covid_1_date = col_date(format="%Y-%m-%d"),
    primary_care_probable_covid_2_date = col_date(format="%Y-%m-%d"),
    primary_care_probable_covid_3_date = col_date(format="%Y-%m-%d"),
    primary_care_probable_covid_4_date = col_date(format="%Y-%m-%d"),
    primary_care_probable_covid_5_date = col_date(format="%Y-%m-%d"),

    primary_care_suspected_covid_1_date = col_date(format="%Y-%m-%d"),
    primary_care_suspected_covid_2_date = col_date(format="%Y-%m-%d"),
    primary_care_suspected_covid_3_date = col_date(format="%Y-%m-%d"),
    primary_care_suspected_covid_4_date = col_date(format="%Y-%m-%d"),
    primary_care_suspected_covid_5_date = col_date(format="%Y-%m-%d"),

    covid_vax_1_date = col_date(format="%Y-%m-%d"),
    covid_vax_2_date = col_date(format="%Y-%m-%d"),
    covid_vax_3_date = col_date(format="%Y-%m-%d"),
    covid_vax_4_date = col_date(format="%Y-%m-%d"),

    covid_vax_pfizer_1_date = col_date(format="%Y-%m-%d"),
    covid_vax_pfizer_2_date = col_date(format="%Y-%m-%d"),
    covid_vax_pfizer_3_date = col_date(format="%Y-%m-%d"),
    covid_vax_pfizer_4_date = col_date(format="%Y-%m-%d"),

    covid_vax_az_1_date = col_date(format="%Y-%m-%d"),
    covid_vax_az_2_date = col_date(format="%Y-%m-%d"),
    covid_vax_az_3_date = col_date(format="%Y-%m-%d"),
    covid_vax_az_4_date = col_date(format="%Y-%m-%d"),

    positive_test_1_date = col_date(format="%Y-%m-%d"),
    positive_test_2_date = col_date(format="%Y-%m-%d"),
    primary_care_covid_case_1_date = col_date(format="%Y-%m-%d"),
    primary_care_covid_case_2_date = col_date(format="%Y-%m-%d"),
    emergency_1_date = col_date(format="%Y-%m-%d"),
    emergency_2_date = col_date(format="%Y-%m-%d"),
    covidadmitted_1_date = col_date(format="%Y-%m-%d"),
    covidadmitted_2_date = col_date(format="%Y-%m-%d"),
    coviddeath_date = col_date(format="%Y-%m-%d"),
    death_date = col_date(format="%Y-%m-%d"),

    # Put into 'category'
    bmi = col_character(),

    chronic_cardiac_disease = col_logical(),
    heart_failure = col_logical(),
    other_heart_disease = col_logical(),

    dialysis = col_logical(),

    diabetes = col_logical(),
    chronic_liver_disease = col_logical(),

    current_copd = col_logical(),
    cystic_fibrosis = col_logical(),
    other_resp_conditions = col_logical(),

    lung_cancer = col_logical(),
    haematological_cancer = col_logical(),
    cancer_excl_lung_and_haem = col_logical(),

    chemo_or_radio = col_logical(),
    solid_organ_transplantation = col_logical(),
    bone_marrow_transplant = col_logical(),
    sickle_cell_disease = col_logical(),
    permanant_immunosuppression = col_logical(),
    temporary_immunosuppression = col_logical(),
    asplenia = col_logical(),
    dmards = col_logical(),

    dementia = col_logical(),
    other_neuro_conditions = col_logical(),
    LD_incl_DS_and_CP = col_logical(),
    psychosis_schiz_bipolar = col_logical(),
    flu_vaccine = col_logical()

  ),

  na = character() # more stable to convert to missing later
)

## Fill in unknown ethnicity from GP records with ethnicity from SUS (secondary care)
data_extract0 <- data_extract0 %>%
  mutate(ethnicity = ifelse(ethnicity == "", ethnicity_6_sus, ethnicity)) %>%
  select(-ethnicity_6_sus)

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



#  SECTION TO SORT OUT BAD DUMMY DATA ----

## Extract all event date columns and collapse into one 'date' and one 'event' column and order by 
## date (for each patient_id)
data_dates_reordered_long <- data_extract %>%
  select(patient_id, matches("^(.*)_(\\d+)_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c("event", "index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, event, date) %>%
  group_by(patient_id, event) %>%
  mutate(
    index = row_number(),
    name = paste0(event,"_",index,"_date")
  ) %>%
  ungroup() %>%
  select(
    patient_id, name, event, index, date
  )

## Convert back to wide format (columns now in date order)
data_dates_reordered_wide <- data_dates_reordered_long %>%
  arrange(name, patient_id) %>%
  pivot_wider(
    id_cols=c(patient_id),
    names_from = name,
    values_from = date
  )

## Merge back non-date columns
data_extract_reordered <- left_join(
  data_extract %>% select(-matches("^(.*)_(\\d+)_date")),
  data_dates_reordered_wide,
  by="patient_id"
)

## Format variables (i.e, set factor levels)
data_processed <- data_extract_reordered %>%
  mutate(

    start_date = as.Date(gbl_vars$start_date), # i.e., this is interpreted later as [midnight at the _end of_ the start date] = [midnight at the _start of_ start date + 1], So that for example deaths on start_date+1 occur at t=1, not t=0.
    end_date = as.Date(gbl_vars$end_date),
    censor_date = pmin(end_date, death_date, na.rm=TRUE),

    sex = fct_case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      #sex == "I" ~ "Inter-sex",
      #sex == "U" ~ "Unknown",
      TRUE ~ NA_character_
    ),

    ageband = cut(
      age,
      breaks=c(-Inf, 18, 50, 60, 70, 80, Inf),
      labels=c("under 18", "18-49", "50s", "60s", "70s", "80+"),
      right=FALSE
    ),

    ethnicity = fct_case_when(
      ethnicity == "1" ~ "White",
      ethnicity == "4" ~ "Black",
      ethnicity == "3" ~ "South Asian",
      ethnicity == "2" ~ "Mixed",
      ethnicity == "5" ~ "Other",
      #TRUE ~ "Unknown",
      TRUE ~ NA_character_

    ),

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

    region = factor(region,
                    levels= c(
                      "East",
                      "East Midlands",
                      "London",
                      "North East",
                      "North West",
                      "South East",
                      "South West",
                      "West Midlands",
                      "Yorkshire and The Humber"
                    )
    ),
    
    stp = as.factor(stp),
    
    care_home_type = as.factor(care_home_type),

    bmi = as.factor(bmi),

    cause_of_death = fct_case_when(
      !is.na(coviddeath_date) ~ "covid-related",
      !is.na(death_date) ~ "not covid-related",
      TRUE ~ NA_character_),

    noncoviddeath_date = if_else(!is.na(death_date) & is.na(coviddeath_date), death_date, as.Date(NA_character_))

  ) %>%
  droplevels() %>%
  filter(
    !is.na(age),
    !is.na(sex),
    !is.na(imd),
    !is.na(ethnicity),
    !is.na(region)
  ) %>%
  ## TEMPORARY STEP TO REDUCE DATASET SIZE -- REMOVE FOR REAL RUN!
  sample_n(tbl=., size=min(c(100000, nrow(.))))


# Create one-row-per-event datasets ----

## For hospitalisation/discharge
data_admissions <- data_processed %>%
    select(patient_id, matches("^admitted\\_\\d+\\_date"), matches("^discharged\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(".value", "index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_drop_na = TRUE
    ) %>%
    select(patient_id, index, admitted_date=admitted, discharged_date = discharged) %>%
    arrange(patient_id, admitted_date)

## For covid in primary care
data_pr_suspected_covid <- data_processed %>%
  select(patient_id, matches("^primary_care_suspected_covid\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "suspected_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

data_pr_probable_covid <- data_processed %>%
  select(patient_id, matches("^primary_care_probable_covid\\_\\d+\\_date")) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "probable_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

## For vaccination
data_vax <- local({

  data_vax_all <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)

  data_vax_pf <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_pfizer\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_pf_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)

  data_vax_az <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_az\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_az_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)

  data_vax_all %>%
    left_join(data_vax_pf, by=c("patient_id", "date")) %>%
    left_join(data_vax_az, by=c("patient_id", "date")) %>%
    mutate(
      vaccine_type = fct_case_when(
        !is.na(vax_az_index) & is.na(vax_pf_index) ~ "Ox-AZ",
        is.na(vax_az_index) & !is.na(vax_pf_index) ~ "Pf-BN",
        is.na(vax_az_index) & is.na(vax_pf_index) ~ "Unknown",
        !is.na(vax_az_index) & !is.na(vax_pf_index) ~ "Both",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(patient_id, date)

})


# Save datasets as .rds files ----
write_rds(data_processed, here::here("output", "data", "data_all.rds"), compress="gz")
write_rds(data_vax, here::here("output", "data", "data_long_vax_dates.rds"), compress="gz")
write_rds(data_admissions, here::here("output", "data", "data_long_admission_dates.rds"), compress="gz")
write_rds(data_pr_probable_covid, here::here("output", "data", "data_long_pr_probable_covid_dates.rds"), compress="gz")
write_rds(data_pr_suspected_covid, here::here("output", "data", "data_long_pr_suspected_covid_dates.rds"), compress="gz")

