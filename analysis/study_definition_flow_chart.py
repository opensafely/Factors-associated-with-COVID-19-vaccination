######################################

# This script extracts all data relating to the study population variables so that inclusion/exclusion
# numbers can be calculated in respect to the study population

######################################

# --- IMPORT STATEMENTS ---

## Import code building blocks from cohort extractor package
from cohortextractor import (
  StudyDefinition,
  patients,
  codelist_from_csv,
  codelist,
  filter_codes_by_category,
  combine_codelists,
)

# Import codelists
from codelists import *
  
  
# --- DEFINE STUDY POPULATION ---
  
## Define study start and end variables explicitly
start_date = "2020-12-07"
end_date = "2021-03-17"

## Define study population and variables
study = StudyDefinition(
  
  # Configure the expectations framework
  default_expectations={
    "date": {"earliest": "1970-01-01", "latest": end_date},
    "rate": "uniform",
    "incidence": 0.2,
  },
  
  # Set index date to start date
  index_date = start_date,
  
  # Define the study population
  population = patients.all(
    # """
       #  NOT has_died
       #  AND
       #  registered
       #  AND
       #  age >= 70
       #  AND
       #  has_follow_up_previous_year
       #  AND
       #  NOT nursing_residential_care
       #  AND
       #  (sex = "M" OR sex = "F")
       #  AND
       #  imd != "0"
    #     """,
  ),
  
  # Outcome
  
  ### Any COVID vaccination (first dose)
  covid_vax_1_date = patients.with_vaccination_record(
    returning = "date",
    tpp = {"target_disease_matches": "SARS-2 CORONAVIRUS",},
    emis = {"procedure_codes": covid_vaccine_EMIS_codes,},
    find_first_match_in_period = True,
    on_or_after = "index_date + 1 day",
    date_format = "YYYY-MM-DD",
    return_expectations = {"date": 
        {"earliest": "2020-12-08",  # first vaccine administered on the 8/12
          "latest": end_date,}
    },
  ),
  
  # Inclusion/exclusion variables
  
  ## Alive
  has_died = patients.died_from_any_cause(
      on_or_before="index_date",
      returning="binary_flag",
  ),
    
  ## Registered
  registered = patients.satisfying(
    "registered_at_start",
    registered_at_start = patients.registered_as_of(start_date),
  ),
  
  ## Age
  age = patients.age_as_of(
    "2020-03-31",
    return_expectations = {
      "rate": "universal",
      "int": {"distribution": "population_ages"},
      "incidence" : 0.001
    },
  ),
  
  ## At least one year of follow-up
  has_follow_up_previous_year = patients.registered_with_one_practice_between(
    start_date = "index_date - 1 year",
    end_date = "index_date",
    return_expectations = {"incidence": 0.95},
  ),
  
  ## Care home
  nursing_residential_care = patients.with_these_clinical_events(
      nursing_residential_care_codes,
      returning = "binary_flag",
      find_last_match_in_period = True,
      on_or_before = "index_date",
  ),
  
  ### Sex
  sex = patients.sex(
    return_expectations = {
      "rate": "universal",
      "category": {"ratios": {"M": 0.49, "F": 0.51}},
    }
  ),
  
  ## Index of multiple deprivation
  imd = patients.categorised_as(
    {"0": "DEFAULT",
      "1": """index_of_multiple_deprivation >=1 AND index_of_multiple_deprivation < 32844*1/5""",
      "2": """index_of_multiple_deprivation >= 32844*1/5 AND index_of_multiple_deprivation < 32844*2/5""",
      "3": """index_of_multiple_deprivation >= 32844*2/5 AND index_of_multiple_deprivation < 32844*3/5""",
      "4": """index_of_multiple_deprivation >= 32844*3/5 AND index_of_multiple_deprivation < 32844*4/5""",
      "5": """index_of_multiple_deprivation >= 32844*4/5 """,
    },
    index_of_multiple_deprivation = patients.address_as_of(
      "index_date",
      returning = "index_of_multiple_deprivation",
      round_to_nearest = 100,
    ),
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "0": 0.01,
          "1": 0.20,
          "2": 0.20,
          "3": 0.20,
          "4": 0.20,
          "5": 0.19,
        }},
    },
  ),
  
  ### Death
  death_date = patients.died_from_any_cause(
    on_or_after = "index_date + 1 day",
    returning = "date_of_death",
    date_format = "YYYY-MM-DD",
  ),
  
)


