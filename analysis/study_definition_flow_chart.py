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
end_date = "2021-04-01"

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
    #     NOT has_died
    #     AND
    #     registered
    #     AND
    #     age >= 80
    #     AND
    #     has_follow_up_previous_year
    #     AND
    #     NOT care_home_type
    #     AND
    #     (sex = "M" OR sex = "F")
    #     AND
    #     imd > 0
    #     AND
    #     region
    #     AND
    #     rural_urban
    #     AND
    #     (ethnicity OR ethnicity_other OR ethnicity_not_given OR ethnicity_not_stated OR ethnicity_no_record)
    #     AND
    #     stp
    #     """,
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
  care_home_type = patients.care_home_status_as_of(
    "index_date",
    categorised_as = {
      "Carehome": """
              IsPotentialCareHome
              AND LocationDoesNotRequireNursing='Y'
              AND LocationRequiresNursing='N'
            """,
      "Nursinghome": """
              IsPotentialCareHome
              AND LocationDoesNotRequireNursing='N'
              AND LocationRequiresNursing='Y'
            """,
      "Mixed": "IsPotentialCareHome",
      "": "DEFAULT",  # use empty string
    },
    return_expectations = {
      "category": {"ratios": {"Carehome": 0.05, "Nursinghome": 0.05, "Mixed": 0.05, "": 0.85, }, },
      "incidence": 1,
    },
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
  
  ## Region - NHS England 9 regions
  region = patients.registered_practice_as_of(
    "index_date",
    returning = "nuts1_region_name",
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "North East": 0.1,
          "North West": 0.1,
          "Yorkshire and The Humber": 0.1,
          "East Midlands": 0.2,
          "West Midlands": 0.1,
          "East": 0.1,
          "London": 0.2,
          "South East": 0.1,},},
    },
  ),
  
  ## STP (regional grouping of practices)
  stp = patients.registered_practice_as_of("index_date",
                                           returning = "stp_code",
                                           return_expectations = {
                                             "rate": "universal",
                                             "category": {
                                               "ratios": {
                                                 "STP1": 0.1,
                                                 "STP2": 0.1,
                                                 "STP3": 0.1,
                                                 "STP4": 0.1,
                                                 "STP5": 0.1,
                                                 "STP6": 0.1,
                                                 "STP7": 0.1,
                                                 "STP8": 0.1,
                                                 "STP9": 0.1,
                                                 "STP10": 0.1,}},
                                           },
  ),
  
  ### Urban vs rural
  rural_urban = patients.address_as_of(
    "index_date",
    returning = "rural_urban_classification",
    return_expectations = {
      "rate": "universal",
      "category": {"ratios": {"rural": 0.1, "urban": 0.9}},
    },
  ),
  
  ### Ethnicity
  ethnicity = patients.with_these_clinical_events(
    ethnicity_codes,
    returning = "category",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    return_expectations = {
      "category": {
        "ratios": {
          "1": 0.0625,
          "2": 0.0625,
          "3": 0.0625,
          "4": 0.0625,
          "5": 0.0625,
          "6": 0.0625,
          "7": 0.0625,
          "8": 0.0625,
          "9": 0.0625,
          "10": 0.0625,
          "11": 0.0625,
          "12": 0.0625,
          "13": 0.0625,
          "14": 0.0625,
          "15": 0.0625,
          "16": 0.0625,
        }
      },
      "rate": "universal",
    },
  ),
  
  ### Any other ethnicity code
  ethnicity_other = patients.with_these_clinical_events(
    ethnicity_other_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### Ethnicity not given - patient refused
  ethnicity_not_given = patients.with_these_clinical_events(
    ethnicity_not_given_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### Ethnicity not stated
  ethnicity_not_stated = patients.with_these_clinical_events(
    ethnicity_not_stated_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### Ethnicity no record
  ethnicity_no_record = patients.with_these_clinical_events(
    ethnicity_no_record_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ), 
  
)


