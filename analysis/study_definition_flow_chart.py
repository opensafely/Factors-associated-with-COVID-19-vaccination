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
    #     NOT nursing_residential_care
    #     AND
    #     (sex = "M" OR sex = "F")
    #     AND
    #     imd > 0
    #     AND
    #     region
    #     AND
    #     rural_urban > 0
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
  
  ### Practice id
  practice_id = patients.registered_practice_as_of(
    "index_date",  # day before vaccine campaign start
    returning = "pseudo_id",
    return_expectations = {
      "int": {"distribution": "normal", "mean": 10, "stddev": 1},
      "incidence": 1,
    },
  ),
  
  ### Practice id at end (to check people moving practices during study)
  practice_id_at_end = patients.registered_practice_as_of(
    end_date,
    returning = "pseudo_id",
    return_expectations = {
      "int": {"distribution": "normal", "mean": 10, "stddev": 1},
      "incidence": 1,
    },
  ),
  
  ### Practice id at death (to check people moving practices during study)
  practice_id_at_death = patients.registered_practice_as_of(
    "death_date",
    returning = "pseudo_id",
    return_expectations = {
      "int": {"distribution": "normal", "mean": 10, "stddev": 1},
      "incidence": 1,
    },
  ),
  
  ### Same practice
  practice_id_same = patients.registered_with_one_practice_between(
      start_date = "index_date",
      end_date = end_date,
      return_expectations = {"incidence": 0.95},
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
      "category": {"ratios": {
        1: 0.125, 
        2: 0.125, 
        3: 0.125, 
        4: 0.125, 
        5: 0.125, 
        6: 0.125, 
        7: 0.125, 
        8: 0.125}
      },
    },
  ),
  
)


