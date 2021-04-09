######################################

# This script provides the formal specification of the study data that will be extracted from 
# the OpenSAFELY database.

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

## Import codelists from codelist.py (which pulls them from the codelist folder)
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
  population = patients.satisfying(
    """
        NOT has_died
        AND
        registered
        AND
        age >= 80
        AND
        has_follow_up_previous_year
        AND
        NOT nursing_residential_care
        AND
        (sex = "M" OR sex = "F")
        AND
        (smoking_status = "E" OR smoking_status = "E" OR smoking_status = "N")
        AND
        imd > 0
        AND
        region
        AND
        rural_urban > 0
        AND
        (ethnicity OR ethnicity_other OR ethnicity_not_given OR ethnicity_not_stated OR ethnicity_no_record)
        AND
        stp
        """,
    
    has_died = patients.died_from_any_cause(
      on_or_before="index_date",
      returning="binary_flag",
    ),
    
    registered = patients.satisfying(
      "registered_at_start",
      registered_at_start = patients.registered_as_of(start_date),
    ),
    
    has_follow_up_previous_year = patients.registered_with_one_practice_between(
      start_date = "index_date - 1 year",
      end_date = "index_date",
      return_expectations = {"incidence": 0.95},
    ),
    
    nursing_residential_care = patients.with_these_clinical_events(
        nursing_residential_care_codes,
        returning = "binary_flag",
        find_last_match_in_period = True,
        on_or_before = "index_date",
    ),
  ),
  
  
  ## OUTCOMES
  
  ### Any COVID vaccination (first dose)
  covid_vax_1_date = patients.with_tpp_vaccination_record(
    target_disease_matches = "SARS-2 CORONAVIRUS",
    on_or_after = "index_date + 1 day",  
    find_first_match_in_period = True,returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {"date": {"earliest": "2020-12-08",  # first vaccine administered on the 8/12
      "latest": end_date,}
    },
  ),
  
  
  ## CENSORING
  
  ### Death
  death_date = patients.died_from_any_cause(
    on_or_after = "index_date + 1 day",
    returning = "date_of_death",
    include_month = True,
    include_day = True,
  ),
  
  ### De-registration
  dereg_date = patients.date_deregistered_from_all_supported_practices(
    on_or_after = "index_date + 1 day",
    date_format = "YYYY-MM-DD",
  ),
  
  
  ## DEMOGRAPHIC INFORMATION
  
  ### Age
  age = patients.age_as_of(
    "2020-03-31",
    return_expectations = {
      "rate": "universal",
      "int": {"distribution": "population_ages"},
      "incidence" : 0.001
    },
  ),
  
  ### Sex
  sex = patients.sex(
    return_expectations = {
      "rate": "universal",
      "category": {"ratios": {"M": 0.49, "F": 0.51}},
    }
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
  
  
  ## CLINICAL MEASUREMENTS & COMORBIDITIES CONSIDERED AS POTENTIAL RISK FACTORS
  
  ### Smoking status
  smoking_status = patients.categorised_as(
    { "S": "most_recent_smoking_code = 'S' OR smoked_last_18_months",
      "E": """ 
               (most_recent_smoking_code = 'E' 
               OR 
               (most_recent_smoking_code = 'N' AND ever_smoked)) 
               AND 
               NOT smoked_last_18_months
               """,
      "N": "most_recent_smoking_code = 'N' AND NOT ever_smoked",
      "M": "DEFAULT",
    },
    
    return_expectations = {"category": {"ratios": {"S": 0.6, "E": 0.1, "N": 0.2, "M": 0.1}},
    "incidence" : 1},
    
    most_recent_smoking_code = patients.with_these_clinical_events(
      clear_smoking_codes,
      find_last_match_in_period = True,
      on_or_before = "index_date",
      returning = "category",
    ),
    
    ever_smoked = patients.with_these_clinical_events(
      filter_codes_by_category(clear_smoking_codes, include = ["S", "E"]),
      on_or_before = "index_date",
    ),
    
    smoked_last_18_months = patients.with_these_clinical_events(
      filter_codes_by_category(clear_smoking_codes, include = ["S"]),
      between = ["index_date - 18 months", "index_date"],
    ),
  ),
  
  ### BMI
  bmi = patients.with_these_clinical_events(
    bmi_codes,
    returning = "numeric_value",
    ignore_missing_values = True,
    find_last_match_in_period = True,
    on_or_before = "index_date",
    return_expectations = {
      "float": {"distribution": "normal", "mean": 25, "stddev": 5},
    },
  ),
  
  ### All BMI coded terms
  bmi_stage_date = patients.with_these_clinical_events(
    bmi_stage_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### Severe Obesity code recorded
  sev_obesity = patients.with_these_clinical_events(
    sev_obesity_codes,
    returning = "date",
    ignore_missing_values = True,
    find_last_match_in_period = True,
    on_or_after = "bmi_stage_date",
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### Chronic heart disease
  chronic_heart_disease = patients.with_these_clinical_events(
    chronic_heart_disease_codes,
    on_or_before = "index_date",
    returning = "binary_flag",
    return_expectations = {"incidence": 0.01},
  ),
  
  ### Diabetes
  diabetes = patients.with_these_clinical_events(
    diabetes_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    return_expectations = {"incidence": 0.01},
  ),
  
  ### Chronic kidney disease diagnostic codes
  chronic_kidney_disease_diagnostic = patients.with_these_clinical_events(
    chronic_kidney_disease_diagnostic_codes,
    on_or_before = "index_date",
    returning = "binary_flag",
    return_expectations = {"incidence": 0.01},
  ),
  
  ### Chronic kidney disease codes - all stages
  chronic_kidney_disease_all_stages = patients.with_these_clinical_events(
    chronic_kidney_disease_codes_all_stages,
    on_or_before = "index_date",
    returning = "binary_flag",
    return_expectations = {"incidence": 0.01},
  ),
  
  ### Chronic kidney disease codes-stages 3 - 5
  chronic_kidney_disease_all_stages_1_5 = patients.with_these_clinical_events(
    chronic_kidney_disease_codes_all_stages_1_5,
    on_or_before = "index_date",
    returning = "binary_flag",
    return_expectations = {"incidence": 0.01},
  ),
  
  ### Severe mental illness
  sev_mental_ill = patients.with_these_clinical_events(
    sev_mental_ill_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### Learning disabilities
  learning_disability = patients.with_these_clinical_events(
    learning_disability_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### Chronic neurological disease (including Significant Learning Disorder)
  chronic_neuro_dis_inc_sig_learn_dis = patients.with_these_clinical_events(
    chronic_neuro_dis_inc_sig_learn_dis_codes,
    returning = "date",
    find_first_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### Stroke
  stroke = patients.with_these_clinical_events(
    stroke_codes,
    on_or_before = "index_date",
    returning = "binary_flag",
    return_expectations = {"incidence": 0.01, },
  ),
  
  ### Asplenia or Dysfunction of the Spleen codes
  asplenia = patients.with_these_clinical_events(
    asplenia_codes,
    on_or_before = "index_date",
    returning = "binary_flag",
    return_expectations = {"incidence": 0.01, },
  ),
  
  ### Chronic liver disease
  chronic_liver_disease = patients.with_these_clinical_events(
    chronis_liver_disease_codes,
    on_or_before = "index_date",
    returning = "binary_flag",
    return_expectations = {"incidence": 0.01},
  ),
  
  ### Chronic respiratory disease
  chronis_respiratory_disease = patients.with_these_clinical_events(
    chronis_respiratory_disease_codes,
    returning = "date",
    find_first_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### Immunosuppression diagnosis
  immunosuppression_diagnosis = patients.with_these_clinical_events(
    immunosuppression_diagnosis_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### Immunosuppression medication
  immunosuppression_medication = patients.with_these_medications(
    immunosuppression_medication_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    on_or_after = "2020-07-01",
    date_format = "YYYY-MM-DD",
  ),
  
  
  ## GEOGRAPHICAL/DEPRIVATION
  
  ### Practice
  practice_id = patients.registered_practice_as_of(
    "index_date",  # day before vaccine campaign start
    returning = "pseudo_id",
    return_expectations = {
      "int": {"distribution": "normal", "mean": 1000, "stddev": 100},
      "incidence": 1,
    },
  ),
  
  ### Index of multiple deprivation
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
  
  ### Region - NHS England 9 regions
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
  
  ### STP (regional grouping of practices)
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
  
  
  ## OTHER FACTORS
  
  ### Flu Vaccine: Last 5 years prior to march 31st 2020
  flu_vaccine = patients.satisfying(
    """
        flu_vaccine_tpp_table>0 OR
        flu_vaccine_med>0 OR
        flu_vaccine_clinical>0
        """,
    
    flu_vaccine_tpp_table = patients.with_tpp_vaccination_record(
      target_disease_matches = "INFLUENZA",
      between = ["2015-04-01", "2020-03-31"], 
      returning = "binary_flag",
    ),
    
    flu_vaccine_med = patients.with_these_medications(
      flu_med_codes,
      between = ["2015-04-01", "2020-03-31"], 
      returning = "binary_flag",
    ),
    flu_vaccine_clinical = patients.with_these_clinical_events(
      flu_clinical_given_codes,
      ignore_days_where_these_codes_occur = flu_clinical_not_given_codes,
      between = ["2015-04-01", "2020-03-31"], 
      returning = "binary_flag",
    ),
    
    return_expectations = {"incidence": 0.5, },
  ),
  
  ### History of covid
  prior_covid_date = patients.with_these_clinical_events(
    combine_codelists(
      covid_primary_care_code,
      covid_primary_care_positive_test,
      covid_primary_care_sequalae,
    ),
    returning = "date",
    date_format = "YYYY-MM-DD",
    on_or_before = "index_date",
    find_first_match_in_period = True,
    return_expectations = {"rate": "exponential_increase"},
  ),
  
  ### Prior covid related hospital admission
  prior_covidadmitted_date = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_diagnoses = covid_codes,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  ### Shielding
  shielding_date = patients.with_these_clinical_events(
    shielding_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### Lower Risk from COVID-19 codes
  nonshield_dat = patients.with_these_clinical_events(
    nonshield_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  ### To represent household contact of shielding individual
  hhld_imdef_dat = patients.with_these_clinical_events(
    hhld_imdef_codes,
    returning = "date",
    find_last_match_in_period = True,
    on_or_before = "index_date",
    date_format = "YYYY-MM-DD",
  ),
  
  
)

