######################################

# Some covariates used in the study are created from codelists of clinical conditions or 
# numerical values available on a patient's records.
# This script fetches all of the codelists identified in codelists.txt from OpenCodelists.

######################################


# --- IMPORT STATEMENTS ---

## Import code building blocks from cohort extractor package
from cohortextractor import (codelist, codelist_from_csv, combine_codelists)


# --- CODELISTS ---

## First COVID vaccination administration in EMIS
covid_vaccine_EMIS_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-covadm1.csv",
    system = "snomed",
    column = "code",
)

## Patients in long-stay nursing and residential care
nursing_residential_care_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-longres.csv",
    system = "snomed",
    column = "code",
)

## Ethnicity
ethnicity_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-eth2001.csv",
    system="snomed",
    column="code",
    category_column="grouping_16_id",
)

## Any other ethnicity code
ethnicity_other_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-non_eth2001.csv",
    system = "snomed",
    column = "code",
)

## Ethnicity not given - patient refused
ethnicity_not_given_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-eth_notgiptref.csv",
    system = "snomed",
    column = "code",
)

## Ethnicity not stated
ethnicity_not_stated_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-eth_notstated.csv",
    system = "snomed",
    column = "code",
)

# Ethnicity no record
ethnicity_no_record_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-eth_norecord.csv",
    system = "snomed",
    column = "code",
)



## BMI
bmi_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-bmi.csv",
    system = "snomed",
    column = "code",
)

# All BMI coded terms
bmi_stage_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-bmi_stage.csv",
    system = "snomed",
    column = "code",
)

# Severe Obesity
sev_obesity_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-sev_obesity.csv",
    system = "snomed",
    column = "code",
)

# Chronic heart disease codes
chronic_heart_disease_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-chd_cov.csv",
  system = "snomed",
  column = "code",
)

## Diabetes
diabetes_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-diab.csv",
  system = "snomed",
  column = "code",
)

## Chronic kidney disease diagnostic codes
chronic_kidney_disease_diagnostic_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ckd_cov.csv",
  system = "snomed",
  column = "code",
)

## Chronic kidney disease codes - all stages
chronic_kidney_disease_all_stages_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ckd15.csv",
  system = "snomed",
  column = "code",
)

# Chronic kidney disease codes-stages 3 - 5
chronic_kidney_disease_all_stages_3_5_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ckd35.csv",
  system = "snomed",
  column = "code",
)

## Severe mental illness
sev_mental_ill_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-sev_mental.csv",
  system = "snomed",
  column = "code",
)

## Learning disabilities
learning_disability_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-learndis.csv",
  system = "snomed",
  column = "code",
)

## Chronic Neurological Disease including Significant Learning Disorder
chronic_neuro_dis_inc_sig_learn_dis_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-cns_cov.csv",
  system = "snomed",
  column = "code",
)


## Asplenia or Dysfunction of the Spleen codes
asplenia_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-spln_cov.csv",
  system = "snomed",
  column = "code",
)

## Chronic Liver disease codes
chronis_liver_disease_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-cld.csv",
  system = "snomed",
  column = "code",
)

## Chronic Respiratory Disease
chronis_respiratory_disease_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-resp_cov.csv",
  system = "snomed",
  column = "code",
)

## Immunosuppression diagnosis codes
immunosuppression_diagnosis_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-immdx_cov.csv",
  system = "snomed",
  column = "code",
)

## Immunosuppression medication codes
immunosuppression_medication_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-immrx.csv",
  system = "snomed",
  column = "code",
)

## Flu Vaccine
flu_med_codes = codelist_from_csv(
  "codelists/opensafely-influenza-vaccination.csv",  
  system = "snomed",
  column = "snomed_id",
)

flu_clinical_given_codes = codelist_from_csv(
  "codelists/nhsd-primary-care-domain-refsets-flu_cod.csv",
  system = "snomed",
  column = "code",
)

flu_clinical_not_given_codes = codelist_from_csv(
  "codelists/nhsd-primary-care-domain-refsets-fludec_cod.csv",
  system = "snomed",
  column = "code",
)

## History of covid
covid_codes = codelist_from_csv(
  "codelists/opensafely-covid-identification.csv",
  system="icd10",
  column="icd10_code",
)

covid_primary_care_code = codelist_from_csv(
  "codelists/user-will-hulme-covid-identification-in-primary-care-probable-covid-clinical-code.csv",
  system = "snomed",
  column = "code",
)

covid_primary_care_positive_test = codelist_from_csv(
  "codelists/user-will-hulme-covid-identification-in-primary-care-probable-covid-positive-test.csv",
  system = "snomed",
  column = "code",
)

covid_primary_care_sequalae = codelist_from_csv(
  "codelists/user-will-hulme-covid-identification-in-primary-care-probable-covid-sequelae.csv",
  system = "snomed",
  column = "code",
)

## Shielding
high_risk_codes = codelist(
    ['1300561000000107'], system="snomed")

## Lower Risk from COVID-19 codes
not_high_risk_codes = codelist(
    ['1300591000000101', '1300571000000100'], system="snomed")
