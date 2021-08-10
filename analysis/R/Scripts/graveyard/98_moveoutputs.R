######################################

# This script copies EMIS/TPP outputs needed for meta analysis
# Models are imported, model summaries are created in a standard format, and then saved in the model folder
# (could change to just copy existing outputs if they were suitably formatted)

######################################


## Import libraries ----
library('here')
library('readr')
library('survival')
library('tidyr')
library('glue')
library('fs')

## Create output directory ----
dir_create(here("output", "model"))

## Custom function ----
tidy_wald <- function(x, conf.int = TRUE, conf.level = .95, exponentiate = TRUE, ...) {
  
  # to use Wald CIs instead of profile CIs.
  ret <- broom::tidy(x, conf.int = FALSE, conf.level = conf.level, exponentiate = TRUE)
  
  if(conf.int){
    ci <- confint.default(x, level = conf.level)
    if(exponentiate){ci = exp(ci)}
    ci <- as_tibble(ci, rownames = "term")
    names(ci) <- c("term", "conf.low", "conf.high")
    
    ret <- dplyr::left_join(ret, ci, by = "term")
  }
  ret
}

# import data and model, and create tidy summary ----
# need to import data too otherwise tidy functions don't work

data_tte <- read_rds(here("output", "data", "data_modelling.rds"))
model <- read_rds(here("output", "model", "mod_strat_coxph_adj.rds"))
tidy_model <- broom.helpers::tidy_plus_plus(model, tidy_fun = tidy_wald, exponentiate = FALSE)

if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  for(backend in c("tpp", "emis")){
    write_csv(tidy_model, here("output", "model", glue("tidy_{backend}.csv")))
  }
} else {
  write_csv(tidy_model, here("output", "model", glue("tidy_{Sys.getenv('OPENSAFELY_BACKEND')}.csv")))
}

