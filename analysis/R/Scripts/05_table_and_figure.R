######################################

# This script:
# - imports cox model
# - saves model summaries (tables and figures)

######################################


# Preliminaries ----

## Import libraries
library('tidyverse')
library('lubridate')
library('survival')
library('gtsummary')
library('gt')
library('survminer')

## Create output directory
dir.create(here::here("output", "model"), showWarnings = FALSE, recursive=TRUE)

## Import processed data
data_tte <- read_rds(here::here("output", "data", "data_modelling.rds"))

## Import model Stratified Cox PH model
mod.strat.coxph.adj <- read_rds(here::here("output", "model", "mod_strat_coxph_adj.rds"))

## Function to plot stratified cox model
forest_from_gt <- function(gt_obj){
  
  # Extract model information from tbl_regression object
  plot_data <- gt_obj %>%
    as_gt() %>%
    .$`_data` %>%
    filter(
      !is.na(term)
    ) %>%
    mutate(
      var_label = if_else(row_type=="label", "", var_label),
      label = if_else(reference_row %in% TRUE, paste0(label, " (ref)"),label),
      estimate = if_else(reference_row %in% TRUE, 1, estimate),
      variable = fct_inorder(variable),
      variable_card = as.numeric(variable)%%2,
    ) %>%
    group_by(variable) %>%
    mutate(
      variable_card = if_else(row_number()!=1, 0, variable_card),
      level = fct_rev(fct_inorder(paste(variable, label, sep="__"))),
      level_label = label
    ) %>%
    ungroup() %>%
    droplevels()
  
  var_lookup <- plot_data$var_label
  var_lookup[36] <- "Clinical Risk Groups"
  var_lookup[42] <- "Other Groups"
  names(var_lookup) <- plot_data$variable
  
  level_lookup <- plot_data$level
  names(level_lookup) <- str_to_title(gsub("_", " ", plot_data$level_label))
  names(level_lookup)[1] <- "70-74 (ref)"
  names(level_lookup)[7] <- "Female (ref)"
  names(level_lookup)[9] <- "White - British (ref)"
  names(level_lookup)[28] <- "1 (most deprived) (ref)"
  names(level_lookup)[32] <- "5 (least deprived)"
  names(level_lookup)[34] <- "Chronic Kidney Disease"
  names(level_lookup)[38] <- "Chronic Neurological Disease (including learning disablilty)"
  names(level_lookup)[41] <- "Sev Mental Illness"
  
  # Plot
  ggplot(plot_data) +
    geom_point(aes(x=estimate, y=level)) +
    geom_linerange(aes(xmin=conf.low, xmax=conf.high, y=level)) +
    geom_vline(aes(xintercept=1), colour='black', alpha=0.8)+
    facet_grid(rows=vars(variable), scales="free_y", switch="y", space="free_y", labeller = labeller(variable = var_lookup))+
    scale_x_log10()+
    scale_y_discrete(breaks=level_lookup, labels=names(level_lookup))+
    geom_rect(aes(alpha = variable_card), xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf, fill='grey', colour="transparent") +
    scale_alpha_continuous(range=c(0,0.3), guide=FALSE)+
    labs(
      y="",
      x="Hazard ratio",
      colour=NULL
    ) +
    theme_minimal() +
    theme(
      strip.placement = "outside",
      strip.background = element_rect(fill="transparent", colour="transparent"),
      strip.text.y.left = element_text(angle = 0, hjust=1),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.spacing = unit(0, "lines")
    )
}


# Counts ----
obs <- data_tte %>%
  select(-patient_id, -covid_vax, - follow_up_time, -practice_id_latest_active_registration) %>%
  tbl_summary()

obs$inputs$data <- NULL

obs <- obs$table_body %>%
  select(group = variable, variable = label, n_obs = stat_0) %>%
  separate(n_obs, c("n_obs","perc"), sep = "([(])") %>%
  mutate(n_obs = gsub(" ", "", n_obs),
         n_obs = as.numeric(gsub(",", "", n_obs))) %>%
  filter(!(is.na(n_obs))) %>%
  select(-perc)

events <- data_tte %>%
  filter(covid_vax == 1) %>%
  select(-patient_id, -covid_vax, - follow_up_time, -practice_id_latest_active_registration) %>%
  tbl_summary()

events$inputs$data <- NULL

events <- events$table_body %>%
  select(group = variable, variable = label, n_events = stat_0) %>%
  separate(n_events, c("n_events","perc"), sep = "([(])") %>%
  mutate(n_events = gsub(" ", "", n_events),
         n_events = as.numeric(gsub(",", "", n_events))) %>%
  filter(!(is.na(n_events))) %>%
  select(-perc)

obs_events <- left_join(obs, events) %>%
  mutate(n_obs = ifelse(n_obs < 5, "[redacted]", n_obs),
         n_events = ifelse(n_events < 5, "[redacted]", n_events))

head(obs_events)

write_csv(obs_events, here::here("output",  "model", "counts_redacted.csv"))


# Output model results ----

## Summary table
#tab_mod1 <- gtsummary::tbl_regression(mod.strat.coxph.adj, exp = TRUE)
#head(tab_mod1$table_body)
#gtsave(tab_mod1 %>% as_gt(), here::here("output", "models", "final", "tab_strat_coxph.html"))
#write_csv(tab_mod1$table_body, here::here("output",  "models", "final", "tab_strat_coxph.csv"))

tab_mod1  <- summary(mod.strat.coxph.adj)$coefficients  %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Variable") %>%
  mutate(LCI = round(exp(coef - 1.96*`se(coef)`), digits = 2),
         UCI = round(exp(coef + 1.96*`se(coef)`), digits = 2),
         HR = round(exp(coef), digits = 2),
         `95% CI` = paste(" (", LCI, " - ", UCI, ")", sep = ""),
         `p-value` = round(`Pr(>|z|)`, digits = 4)) %>%
  select(Variable, HR, `95% CI`, `p-value`)

write_csv(tab_mod1, here::here("output",  "model", "tab_strat_coxph.csv"))


tidy_coxph <- function(x, conf.int = TRUE, conf.level = .95, exponentiate = TRUE, ...) {
  
  # to use Wald CIs instead of profile CIs.
  ret <- broom::tidy(x, conf.int = FALSE, conf.level = conf.level, exponentiate = exponentiate)
  
  if(conf.int){
    ci <- confint.default(x, level = conf.level)
    if(exponentiate){ci = exp(ci)}
    ci <- as_tibble(ci, rownames = "term")
    names(ci) <- c("term", "conf.low", "conf.high")
    
    ret <- dplyr::left_join(ret, ci, by = "term")
  }
  ret
}


## Summary table for plot
tbl_summary <- tbl_regression(
  x = mod.strat.coxph.adj,
  pvalue_fun = ~style_pvalue(.x, digits=3),
  tidy_fun = tidy_coxph,
  exponentiate= TRUE,
  label = list(ageband = "Age Band", sex = "Sex", ethnicity = "Ethnicity", 
               imd = "IMD")
)

tbl_summary$table_body$variable <- str_to_title(gsub("_", " ", tbl_summary$table_body$variable))

## Forest plot
plot_coxph <- forest_from_gt(tbl_summary)
ggsave(
  here::here("output", "model", "plot_strat_coxph.svg"),
  plot_coxph,
  units = "cm", width = 40, height = 20
)






