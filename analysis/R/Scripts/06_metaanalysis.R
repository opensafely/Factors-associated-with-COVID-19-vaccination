######################################

# This script:
# - imports outputted model summaries from EMIS and TPP EHR backends
# - combined the model coefficients using inverse-variance-weighted averages
# - saves and plots outputs
# IMPORTANT: this script should only be run OFFLINE

######################################



# Preliminaries ----

## Import libraries
library('here')
library('tidyverse')
library('lubridate')
library('survival')
library('fs')

## Create output directory
dir_create(here("released_outputs", "combined"))

# Import models
dummy_data <- TRUE

if(dummy_data){
  tidy_tpp <- read_csv(fs::path(here("output", "model"), "tidy_tpp.csv"))
  tidy_emis <- read_csv(fs::path(here("output", "model"), "tidy_emis.csv"))
} else{
  tidy_tpp <- read_csv(fs::path(here("released_output", "tpp", "model"), "tidy_tpp.csv"))
  tidy_emis <- read_csv(fs::path(here("released_output", "emis", "model"), "tidy_emis.csv"))
}


# Combine estimates ----
tidy_stack <- 
  bind_rows(
    tidy_tpp, 
    tidy_emis
  ) 

tidy_combined <- tidy_stack %>%
  group_by(term, variable, var_label, var_class, var_type, var_nlevels, contrasts, contrasts_type, reference_row, label) %>% 
  summarise(
    n_obs = sum(n_obs),
    n_event = sum(n_event),
    exposure = sum(exposure),
    estimate = weighted.mean(estimate, std.error^-2),
    std.error = sqrt(1/sum(std.error^-2)),
    statistic = estimate/std.error,
    p.value = 2 * pmin(pnorm(statistic), pnorm(-statistic)),
    conf.low = estimate + qnorm(0.025)*std.error,
    conf.high = estimate + qnorm(0.975)*std.error,
    
  ) %>% ungroup()

write_csv(tidy_combined, here("released_outputs", "combined", "meta_estimates.csv"))


# Results ----

## Table
tab_mod1  <- tidy_combined  %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Variable") %>%
  mutate(LCI = round(exp(conf.low), digits = 2),
         UCI = round(exp(conf.high), digits = 2),
         HR = round(exp(estimate), digits = 2),
         `95% CI` = paste(" (", LCI, " - ", UCI, ")", sep = ""),
         `p-value` = round(`p.value`, digits = 4)) %>%
  select(Variable = term, HR, `95% CI`, `p-value`)

write_csv(tab_mod1, here::here("released_outputs",  "combined", "tab_strat_coxph.csv"))

## Forest plot
plot_data <- tidy_combined %>%
  mutate(order = factor(variable, levels = c("ageband", "sex", "ethnicity", "imd", "immunosuppression", "ckd",
                                             "chronic_respiratory_disease", "diabetes", "chronic_liver_disease",
                                             "chronic_neuro_dis_inc_sig_learn_dis", "chronic_heart_disease", "asplenia",
                                             "sev_mental_ill", "morbid_obesity"))) %>%
  filter(!is.na(term)) %>%
  arrange(order) %>%
  mutate(
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

var_lookup <- c("Age Band", "Age Band", "Age Band", "Age Band", "Age Band", "Age Band", 
                "Sex", "Sex", 
                "Ethnicity", "Ethnicity", "Ethnicity", "Ethnicity", "Ethnicity", "Ethnicity", "Ethnicity", 
                "Ethnicity", "Ethnicity", "Ethnicity", "Ethnicity", "Ethnicity", "Ethnicity", "Ethnicity",
                "Ethnicity", "Ethnicity", "Ethnicity", "Ethnicity", "Ethnicity",
                "IMD", "IMD", "IMD",  "IMD", "IMD", "Clinical Risk Groups", "", "", "", "", "", "", "", "",      
                "Other Groups")
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
plot_coxph <- ggplot(plot_data) +
  geom_point(aes(x=estimate, y=level)) +
  geom_linerange(aes(xmin=conf.low, xmax=conf.high, y=level)) +
  geom_vline(aes(xintercept=1), colour='black', alpha=0.8)+
  facet_grid(rows=vars(variable), scales="free_y", switch="y", space="free_y", labeller = labeller(variable = var_lookup))+
  scale_x_log10()+
  scale_y_discrete(breaks=level_lookup, labels=names(level_lookup)) +
  geom_rect(aes(alpha = variable_card), xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf, fill='grey', colour="transparent") +
  scale_alpha_continuous(range=c(0,0.3), guide="none")+
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

ggsave(
  here::here("released_outputs", "combined", "plot_strat_coxph.svg"),
  plot_coxph,
  units = "cm", width = 40, height = 20
)




