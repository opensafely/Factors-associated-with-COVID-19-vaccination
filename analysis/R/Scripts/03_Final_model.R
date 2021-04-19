######################################

# This script:
# - imports processed data
# - fits a mixed effects cox model using the coxme package
# - saves model summaries (tables and figures)

######################################


# Preliminaries ----

## Import libraries
library('tidyverse')
library('lubridate')
library('survival')
library('coxme')
library('gtsummary')
library('gt')
library('survminer')
#library('ehahelper')

## Create output directory
dir.create(here::here("output", "models"), showWarnings = FALSE, recursive=TRUE)

## Function to plot stratified cox model
ggforest2 <-  function (model, data = NULL, main = "Hazard ratio", cpositions = c(0.02, 
                                                                                  0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2) 
{
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- survminer:::.get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(broom::tidy(model, conf.int = TRUE))
  gmodel <- broom::glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if(var %in% colnames(data)) {
      if (terms[i] %in% c("factor", "character")) {
        adf <- as.data.frame(table(data[, var]))
        cbind(var = var, adf, pos = 1:nrow(adf))
      }
      else if (terms[i] == "numeric") {
        data.frame(var = var, Var1 = "", Freq = nrow(data), 
                   pos = 1)
      }
      else {
        vars = grep(paste0("^", var, "*."), coef$term, 
                    value = TRUE)
        data.frame(var = vars, Var1 = "", Freq = nrow(data), 
                   pos = seq_along(vars))
      }
    } else {
      message(var, "is not found in data columns, and will be skipped.")
    }    
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds, ])[, c("var", "level", 
                                                "N", "p.value", "estimate", "conf.low", "conf.high", 
                                                "pos")]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 
                                                              4], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, 
                                       noDigits + 1), " ", ifelse(toShowExpClean$p.value < 0.05, 
                                                                  "*", ""), ifelse(toShowExpClean$p.value < 0.01, "*", 
                                                                                   ""), ifelse(toShowExpClean$p.value < 0.001, "*", ""))
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"], 
                              " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, 
  ]
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, 
                  na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(grid::convertX(unit(theme_get()$text$size, 
                                                             "pt"), "mm"))
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) + 
    geom_rect(aes(xmin = seq_along(var) - 0.5, xmax = seq_along(var) + 
                    0.5, ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]), 
                  fill = ordered(seq_along(var)%%2 + 1))) + scale_fill_manual(values = c("#FFFFFF33", 
                                                                                         "#00000033"), guide = "none") + geom_point(pch = 15, 
                                                                                                                                    size = 4) + geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), 
                                                                                                                                                              width = 0.15) + geom_hline(yintercept = 1, linetype = 3) + 
    coord_flip(ylim = exp(rangeplot)) + ggtitle(main) + scale_y_log10(name = "", 
                                                                      labels = sprintf("%g", breaks), expand = c(0.02, 0.02), 
                                                                      breaks = breaks) + theme_light() + theme(panel.grid.minor.y = element_blank(), 
                                                                                                               panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), 
                                                                                                               legend.position = "none", panel.border = element_blank(), 
                                                                                                               axis.title.y = element_blank(), axis.text.y = element_blank(), 
                                                                                                               axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
    xlab("") + annotate(geom = "text", x = x_annotate, y = exp(y_variable), 
                        label = toShowExpClean$var, fontface = "bold", hjust = 0, 
                        size = annot_size_mm) + annotate(geom = "text", x = x_annotate, 
                                                         y = exp(y_nlevel), hjust = 0, label = toShowExpClean$level, 
                                                         vjust = -0.1, size = annot_size_mm) + annotate(geom = "text", 
                                                                                                        x = x_annotate, y = exp(y_nlevel), label = toShowExpClean$N, 
                                                                                                        fontface = "italic", hjust = 0, vjust = ifelse(toShowExpClean$level == 
                                                                                                                                                         "", 0.5, 1.1), size = annot_size_mm) + annotate(geom = "text", 
                                                                                                                                                                                                         x = x_annotate, y = exp(y_cistring), label = toShowExpClean$estimate.1, 
                                                                                                                                                                                                         size = annot_size_mm, vjust = ifelse(toShowExpClean$estimate.1 == 
                                                                                                                                                                                                                                                "reference", 0.5, -0.1)) + annotate(geom = "text", 
                                                                                                                                                                                                                                                                                    x = x_annotate, y = exp(y_cistring), label = toShowExpClean$ci, 
                                                                                                                                                                                                                                                                                    size = annot_size_mm, vjust = 1.1, fontface = "italic") + 
    annotate(geom = "text", x = x_annotate, y = exp(y_stars), 
             label = toShowExpClean$stars, size = annot_size_mm, 
             hjust = -0.2, fontface = "italic") + annotate(geom = "text", 
                                                           x = 0.5, y = exp(y_variable), label = paste0("# Events: ", 
                                                                                                        gmodel$nevent, "; Global p-value (Log-Rank): ", format.pval(gmodel$p.value.log, 
                                                                                                                                                                    eps = ".001"), " \nAIC: ", round(gmodel$AIC, 
                                                                                                                                                                                                     2), "; Concordance Index: ", round(gmodel$concordance, 
                                                                                                                                                                                                                                        2)), size = annot_size_mm, hjust = 0, vjust = 1.2, 
                                                           fontface = "italic")
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}

## Import processed data
data_tte <- read_rds(here::here("output", "data", "data_all.rds"))

## Converts logical to integer so that model coefficients print nicely in gtsummary methods
data_cox <- data_tte %>%
  mutate(
    across(
      where(is.logical),
      ~.x*1L
    )
  )

## Exclude practices with less than 100 registered patients

### Registered patients counts
practice_counts <- data_cox %>% 
  group_by(practice_id) %>%
  summarise(`Number_of_registered_patients` = n())

## Exclude 
data_cox_stratification <- data_cox %>%
  filter(practice_id %in% subset(practice_counts, Number_of_registered_patients >= 100)$practice_id)



# MODELS ----

## Cox PH model - unadjusted
mod.coxph.unadj <- coxph(Surv(follow_up_time, covid_vax) ~ 1, data = data_cox_stratification)

write_rds(mod.coxph.unadj, here::here("output", "models", "mod_coxph_unadj.rds"), compress="gz")

## Cox PH model - adjusted; baseline demographics, comorbs, geographical, flu and shielding 
mod.coxph.adj <- coxph(Surv(follow_up_time, covid_vax) ~ 
                         ageband + sex + ethnicity + morbid_obesity +
                         chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                         chronic_kidney_disease_all_stages_1_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                         asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                         immunosuppression_medication + imd + stp + region + rural_urban + flu_vaccine + shielded +
                         shielded_since_feb_15,
                       data = data_cox_stratification)

write_rds(mod.coxph.adj, here::here("output", "models", "mod_coxph_adj.rds"), compress="gz")

# Cox model - adjusted; baseline demographics, comorbs, geographical, flu, shielding & practice id
mod.coxph.adjb <- coxph(Surv(follow_up_time, covid_vax) ~
                          ageband + sex + ethnicity + morbid_obesity +
                          chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                          chronic_kidney_disease_all_stages_1_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                          asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                          immunosuppression_medication + imd + region + rural_urban + flu_vaccine + shielded +
                          shielded_since_feb_15 + strata(practice_id),
                        data = data_cox_stratification)

write_rds(mod.coxph.adjb, here::here("output", "models", "mod_coxph_adjb.rds"), compress="gz")

# Cox model with RE for practice - adjusted; baseline demographics, comorbs, geographical, flu, shielding & practice id
mod.coxph.adjc <- coxph(Surv(follow_up_time, covid_vax) ~
                          ageband + sex + ethnicity + morbid_obesity +
                          chronic_heart_disease + diabetes + chronic_kidney_disease_diagnostic + chronic_kidney_disease_all_stages +
                          chronic_kidney_disease_all_stages_1_5 + sev_mental_ill + learning_disability + chronic_neuro_dis_inc_sig_learn_dis +
                          asplenia + chronic_liver_disease + chronis_respiratory_disease + immunosuppression_diagnosis +
                          immunosuppression_medication + imd + region + rural_urban + flu_vaccine + shielded +
                          shielded_since_feb_15 + frailty(practice_id),
                        data = data_cox_stratification)

write_rds(mod.coxph.adjc, here::here("output", "models", "mod_coxph_adjc.rds"), compress="gz")



# Output model coefficients ----

## Summary tables

### Cox PH model - adjusted
tab_mod1 <- gtsummary::tbl_regression(mod.coxph.adj, exp = TRUE)
gtsave(tab_mod1 %>% as_gt(), here::here("output", "models", "tab_coxph.html"))
write_csv(tab_mod1$table_body, here::here("output",  "models", "tab_coxph.csv"))

### Stratified Cox PH model - adjusted
tab_mod2 <- gtsummary::tbl_regression(mod.coxph.adjb, exp = TRUE)
gtsave(tab_mod2 %>% as_gt(), here::here("output", "models", "tab_coxphb.html"))
write_csv(tab_mod2$table_body, here::here("output",  "models", "tab_coxphb.csv"))

### Cox PH model with REs- adjusted
tab_mod3 <- gtsummary::tbl_regression(mod.coxph.adjc, exp = TRUE)
gtsave(tab_mod3 %>% as_gt(), here::here("output", "models", "tab_coxphc.html"))
write_csv(tab_mod3$table_body, here::here("output",  "models", "tab_coxphc.csv"))


# ## Mixed effects Cox model - adjusted
# tab_mod2 <- gtsummary::tbl_regression(mod.coxme.adj, exp = TRUE)
# gtsave(tab_mod2 %>% as_gt(), here::here("output", "models", "tab_coxme.html"))
# write_csv(tab_mod2$table_body, here::here("output",  "models", "tab_coxme.csv"))

## Forest plots

### Cox PH model - adjusted
plot_coxph <- ggforest(mod.coxph.adj, data = data_cox_stratification)
ggsave(
  here::here("output", "models", "plot_coxph.svg"),
  plot_coxph,
  units = "cm", width = 20, height = 30
)

### Stratified Cox PH model - adjusted
plot_coxph <- ggforest2(mod.coxph.adjb, data = data_cox_stratification)
ggsave(
  here::here("output", "models", "plot_coxphb.svg"),
  plot_coxph,
  units = "cm", width = 20, height = 30
)

### Cox PH model with REs - adjusted
# plot_coxph <- ggforest2(mod.coxph.adjc, data = data_cox)
# ggsave(
#   here::here("output", "models", "plot_coxphc.svg"),
#   plot_coxph,
#   units = "cm", width = 20, height = 30
# )

### Mixed effects Cox model - adjusted
# plot_coxme <- survminer::ggforest(mod.coxme.adj, data = data_cox)
# ggsave(
#   here::here("output", "models", "plot_coxph.svg"),
#   plot_coxph,
#   units = "cm", width = 20, height = 25
# )




