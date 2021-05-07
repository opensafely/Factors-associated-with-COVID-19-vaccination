######################################

# This script:
# - imports fitted stratified Cox model
# - reports cumulative baseline hazard for each strata

######################################


## Import libraries
library('tidyverse')
library('survival')


## Create output directory
dir.create(here::here("output", "models", "final"), showWarnings = FALSE, recursive=TRUE)


# function to get survival estimates for each strata over time
# basically like broom::tidy.coxph() but with a bit more stuff
tidy_surv <-
  function(
    survfit, # fitted survival model using coxph()
    times = NULL, # return estimates for these times only. If NULL then uses all event times
    addtimezero=FALSE # add estimates at time zero (ie where survival=1 and hazard=0)
  ) {
    
    # tidy post-fit survival dataset, with extra estimates than provided by broom::tidy.coxph
    
    strata <- names(survfit$strata)
    mintime <- min(survfit$time)
    timezero <- min(0, mintime-1)
    
    output0 <- broom::tidy(survfit)
    
    if (is.null(strata)){
      output0$strata = "1"
    }
    
    
    output <- output0 %>%
      group_by(strata) %>%
      transmute(
        strata,
        time,
        leadtime = lead(time),
        interval = leadtime - time,
        n.risk,
        n.event,
        n.censor,
        estimate,
        std.error,
        conf.low, conf.high,
        haz = -(estimate-lag(estimate, n=1L, default=1))/lag(estimate, n=1L, default=1),
        cml.haz = cumsum(haz),
        cml.haz2 = -log(estimate),
        haz2 = cml.haz2-lag(cml.haz2, n=1L, default=0),
        # log(-log()) scale
        llsurv = log(-log(estimate)),
        llsurv.low = log(-log(conf.low)),
        llsurv.high = log(-log(conf.high)),
      ) %>%
      ungroup()
    
    if(!is.null(times)){
      # this fills in (with constant interpolation) estimates for times where no events occurred
      # it does NOT re-estimate hazard over the smaller time period
      output <-
        output %>%
        group_by(strata) %>%
        complete(
          time = times,
          fill = list(n.event = 0, n.censor = 0)
        ) %>%
        fill(n.risk, .direction = c("up")) %>%
        fill(
          estimate, std.error, conf.low, conf.high,
          haz, cml.haz, haz2, cml.haz2,
          llsurv, llsurv.low, llsurv.high
        ) %>%
        ungroup()
    }
    
    if(addtimezero){
      
      time0row <- tibble(
        time = timezero,
        leadtime = mintime,
        interval = leadtime-time,
        
        estimate=1,
        conf.low=1,
        conf.high=1,
        std.error=0,
        
        haz=0,
        cml.haz=0,
      ) %>%
        right_join(
          tibble(strata=strata, time=timezero),
          by="time"
        )
      
      
      output <- bind_rows(
        output,
        time0row
      ) %>%
        arrange(strata, time) %>%
        fill(estimate, std.error, conf.low, conf.high, haz, cml.haz)
    }
    
    return(output)
  }



## Import processed data
data_tte <- read_rds(here::here("output", "data", "data_modelling.rds"))

## Converts logical to integer so that model coefficients print nicely in gtsummary methods
data_cox <- data_tte %>%
  mutate(
    across(
      where(is.logical),
      ~.x*1L
    )
  )


mod.strat.coxph.adj <- read_rds(here::here("output", "models", "final", "mod_strat_coxph_adj.rds"))


# get strata-specific estimates based on mean-centered covariates 
# mean-centered doesn't make much sense but strata-specific cumulative hazard is multiplicative wrt covariates,
# so relative differences between strata are invariant
fit <- survfit(mod.strat.coxph.adj, conf.type="log-log") 

# tidy output, with one estimate per time point
strata_estimates <- tidy_surv(fit, times=1:max(fit$time), addtimezero=TRUE)


plot_strata_cmlhaz <- ggplot(strata_estimates)+
  geom_step(aes(x=time, y=cml.haz, group=strata), alpha=0.2)+
  #geom_ribbon(aes(x=time, ymin=conf.low, ymax=conf.high, group=strata), colour="transparent", alpha=0.2)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


ggsave(
  here::here("output", "models", "final", "plot_strata_cmlhaz.svg"),
  plot_strata_cmlhaz,
  units = "cm", width = 20, height = 20
)

plot_strata_llsurv <- ggplot(strata_estimates)+
  geom_step(aes(x=time, y=llsurv, group=strata), alpha=0.2)+
  #geom_ribbon(aes(x=time, ymin=llsurv.low, ymax=llsurv.high, group=strata), colour="transparent", alpha=0.2)+
  scale_x_log10()+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  here::here("output", "models", "final", "plot_strata_llsurv.svg"),
  plot_strata_llsurv,
  units = "cm", width = 20, height = 20
)




plot_strata_cmlhaz28 <- ggplot(strata_estimates %>% filter(time==28))+
  geom_histogram(aes(x=cml.haz), alpha=0.2)+
  theme_bw()

ggsave(
  here::here("output", "models", "final", "plot_strata_cmlhaz28.svg"),
  plot_strata_cmlhaz28,
  units = "cm", width = 20, height = 20
)


plot_strata_cmlhaz56 <- ggplot(strata_estimates %>% filter(time==56))+
  geom_histogram(aes(x=cml.haz), alpha=0.2)+
  theme_bw()

ggsave(
  here::here("output", "models", "final", "plot_strata_cmlhaz56.svg"),
  plot_strata_cmlhaz56,
  units = "cm", width = 20, height = 20
)







Qdata <- function(data, groupvar, contvar){
  
  data %>%
    group_by({{groupvar}}) %>%
    summarise(
      n=n(),
      n_valid=sum(!is.na({{contvar}})),
      pct_valid = n_valid/n,
      Q10=quantile({{contvar}}, 0.1, na.rm=TRUE),
      Q25=quantile({{contvar}}, 0.25, na.rm=TRUE),
      Q50=quantile({{contvar}}, 0.5, na.rm=TRUE),
      Q75=quantile({{contvar}}, 0.75, na.rm=TRUE),
      Q90=quantile({{contvar}}, 0.9, na.rm=TRUE),
      #mean=mean({{contvar}}, na.rm=TRUE),
      #bsci = list(Hmisc::smean.cl.boot({{contvar}}, conf.int=0.95, B=1000, reps=FALSE)),
      #mean.ll = map_dbl(bsci, ~.[2]),
      #mean.ul = map_dbl(bsci, ~.[3])
    ) %>% 
    #select(-bsci) %>%
    ungroup() 
}



strata_estimates_snapshots <- strata_estimates %>% 
  filter(time %in% c(7,28,56)) %>%
  mutate(time_cat = fct_reorder(paste0("Day ",time), time))

plotfacethistdata <- function(data, catvar, contvar, filt=TRUE){
  #  function to get data for simple faceted histogram
  plotdata <-
    data %>% 
    filter(filt) %>%
    select(all_of(c(catvar, contvar))) %>%
    rename(variable=all_of(catvar), contvariable=all_of(contvar)) %>%
    mutate(
      variable_explicit_na=(fct_explicit_na(variable, na_level="(Missing)")),
    )
  
  Qdata(plotdata, variable_explicit_na, contvariable)
}

plotfacethist <- function(data, catvar, contvar, contname, subtitle=NULL, 
                          ywrapwidth=Inf, breakint=50, titlewrapwidth=40, ylim_upper=NULL){
  #  function to plot faceted histogram, taking data from plotfacethistdata
  
  contdata <- 
    data %>% 
    select(all_of(c(catvar, contvar))) %>%
    rename(variable=all_of(catvar), contvariable=all_of(contvar)) %>%
    mutate(
      variable_explicit_na=(fct_explicit_na(variable, na_level="(Missing)")),
    )
  
  plotdata <- plotfacethistdata(contdata, "variable_explicit_na", "contvariable")
  
    ggplot()+
    geom_histogram(
      data = contdata %>% filter(!is.na(contvariable)),
      aes(x=contvariable), colour="white", fill="darkgrey", size=0.5, binwidth=0.05, boundary = 0.5, closed = "left"
    ) +
    geom_point(data=plotdata, aes(y=-breakint/3, x=Q50), colour='darkgrey', size=2, alpha=0.5)+
    geom_linerange(data=plotdata, aes(y=-breakint/3, xmin=Q25, xmax=Q75), colour='darkgrey', size=1.5, alpha=0.5)+
    geom_linerange(data=plotdata, aes(y=-breakint/3, xmin=Q10, xmax=Q90), colour='darkgrey', size=0.5, alpha=0.5)+
    geom_hline(yintercept=0)+
    facet_wrap(vars(variable_explicit_na), ncol=1, strip.position="top")+#, space='free_y', scales="free_y")+
    scale_y_continuous(breaks = seq(0, ylim_upper, breakint), limits = c(-(breakint/1.8), NA))+
    labs(
      x=contname, y=NULL)+
    theme_bw(base_size = 10) +
    theme(
      panel.border = element_blank(),
      #axis.line.x = element_line(colour = "black"),
      panel.grid = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "lightgrey"),
      strip.background = element_blank(),
      
      plot.title = element_text(hjust = 0),
      plot.title.position = "plot",
      plot.caption.position =  "plot",
      plot.caption = element_text(hjust = 0, face = "italic"),
      strip.text = element_text(angle = 0, hjust = 1),
    )+
    NULL
  
}


CHsnapshots <- plotfacethist(strata_estimates_snapshots, "time_cat", "cml.haz", "Cumul. hazard", breakint=50, ylim_upper=1000)

ggsave(
  here::here("output", "models", "final", "plot_strata_snapshots.svg"),
  CHsnapshots,
  units = "cm", width = 20, height = 20
)


