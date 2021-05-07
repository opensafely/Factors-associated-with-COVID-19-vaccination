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




plot_strata_cmlhaz30 <- ggplot(strata_estimates %>% filter(time==30))+
  geom_histogram(aes(x=cml.haz), alpha=0.2)+
  theme_bw()

ggsave(
  here::here("output", "models", "final", "plot_strata_cmlhaz30.svg"),
  plot_strata_cmlhaz30,
  units = "cm", width = 20, height = 20
)