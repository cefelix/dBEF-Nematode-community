####data and packages####
library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(hexbin)
library(GGally)

# a seed:
SEED = 22061996

####Hill_q0 11: Hill_q0 ~ sowndivLog*treatment + (1|block/plot), fam=negBinom####
dBEF_nem21 %>% filter(is.na(Hill_q0)==FALSE) %>% pull(Hill_q0) %>% density %>% plot()
#as it is a count, lets try negbinomial:
SEED = 22061996
m.Hill_q0.11 <- brm(Hill_q0 ~ sowndivLog*treatment + (1|block/plot), 
                            data = dBEF_nem21, family = "negbinomial",
                            chains = 3,
                            cores = 3,
                            iter = 2000, warmup = 1000,
                            seed = SEED,
                            control = list(adapt_delta = 0.9) ) 

m.Hill_q0.12 <- update(m.Hill_q0.11,
                       control = list(adapt_delta = 0.99))

pp_check(m.Hill_q0.12, ndraws=100) # thats off

####H1 Shannon ####
dBEF_nem21 %>% filter(is.na(Shannon_H)==FALSE) %>% pull(Shannon_H) %>% density %>% plot()
SEED = 22061996
m.Shannon_H.11 <- brm(Shannon_H ~ sowndivLog*treatment + (1|block/plot), 
                    data = dBEF_nem21, family = "gaussian",
                    chains = 3,
                    cores = 3,
                    iter = 2000, warmup = 1000,
                    seed = SEED,
                    control = list(adapt_delta = 0.9) ) 
  #10 divergent transitions

m.Shannon_H.12 <- update(m.Shannon_H.11,
                         control=list(adapt_delta=0.99))
  #7 divergent transitions
m.Shannon_H.13 <- update(m.Shannon_H.12,
                         control=list(adapt_delta=0.999))

pp_check(m.Shannon_H.13, ndraws=100) # a little overpredictive in quartile below mean

#### hill q1 ####
dBEF_nem21 %>% filter(is.na(Hill_q1)==FALSE) %>% pull(Hill_q1) %>% density %>% plot()
SEED = 22061996
m.Hill_q1.11 <- brm(Hill_q1 ~ sowndivLog*treatment + (1|block/plot), 
                      data = dBEF_nem21, family = "gaussian",
                      chains = 3,
                      cores = 3,
                      iter = 2000, warmup = 1000,
                      seed = SEED,
                      control = list(adapt_delta = 0.9) ) 

m.Hill_q1.12 <- update(m.Hill_q1.11,
                       control = list(adapt_delta = 0.99) )

pp_check(m.Hill_q1.12, ndraws=100)
mcmc_plot(m.Hill_q1.12 )

#### Hill q1 without 60 species ####
SEED = 22061996
dat <- subset(dBEF_nem21, sowndiv != 60)
dat %>% filter(is.na(Hill_q1.Fu)==FALSE) %>% pull(Hill_q1.Fu) %>% density %>% plot()

  m.41.Hill_q1 <- brm(Hill_q1 ~ sowndivLog*treatment + (1|block/plot), 
                      data = dat, family = "gaussian",
                      chains = 3,
                      cores = 3,
                      iter = 2000, warmup = 1000,
                      seed = SEED,
                      control = list(adapt_delta = 0.99) ) 
  
 
  
  pp_check(m.41.Hill_q1, ndraws=100)
  
  m.41b.Hill_q1 <- brm(Hill_q1 ~ sowndivLog*treatment + (1|block/plot), 
                      data = dat, family = "gamma",
                      chains = 3,
                      cores = 3,
                      iter = 2000, warmup = 1000,
                      seed = SEED,
                      control = list(adapt_delta = 0.99) ) 
  pp_check(m.41b.Hill_q1, ndraws=100)
  
#for Fu  
  dat %>% filter(is.na(Hill_q1.Fu)==FALSE) %>% pull(Hill_q1.Fu) %>% density %>% plot()
  m.41.Hill_q1.Fu <- brm(Hill_q1.Fu ~ sowndivLog*treatment + (1|block/plot), 
                      data = dat, family = "gaussian",
                      chains = 3,
                      cores = 3,
                      iter = 2000, warmup = 1000,
                      seed = SEED,
                      control = list(adapt_delta = 0.99) )
  pp_check(m.41.Hill_q1.Fu, ndraws=100) #thats off
  
  m.41b.Hill_q1.Fu <- brm(Hill_q1.Fu ~ sowndivLog*treatment + (1|block/plot), 
                         data = dat, family = "gamma",
                         chains = 3,
                         cores = 3,
                         iter = 2000, warmup = 1000,
                         seed = SEED,
                         control = list(adapt_delta = 0.99) ) #all good
  pp_check(m.41b.Hill_q1.Fu, ndraws=100) #better than the gaussian
  conditional_effects(m.41b.Hill_q1.Fu)
  
  loo(m.41.Hill_q1.Fu, m.41b.Hill_q1.Fu)
  
  #for Ba 
  dat %>% filter(is.na(Hill_q1.Ba)==FALSE) %>% pull(Hill_q1.Ba) %>% density %>% plot()
  m.41.Hill_q1.Ba <- brm(Hill_q1.Ba ~ sowndivLog*treatment + (1|block/plot), 
                         data = dat, family = "gaussian",
                         chains = 3,
                         cores = 3,
                         iter = 2000, warmup = 1000,
                         seed = SEED,
                         control = list(adapt_delta = 0.99) )
  pp_check(m.41.Hill_q1.Ba, ndraws=100) #all good
  
  m.41b.Hill_q1.Ba <- brm(Hill_q1.Ba ~ sowndivLog*treatment + (1|block/plot), 
                         data = dat, family = "gamma",
                         chains = 3,
                         cores = 3,
                         iter = 2000, warmup = 1000,
                         seed = SEED,
                         control = list(adapt_delta = 0.99) )
  pp_check(m.41b.Hill_q1.Ba, ndraws=100) #looks better than gaussian
  loo(m.41.Hill_q1.Ba, m.41b.Hill_q1.Ba)
  conditional_effects(m.41b.Hill_q1.Ba)
  
  #for Pl  
  dat %>% filter(is.na(Hill_q1.Pl)==FALSE) %>% pull(Hill_q1.Pl) %>% density %>% plot()
  m.41.Hill_q1.Pl <- brm(Hill_q1.Pl ~ sowndivLog*treatment + (1|block/plot), 
                         data = dat, family = "gaussian",
                         chains = 3,
                         cores = 3,
                         iter = 2000, warmup = 1000,
                         seed = SEED,
                         control = list(adapt_delta = 0.99) ) 
  pp_check(m.41.Hill_q1.Pl, ndraws=100) #all good
  
  m.41b.Hill_q1.Pl <- brm(Hill_q1.Pl ~ sowndivLog*treatment + (1|block/plot), 
                         data = dat, family = "gamma",
                         chains = 3,
                         cores = 3,
                         iter = 2000, warmup = 1000,
                         seed = SEED,
                         control = list(adapt_delta = 0.99) ) #all good
  pp_check(m.41b.Hill_q1.Pl, ndraws=100) #okayish
  conditional_effects( m.41b.Hill_q1.Pl)
  
  #for Pr  
  dat %>% filter(is.na(Hill_q1.Pr)==FALSE) %>% pull(Hill_q1.Pr) %>% density %>% plot()
  m.41.Hill_q1.Pr <- brm(Hill_q1.Pr ~ sowndivLog*treatment + (1|block/plot), 
                         data = dat, family = "gaussian",
                         chains = 3,
                         cores = 3,
                         iter = 2000, warmup = 1000,
                         seed = SEED,
                         control = list(adapt_delta = 0.99) ) #2 divergent transitions
  m.42.Hill_q1.Pr <- update(m.41.Hill_q1.Pr, control = list(adapt_delta = 0.999)  ) #all good
  pp_check(m.42.Hill_q1.Pr, ndraws=100)
  
  m.41b.Hill_q1.Pr <- brm(Hill_q1.Pr ~ sowndivLog*treatment + (1|block/plot), 
                         data = dat, family = "gamma",
                         chains = 3,
                         cores = 3,
                         iter = 2000, warmup = 1000,
                         seed = SEED,
                         control = list(adapt_delta = 0.99) ) #all good
  pp_check(m.41b.Hill_q1.Pr, ndraws=100) #better than gaussian
  conditional_effects( m.41b.Hill_q1.Pr)
  
  #for Om  
  dat %>% filter(is.na(Hill_q1.Om)==FALSE) %>% pull(Hill_q1.Om) %>% density %>% plot()
  sum(dat$Hill_q1.Om == TRUE )
  m.41.Hill_q1.Om <- brm(Hill_q1.Om ~ sowndivLog*treatment + (1|block/plot), 
                         data = dat, family = "gaussian",
                         chains = 3,
                         cores = 3,
                         iter = 2000, warmup = 1000,
                         seed = SEED,
                         control = list(adapt_delta = 0.99) ) #1 divergent transition
  
  m.42.Hill_q1.Om <- update(m.41.Hill_q1.Om ,
                            control = list(adapt_delta = 0.999))

  pp_check(m.42.Hill_q1.Om, ndraws = 100) #completely off
  
  m.41b.Hill_q1.Om <- brm(Hill_q1.Om ~ sowndivLog*treatment + (1|block/plot), 
                         data = dat, family = "gamma",
                         chains = 3,
                         cores = 3,
                         iter = 2000, warmup = 1000,
                         seed = SEED,
                         control = list(adapt_delta = 0.99) )
  pp_check( m.41b.Hill_q1.Om, ndraws=100) #still completely off
  conditional_effects( m.41b.Hill_q1.Om)

####save ####
save( m.41.Hill_q1,
        #m.41b.Hill_q1, 
      m.41.Hill_q1.Fu,
        m.41b.Hill_q1.Fu,
      m.41.Hill_q1.Ba,
        m.41b.Hill_q1.Ba,
      m.41b.Hill_q1.Pl,
        m.41.Hill_q1.Pl,
      m.41.Hill_q1.Pr, m.42.Hill_q1.Pr,
        m.41b.Hill_q1.Pr,
      m.41.Hill_q1.Om, m.42.Hill_q1.Om,
        m.41b.Hill_q1.Om,
     file = "./statistics/brms/231127_Hill.RData" )

load(file = "./statistics/brms/231124_Hill.RData" )
conditional_effects(m.41.Hill_q1 , prob=0.89)
