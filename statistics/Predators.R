####data and packages####
library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(hexbin)
library(GGally)#
library(emmeans)

# a seed:
SEED = 19111996



#### 11a hurdle: Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot), fam=hurdle_lognormal ####

m.Pr.hurdle11a <- brm(
  bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) # 1 divergent transition


pp_check(m.Pr.hurdle11a, ndraws=100)+
  xlim(0,300)


#### 21a hurdle: Pr_per100g ~ sowndivLog*treatment + (1|block/plot), fam=hurdle_lognormal ####

m.Pr.hurdle21a <- brm(
  bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle21a, ndraws=100)+
  xlim(0,300)

####21b hurdle: Pr_per100g ~ sowndivLog*treatment + (1|block/plot), hu~term, fam=hurdle_lognormal ####

m.Pr.hurdle21b <- brm(
  bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~  sowndivLog*treatment + (1|block/plot)),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle21b, ndraws=100)+
  xlim(0,300) 

#model predictions:
  predictions <- conditional_effects(m.Pr.hurdle21b)[[3]]
  predictions$estimate__

#plot results
  ggplot(dBEF_nem21, aes(x=sowndivLog, y=Pr_per100g, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=estimate__, col=treatment))
  
  summary(m.Pr.hurdle21a)
  
  loo(m.Pr.hurdle11a, m.Pr.hurdle21a, m.Pr.hurdle21b)
  #looic 11a: 1867.8 // 21a: 1867.6 // 21b: 1859.0
    
  
  
#### 31a hurdle: Pr_per100g ~ sowndivLog*treatment + (1|block/plot), fam=hurdle_lognormal, no 60 sp. plots ####
  dat = dBEF_nem21 %>% filter(sowndiv != 60)
  
  m.Pr.hurdle31a <- brm(
    bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
       hu ~  1),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99))
  
  pp_check(m.Pr.hurdle31a, ndraws = 100)
  

#### 31b hurdle: Pr_per100g ~ sowndivLog*treatment + (1|block/plot), hu~term, fam=hurdle_lognormal, no 60 sp. plots ####
  dat = dBEF_nem21 %>% filter(sowndiv != 60)

  m.Pr.hurdle31b <- brm(
    bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
       hu ~  sowndivLog*treatment + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99))
  
  pp_check(m.Pr.hurdle31b, ndraws = 100)
  
  loo(m.Pr.hurdle31a, m.Pr.hurdle31b) #looic 31a: 1769.9 // 31b: 1760.4
    #31a has one bad pareto k value!
  
  
  
#model predictions:
  predictions <- conditional_effects(m.Pr.hurdle31a)[[3]]
  predictions$estimate__
  
  #plot results
  ggplot(dat, aes(x=sowndivLog, y=Pr_per100g, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=estimate__, col=treatment))
  
  summary(m.Pr.hurdle31a)
  
  
#### save hurdle models ####
save(m.Pr.hurdle11a,
     m.Pr.hurdle21a,
     m.Pr.hurdle21b,
     m.Pr.hurdle31a,
     m.Pr.hurdle31b,
     file = "./statistics/brms/231117_Pr.RData")

load(file = "./statistics/brms/231108_Pr_hurdle.RData")

