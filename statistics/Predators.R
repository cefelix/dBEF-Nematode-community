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


#### 21a hurdle: Pr_per100g ~ sowndivLog*treatment + (1|block/plot), hu~1, fam=hurdle_lognormal ####

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
    
  
  
#### 31a hurdle: Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~1, fam=hurdle_lognormal, no 60 sp. plots ####
  SEED = 19111996
  dat = dBEF_nem21 %>% filter(sowndiv != 60)
  #standardize:  
  dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                        .after = sowndivLog)
  (dat$Pr_per100g !=0) %>% sum() #171
  
  
  m.Pr.hurdle31a <- brm(
    bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~  1),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) # all good
  
  pp_check(m.Pr.hurdle31a, ndraws = 100)
  
#### 31a_Pr3-5 hurdle: Pr3-5 ~ sowndivLogStd*treatment + (1|block/plot), hu~1, fam=hurdle_lognormal, no 60 sp. plots ####
  SEED = 19111996
  dat = dBEF_nem21 %>% filter(sowndiv != 60)
  dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                        .after = sowndivLog)
  (dat$Pr3!=0) %>% sum() #0
  (dat$Pr4!=0) %>% sum() #166
  (dat$Pr5!=0) %>% sum() #23
  
  m.Pr4.hurdle31a <- brm(
    bf(Pr4 ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~  1),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #2 divergent transitions
  
  m.Pr4.hurdle32a <- update(m.Pr4.hurdle31a,
                            control = list(adapt_delta=0.999)) #1 divergent transition

  m.Pr5.hurdle31a <- brm(
    bf(Pr5 ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~  1),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) 
  #1 divergent transition
  
  m.Pr5.hurdle32a <- update(m.Pr5.hurdle31a,
                            iter = 4000, warmup = 2000,
                            control = list(adapt_delta=0.9999))
  #3 divergent transitions
  

#### 31b hurdle: Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~term, fam=hurdle_lognormal, no 60 sp. plots ####
  dat = dBEF_nem21 %>% filter(sowndiv != 60)

  m.Pr.hurdle31b <- brm(
    bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
       hu ~  sowndivLogStd*treatment + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #all good
  
  pp_check(m.Pr.hurdle31b, ndraws = 100)
  
  loo(m.Pr.hurdle32a, m.Pr.hurdle31b) 

  
  #### 31b_Pr3-5 hurdle: Pr3-5 ~ sowndivLogStd*treatment + (1|block/plot), hu~1, fam=hurdle_lognormal, no 60 sp. plots ####
  #leaving out Pr3, as there are none:
  SEED = 19111996
  (dat$Pr3!=0) %>% sum()
  
  
  m.Pr4.hurdle31b <- brm(
    bf(Pr4 ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ sowndivLogStd*treatment + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #all good
  
  pp_check(m.Pr4.hurdle31b, ndraws=100)+
    xlim(0,200)
  
  m.Pr5.hurdle31b <- brm(
    bf(Pr5 ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ sowndivLogStd*treatment + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) 
  #2 divergent transitions, 1000 after warmup exceeded max_treedepth, bulk/tail ESS too low
  
  
#### 41a hurdle: Pr_per100g ~ sowndivLog*treatment + (1|block/plot), hu~1, fam=hurdle_lognormal, no 60 sp. plots ####
  dat = dBEF_nem21 %>% filter(sowndiv != 60)
  
  m.Pr.hurdle41a <- brm(
    bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
       hu ~  1),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99))
  
  pp_check(m.Pr.hurdle41a, ndraws = 100)
  
  
#### 41b hurdle: Pr_per100g ~ sowndivLog*treatment + (1|block/plot), hu~term, fam=hurdle_lognormal, no 60 sp. plots ####
  dat = dBEF_nem21 %>% filter(sowndiv != 60)
  
  m.Pr.hurdle41b <- brm(
    bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
       hu ~  sowndivLog*treatment + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99))
  
  pp_check(m.Pr.hurdle41b, ndraws = 100)
  
  
  
#model predictions:
  predictions <- conditional_effects(m.Pr.hurdle31b)[[3]]
  predictions$estimate__
  
  dat <- subset(dBEF_nem21, sowndiv != 60)
  #plot results
  ggplot(dat, aes(x=sowndivLog, y=Pr_per100g, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=estimate__, col=treatment))
  
  summary(m.Pr.hurdle31a)

  

  
  
#### save hurdle models ####
save(#m.Pr.hurdle11a,
     #m.Pr.hurdle21a,
     #m.Pr.hurdle21b,
     m.Pr.hurdle31a, #m.Pr.hurdle32a,
       m.Pr4.hurdle31a, m.Pr4.hurdle32a,
       m.Pr5.hurdle31a, m.Pr5.hurdle32a,
     m.Pr.hurdle31b,
      m.Pr4.hurdle31b,
      m.Pr5.hurdle31b,  
     #m.Pr.hurdle41a,
     #m.Pr.hurdle41b,
     file = "./statistics/brms/231127_Pr.RData")

load(file = "./statistics/brms/231122_Pr.RData")

