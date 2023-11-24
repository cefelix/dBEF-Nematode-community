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
SEED = 19111996



#### 11a hurdle: Om_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~1, family = hurdle_lognormal ####
SEED = 19111996

m.Om.hurdle11a <- brm(
  bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
  #all fine

pp_check(m.Om.hurdle12a, ndraws=100)+ 
  xlim(0,300)

#### 11b hurdle: Om_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~term, family = hurdle_lognormal ####
SEED = 19111996

m.Om.hurdle11b <- brm(
  bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
  #5 divergent transitions, tail EES too low

m.Om.hurdle12b <- update(m.Om.hurdle11b,
                         iter = 3000, warmup = 1500,
                         control = list(adapt_delta = 0.999)) #all fine



#### 21a hurdle: Om_per100g ~ sowndivLog*treatment + (1|block/plot), hu~1, family = hurdle_lognormal ####
SEED = 19111996

m.Om.hurdle21a <- brm( #former 51
  bf(Om_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

m.Om.hurdle22a <- update(m.Om.hurdle51,
                        iter = 3000, warmup = 1500,
                        control = list(adapt_delta = 0.999))

pp_check(m.Om.hurdle22a, ndraws=100)

#### 21b hurdle: Om_per100g ~ sowndivLog*treatment + (1|block/plot), hu~term, family = hurdle_lognormal ####
SEED = 19111996

m.Om.hurdle21b <- brm( #former 51
  bf(Om_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ sowndivLog*treatment + (1|block/plot)),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
  #39 divergent transitions, bulk EES too low

m.Om.hurdle22b <- update(m.Om.hurdle51,
                         iter = 3000, warmup = 1500,
                         control = list(adapt_delta = 0.9999))

pp_check(m.Om.hurdle22b, ndraws=100)

#### 31a hurdle: Om_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~1, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 19111996
dat <- subset(dBEF_nem21, sowndiv != 60)

m.Om.hurdle31a <- brm(
  bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
  #all fine


pp_check(m.Om.hurdle31a, ndraws=100)+ 
  xlim(0,300)


####31a_Om3-5 hurdle: Om2-5 ~ sowndivLogStd*treatment + (1|block/plot), hu~1, family = hurdle_lognormal, no 60 sp plots ####
SEED = 22061996
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60)

sum(dat$Om3!=0) #0
sum(dat$Om4!=0) #97
sum(dat$Om5!=0) #1

m.Om4.hurdle31a <- brm(
  bf(Om4 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #4 divergent transitions

m.Om4.hurdle32a <- update(m.Om4.hurdle31a,
                          control=list(adapt_delta=0.99999))
                          #less adapt_delta doesnt do it

m.Om5.hurdle31a <- brm(
  bf(Om5 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #don't do it
  #2999 transitions exceeded max_treedepth, max Rhat=2.78, bulk/tail ESS too low

#### 31b hurdle: Om_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~term, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 19111996
dat <- subset(dBEF_nem21, sowndiv != 60)


m.Om.hurdle31b <- brm(
  bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
  #1 divergent transition 

m.Om.hurdle32b <- update(m.Om.hurdle31b,
                         #iter = 3000, warmup = 1500,
                         control = list(adapt_delta = 0.999)) #all fine



####31b_Om2-5 hurdle: Om2-5 ~ sowndivLogStd*treatment + (1|block/plot), hu~term, family = hurdle_lognormal, no 60 sp plots ####
SEED = 22061996
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60)

sum(dat$Om3!=0) #0
sum(dat$Om4!=0) #97
sum(dat$Om5!=0) #1

m.Om4.hurdle31b <- brm(
  bf(Om4 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good

m.Om5.hurdle31b <- brm(
  bf(Om5 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #pointless, 1 sample !=0

#### 41a hurdle: Om_per100g ~ sowndivLog*treatment + (1|block/plot), hu~1, fam=hurdle_lognormal, no 60 sp. plots ####
dat = dBEF_nem21 %>% filter(sowndiv != 60)
SEED=19111996

m.Om.hurdle41a <- brm(
  bf(Om_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~  1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #4 divergent transitions

m.Om.hurdle42a <- update(m.Om.hurdle41a,
                         control=list(adapt_delta=0.9999))

pp_check(m.Om.hurdle41a, ndraws = 100)

#### 41b hurdle: Om_per100g ~ sowndivLog*treatment + (1|block/plot), hu~term, fam=hurdle_lognormal, no 60 sp. plots ####
dat = dBEF_nem21 %>% filter(sowndiv != 60)

m.Om.hurdle41b <- brm(
  bf(Om_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~  sowndivLog*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) # 1 divergent transitions

m.Om.hurdle42b <- update(m.Om.hurdle41b,
                         control=list(adapt_delta=0.9999)) 

pp_check(m.Om.hurdle41b, ndraws = 100)


#### saving modles ####

save(m.Om.hurdle11a, 
     m.Om.hurdle11b, m.Om.hurdle12b,
     m.Om.hurdle21a, 
     m.Om.hurdle21b, 
     m.Om.hurdle31a,
       m.Om4.hurdle31a, m.Om4.hurdle32a, 
       m.Om5.hurdle31a, #pointless as only 1 samples !=0
     m.Om.hurdle31b, m.Om.hurdle32b,
       m.Om4.hurdle31b,
       m.Om5.hurdle31b, #pointless as only 1 samples !=0
     
     m.Om.hurdle41a, m.Om.hurdle42a,
     m.Om.hurdle41b, m.Om.hurdle42b,
     
     file = "./statistics/brms/231124_Om.RData")

load(file = "./statistics/brms/231122_Om.RData")

conditional_effects(m.Om.hurdle21)

