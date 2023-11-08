library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)
SEED=22061996

#### 11 hurdle: Ba_per100gLog ~ sowndivLog*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Ba17.hurdle11 <- brm(
  bf(Ba_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Ba17.hurdle11, ndraws=100)

#### 21 hurdle: Ba_per100g ~ sowndivLog*treatment + (1|block/plot), family = hurdle_lognormal ####
m.Ba.hurdle21 <- brm(
  bf(Ba_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Ba.hurdle21, ndraws=100)+
  xlim(0,1000)


#### 31 hurdle: Ba_per100gLog ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Ba.hurdle31 <- brm(
  bf(Ba_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Ba.hurdle31, ndraws=100)

#### 41 hurdle: Ba_per100gLog ~ sowndiv*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Ba.hurdle41 <- brm(
  bf(Ba_per100gLog.hurdle ~ sowndiv*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Ba.hurdle41, ndraws=100)

####51 Ba_per100gLog ~ sowndivLog*treatment + (1|block/plot), family = gaussian ####
m.Ba17.51 <- brm(
  Ba_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot),
  data = dBEF_nem17, 
  family = gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

m.Ba17.52 <- update(m.Ba17.51, 
                    control = list(adapt_delta=0.999))

pp_check(m.Ba17.52, ndraws=100)

#### save hurdle models ####
save(m.Ba17.hurdle11, #log(Ba) ~ sowndivLog, fam=hurdle_gaussian
     #m.Ba.hurdle21, #Ba  ~ sowndivLog, fam=hurdle_lognormal
     #m.Ba.hurdle31, #log(Ba) ~ sowndivLogStd, fam=hurdle_gaussian
     #m.Ba.hurdle41, #log(Ba) ~ sowndiv, fam=hurdle_gaussian
     m.Ba17.51, m.Ba17.52,
     file="./statistics/17_remodelling/brms/231108_Ba17_hurdle.RData")

load(file="./statistics/17_remodelling/brms/231107_Ba17_hurdle.RData")
conditional_effects(m.Ba17.52)
