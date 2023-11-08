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

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
dBEF_nemSH1 <- subset(dBEF_nem, SH == 1)
dBEF_nemSH5 <- subset(dBEF_nem, SH == 5)
dBEF_nemSH15 <- subset(dBEF_nem, SH == 15)
dBEF_nemSH19 <- subset(dBEF_nem, SH == 19)

dBEF_nem21 <- subset(dBEF_nem, year==2021)

#### 11 hurdle: Om_per100gLog ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Om.hurdle11 <- brm(
  bf(Om_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

m.Om.hurdle12 <- update(m.Om.hurdle11, 
                        seed = SEED,
                        control = list(adapt_delta = 0.999))

pp_check(m.Om.hurdle32, ndraws=100)

#### 21 hurdle: Om_per100gLog ~ sowndivLog*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Om.hurdle21 <- brm(
  bf(Om_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Om.hurdle21, ndraws=100)


#### 31 hurdle: Om_per100gLog ~ sowndiv*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Om.hurdle31 <- brm(
  bf(Om_per100gLog.hurdle ~ sowndiv*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

m.Om.hurdle32 <- update(m.Om.hurdle31, 
                        seed = SEED,
                        control = list(adapt_delta = 0.999))

pp_check(m.Om.hurdle32, ndraws=100)


#### 41 hurdle: Om_per100g ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_lognormal ####
m.Om.hurdle41 <- brm(
  bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

m.Om.hurdle42 <- update(m.Om.hurdle41,
                        iter = 3000, warmup = 1500,
                        control = list(adapt_delta = 0.999))
  #2 divergent transitions
m.Om.hurdle43 <- update(m.Om.hurdle42,
                        control = list(adapt_delta = 0.9999))
  #1 divergent transition
m.Om.hurdle44 <- update(m.Om.hurdle42,
                        control = list(adapt_delta = 0.99999))
  #2 divergent transitions, screw it

pp_check(m.Om.hurdle43, ndraws=100)+ 
  xlim(0,300)

#### 51 hurdle: Om_per100g ~ sowndivLog*treatment + (1|block/plot), family = hurdle_lognormal ####
m.Om.hurdle51 <- brm(
  bf(Om_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

m.Om.hurdle52 <- update(m.Om.hurdle51,
                        iter = 3000, warmup = 1500,
                        control = list(adapt_delta = 0.999))

pp_check(m.Om.hurdle51, ndraws=100)

#### saving the hurdle modles ####

save(m.Om.hurdle12,
     m.Om.hurdle21,
     m.Om.hurdle31, m.Om.hurdle32,
     m.Om.hurdle41, m.Om.hurdle42,
     m.Om.hurdle43,m.Om.hurdle44,
     file = "./statistics/brms/231108_Om_hurdle.RData")

load(file = "./statistics/brms/231105_Om_hurdle.RData")

