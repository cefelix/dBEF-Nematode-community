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

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
dBEF_nemSH1 <- subset(dBEF_nem, SH == 1)
dBEF_nemSH5 <- subset(dBEF_nem, SH == 5)
dBEF_nemSH15 <- subset(dBEF_nem, SH == 15)
dBEF_nemSH19 <- subset(dBEF_nem, SH == 19)

dBEF_nem21 <- subset(dBEF_nem, year==2021)



#### 11a hurdle: Om_per100g ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_lognormal ####
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


m.Om.hurdle12a <- update(m.Om.hurdle11a,
                         iter = 3000, warmup = 1500,
                         control = list(adapt_delta = 0.9999))
 

pp_check(m.Om.hurdle12a, ndraws=100)+ 
  xlim(0,300)

#### 21a hurdle: Om_per100g ~ sowndivLog*treatment + (1|block/plot), family = hurdle_lognormal ####
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

m.Om.hurdle22b <- update(m.Om.hurdle51,
                         iter = 3000, warmup = 1500,
                         control = list(adapt_delta = 0.999))

pp_check(m.Om.hurdle22b, ndraws=100)


#### saving modles ####

save(m.Om.hurdle11a, m.Om.hurdle12a,
     m.Om.hurdle21a, m.Om.hurdle22a,
     m.Om.hurdle21b, m.Om.hurdle22b,
     file = "./statistics/brms/231117_Om.RData")

load(file = "./statistics/brms/231108_Om_hurdle.RData")

conditional_effects(m.Om.hurdle21)

