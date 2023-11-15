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
#use full history as reference!


hurdle_gamma()


# a seed:
SEED = 22061996

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
dBEF_nemSH1 <- subset(dBEF_nem, SH == 1)
dBEF_nemSH5 <- subset(dBEF_nem, SH == 5)
dBEF_nemSH15 <- subset(dBEF_nem, SH == 15)
dBEF_nemSH19 <- subset(dBEF_nem, SH == 19)

dBEF_nem21 <- subset(dBEF_nem, year==2021)

#### hurdle: Om_per100gLog ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Pr.hurdle11 <- brm(
  bf(Pr_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle11, ndraws=100)


#### hurdle: Pr_per100gLog ~ sowndiv*treatment + (1|block/plot), fam=hurdle_gaussian ####
SEED = 19111996
m.Pr.hurdle21 <- brm(
  bf(Pr_per100gLog.hurdle ~ sowndiv*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

m.Pr.hurdle22 <- update(m.Pr.hurdle21,
                        control = list(adapt_delta = 0.999))

pp_check(m.Pr.hurdle22, ndraws=100)

#### hurdle: Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot), fam=hurdle_lognormal ####
SEED = 19111996
m.Pr.hurdle31 <- brm(
  bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle31, ndraws=100)+
  xlim(0,300)

#### hurdle: Pr_per100g ~ sowndiv*treatment + (1|block/plot), fam=hurdle_lognormal ####
SEED = 19111996
m.Pr.hurdle41 <- brm(
  bf(Pr_per100g ~ sowndiv*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

m.Pr.hurdle42 <- update(m.Pr.hurdle41,
  control=list(adapt_delta=0.999)
)

pp_check(m.Pr.hurdle42, ndraws=100)+
  xlim(0,300)

#### 51 hurdle: Pr_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot), fam=hurdle_gaussian ####
SEED = 19111996
m.Pr.hurdle51b <- brm(
  bf(Pr_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot),
     hu ~ sowndivLog*treatment + (1|block/plot)), #should predict 
  data = dBEF_nem21, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle51, ndraws=100)


emt = emtrends(m.Pr.hurdle51, "treatment", var="sowndivLog")
summary(emt, point.est=mean)
summary(emt, point.est=mean, level = .9) #get slopes

emt.pairs <- pairs(emt)
summary(emt.pairs, point.est=mean, level = .95) #get differences in slopes betweem treatments
bayestestR::p_direction(emt.pairs) #probability of difference between the treatments in the output to be negative
bayestestR::p_direction(m.Pr.hurdle51)
#### 61 hurdle: Pr_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot), fam=hurdle_gaussian ####
SEED = 19111996
m.Pr.hurdle61 <- brm(
  bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle61, ndraws=100)+
  xlim(0,300)



####comparison - waic ####

waic(m.Pr.hurdle31, m.Pr.hurdle42)
loo(m.Pr.hurdle31, m.Pr.hurdle42)

waic(m.Pr.hurdle31)
loo(m.Pr.hurdle31)
#### save hurdle models ####
save(m.Pr.hurdle11, #Pr.Log ~ sowndivLogStd, fam=hurdle_gaussian
     m.Pr.hurdle21, m.Pr.hurdle22, #Pr.Log ~ sowndiv, fam=hurdle_gaussian
     m.Pr.hurdle31, #Pr ~ sowndivLogStd, fam=hurdle_lognormal
     m.Pr.hurdle42, #Pr ~ sowndiv, fam=hurdle_lognormal
     m.Pr.hurdle51, #Pr.Log ~ sowndivLog, fam=hurdle_gaussian
     m.Pr.hurdle61, #Pr ~ sowndivLog, fam=hurdle_lognormal
     
     
     file = "./statistics/brms/231108_Pr_hurdle.RData")

load(file = "./statistics/brms/231108_Pr_hurdle.RData")
conditional_effects(m.Pr.hurdle51)

