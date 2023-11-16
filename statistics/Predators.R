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

#### 11a hurdle: Pr_per100g ~ sowndiv*treatment + (1|block/plot), fam=hurdle_lognormal ####
SEED = 19111996
m.Pr.hurdle11a <- brm(
  bf(Pr_per100g ~ sowndiv*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

m.Pr.hurdle12a <- update(m.Pr.hurdle11a,
                        control=list(adapt_delta=0.999)
)

pp_check(m.Pr.hurdle12a, ndraws=100)+
  xlim(0,300)



#### 11b hurdle: Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot), fam=hurdle_lognormal ####
SEED = 19111996
m.Pr.hurdle11b <- brm(
  bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle11b, ndraws=100)+
  xlim(0,300)


#### 11c hurdle: Pr_per100g ~ sowndivLog*treatment + (1|block/plot), fam=hurdle_lognormal ####
SEED = 19111996
m.Pr.hurdle11c <- brm(
  bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle11c, ndraws=100)+
  xlim(0,300)

####21 hurdle: Pr_per100g ~ sowndivLog*treatment + (1|block/plot), hu~term, fam=hurdle_lognormal ####
SEED = 19111996
m.Pr.hurdle21 <- brm(
  bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~  sowndivLog*treatment + (1|block/plot)),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle71, ndraws=100)+
  xlim(0,300) 

#model predictions:
  predictions <- conditional_effects(m.Pr.hurdle71)[[3]]
  predictions$estimate__

#plot results
  ggplot(dBEF_nem21, aes(x=sowndivLog, y=Pr_per100g, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=estimate__, col=treatment))

#### 31 hurdle: Pr_per100g ~ sowndivLog*treatment + (1|block/plot), hu~term, fam=hurdle_lognormal, no 60 sp. plots ####
  
  
  dat = dBEF_nem21 %>% filter(sowndiv != 60)
  SEED = 19111996
  m.Pr.hurdle81 <- brm(
    bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
       hu ~  sowndivLog*treatment + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99))
  
  pp_check(m.Pr.hurdle81, ndraws = 100)
  
#model predictions:
  predictions <- conditional_effects(m.Pr.hurdle81)[[3]]
  predictions$estimate__
  
  #plot results
  ggplot(dat, aes(x=sowndivLog, y=Pr_per100g, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=estimate__, col=treatment))
  
  summary(m.Pr.hurdle81)
  

####comparison - waic ####

waic(m.Pr.hurdle31, m.Pr.hurdle42)
loo(m.Pr.hurdle31, m.Pr.hurdle42)

waic(m.Pr.hurdle31)
loo(m.Pr.hurdle31)
#### save hurdle models ####
save(#m.Pr.hurdle11a, m.Pr.hurdle11b, m.Pr.hurdle11c, #Pr.Log ~ sowndivLogStd, fam=hurdle_gaussian
     m.Pr.hurdle21,
     m.Pr.hurdle31,
     file = "./statistics/brms/231116_Pr_hurdle.RData")

load(file = "./statistics/brms/231108_Pr_hurdle.RData")

