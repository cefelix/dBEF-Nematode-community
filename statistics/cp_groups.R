####data and packages####
#problem: we have quite a few samples where we have zero nematodes of a certain group
#thus, we need to use distributions which allow to contain zeros: 
#https://www.andrewheiss.com/blog/2022/05/09/hurdle-lognormal-gaussian-brms/


library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
dBEF_nem %>%
  str()

dBEF_nem21 <- subset(dBEF_nem, year==2021)
dBEF_nem17 <- subset(dBEF_nem, year==2017)

SEED = 22061996



#### 2021's data####


#explore the data:
    par(mfrow = c(1,1))
    hist(dBEF_nem$cp1_per100g, breaks = seq(min(dBEF_nem$cp1_per100g), max(dBEF_nem$cp1_per100g), length.out=30))
    
    #how many samples have zero cp1 nematodes: 
    sum(dBEF_nem21$cp1_per100g == 0) #116 of 240
    
    ggplot(data = dBEF_nem21, aes(x = sowndivLog , y = cp1_per100g, col=treatment))+
      geom_point()+
      facet_wrap(~treatment)



#### hurdle: cp1.Log ~ sowndivLogStd*treatment + (1|block/plot), fam = hurdle_gaussian ####

m.cp1.hurdle11 <- brm( 
  bf(cp1_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
                    hu ~ 1),
  data = dBEF_nem21,
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.9))
    
summary(m.cp1.hurdle11)
pp_check(m.cp1.hurdle11, ndraws = 100)
  



#### hurdle: cp2.Log ~ sowndivLogStd*treatment + (1|block/plot), fam = hurdle_gaussian  ####
sum(dBEF_nem$cp2_per100g == 0) # 1 sample with zero cp2 nematodes
    
m.cp2.hurdle11 <- brm( 
  bf(cp2_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21,
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.9))

summary(m.cp2.hurdle11)
pp_check(m.cp2.hurdle11, ndraws = 100) #overestimating probability mass in the zero's,
  # --> maybe compare against a non-hurdle model (but loo is not possible due to custom family!)


#### hurdle: cp3.Log ~ sowndivLogStd*treatment + (1|block/plot), fam = hurdle_gaussian #### 
sum(dBEF_nem21$cp3_per100g == 0) # 5 sample2 with zero cp3 nematodes

m.cp3.hurdle11 <- brm( 
  bf(cp3_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21,
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.9))

summary(m.cp3.hurdle11)
pp_check(m.cp3.hurdle11, ndraws=100) #weird shape of y, maybe there is a better family?

#### hurdle: cp4.Log ~ sowndivLogStd*treatment + (1|block/plot), fam = hurdle_gaussian####
hist(dBEF_nem21$cp4_per100g, breaks = seq(min(dBEF_nem21$cp4_per100g), max(dBEF_nem21$cp4_per100g), length.out=30))
sum(dBEF_nem$cp4_per100g == 0) # 19 samples with zero cp4 nematodes

m.cp4.hurdle11 <- brm( 
  bf(cp4_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21,
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.9))

summary(m.cp4.hurdle11)
pp_check(m.cp4.hurdle11, ndraws=100)




#### hurdle: cp5.Log ~ sowndivLogStd*treatment + (1|block/plot), fam = hurdle_gaussian  ####
hist(dBEF_nem21$cp5_per100g, breaks = seq(min(dBEF_nem21$cp5_per100g), max(dBEF_nem21$cp5_per100g), length.out=30))
  sum(dBEF_nem21$cp5_per100g == 0) # 213 sample with zero cp5 nematodes
  #how many of the 19 zero cp4 samples also have zero cp5 individuals?
  sum(dBEF_nem21$cp4_per100g == 0 & dBEF_nem21$cp5_per100g == 0) #15 out of 19 
    #if a sample has zero cp4 individuals, its quite likely that it has zero cp5 individuals
    #on the other hand, if a sample has zero cp5 individuals, it does not mean that it is likely to have zero cp4 individuals!
  
m.cp5.hurdle11 <- brm( 
  bf(cp5_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21,
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.9))
  
  summary(m.cp5.hurdle11)
  pp_check(m.cp5.hurdle11, ndraws=100) # not a good fit
