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


#### exploration ####
  dBEF_nemSH1 <- subset(dBEF_nem, SH == 1)
  dBEF_nemSH5 <- subset(dBEF_nem, SH == 5)
  dBEF_nemSH15 <- subset(dBEF_nem, SH == 15)
  dBEF_nemSH19 <- subset(dBEF_nem, SH == 19)
  
p.1 <- ggplot(dBEF_nemSH1, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  labs(title = "SH 1")+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")

p.5 <- ggplot(dBEF_nemSH5, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")+
  labs(title = "SH 5")

p.15 <- ggplot(dBEF_nemSH15, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 15")+
  geom_smooth(method="lm")

p.19 <- ggplot(dBEF_nemSH19, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 19")+
  geom_smooth(method="lm")

grid.arrange(p.1, p.5, p.15, p.19)
rm(p.1, p.5, p.15, p.19,
   dBEF_nemSH19, dBEF_nemSH15, dBEF_nemSH5, dBEF_nemSH1)



#### 11a hurdle: Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~1, family = hurdle_lognormal ####
SEED = 22061996
m.Ba.hurdle11a <- brm(
  bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Ba.hurdle11a, ndraws=100)+
  xlim(0,1000)

#### 11b hurdle: Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~term, family = hurdle_lognormal ####
SEED = 22061996
m.Ba.hurdle11b <- brm(
  bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))


#### 31a hurdle: Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~1, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 22061996
#exclude 60 sp.:
  dat <- subset(dBEF_nem21, sowndiv != 60) 
  #standardize:  
  dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                 .after = sowndivLog)

m.Ba.hurdle31a <- brm(
  bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good

pp_check(m.Ba.hurdle31a, ndraws=100)+
  xlim(0,1000)

#### 31a_Ba1-4 hurdle: Ba1-4 ~ sowndivLogStd*treatment + (1|block/plot), hu~1, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 22061996
dat <- subset(dBEF_nem21, sowndiv != 60)
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog)

  (dat$Ba1!=0) %>% sum() # 119 of 228
  (dat$Ba2!=0) %>% sum() # 206 of 228 
  (dat$Ba3!=0) %>% sum() #  10 of 228 probably wont work
  (dat$Ba4!=0) %>% sum() #  89 of 228
  (dat$Ba5!=0) %>% sum() #   0 of 228  

#the models for Bacterivores in cp1-cp4
m.Ba1.hurdle31a <- brm(
  bf(Ba1 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~  1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) 
  #all good

  pp_check(m.Ba1.hurdle31a, ndraws = 100)+
    xlim(0,200) #thats okay

m.Ba2.hurdle31a <- brm(
  bf(Ba2 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~  1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) 

pp_check(m.Ba2.hurdle31a, ndraws = 100)+
  xlim(0,600) #okay

m.Ba3.hurdle31a <- brm(
  bf(Ba3 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~  1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) 
  #62 divergent transitions, 1592 exceedeed max_treedepth, max Rhat=1.37

dat$Ba3 %>% density() %>% plot()
m.Ba3.hurdle32a <- update(m.Ba3.hurdle31a,
                          iter = 4000, warmup = 2000,
                          control = list(adapt_delta=0.999,
                                         max_treedepth=12))
 #72 divergent transizions

pp_check(m.Ba3.hurdle31a, ndraws = 100)+xlim(0,10) #that's completely off


m.Ba4.hurdle31a <- brm(
  bf(Ba4 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~  1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))  #3 divergent transitions

m.Ba4.hurdle32a <- update(m.Ba4.hurdle31a, 
                          control=list(adapt_delta=0.999))

pp_check(m.Ba4.hurdle32a, ndraws = 100)+
  xlim(0,50) #bearable

#### 31b hurdle: Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~term, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 22061996
dat <- subset(dBEF_nem21, sowndiv != 60)
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog)

m.Ba.hurdle31b <- brm(
  bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all fine

pp_check(m.Ba.hurdle31b, ndraws=100)+
  xlim(0,1000)

#### 31b_Ba1-4 hurdle: Ba1-4 ~ sowndivLogStd*treatment + (1|block/plot), hu~term, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 22061996
dat <- subset(dBEF_nem21, sowndiv != 60)
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog)

m.Ba1.hurdle31b <- brm(
  bf(Ba1 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~  sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all fine

pp_check(m.Ba1.hurdle31b, ndraws=100)+
  xlim(0,200) 

m.Ba2.hurdle31b <- brm(
  bf(Ba2 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~  sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good

pp_check(m.Ba2.hurdle31b, ndraws=100)+
  xlim(0,700) 

m.Ba3.hurdle31b <- brm(
  bf(Ba3 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~  sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) 
#86 divergent transitions, 1934 exceeded max_treedepth, bulk ESS too low

  m.Ba3.hurdle32b <- update(m.Ba3.hurdle31b,
                            iter = 4000, warmup = 2000,
                            control=list(adapt_delta=0.999, 
                                         max_treedepth=15)) 
  
  pp_check(m.Ba3.hurdle31b, ndraws=100)+
    xlim(0,40) 
  

m.Ba4.hurdle31b <- brm(
  bf(Ba4 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~  sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) 
  #all good

pp_check(m.Ba4.hurdle31b, ndraws=100)+
  xlim(0,100) 


###31ab Ba model comparisons####
m.Ba.hurdle31a <- add_criterion(m.Ba.hurdle31a, "loo")
m.Ba.hurdle31b <- add_criterion(m.Ba.hurdle31b, "loo")

m.Ba.hurdle31a$criteria  #looic: 2502.1 (for looic: smaller = better)
m.Ba.hurdle31b$criteria  #looic: 2500

#high pareto-k values: look whether p_loo is bigger/smaller than true number of parameters
#as described here: https://mc-stan.org/loo/reference/loo-glossary.html
#to get the true number of paremeters, follow this:
#https://discourse.mc-stan.org/t/determine-number-of-parameters-in-brms-gamm-to-compare-to-p-loo-value/23804/3

#for this you have to re-compile the model: https://github.com/quentingronau/bridgesampling/issues/7
rstan::get_num_upars(m.Ba.hurdle31b$fit) #177, thus p > N/5
m.Ba.hurdle31b$criteria  #p_loo = 38.1
  
pp_check(m.Ba.hurdle31b, ndraws=100)+
  xlim(0,1200)

#for elpd: higher = better  
loo_compare(m.Ba.hurdle31a, m.Ba.hurdle31b) 


#### 41a hurdle: Ba_per100g ~ sowndivLog*treatment + (1|block/plot), hu~1, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 22061996
dat <- subset(dBEF_nem21, sowndiv != 60)
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog)

m.Ba.hurdle41a <- brm(
  bf(Ba_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Ba.hurdle41a, ndraws=100)+
  xlim(0,1000)

#### 41b hurdle: Ba_per100g ~ sowndivLog*treatment + (1|block/plot), hu~term, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 22061996
dat <- subset(dBEF_nem21, sowndiv != 60)
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog)

m.Ba.hurdle41b <- brm(
  bf(Ba_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ sowndivLog*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Ba.hurdle41b, ndraws = 100)+
  xlim(0,1000)



#### save hurdle models ####
save(#m.Ba.hurdle11a,
     #m.Ba.hurdle11b, 
     m.Ba.hurdle31a,
      m.Ba1.hurdle31a,
      m.Ba2.hurdle31a,
      m.Ba3.hurdle31a, 
      m.Ba4.hurdle31a, m.Ba4.hurdle32a,
     m.Ba.hurdle31b,
       m.Ba1.hurdle31b, 
       m.Ba2.hurdle31b,
       m.Ba3.hurdle31b, 
       m.Ba4.hurdle31b, 
     
     #m.Ba.hurdle41a,
     #m.Ba.hurdle41b,
     
     file="./statistics/brms/231127_Ba.RData")

load(file="./statistics/brms/231127_Ba.RData")
conditional_effects(m.Ba.hurdle11)
