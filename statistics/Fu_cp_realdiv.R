library(brms)
library(rstan)
library(ggplot2)

# fitting trophic group densities ~ funcdiv
#data:
#exclude 60 sp.:
dat <- subset(dBEF_nem21, sowndiv != 60) 
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog) %>%
  mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
         .after = realdivLog) %>%
  mutate(funcdivStd = ( (funcdiv - mean(funcdiv)) / sd(funcdiv) ),
         .after = funcdiv)


datW1 <- subset(dat, week=="W1")
datW2 <- subset(dat, week=="W2")

#priors    
beta_coeff_priors <- prior(normal(0,10), class = "b")  

#### Fu2 ~ realdiv ####
#selecting models based on looic,
#see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff

beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(dat$Fu2==0) #12 -> hu ~ 1

m.Fu2_realdiv_p <- brm(
  bf(Fu2 ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ 1),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #10 div

pp_check(m.Fu2_realdiv_p, ndraws=100)+
  xlim(0,300)
summary(m.Fu2_realdiv_p, prob =0.9)

m.Fu2_realdiv_p2 <- update(m.Fu2_realdiv_p,
                           bf(Fu2 ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                                (1|block/plot),
                              hu ~ 1 ),
                           seed=SEED)
m.Fu2_realdiv_p2 %>% pp_check(ndraws=100)

m.Fu2_realdiv_p31 <- update(m.Fu2_realdiv_p,
                            bf(Fu2 ~ realdivLogStd*treatment + treatment*week + 
                                 (1|block/plot),
                               hu ~ 1 ),
                            seed=SEED)
m.Fu2_realdiv_p31 %>% pp_check(ndraws=100)

m.Fu2_realdiv_p32 <- update(m.Fu2_realdiv_p,
                            bf(Fu2 ~ realdivLogStd*treatment + realdivLogStd*week  + 
                                 (1|block/plot),
                               hu ~ 1 ),
                            seed=SEED)
m.Fu2_realdiv_p32 %>% pp_check(ndraws=100)

m.Fu2_realdiv_p4 <- update(m.Fu2_realdiv_p,
                           bf(Fu2 ~ realdivLogStd*treatment  + week + (1|block/plot),
                              hu ~ 1 ),
                           seed=SEED)
m.Fu2_realdiv_p4 %>% pp_check(ndraws=100)

m.Fu2_realdiv_p5 <- update(m.Fu2_realdiv_p,
                           bf(Fu2 ~ realdivLogStd*treatment + (1|block/plot),
                              hu ~ 1 ),
                           seed=SEED)
m.Fu2_realdiv_p5 %>% pp_check(ndraws=100)

loo.Fu2 <- loo(m.Fu2_realdiv_p, m.Fu2_realdiv_p2, m.Fu2_realdiv_p31, 
               m.Fu2_realdiv_p32, m.Fu2_realdiv_p4, m.Fu2_realdiv_p5)
loo.Fu2

save(m.Fu2_realdiv_p, m.Fu2_realdiv_p2, m.Fu2_realdiv_p31, 
     m.Fu2_realdiv_p32, m.Fu2_realdiv_p4, m.Fu2_realdiv_p5,
     file = "./statistics/brms/240205_Fu2_realdiv.RData")

rm(m.Fu2_realdiv_p, m.Fu2_realdiv_p2, m.Fu2_realdiv_p31, 
   m.Fu2_realdiv_p32, m.Fu2_realdiv_p4, m.Fu2_realdiv_p5)

#### Fu3 ~ realdiv ####
#selecting models based on looic,
#see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff

beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(dat$Fu3==0) #59 -> hu ~ realdivLogStd*treatment + (1|block/plot)

m.Fu3_realdiv_p <- brm(
  bf(Fu3 ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd*treatment + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #10 div

pp_check(m.Fu3_realdiv_p, ndraws=100)+
  xlim(0,300)
summary(m.Fu3_realdiv_p, prob =0.9)

m.Fu3_realdiv_p2 <- update(m.Fu3_realdiv_p,
                           bf(Fu3 ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                                (1|block/plot),
                              hu ~ realdivLogStd*treatment + (1|block/plot) ),
                           seed=SEED)
m.Fu3_realdiv_p2 %>% pp_check(ndraws=100)

m.Fu3_realdiv_p31 <- update(m.Fu3_realdiv_p,
                            bf(Fu3 ~ realdivLogStd*treatment + treatment*week + 
                                 (1|block/plot),
                               hu ~ 1 ),
                            seed=SEED)
m.Fu3_realdiv_p31 %>% pp_check(ndraws=100)

m.Fu3_realdiv_p32 <- update(m.Fu3_realdiv_p,
                            bf(Fu3 ~ realdivLogStd*treatment + realdivLogStd*week  + 
                                 (1|block/plot),
                               hu ~ realdivLogStd*treatment + (1|block/plot) ),
                            seed=SEED)
m.Fu3_realdiv_p32 %>% pp_check(ndraws=100)

m.Fu3_realdiv_p4 <- update(m.Fu3_realdiv_p,
                           bf(Fu3 ~ realdivLogStd*treatment  + week + (1|block/plot),
                              hu ~ realdivLogStd*treatment + (1|block/plot) ),
                           seed=SEED)
m.Fu3_realdiv_p4 %>% pp_check(ndraws=100)

m.Fu3_realdiv_p5 <- update(m.Fu3_realdiv_p,
                           bf(Fu3 ~ realdivLogStd*treatment + (1|block/plot),
                              hu ~ realdivLogStd*treatment + (1|block/plot) ),
                           seed=SEED)
m.Fu3_realdiv_p5 %>% pp_check(ndraws=100)

loo.Fu3 <- loo(m.Fu3_realdiv_p, m.Fu3_realdiv_p2, m.Fu3_realdiv_p31, 
               m.Fu3_realdiv_p32, m.Fu3_realdiv_p4, m.Fu3_realdiv_p5)
loo.Fu3

save(m.Fu3_realdiv_p, m.Fu3_realdiv_p2, m.Fu3_realdiv_p31, 
     m.Fu3_realdiv_p32, m.Fu3_realdiv_p4, m.Fu3_realdiv_p5,
     file = "./statistics/brms/240205_Fu3_realdiv.RData")

rm(m.Fu3_realdiv_p, m.Fu3_realdiv_p2, m.Fu3_realdiv_p31, 
   m.Fu3_realdiv_p32, m.Fu3_realdiv_p4, m.Fu3_realdiv_p5)

#### Fu4 ~ realdiv ####
#selecting models based on looic,
#see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff

beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(dat$Fu4==0) #68/228 -> hu ~ realdivLogStd*treatment + (1|block/plot)

m.Fu4_realdiv_p <- brm(
  bf(Fu4 ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd*treatment + (1|block/plot) ),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #10 div

pp_check(m.Fu4_realdiv_p, ndraws=100)+
  xlim(0,300)
summary(m.Fu4_realdiv_p, prob =0.9)

m.Fu4_realdiv_p2 <- update(m.Fu4_realdiv_p,
                           bf(Fu4 ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                                (1|block/plot),
                              hu ~ realdivLogStd*treatment + (1|block/plot) ),
                           seed=SEED)
m.Fu4_realdiv_p2 %>% pp_check(ndraws=100)

m.Fu4_realdiv_p31 <- update(m.Fu4_realdiv_p,
                            bf(Fu4 ~ realdivLogStd*treatment + treatment*week + 
                                 (1|block/plot),
                               hu ~ realdivLogStd*treatment + (1|block/plot) ),
                            seed=SEED)
m.Fu4_realdiv_p31 %>% pp_check(ndraws=100)

m.Fu4_realdiv_p32 <- update(m.Fu4_realdiv_p,
                            bf(Fu4 ~ realdivLogStd*treatment + realdivLogStd*week  + 
                                 (1|block/plot),
                               hu ~ realdivLogStd*treatment + (1|block/plot) ),
                            seed=SEED)
m.Fu4_realdiv_p32 %>% pp_check(ndraws=100)

m.Fu4_realdiv_p4 <- update(m.Fu4_realdiv_p,
                           bf(Fu4 ~ realdivLogStd*treatment  + week + (1|block/plot),
                              hu ~ realdivLogStd*treatment + (1|block/plot) ),
                           seed=SEED)
m.Fu4_realdiv_p4 %>% pp_check(ndraws=100)

m.Fu4_realdiv_p5 <- update(m.Fu4_realdiv_p,
                           bf(Fu4 ~ realdivLogStd*treatment + (1|block/plot),
                              hu ~ realdivLogStd*treatment + (1|block/plot) ),
                           seed=SEED)
m.Fu4_realdiv_p5 %>% pp_check(ndraws=100)

loo.Fu4 <- loo(m.Fu4_realdiv_p, m.Fu4_realdiv_p2, m.Fu4_realdiv_p31, 
               m.Fu4_realdiv_p32, m.Fu4_realdiv_p4, m.Fu4_realdiv_p5)
loo.Fu4

save(m.Fu4_realdiv_p, m.Fu4_realdiv_p2, m.Fu4_realdiv_p31, 
     m.Fu4_realdiv_p32, m.Fu4_realdiv_p4, m.Fu4_realdiv_p5,
     file = "./statistics/brms/240205_Fu4_realdiv.RData")

rm(m.Fu4_realdiv_p, m.Fu4_realdiv_p2, m.Fu4_realdiv_p31, 
   m.Fu4_realdiv_p32, m.Fu4_realdiv_p4, m.Fu4_realdiv_p5)

#### Fu5 ~ realdiv ####
#selecting models based on looic,
#see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff

beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(dat$Fu5==0) #228/228 -> no modelling possible

#### model selection####
load("./statistics/brms/240205_Fu2_realdiv.RData")
load("./statistics/brms/240205_Fu3_realdiv.RData")
load("./statistics/brms/240205_Fu4_realdiv.RData")

loo.Fu2 <- loo(m.Fu2_realdiv_p, m.Fu2_realdiv_p2, m.Fu2_realdiv_p31, 
               m.Fu2_realdiv_p32, m.Fu2_realdiv_p4, m.Fu2_realdiv_p5)
loo.Fu2 #p5

loo.Fu3 <- loo(m.Fu3_realdiv_p, m.Fu3_realdiv_p2, m.Fu3_realdiv_p31, 
               m.Fu3_realdiv_p32, m.Fu3_realdiv_p4, m.Fu3_realdiv_p5)
loo.Fu3 #p5

loo.Fu4 <- loo(m.Fu4_realdiv_p, m.Fu4_realdiv_p2, m.Fu4_realdiv_p31, 
               m.Fu4_realdiv_p32, m.Fu4_realdiv_p4, m.Fu4_realdiv_p5)
loo.Fu4 #p5

#posterior predictive checks on selected models:
pp_check(m.Fu2_realdiv_p5, ndraws = 100) +
  xlim(0,2000)
pp_check(m.Fu3_realdiv_p5, ndraws = 100) +
  xlim(0,500)
pp_check(m.Fu4_realdiv_p5, ndraws = 100)

#conditional effects:
conditional_effects(m.Fu2_realdiv_p5) #slightly positive trend in t1/t2, slightly neg in t3
conditional_effects(m.Fu3_realdiv_p5) #positive in all treatments
conditional_effects(m.Fu4_realdiv_p5) #slightly positive in t1/t3, neutral in t2


#save the selected models
save(m.Fu2_realdiv_p5, m.Fu3_realdiv_p5, m.Fu4_realdiv_p5,
     file = "./statistics/brms/240205_Fu_cp_realdiv_mselect.RData")
