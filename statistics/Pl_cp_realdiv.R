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

#### Pl2 ~ realdiv ####
#selecting models based on looic,
#see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff

beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(dat$Pl2==0) #13 -> hu ~ 1

m.Pl2_realdiv_p <- brm(
  bf(Pl2 ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ 1),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #10 div

pp_check(m.Pl2_realdiv_p, ndraws=100)+
  xlim(0,300)
summary(m.Pl2_realdiv_p, prob =0.9)

m.Pl2_realdiv_p2 <- update(m.Pl2_realdiv_p,
                           bf(Pl2 ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                                (1|block/plot),
                              hu ~ 1),
                           seed=SEED)
m.Pl2_realdiv_p2 %>% pp_check(ndraws=100)

m.Pl2_realdiv_p31 <- update(m.Pl2_realdiv_p,
                            bf(Pl2 ~ realdivLogStd*treatment + treatment*week + 
                                 (1|block/plot),
                               hu ~ 1),
                            seed=SEED)
m.Pl2_realdiv_p31 %>% pp_check(ndraws=100)

m.Pl2_realdiv_p32 <- update(m.Pl2_realdiv_p,
                            bf(Pl2 ~ realdivLogStd*treatment + realdivLogStd*week  + 
                                 (1|block/plot),
                               hu ~ 1),
                            seed=SEED)
m.Pl2_realdiv_p32 %>% pp_check(ndraws=100)

m.Pl2_realdiv_p4 <- update(m.Pl2_realdiv_p,
                           bf(Pl2 ~ realdivLogStd*treatment  + week + (1|block/plot),
                              hu ~ 1),
                           seed=SEED)
m.Pl2_realdiv_p4 %>% pp_check(ndraws=100)

m.Pl2_realdiv_p5 <- update(m.Pl2_realdiv_p,
                           bf(Pl2 ~ realdivLogStd*treatment + (1|block/plot),
                              hu ~ 1),
                           seed=SEED)
m.Pl2_realdiv_p5 %>% pp_check(ndraws=100)

loo.Pl2 <- loo(m.Pl2_realdiv_p, m.Pl2_realdiv_p2, m.Pl2_realdiv_p31, 
               m.Pl2_realdiv_p32, m.Pl2_realdiv_p4, m.Pl2_realdiv_p5)
loo.Pl2

save(m.Pl2_realdiv_p, m.Pl2_realdiv_p2, m.Pl2_realdiv_p31, 
     m.Pl2_realdiv_p32, m.Pl2_realdiv_p4, m.Pl2_realdiv_p5,
     file = "./statistics/brms/240205_Pl2_realdiv.RData")

rm(m.Pl2_realdiv_p, m.Pl2_realdiv_p2, m.Pl2_realdiv_p31, 
   m.Pl2_realdiv_p32, m.Pl2_realdiv_p4, m.Pl2_realdiv_p5)

#### Pl3 ~ realdiv ####
#selecting models based on looic,
#see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff

beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(dat$Pl3==0) #10 -> hu ~ 1

m.Pl3_realdiv_p <- brm(
  bf(Pl3 ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ 1),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #10 div

pp_check(m.Pl3_realdiv_p, ndraws=100)+
  xlim(0,300)
summary(m.Pl3_realdiv_p, prob =0.9)

m.Pl3_realdiv_p2 <- update(m.Pl3_realdiv_p,
                           bf(Pl3 ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                                (1|block/plot),
                              hu ~ 1),
                           seed=SEED)
m.Pl3_realdiv_p2 %>% pp_check(ndraws=100)

m.Pl3_realdiv_p31 <- update(m.Pl3_realdiv_p,
                            bf(Pl3 ~ realdivLogStd*treatment + treatment*week + 
                                 (1|block/plot),
                               hu ~ 1),
                            seed=SEED)
m.Pl3_realdiv_p31 %>% pp_check(ndraws=100)

m.Pl3_realdiv_p32 <- update(m.Pl3_realdiv_p,
                            bf(Pl3 ~ realdivLogStd*treatment + realdivLogStd*week  + 
                                 (1|block/plot),
                               hu ~ 1),
                            seed=SEED)
m.Pl3_realdiv_p32 %>% pp_check(ndraws=100)

m.Pl3_realdiv_p4 <- update(m.Pl3_realdiv_p,
                           bf(Pl3 ~ realdivLogStd*treatment  + week + (1|block/plot),
                              hu ~ 1),
                           seed=SEED)
m.Pl3_realdiv_p4 %>% pp_check(ndraws=100)

m.Pl3_realdiv_p5 <- update(m.Pl3_realdiv_p,
                           bf(Pl3 ~ realdivLogStd*treatment + (1|block/plot),
                              hu ~ 1),
                           seed=SEED)
m.Pl3_realdiv_p5 %>% pp_check(ndraws=100)

loo.Pl3 <- loo(m.Pl3_realdiv_p, m.Pl3_realdiv_p2, m.Pl3_realdiv_p31, 
               m.Pl3_realdiv_p32, m.Pl3_realdiv_p4, m.Pl3_realdiv_p5)
loo.Pl3

save(m.Pl3_realdiv_p, m.Pl3_realdiv_p2, m.Pl3_realdiv_p31, 
     m.Pl3_realdiv_p32, m.Pl3_realdiv_p4, m.Pl3_realdiv_p5,
     file = "./statistics/brms/240205_Pl3_realdiv.RData")

rm(m.Pl3_realdiv_p, m.Pl3_realdiv_p2, m.Pl3_realdiv_p31, 
   m.Pl3_realdiv_p32, m.Pl3_realdiv_p4, m.Pl3_realdiv_p5)

#### Pl4 ~ realdiv ####
#selecting models based on looic,
#see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff

beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(dat$Pl4==0) #110/228 -> hu ~ realdivLogStd*treatment + (1|block/plot)

m.Pl4_realdiv_p <- brm(
  bf(Pl4 ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd*treatment  + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #10 div

pp_check(m.Pl4_realdiv_p, ndraws=100)+
  xlim(0,300)
summary(m.Pl4_realdiv_p, prob =0.9)

m.Pl4_realdiv_p2 <- update(m.Pl4_realdiv_p,
                           bf(Pl4 ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                                (1|block/plot),
                              hu ~ realdivLogStd*treatment  + (1|block/plot)),
                           seed=SEED)
m.Pl4_realdiv_p2 %>% pp_check(ndraws=100)

m.Pl4_realdiv_p31 <- update(m.Pl4_realdiv_p,
                            bf(Pl4 ~ realdivLogStd*treatment + treatment*week + 
                                 (1|block/plot),
                               hu ~ realdivLogStd*treatment  + (1|block/plot)),
                            seed=SEED)
m.Pl4_realdiv_p31 %>% pp_check(ndraws=100)

m.Pl4_realdiv_p32 <- update(m.Pl4_realdiv_p,
                            bf(Pl4 ~ realdivLogStd*treatment + realdivLogStd*week  + 
                                 (1|block/plot),
                               hu ~ realdivLogStd*treatment  + (1|block/plot)),
                            seed=SEED)
m.Pl4_realdiv_p32 %>% pp_check(ndraws=100)

m.Pl4_realdiv_p4 <- update(m.Pl4_realdiv_p,
                           bf(Pl4 ~ realdivLogStd*treatment  + week + (1|block/plot),
                              hu ~ realdivLogStd*treatment  + (1|block/plot)),
                           seed=SEED)
m.Pl4_realdiv_p4 %>% pp_check(ndraws=100)

m.Pl4_realdiv_p5 <- update(m.Pl4_realdiv_p,
                           bf(Pl4 ~ realdivLogStd*treatment + (1|block/plot),
                              hu ~ realdivLogStd*treatment  + (1|block/plot)),
                           seed=SEED)
m.Pl4_realdiv_p5 %>% pp_check(ndraws=100)

loo.Pl4 <- loo(m.Pl4_realdiv_p, m.Pl4_realdiv_p2, m.Pl4_realdiv_p31, 
               m.Pl4_realdiv_p32, m.Pl4_realdiv_p4, m.Pl4_realdiv_p5)
loo.Pl4

save(m.Pl4_realdiv_p, m.Pl4_realdiv_p2, m.Pl4_realdiv_p31, 
     m.Pl4_realdiv_p32, m.Pl4_realdiv_p4, m.Pl4_realdiv_p5,
     file = "./statistics/brms/240205_Pl4_realdiv.RData")

rm(m.Pl4_realdiv_p, m.Pl4_realdiv_p2, m.Pl4_realdiv_p31, 
   m.Pl4_realdiv_p32, m.Pl4_realdiv_p4, m.Pl4_realdiv_p5)

####Pl5 ~ realdiv ####
sum(dat$Pl5==0) #224/228 -> no model possible

####select models ####
load("./statistics/brms/240205_Pl2_realdiv.RData")
load("./statistics/brms/240205_Pl3_realdiv.RData")
load("./statistics/brms/240205_Pl4_realdiv.RData")

loo.Pl2 <- loo(m.Pl2_realdiv_p, m.Pl2_realdiv_p2, m.Pl2_realdiv_p31, 
               m.Pl2_realdiv_p32, m.Pl2_realdiv_p4, m.Pl2_realdiv_p5)
loo.Pl2 #p5

loo.Pl3 <- loo(m.Pl3_realdiv_p, m.Pl3_realdiv_p2, m.Pl3_realdiv_p31, 
               m.Pl3_realdiv_p32, m.Pl3_realdiv_p4, m.Pl3_realdiv_p5)
loo.Pl3 #p5

loo.Pl4 <- loo(m.Pl4_realdiv_p, m.Pl4_realdiv_p2, m.Pl4_realdiv_p31, 
               m.Pl4_realdiv_p32, m.Pl4_realdiv_p4, m.Pl4_realdiv_p5)
loo.Pl4 #p5

#posterior predictive checks on selected models:
pp_check(m.Pl2_realdiv_p5, ndraws = 100) +
  xlim(0,1000)
pp_check(m.Pl3_realdiv_p5, ndraws = 100) +
  xlim(0,1000)
pp_check(m.Pl4_realdiv_p5, ndraws = 100) +
  xlim(0,150)

#conditional effects
conditional_effects(m.Pl2_realdiv_p5) #positive trend in t1, slightly pos. in t2, slightly neg. in t3
conditional_effects(m.Pl3_realdiv_p5) #slightly pos in t1/t2/t2 at same level
conditional_effects(m.Pl4_realdiv_p5) #neutral at same level


save(m.Pl2_realdiv_p5, m.Pl3_realdiv_p5, m.Pl4_realdiv_p5,
     file = "./statistics/brms/240205_Pl_cp_realdiv_mselect.RData")
