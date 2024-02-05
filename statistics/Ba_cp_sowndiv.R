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


#### Ba1 ~ sowndiv ####
#selecting models based on looic,
#see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff

beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(dat$Ba1==0) #109 -> hu ~ treatment*sowndivLogStd

m.Ba1_sowndiv_p <- brm(
  bf(Ba1 ~ sowndivLogStd*treatment*week + (1|block/plot),
     hu ~ sowndivLogStd*treatment  + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #10 div

pp_check(m.Ba1_sowndiv_p, ndraws=100)+
  xlim(0,300)
summary(m.Ba1_sowndiv_p, prob =0.9)
conditional_effects(m.Ba1_sowndiv_p)

m.Ba1_sowndiv_p2 <- update(m.Ba1_sowndiv_p,
                           bf(Ba1 ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                (1|block/plot),
                              hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                           seed=SEED)
  m.Ba1_sowndiv_p2 %>% pp_check(ndraws=100)

m.Ba1_sowndiv_p31 <- update(m.Ba1_sowndiv_p,
                           bf(Ba1 ~ sowndivLogStd*treatment + treatment*week + 
                                (1|block/plot),
                              hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                           seed=SEED)
  m.Ba1_sowndiv_p31 %>% pp_check(ndraws=100)

m.Ba1_sowndiv_p32 <- update(m.Ba1_sowndiv_p,
                            bf(Ba1 ~ sowndivLogStd*treatment + sowndivLogStd*week  + 
                                 (1|block/plot),
                               hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                            seed=SEED)
  m.Ba1_sowndiv_p32 %>% pp_check(ndraws=100)

m.Ba1_sowndiv_p4 <- update(m.Ba1_sowndiv_p,
                            bf(Ba1 ~ sowndivLogStd*treatment  + week + (1|block/plot),
                               hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                            seed=SEED)
  m.Ba1_sowndiv_p4 %>% pp_check(ndraws=100)

m.Ba1_sowndiv_p5 <- update(m.Ba1_sowndiv_p,
                           bf(Ba1 ~ sowndivLogStd*treatment + (1|block/plot),
                              hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                           seed=SEED)
  m.Ba1_sowndiv_p5 %>% pp_check(ndraws=100)

loo.Ba1 <- loo(m.Ba1_sowndiv_p, m.Ba1_sowndiv_p2, m.Ba1_sowndiv_p31, 
               m.Ba1_sowndiv_p32, m.Ba1_sowndiv_p4, m.Ba1_sowndiv_p5)
loo.Ba1

#### Ba2 ~ sowndiv ####
  #selecting models based on looic,
  #see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff
  
  beta_coeff_priors <- prior(normal(0,10), class = "b")  
  SEED = 22061996
  
  sum(dat$Ba2==0) #22 -> hu ~ 1
  
  m.Ba2_sowndiv_p <- brm(
    bf(Ba2 ~ sowndivLogStd*treatment*week + (1|block/plot),
       hu ~ 1 ),
    data = dat, 
    prior = beta_coeff_priors,
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #10 div
  
  pp_check(m.Ba2_sowndiv_p, ndraws=100)+
    xlim(0,300)
  summary(m.Ba2_sowndiv_p, prob =0.9)
  conditional_effects(m.Ba2_sowndiv_p)
  
  m.Ba2_sowndiv_p2 <- update(m.Ba2_sowndiv_p,
                             bf(Ba2 ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                  (1|block/plot),
                                hu ~ 1 ),
                             seed=SEED)
  m.Ba2_sowndiv_p2 %>% pp_check(ndraws=100)
  
  m.Ba2_sowndiv_p31 <- update(m.Ba2_sowndiv_p,
                              bf(Ba2 ~ sowndivLogStd*treatment + treatment*week + 
                                   (1|block/plot),
                                 hu ~ 1 ),
                              seed=SEED)
  m.Ba2_sowndiv_p31 %>% pp_check(ndraws=100)
  
  m.Ba2_sowndiv_p32 <- update(m.Ba2_sowndiv_p,
                              bf(Ba2 ~ sowndivLogStd*treatment + sowndivLogStd*week  + 
                                   (1|block/plot),
                                 hu ~ 1 ),
                              seed=SEED)
  m.Ba2_sowndiv_p32 %>% pp_check(ndraws=100)
  
  m.Ba2_sowndiv_p4 <- update(m.Ba2_sowndiv_p,
                             bf(Ba2 ~ sowndivLogStd*treatment  + week + (1|block/plot),
                                hu ~ 1 ),
                             seed=SEED)
  m.Ba2_sowndiv_p4 %>% pp_check(ndraws=100)
  
  m.Ba2_sowndiv_p5 <- update(m.Ba2_sowndiv_p,
                             bf(Ba2 ~ sowndivLogStd*treatment + (1|block/plot),
                                hu ~ 1 ),
                             seed=SEED)
  m.Ba2_sowndiv_p5 %>% pp_check(ndraws=100)
  
  loo.Ba2 <- loo(m.Ba2_sowndiv_p, m.Ba2_sowndiv_p2, m.Ba2_sowndiv_p31, 
                 m.Ba2_sowndiv_p32, m.Ba2_sowndiv_p4, m.Ba2_sowndiv_p5)
  loo.Ba2
  
#### Ba3 ~ sowndiv ####
  #selecting models based on looic,
  #see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff
  
  beta_coeff_priors <- prior(normal(0,10), class = "b")  
  SEED = 22061996
  
  sum(dat$Ba3==0) #218/228 -> not worth fitting models
  
#### Ba4 ~ sowndiv ####
  #selecting models based on looic,
  #see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff
  
  beta_coeff_priors <- prior(normal(0,10), class = "b")  
  SEED = 22061996
  
  sum(dat$Ba4==0) #139 -> hu ~ treatment*sowndivLogStd
  
  m.Ba4_sowndiv_p <- brm(
    bf(Ba4 ~ sowndivLogStd*treatment*week + (1|block/plot),
       hu ~ sowndivLogStd*treatment  + (1|block/plot)),
    data = dat, 
    prior = beta_coeff_priors,
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #10 div
  
  pp_check(m.Ba4_sowndiv_p, ndraws=100)+
    xlim(0,300)
  summary(m.Ba4_sowndiv_p, prob =0.9)
  conditional_effects(m.Ba4_sowndiv_p)
  
  m.Ba4_sowndiv_p2 <- update(m.Ba4_sowndiv_p,
                             bf(Ba4 ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                  (1|block/plot),
                                hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                             seed=SEED)
  m.Ba4_sowndiv_p2 %>% pp_check(ndraws=100)
  
  m.Ba4_sowndiv_p31 <- update(m.Ba4_sowndiv_p,
                              bf(Ba4 ~ sowndivLogStd*treatment + treatment*week + 
                                   (1|block/plot),
                                 hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                              seed=SEED)
  m.Ba4_sowndiv_p31 %>% pp_check(ndraws=100)
  
  m.Ba4_sowndiv_p32 <- update(m.Ba4_sowndiv_p,
                              bf(Ba4 ~ sowndivLogStd*treatment + sowndivLogStd*week  + 
                                   (1|block/plot),
                                 hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                              seed=SEED)
  m.Ba4_sowndiv_p32 %>% pp_check(ndraws=100)
  
  m.Ba4_sowndiv_p4 <- update(m.Ba4_sowndiv_p,
                             bf(Ba4 ~ sowndivLogStd*treatment  + week + (1|block/plot),
                                hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                             seed=SEED)
  m.Ba4_sowndiv_p4 %>% pp_check(ndraws=100)
  
  m.Ba4_sowndiv_p5 <- update(m.Ba4_sowndiv_p,
                             bf(Ba4 ~ sowndivLogStd*treatment + (1|block/plot),
                                hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                             seed=SEED)
  m.Ba4_sowndiv_p5 %>% pp_check(ndraws=100)
  
  loo.Ba4 <- loo(m.Ba4_sowndiv_p, m.Ba4_sowndiv_p2, m.Ba4_sowndiv_p31, 
                 m.Ba4_sowndiv_p32, m.Ba4_sowndiv_p4, m.Ba4_sowndiv_p5)
  loo.Ba4  
  
#### Ba5 ~ sowndiv ####
  #selecting models based on looic,
  #see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff
  
  beta_coeff_priors <- prior(normal(0,10), class = "b")  
  SEED = 22061996
  
  sum(dat$Ba5==0) #228/228 -> no model fitting possible
  
#### save best fit models ####
  save(m.Ba1_sowndiv_p, m.Ba2_sowndiv_p, m.Ba4_sowndiv_p,
       file = "./statistics/brms/240205_Ba_cp_sowndiv_mselect.RData")
  
 