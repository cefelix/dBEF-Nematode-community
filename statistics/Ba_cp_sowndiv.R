library(brms)
library(rstan)
library(ggplot2)
library(dplyr)

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

save(m.Ba1_sowndiv_p, m.Ba1_sowndiv_p2, m.Ba1_sowndiv_p31, 
     m.Ba1_sowndiv_p32, m.Ba1_sowndiv_p4, m.Ba1_sowndiv_p5,
     file = "./statistics/brms/240205_Ba1_sowndiv.RData")

rm(m.Ba1_sowndiv_p, m.Ba1_sowndiv_p2, m.Ba1_sowndiv_p31, 
   m.Ba1_sowndiv_p32, m.Ba1_sowndiv_p4, m.Ba1_sowndiv_p5)

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
  
  save(m.Ba2_sowndiv_p, m.Ba2_sowndiv_p2, m.Ba2_sowndiv_p31, 
       m.Ba2_sowndiv_p32, m.Ba2_sowndiv_p4, m.Ba2_sowndiv_p5,
       file = "./statistics/brms/240205_Ba2_sowndiv.RData")
  
  rm(m.Ba2_sowndiv_p, m.Ba2_sowndiv_p2, m.Ba2_sowndiv_p31, 
     m.Ba2_sowndiv_p32, m.Ba2_sowndiv_p4, m.Ba2_sowndiv_p5)
  
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
  
  save(m.Ba4_sowndiv_p, m.Ba4_sowndiv_p2, m.Ba4_sowndiv_p31, 
       m.Ba4_sowndiv_p32, m.Ba4_sowndiv_p4, m.Ba4_sowndiv_p5,
       file = "./statistics/brms/240205_Ba4_sowndiv.RData")
  
  rm(m.Ba4_sowndiv_p, m.Ba4_sowndiv_p2, m.Ba4_sowndiv_p31, 
     m.Ba4_sowndiv_p32, m.Ba4_sowndiv_p4, m.Ba4_sowndiv_p5)
  
#### Ba5 ~ sowndiv ####
  #selecting models based on looic,
  #see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff
  
  beta_coeff_priors <- prior(normal(0,10), class = "b")  
  SEED = 22061996
  
  sum(dat$Ba5==0) #228/228 -> no model fitting possible
  
#### model selection ####
  load("./statistics/brms/240205_Ba1_sowndiv.RData")
  load("./statistics/brms/240205_Ba2_sowndiv.RData")
  load("./statistics/brms/240205_Ba4_sowndiv.RData")
  
  loo.Ba1 <- loo(m.Ba1_sowndiv_p, m.Ba1_sowndiv_p2, m.Ba1_sowndiv_p31, 
                 m.Ba1_sowndiv_p32, m.Ba1_sowndiv_p4, m.Ba1_sowndiv_p5)
  loo.Ba1 #p5
  
  loo.Ba2 <- loo(m.Ba2_sowndiv_p, m.Ba2_sowndiv_p2, m.Ba2_sowndiv_p31, 
                 m.Ba2_sowndiv_p32, m.Ba2_sowndiv_p4, m.Ba2_sowndiv_p5)
  loo.Ba2 #p5
  
  loo.Ba4 <- loo(m.Ba4_sowndiv_p, m.Ba4_sowndiv_p2, m.Ba4_sowndiv_p31, 
                 m.Ba4_sowndiv_p32, m.Ba4_sowndiv_p4, m.Ba4_sowndiv_p5)
  loo.Ba4 #p5
  
  #posterior predictive checks on selected models:
  pp_check(m.Ba1_sowndiv_p5, ndraws = 100) +
    xlim(0,200)
  pp_check(m.Ba2_sowndiv_p5, ndraws = 100) +
    xlim(0,500)
  pp_check(m.Ba4_sowndiv_p5, ndraws = 100) +
    xlim(0,150)
  
  #conditional effects
  conditional_effects(m.Ba1_sowndiv_p5) #t1 substantially higher than t2, trends probably non-significant
  conditional_effects(m.Ba2_sowndiv_p5) #slightly positive trend in t1, slightly negativ in t2/t3
  conditional_effects(m.Ba4_sowndiv_p5) #strong negative trend in t3, slightly positive in t2/t3. t3 might be driven mainly by 2 outliers at sowndiv=1
    
    ggplot(dat, aes(x=sowndivLog, y=Ba4, col=treatment))+
      geom_jitter(width=.2) #check out the blue outliers in treatment 3:
    dat %>% filter(treatment == 3 & sowndivLog == 0 & Ba4 > 75) %>%
      select(Sample) #B4A13D1 
    dat %>% filter(treatment == 3 & sowndivLog == 1 & Ba4 > 50) %>%
      select(Sample) #B3A08D1
    dat %>% filter(treatment == 3 & sowndivLog == 2 & Ba4 > 55) %>%
      select(Sample) #B4A07D1
    #the 3 high B4 densities are from different plots
  
  
  save(m.Ba1_sowndiv_p5, m.Ba2_sowndiv_p5, m.Ba4_sowndiv_p5,
       file = "./statistics/brms/240205_Ba_cp_sowndiv_mselect.RData")
  
 