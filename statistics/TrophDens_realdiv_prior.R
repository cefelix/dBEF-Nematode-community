library(brms)
library(rstan)
library(ggplot2)

# fitting trophic group densities ~ realdiv

#data:
#exclude 60 sp.:
dat <- subset(dBEF_nem21, sowndiv != 60) 
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog) %>%
  mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
         .after = realdivLog)


datW1 <- subset(dat, week=="W1")
datW2 <- subset(dat, week=="W2")

#priors    
beta_coeff_priors <- prior(normal(0,10), class = "b")  

####Ba ~ realdiv, both weeks:  ####
#selecting models based on looic,
#see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff

beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Ba_realdiv_p <- brm(
  bf(Ba_per100g ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd + treatment + week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good 

summary(m.Ba_realdiv_p)
pp_check(m.Ba_realdiv_p, ndraws=100)+
  xlim(0,2000)

emtrends(m.Ba_realdiv_p, specs = c("treatment", "week"), var="realdivLogStd") %>%
  summary() #CIs overlap

#remove 3way interaction:
m.Ba_realdiv_p2 <- update(m.Ba_realdiv_p, 
                          bf(Ba_per100g ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
emtrends(m.Ba_realdiv_p2, specs = c("treatment", "week"), var="realdivLogStd") %>%
  summary() #CIs overlap

#remove realdivLogStd*week
m.Ba_realdiv_p31 <- update(m.Ba_realdiv_p, 
                           bf(Ba_per100g ~ realdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

summary(m.Ba_realdiv_p31, prob=0.9)
emtrends(m.Ba_realdiv_p31, specs = c("treatment", "week"), var="realdivLogStd") %>%
  summary()

#remove treatment*week 
m.Ba_realdiv_p32 <- update(m.Ba_realdiv_p, 
                           bf(Ba_per100g ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot),
                              hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)
emtrends(m.Ba_realdiv_p32, specs = c("treatment", "week"), var="realdivLogStd") %>%
  summary()

#remove realdivLogStd*week and treatment*week 
m.Ba_realdiv_p4 <- update(m.Ba_realdiv_p, 
                          bf(Ba_per100g ~ realdivLogStd*treatment + week + (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Ba_realdiv_p4, prob=0.9)

#remove week
m.Ba_realdiv_p5 <- update(m.Ba_realdiv_p, 
                          bf(Ba_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Ba_realdiv_p5, prob=0.9)


#remove hu~treatment
m.Ba_realdiv_p6 <- update(m.Ba_realdiv_p, 
                          bf(Ba_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd + week + (1|block/plot)), 
                          seed = SEED)
#remove hu~realdiv
m.Ba_realdiv_p7 <- update(m.Ba_realdiv_p, 
                          bf(Ba_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ week + (1|block/plot)), 
                          seed = SEED)
summary(m.Ba_realdiv_p7, prob=0.9)

#remove hu~week
m.Ba_realdiv_p8 <- update(m.Ba_realdiv_p, 
                          bf(Ba_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~1), 
                          seed = SEED)
summary(m.Ba_realdiv_p8, prob=0.9)


loo.Ba <- loo(m.Ba_realdiv_p, m.Ba_realdiv_p2, m.Ba_realdiv_p31, m.Ba_realdiv_p32, m.Ba_realdiv_p4,
              m.Ba_realdiv_p5, m.Ba_realdiv_p6, m.Ba_realdiv_p7, m.Ba_realdiv_p8 )
loo.Ba

save(m.Ba_realdiv_p, m.Ba_realdiv_p2, m.Ba_realdiv_p31, m.Ba_realdiv_p32, m.Ba_realdiv_p4,
     m.Ba_realdiv_p5, m.Ba_realdiv_p6, m.Ba_realdiv_p7, m.Ba_realdiv_p8,
     file="./statistics/brms/231219_Ba_realdiv_priors.RData")

rm(m.Ba_realdiv_p, m.Ba_realdiv_p2, m.Ba_realdiv_p31, m.Ba_realdiv_p32, m.Ba_realdiv_p4,
   m.Ba_realdiv_p5, m.Ba_realdiv_p6, m.Ba_realdiv_p7, m.Ba_realdiv_p8)

####Fu ~ realdiv, both weeks:  ####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Fu_realdiv_p <- brm(
  bf(Fu_per100g ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd + treatment + week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good 

summary(m.Fu_realdiv_p)
pp_check(m.Fu_realdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Fu_realdiv_p2 <- update(m.Fu_realdiv_p, 
                          bf(Fu_per100g ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)

#remove realdivLogStd*week
m.Fu_realdiv_p31 <- update(m.Fu_realdiv_p, 
                           bf(Fu_per100g ~ realdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove treatment*week 
m.Fu_realdiv_p32 <- update(m.Fu_realdiv_p, 
                           bf(Fu_per100g ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot),
                              hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove realdivLogStd*week and treatment*week 
m.Fu_realdiv_p4 <- update(m.Fu_realdiv_p, 
                          bf(Fu_per100g ~ realdivLogStd*treatment + week + (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Fu_realdiv_p4, prob=0.9)

#remove week
m.Fu_realdiv_p5 <- update(m.Fu_realdiv_p, 
                          bf(Fu_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Fu_realdiv_p5, prob=0.9)


#remove hu~treatment
m.Fu_realdiv_p6 <- update(m.Fu_realdiv_p, 
                          bf(Fu_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd + week + (1|block/plot)), 
                          seed = SEED)
#remove hu~realdiv
m.Fu_realdiv_p7 <- update(m.Fu_realdiv_p, 
                          bf(Fu_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ week + (1|block/plot)), 
                          seed = SEED)
summary(m.Fu_realdiv_p7, prob=0.9)

#remove hu~week
m.Fu_realdiv_p8 <- update(m.Fu_realdiv_p, 
                          bf(Fu_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~1), 
                          seed = SEED)
summary(m.Fu_realdiv_p8, prob=0.9)


loo.Fu <- loo(m.Fu_realdiv_p, m.Fu_realdiv_p2, m.Fu_realdiv_p31, m.Fu_realdiv_p32, m.Fu_realdiv_p4,
              m.Fu_realdiv_p5, m.Fu_realdiv_p6, m.Fu_realdiv_p7, m.Fu_realdiv_p8 )

loo.Fu

save(m.Fu_realdiv_p, m.Fu_realdiv_p2, m.Fu_realdiv_p31, m.Fu_realdiv_p32, m.Fu_realdiv_p4,
     m.Fu_realdiv_p5, m.Fu_realdiv_p6, m.Fu_realdiv_p7, m.Fu_realdiv_p8,
     file="./statistics/brms/231219_Fu_realdiv_priors.RData")  

rm(m.Fu_realdiv_p, m.Fu_realdiv_p2, m.Fu_realdiv_p31, m.Fu_realdiv_p32, m.Fu_realdiv_p4,
   m.Fu_realdiv_p5, m.Fu_realdiv_p6, m.Fu_realdiv_p7, m.Fu_realdiv_p8)

####Pl ~ realdiv, both weeks:  ####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Pl_realdiv_p <- brm(
  bf(Pl_per100g ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd + treatment + week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good 

summary(m.Pl_realdiv_p)
pp_check(m.Pl_realdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Pl_realdiv_p2 <- update(m.Pl_realdiv_p, 
                          bf(Pl_per100g ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)

#remove realdivLogStd*week
m.Pl_realdiv_p31 <- update(m.Pl_realdiv_p, 
                           bf(Pl_per100g ~ realdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove treatment*week 
m.Pl_realdiv_p32 <- update(m.Pl_realdiv_p, 
                           bf(Pl_per100g ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot),
                              hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove realdivLogStd*week and treatment*week 
m.Pl_realdiv_p4 <- update(m.Pl_realdiv_p, 
                          bf(Pl_per100g ~ realdivLogStd*treatment + week + (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Pl_realdiv_p4, prob=0.9)

#remove week
m.Pl_realdiv_p5 <- update(m.Pl_realdiv_p, 
                          bf(Pl_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Pl_realdiv_p5, prob=0.9)


#remove hu~treatment
m.Pl_realdiv_p6 <- update(m.Pl_realdiv_p, 
                          bf(Pl_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd + week + (1|block/plot)), 
                          seed = SEED)
#remove hu~realdiv
m.Pl_realdiv_p7 <- update(m.Pl_realdiv_p, 
                          bf(Pl_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ week + (1|block/plot)), 
                          seed = SEED)
summary(m.Pl_realdiv_p7, prob=0.9)

#remove hu~week
m.Pl_realdiv_p8 <- update(m.Pl_realdiv_p, 
                          bf(Pl_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~1), 
                          seed = SEED)
summary(m.Pl_realdiv_p8, prob=0.9)


loo.Pl <- loo(m.Pl_realdiv_p, m.Pl_realdiv_p2, m.Pl_realdiv_p31, m.Pl_realdiv_p32, m.Pl_realdiv_p4,
              m.Pl_realdiv_p5, m.Pl_realdiv_p6, m.Pl_realdiv_p7, m.Pl_realdiv_p8 )
loo.Pl

save(m.Pl_realdiv_p, m.Pl_realdiv_p2, m.Pl_realdiv_p31, m.Pl_realdiv_p32, m.Pl_realdiv_p4,
     m.Pl_realdiv_p5, m.Pl_realdiv_p6, m.Pl_realdiv_p7, m.Pl_realdiv_p8,
     file="./statistics/brms/231219_Pl_realdiv_priors.RData")  

rm(m.Pl_realdiv_p, m.Pl_realdiv_p2, m.Pl_realdiv_p31, m.Pl_realdiv_p32, m.Pl_realdiv_p4,
   m.Pl_realdiv_p5, m.Pl_realdiv_p6, m.Pl_realdiv_p7, m.Pl_realdiv_p8)

####Pr ~ realdiv, both weeks:  ####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Pr_realdiv_p <- brm(
  bf(Pr_per100g ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd + treatment + week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good 

summary(m.Pr_realdiv_p)
pp_check(m.Pr_realdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Pr_realdiv_p2 <- update(m.Pr_realdiv_p, 
                          bf(Pr_per100g ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)

#remove realdivLogStd*week
m.Pr_realdiv_p31 <- update(m.Pr_realdiv_p, 
                           bf(Pr_per100g ~ realdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove treatment*week 
m.Pr_realdiv_p32 <- update(m.Pr_realdiv_p, 
                           bf(Pr_per100g ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot),
                              hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove realdivLogStd*week and treatment*week 
m.Pr_realdiv_p4 <- update(m.Pr_realdiv_p, 
                          bf(Pr_per100g ~ realdivLogStd*treatment + week + (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Pr_realdiv_p4, prob=0.9)

#remove week
m.Pr_realdiv_p5 <- update(m.Pr_realdiv_p, 
                          bf(Pr_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Pr_realdiv_p5, prob=0.9)


#remove hu~treatment
m.Pr_realdiv_p6 <- update(m.Pr_realdiv_p, 
                          bf(Pr_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd + week + (1|block/plot)), 
                          seed = SEED)
#remove hu~realdiv
m.Pr_realdiv_p7 <- update(m.Pr_realdiv_p, 
                          bf(Pr_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ week + (1|block/plot)), 
                          seed = SEED)
summary(m.Pr_realdiv_p7, prob=0.9)

#remove hu~week
m.Pr_realdiv_p8 <- update(m.Pr_realdiv_p, 
                          bf(Pr_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~1), 
                          seed = SEED)
summary(m.Pr_realdiv_p8, prob=0.9)


loo.Pr <- loo(m.Pr_realdiv_p, m.Pr_realdiv_p2, m.Pr_realdiv_p31, m.Pr_realdiv_p32, m.Pr_realdiv_p4,
              m.Pr_realdiv_p5, m.Pr_realdiv_p6, m.Pr_realdiv_p7, m.Pr_realdiv_p8 )

loo.Pr

save(m.Pr_realdiv_p, m.Pr_realdiv_p2, m.Pr_realdiv_p31, m.Pr_realdiv_p32, m.Pr_realdiv_p4,
     m.Pr_realdiv_p5, m.Pr_realdiv_p6, m.Pr_realdiv_p7, m.Pr_realdiv_p8,
     file="./statistics/brms/231219_Pr_realdiv_priors.RData")  

rm(m.Pr_realdiv_p, m.Pr_realdiv_p2, m.Pr_realdiv_p31, m.Pr_realdiv_p32, m.Pr_realdiv_p4,
   m.Pr_realdiv_p5, m.Pr_realdiv_p6, m.Pr_realdiv_p7, m.Pr_realdiv_p8)

####Om ~ realdiv, both weeks:  ####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Om_realdiv_p <- brm(
  bf(Om_per100g ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd + treatment + week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good 

summary(m.Om_realdiv_p)
pp_check(m.Om_realdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Om_realdiv_p2 <- update(m.Om_realdiv_p, 
                          bf(Om_per100g ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)

#remove realdivLogStd*week
m.Om_realdiv_p31 <- update(m.Om_realdiv_p, 
                           bf(Om_per100g ~ realdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove treatment*week 
m.Om_realdiv_p32 <- update(m.Om_realdiv_p, 
                           bf(Om_per100g ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot),
                              hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove realdivLogStd*week and treatment*week 
m.Om_realdiv_p4 <- update(m.Om_realdiv_p, 
                          bf(Om_per100g ~ realdivLogStd*treatment + week + (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Om_realdiv_p4, prob=0.9)

#remove week
m.Om_realdiv_p5 <- update(m.Om_realdiv_p, 
                          bf(Om_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Om_realdiv_p5, prob=0.9)


#remove hu~treatment
m.Om_realdiv_p6 <- update(m.Om_realdiv_p, 
                          bf(Om_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd + week + (1|block/plot)), 
                          seed = SEED)
#remove hu~realdiv
m.Om_realdiv_p7 <- update(m.Om_realdiv_p, 
                          bf(Om_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ week + (1|block/plot)), 
                          seed = SEED)
summary(m.Om_realdiv_p7, prob=0.9)

#remove hu~week
m.Om_realdiv_p8 <- update(m.Om_realdiv_p, 
                          bf(Om_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~1), 
                          seed = SEED)
summary(m.Om_realdiv_p8, prob=0.9)


loo.Om <- loo(m.Om_realdiv_p, m.Om_realdiv_p2, m.Om_realdiv_p31, m.Om_realdiv_p32, m.Om_realdiv_p4,
              m.Om_realdiv_p5, m.Om_realdiv_p6, m.Om_realdiv_p7, m.Om_realdiv_p8 )

loo.Om

save(m.Om_realdiv_p, m.Om_realdiv_p2, m.Om_realdiv_p31, m.Om_realdiv_p32, m.Om_realdiv_p4,
     m.Om_realdiv_p5, m.Om_realdiv_p6, m.Om_realdiv_p7, m.Om_realdiv_p8,
     file="./statistics/brms/231219_Om_realdiv_priors.RData")  

rm(m.Om_realdiv_p, m.Om_realdiv_p2, m.Om_realdiv_p31, m.Om_realdiv_p32, m.Om_realdiv_p4,
   m.Om_realdiv_p5, m.Om_realdiv_p6, m.Om_realdiv_p7, m.Om_realdiv_p8)




