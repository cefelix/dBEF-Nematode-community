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
  mutate(funcdivStd = funcdiv - mean(funcdiv) / sd(funcdiv),
         .after = funcdiv)


datW1 <- subset(dat, week=="W1")
datW2 <- subset(dat, week=="W2")

#priors    
beta_coeff_priors <- prior(normal(0,10), class = "b")  

####Ba ~ funcdiv, both weeks ####
#selecting models based on looic,
#see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff

beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Ba_funcdiv_p <- brm(
  bf(Ba_per100g ~ funcdivLogStd*treatment*week + (1|block/plot),
     hu ~ 1 ),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #1 div

pp_check(m.Ba_funcdiv_p, ndraws=100)+
  xlim(0,2000)
summary(m.Ba_funcdiv_p, prob =0.9)

#as different orientation of funcdivLogStd:treatment2:weekW2 and funcdivLogStd:treatment3:weekW2 -->
#check pairwise with emmeans:
emtrends(m.Ba_funcdiv_p, specs = c("treatment", "week"), var="funcdivLogStd") %>%
  summary() #CIs overlap

#remove 3way interaction:
m.Ba_funcdiv_p2 <- update(m.Ba_funcdiv_p, 
                          bf(Ba_per100g ~ funcdivLogStd*treatment + funcdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ 1), 
                          seed = SEED) #2 div
summary(m.Ba_funcdiv_p2, prob =0.9)

emtrends(m.Ba_funcdiv_p2, specs = c("treatment", "week"), var="funcdivLogStd") %>%
  summary() #CIs overlap

#remove funcdivLogStd*week
m.Ba_funcdiv_p31 <- update(m.Ba_funcdiv_p, 
                           bf(Ba_per100g ~ funcdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ 1 + (1|block/plot)), 
                           seed = SEED) #4 div

summary(m.Ba_funcdiv_p31, prob=0.9)
emtrends(m.Ba_funcdiv_p31, specs = c("treatment", "week"), var="funcdivLogStd") %>%
  summary()

#remove treatment*week 
m.Ba_funcdiv_p32 <- update(m.Ba_funcdiv_p, 
                           bf(Ba_per100g ~ funcdivLogStd*treatment + funcdivLogStd*week + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED)
summary(m.Ba_funcdiv_p32, prob=0.9)

emtrends(m.Ba_funcdiv_p32, specs = c("treatment", "week"), var="funcdivLogStd") %>%
  summary()

#remove funcdivLogStd:week
m.Ba_funcdiv_p4 <- update(m.Ba_funcdiv_p, 
                          bf(Ba_per100g ~ funcdivLogStd*treatment + week + (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)
summary(m.Ba_funcdiv_p4, prob=0.9)
pp_check(m.Ba_funcdiv_p4, ndraws=100)+xlim(0,500)

#remove week
m.Ba_funcdiv_p5 <- update(m.Ba_funcdiv_p, 
                          bf(Ba_per100g ~ funcdivLogStd*treatment + (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)
summary(m.Ba_funcdiv_p5, prob=0.9)
pp_check(m.Ba_funcdiv_p5, ndraws=100)+xlim(0,500)

save(m.Ba_funcdiv_p, m.Ba_funcdiv_p2, m.Ba_funcdiv_p31, m.Ba_funcdiv_p32, m.Ba_funcdiv_p4,
     m.Ba_funcdiv_p5,
     file="./statistics/brms/240108_Ba_funcdiv_priors.RData")

rm(m.Ba_funcdiv_p, m.Ba_funcdiv_p2, m.Ba_funcdiv_p31, m.Ba_funcdiv_p32, m.Ba_funcdiv_p4,
   m.Ba_funcdiv_p5)  

####Fu ~ funcdiv, both weeks:  ####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Fu_funcdiv_p <- brm(
  bf(Fu_per100g ~ funcdivLogStd*treatment*week + (1|block/plot),
     hu ~ 1),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good 

summary(m.Fu_funcdiv_p)
pp_check(m.Fu_funcdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Fu_funcdiv_p2 <- update(m.Fu_funcdiv_p, 
                          bf(Fu_per100g ~ funcdivLogStd*treatment + funcdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)

#remove funcdivLogStd*week
m.Fu_funcdiv_p31 <- update(m.Fu_funcdiv_p, 
                           bf(Fu_per100g ~ funcdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED)

#remove treatment*week 
m.Fu_funcdiv_p32 <- update(m.Fu_funcdiv_p, 
                           bf(Fu_per100g ~ funcdivLogStd*treatment + funcdivLogStd*week + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED)

#remove funcdivLogStd*week and treatment*week 
m.Fu_funcdiv_p4 <- update(m.Fu_funcdiv_p, 
                          bf(Fu_per100g ~ funcdivLogStd*treatment + week + (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)
summary(m.Fu_funcdiv_p4, prob=0.9)

#remove week
m.Fu_funcdiv_p5 <- update(m.Fu_funcdiv_p, 
                          bf(Fu_per100g ~ funcdivLogStd*treatment + (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)
summary(m.Fu_funcdiv_p5, prob=0.9)


save(m.Fu_funcdiv_p, m.Fu_funcdiv_p2, m.Fu_funcdiv_p31, m.Fu_funcdiv_p32, m.Fu_funcdiv_p4,
     m.Fu_funcdiv_p5,
     file="./statistics/brms/240109_Fu_funcdiv_priors.RData")  

rm(m.Fu_funcdiv_p, m.Fu_funcdiv_p2, m.Fu_funcdiv_p31, m.Fu_funcdiv_p32, m.Fu_funcdiv_p4,
   m.Fu_funcdiv_p5)

####Pl ~ funcdiv, both weeks:  ####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Pl_funcdiv_p <- brm(
  bf(Pl_per100g ~ funcdivLogStd*treatment*week + (1|block/plot),
     hu ~ 1),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good 

summary(m.Pl_funcdiv_p)
pp_check(m.Pl_funcdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Pl_funcdiv_p2 <- update(m.Pl_funcdiv_p, 
                          bf(Pl_per100g ~ funcdivLogStd*treatment + funcdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)

#remove funcdivLogStd*week
m.Pl_funcdiv_p31 <- update(m.Pl_funcdiv_p, 
                           bf(Pl_per100g ~ funcdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED)

#remove treatment*week 
m.Pl_funcdiv_p32 <- update(m.Pl_funcdiv_p, 
                           bf(Pl_per100g ~ funcdivLogStd*treatment + funcdivLogStd*week + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED)

#remove funcdivLogStd*week and treatment*week 
m.Pl_funcdiv_p4 <- update(m.Pl_funcdiv_p, 
                          bf(Pl_per100g ~ funcdivLogStd*treatment + week + (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)
summary(m.Pl_funcdiv_p4, prob=0.9)

#remove week
m.Pl_funcdiv_p5 <- update(m.Pl_funcdiv_p, 
                          bf(Pl_per100g ~ funcdivLogStd*treatment + (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)
summary(m.Pl_funcdiv_p5, prob=0.9)



save(m.Pl_funcdiv_p, m.Pl_funcdiv_p2, m.Pl_funcdiv_p31, m.Pl_funcdiv_p32, m.Pl_funcdiv_p4,
     m.Pl_funcdiv_p5, file="./statistics/brms/240109_Pl_funcdiv_priors.RData")  

rm(m.Pl_funcdiv_p, m.Pl_funcdiv_p2, m.Pl_funcdiv_p31, m.Pl_funcdiv_p32, m.Pl_funcdiv_p4,
   m.Pl_funcdiv_p5)

####Pr ~ funcdiv, both weeks:  ####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Pr_funcdiv_p <- brm(
  bf(Pr_per100g ~ funcdivLogStd*treatment*week + (1|block/plot),
     hu ~ funcdivLogStd + treatment + week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good 

summary(m.Pr_funcdiv_p)
pp_check(m.Pr_funcdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Pr_funcdiv_p2 <- update(m.Pr_funcdiv_p, 
                          bf(Pr_per100g ~ funcdivLogStd*treatment + funcdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ funcdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)

#remove funcdivLogStd*week
m.Pr_funcdiv_p31 <- update(m.Pr_funcdiv_p, 
                           bf(Pr_per100g ~ funcdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ funcdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove treatment*week 
m.Pr_funcdiv_p32 <- update(m.Pr_funcdiv_p, 
                           bf(Pr_per100g ~ funcdivLogStd*treatment + funcdivLogStd*week + (1|block/plot),
                              hu ~ funcdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove funcdivLogStd*week and treatment*week 
m.Pr_funcdiv_p4 <- update(m.Pr_funcdiv_p, 
                          bf(Pr_per100g ~ funcdivLogStd*treatment + week + (1|block/plot),
                             hu ~ funcdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Pr_funcdiv_p4, prob=0.9)

#remove week
m.Pr_funcdiv_p5 <- update(m.Pr_funcdiv_p, 
                          bf(Pr_per100g ~ funcdivLogStd*treatment + (1|block/plot),
                             hu ~ funcdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Pr_funcdiv_p5, prob=0.9)


#remove hu~treatment
m.Pr_funcdiv_p6 <- update(m.Pr_funcdiv_p, 
                          bf(Pr_per100g ~ funcdivLogStd*treatment + (1|block/plot),
                             hu ~ funcdivLogStd + week + (1|block/plot)), 
                          seed = SEED)
#remove hu~funcdiv
m.Pr_funcdiv_p7 <- update(m.Pr_funcdiv_p, 
                          bf(Pr_per100g ~ funcdivLogStd*treatment + (1|block/plot),
                             hu ~ week + (1|block/plot)), 
                          seed = SEED)
summary(m.Pr_funcdiv_p7, prob=0.9)

#remove hu~week
m.Pr_funcdiv_p8 <- update(m.Pr_funcdiv_p, 
                          bf(Pr_per100g ~ funcdivLogStd*treatment + (1|block/plot),
                             hu ~1), 
                          seed = SEED)
summary(m.Pr_funcdiv_p8, prob=0.9)


loo.Pr <- loo(m.Pr_funcdiv_p, m.Pr_funcdiv_p2, m.Pr_funcdiv_p31, m.Pr_funcdiv_p32, m.Pr_funcdiv_p4,
              m.Pr_funcdiv_p5, m.Pr_funcdiv_p6, m.Pr_funcdiv_p7, m.Pr_funcdiv_p8 )

loo.Pr

save(m.Pr_funcdiv_p, m.Pr_funcdiv_p2, m.Pr_funcdiv_p31, m.Pr_funcdiv_p32, m.Pr_funcdiv_p4,
     m.Pr_funcdiv_p5, m.Pr_funcdiv_p6, m.Pr_funcdiv_p7, m.Pr_funcdiv_p8,
     file="./statistics/brms/240109_Pr_funcdiv_priors.RData")  

rm(m.Pr_funcdiv_p, m.Pr_funcdiv_p2, m.Pr_funcdiv_p31, m.Pr_funcdiv_p32, m.Pr_funcdiv_p4,
   m.Pr_funcdiv_p5, m.Pr_funcdiv_p6, m.Pr_funcdiv_p7, m.Pr_funcdiv_p8)

####Om ~ funcdiv, both weeks:  ####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Om_funcdiv_p <- brm(
  bf(Om_per100g ~ funcdivLogStd*treatment*week + (1|block/plot),
     hu ~ funcdivLogStd + treatment + week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good 

summary(m.Om_funcdiv_p)
pp_check(m.Om_funcdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Om_funcdiv_p2 <- update(m.Om_funcdiv_p, 
                          bf(Om_per100g ~ funcdivLogStd*treatment + funcdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ funcdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)

#remove funcdivLogStd*week
m.Om_funcdiv_p31 <- update(m.Om_funcdiv_p, 
                           bf(Om_per100g ~ funcdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ funcdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove treatment*week 
m.Om_funcdiv_p32 <- update(m.Om_funcdiv_p, 
                           bf(Om_per100g ~ funcdivLogStd*treatment + funcdivLogStd*week + (1|block/plot),
                              hu ~ funcdivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)

#remove funcdivLogStd*week and treatment*week 
m.Om_funcdiv_p4 <- update(m.Om_funcdiv_p, 
                          bf(Om_per100g ~ funcdivLogStd*treatment + week + (1|block/plot),
                             hu ~ funcdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Om_funcdiv_p4, prob=0.9)

#remove week
m.Om_funcdiv_p5 <- update(m.Om_funcdiv_p, 
                          bf(Om_per100g ~ funcdivLogStd*treatment + (1|block/plot),
                             hu ~ funcdivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED)
summary(m.Om_funcdiv_p5, prob=0.9)


#remove hu~treatment
m.Om_funcdiv_p6 <- update(m.Om_funcdiv_p, 
                          bf(Om_per100g ~ funcdivLogStd*treatment + (1|block/plot),
                             hu ~ funcdivLogStd + week + (1|block/plot)), 
                          seed = SEED)
#remove hu~funcdiv
m.Om_funcdiv_p7 <- update(m.Om_funcdiv_p, 
                          bf(Om_per100g ~ funcdivLogStd*treatment + (1|block/plot),
                             hu ~ week + (1|block/plot)), 
                          seed = SEED)
summary(m.Om_funcdiv_p7, prob=0.9)

#remove hu~week
m.Om_funcdiv_p8 <- update(m.Om_funcdiv_p, 
                          bf(Om_per100g ~ funcdivLogStd*treatment + (1|block/plot),
                             hu ~1), 
                          seed = SEED)
summary(m.Om_funcdiv_p8, prob=0.9)

save(m.Om_funcdiv_p, m.Om_funcdiv_p2, m.Om_funcdiv_p31, m.Om_funcdiv_p32, m.Om_funcdiv_p4,
     m.Om_funcdiv_p5, m.Om_funcdiv_p6, m.Om_funcdiv_p7, m.Om_funcdiv_p8,
     file="./statistics/brms/240109_Om_funcdiv_priors.RData")  

rm(m.Om_funcdiv_p, m.Om_funcdiv_p2, m.Om_funcdiv_p31, m.Om_funcdiv_p32, m.Om_funcdiv_p4,
   m.Om_funcdiv_p5, m.Om_funcdiv_p6, m.Om_funcdiv_p7, m.Om_funcdiv_p8)


#### model selection ####
library(brms)
library(ggplot2)

#selection criteria: picking the most parsimonious model which has an elpd_diff > -4 to the model with the best fit
#saving the selected models in a seperate file:

#Ba ~ funcdiv
load("./statistics/brms/240108_Ba_funcdiv_priors.RData")
loo.Ba <- loo(m.Ba_funcdiv_p, m.Ba_funcdiv_p2, m.Ba_funcdiv_p31, m.Ba_funcdiv_p32, m.Ba_funcdiv_p4,
              m.Ba_funcdiv_p5)
loo.Ba

rm(m.Ba_funcdiv_p, m.Ba_funcdiv_p2, m.Ba_funcdiv_p31, m.Ba_funcdiv_p32, m.Ba_funcdiv_p4) 

#Fu ~ realdviv
load("./statistics/brms/231219_Fu_funcdiv_priors.RData")  
loo.Fu <- loo(m.Fu_funcdiv_p, m.Fu_funcdiv_p2, m.Fu_funcdiv_p31, m.Fu_funcdiv_p32, m.Fu_funcdiv_p4,
              m.Fu_funcdiv_p5)
loo.Fu

rm(m.Fu_funcdiv_p, m.Fu_funcdiv_p2, m.Fu_funcdiv_p31, m.Fu_funcdiv_p32, m.Fu_funcdiv_p4)

#Pl ~ funcdiv
load("./statistics/brms/231219_Pl_funcdiv_priors.RData")  
loo.Pl <- loo(m.Pl_funcdiv_p, m.Pl_funcdiv_p2, m.Pl_funcdiv_p31, m.Pl_funcdiv_p32, m.Pl_funcdiv_p4,
              m.Pl_funcdiv_p5 )
loo.Pl

rm(m.Pl_funcdiv_p, m.Pl_funcdiv_p2, m.Pl_funcdiv_p31, m.Pl_funcdiv_p32, m.Pl_funcdiv_p4)

#Pr ~ funcdiv
load("./statistics/brms/231219_Pr_funcdiv_priors.RData")  
loo.Pr <- loo(m.Pr_funcdiv_p, m.Pr_funcdiv_p2, m.Pr_funcdiv_p31, m.Pr_funcdiv_p32, m.Pr_funcdiv_p4,
              m.Pr_funcdiv_p5, m.Pr_funcdiv_p6, m.Pr_funcdiv_p7, m.Pr_funcdiv_p8 )
loo.Pr

rm(m.Pr_funcdiv_p, m.Pr_funcdiv_p2, m.Pr_funcdiv_p31, m.Pr_funcdiv_p32, m.Pr_funcdiv_p4,
   m.Pr_funcdiv_p5, m.Pr_funcdiv_p6, m.Pr_funcdiv_p8)

#Om ~ funcdiv 
load("./statistics/brms/231219_Om_funcdiv_priors.RData")  

loo.Om <- loo(m.Om_funcdiv_p, m.Om_funcdiv_p2, m.Om_funcdiv_p31, m.Om_funcdiv_p32, m.Om_funcdiv_p4,
              m.Om_funcdiv_p5, m.Om_funcdiv_p6, m.Om_funcdiv_p7, m.Om_funcdiv_p8 )
loo.Om

rm(m.Om_funcdiv_p, m.Om_funcdiv_p2, m.Om_funcdiv_p31, m.Om_funcdiv_p32, m.Om_funcdiv_p4,
   m.Om_funcdiv_p5, m.Om_funcdiv_p6, m.Om_funcdiv_p8)

#save the best fit models:
save(m.Ba_funcdiv_p5, m.Fu_funcdiv_p5, m.Pl_funcdiv_p5, m.Om_funcdiv_p7, m.Pr_funcdiv_p7,
     file = "./statistics/brms/240109_TrophDens_funcdiv_mselect.RData")
