library(brms)
library(rstan)
library(ggplot2)
library(emmeans)
library(bayestestR)
library(dplyr)

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

####Ba ~ realdiv, both weeks: ####
#stepwise elimination: p6, as realdiv:week is marginally significant (90%CI)
#looic: p5, as it the most parsimonous with a elpd_diff less 4 to the best elpd model 

#see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff
load("./statistics/brms/240108_Ba_realdiv_priors.RData")

beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996
sum(subset(dat, week=="W1")$Ba_per100g == 0) #9
sum(subset(dat, week=="W2")$Ba_per100g == 0) #2

#for both weeks  
m.Ba_realdiv_p <- brm(
  bf(Ba_per100g ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ 1 ),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #1 div

pp_check(m.Ba_realdiv_p, ndraws=100)+
  xlim(0,2000)
summary(m.Ba_realdiv_p, prob =0.9)

#as different orientation of realdivLogStd:treatment2:weekW2 and realdivLogStd:treatment3:weekW2 -->
#check pairwise with emmeans:
emtrends(m.Ba_realdiv_p, specs = c("treatment", "week"), var="realdivLogStd") %>%
  summary() #CIs overlap

#remove 3way interaction:
m.Ba_realdiv_p2 <- update(m.Ba_realdiv_p, 
                          bf(Ba_per100g ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ 1), 
                          seed = SEED) #2 div
summary(m.Ba_realdiv_p2, prob =0.9)

emtrends(m.Ba_realdiv_p2, specs = c("treatment", "week"), var="realdivLogStd") %>%
  summary() #CIs overlap

#remove realdivLogStd*week
m.Ba_realdiv_p31 <- update(m.Ba_realdiv_p, 
                           bf(Ba_per100g ~ realdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ 1 + (1|block/plot)), 
                           seed = SEED) #4 div

summary(m.Ba_realdiv_p31, prob=0.9)
emtrends(m.Ba_realdiv_p31, specs = c("treatment", "week"), var="realdivLogStd") %>%
  summary()

#remove treatment*week 
m.Ba_realdiv_p32 <- update(m.Ba_realdiv_p, 
                           bf(Ba_per100g ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED)
summary(m.Ba_realdiv_p32, prob=0.9)

emtrends(m.Ba_realdiv_p32, specs = c("treatment", "week"), var="realdivLogStd") %>%
  summary()

#remove realdivLogStd:week
m.Ba_realdiv_p4 <- update(m.Ba_realdiv_p, 
                          bf(Ba_per100g ~ realdivLogStd*treatment + week + (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)
summary(m.Ba_realdiv_p4, prob=0.9)
pp_check(m.Ba_realdiv_p4, ndraws=100)+xlim(0,500)

#remove week
m.Ba_realdiv_p5 <- update(m.Ba_realdiv_p, 
                           bf(Ba_per100g ~ realdivLogStd*treatment + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED)
summary(m.Ba_realdiv_p5, prob=0.9)
pp_check(m.Ba_realdiv_p5, ndraws=100)+xlim(0,500)

#compare models
save(m.Ba_realdiv_p, m.Ba_realdiv_p2, m.Ba_realdiv_p31, m.Ba_realdiv_p32, m.Ba_realdiv_p4,
     m.Ba_realdiv_p5,
     file="./statistics/brms/240108_Ba_realdiv_priors.RData")

rm(m.Ba_realdiv_p, m.Ba_realdiv_p2, m.Ba_realdiv_p31, m.Ba_realdiv_p32, m.Ba_realdiv_p4,
   m.Ba_realdiv_p5)

####Fu ~ realdiv, both weeks: p4 (looic p5) ####
#stepwise elimination of non-significant terms: p4, as week is marginally significant
#looic: p5, as it is the most parsimonious model and all elpd_diff lay in range of 2 SE_diff

load("./statistics/brms/231219_Fu_realdiv_priors.RData")
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(subset(dat, week=="W1")$Fu_per100g == 0) #3
sum(subset(dat, week=="W2")$Fu_per100g == 0) #1
#use hu~1, as too little zeros to estimate anything

#for both weeks  
m.Fu_realdiv_p <- brm(
  bf(Fu_per100g ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ 1),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #1 div

summary(m.Fu_realdiv_p, prob =0.9)
pp_check(m.Fu_realdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Fu_realdiv_p2 <- update(m.Fu_realdiv_p, 
                          bf(Fu_per100g ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ 1), 
                          seed = SEED) #4 div
summary(m.Fu_realdiv_p2, prob =0.9)

#remove realdivLogStd*week
m.Fu_realdiv_p31 <- update(m.Fu_realdiv_p, 
                           bf(Fu_per100g ~ realdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED) # 1 div
summary(m.Fu_realdiv_p31, prob =0.9)

#remove treatment*week 
m.Fu_realdiv_p32 <- update(m.Fu_realdiv_p, 
                           bf(Fu_per100g ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED)
summary(m.Fu_realdiv_p32, prob =0.9)


#remove realdivLogStd*week and treatment*week 
m.Fu_realdiv_p4 <- update(m.Fu_realdiv_p, 
                          bf(Fu_per100g ~ realdivLogStd*treatment + week + (1|block/plot),
                             hu ~ 1), 
                          seed = SEED) #3 div
summary(m.Fu_realdiv_p4, prob=0.9) #stepwise: keep week, as it is marginally significant at 90% CI

#remove week
m.Fu_realdiv_p5 <- update(m.Fu_realdiv_p, 
                          bf(Fu_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)
summary(m.Fu_realdiv_p5, prob=0.9)

save(m.Fu_realdiv_p, m.Fu_realdiv_p2, m.Fu_realdiv_p31, m.Fu_realdiv_p32, m.Fu_realdiv_p4,
     m.Fu_realdiv_p5,
     file="./statistics/brms/231219_Fu_realdiv_priors.RData")  

rm(m.Fu_realdiv_p, m.Fu_realdiv_p2, m.Fu_realdiv_p31, m.Fu_realdiv_p32,
   m.Fu_realdiv_p4, m.Fu_realdiv_p5)

####Pl ~ realdiv, both weeks: p4 (looic p4) ####
#stepwise elimination of non-significant terms: p4, as week is significant
#looic: p4, as it is the most parsimonous model with a elpd_se_diff of less than 2 SE
sum(subset(dat, week=="W1")$Pl_per100g == 0) #2
sum(subset(dat, week=="W2")$Pl_per100g == 0) #0

load("./statistics/brms/231219_Pl_realdiv_priors.RData")
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Pl_realdiv_p <- brm(
  bf(Pl_per100g ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ 1),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #2 div
summary(m.Pl_realdiv_p, prob =0.9)

pp_check(m.Pl_realdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Pl_realdiv_p2 <- update(m.Pl_realdiv_p, 
                          bf(Pl_per100g ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)
summary(m.Pl_realdiv_p2, prob =0.9)


#remove realdivLogStd*week
m.Pl_realdiv_p31 <- update(m.Pl_realdiv_p, 
                           bf(Pl_per100g ~ realdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED)
summary(m.Pl_realdiv_p31, prob =0.9)

#remove treatment*week 
m.Pl_realdiv_p32 <- update(m.Pl_realdiv_p, 
                           bf(Pl_per100g ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED)
summary(m.Pl_realdiv_p32, prob =0.9)

#remove realdivLogStd*week and treatment*week 
m.Pl_realdiv_p4 <- update(m.Pl_realdiv_p, 
                          bf(Pl_per100g ~ realdivLogStd*treatment + week + (1|block/plot),
                             hu ~ 1), 
                          seed = SEED) #26 div
summary(m.Pl_realdiv_p4, prob=0.95)
#stepwise simplification; keep week as it is significant when considering a 95% CI

#remove week
m.Pl_realdiv_p5 <- update(m.Pl_realdiv_p, 
                          bf(Pl_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ 1), 
                          seed = SEED)
summary(m.Pl_realdiv_p5, prob=0.9)


save(m.Pl_realdiv_p, m.Pl_realdiv_p2, m.Pl_realdiv_p31, m.Pl_realdiv_p32, m.Pl_realdiv_p4,
     m.Pl_realdiv_p5,
     file="./statistics/brms/231219_Pl_realdiv_priors.RData")  

rm(m.Pl_realdiv_p, m.Pl_realdiv_p2, m.Pl_realdiv_p31, m.Pl_realdiv_p32, 
   m.Pl_realdiv_p4, m.Pl_realdiv_p5)

####Pr ~ realdiv, both weeks: p6 (looic p8), simplification CHECKED####
#looic p8, as most parsimonious with elpd not differing significantly (less than 2 SE)
#stepwise removal p6, as hu~term is marginally significant 
sum(subset(dat, week=="W1")$Pr_per100g == 0) #44
sum(subset(dat, week=="W2")$Pr_per100g == 0) #13

load("./statistics/brms/231219_Pr_realdiv_priors.RData")
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996


#for both weeks  
m.Pr_realdiv_p <- brm(
  bf(Pr_per100g ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd*treatment + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good 

summary(m.Pr_realdiv_p, prob=0.9)
pp_check(m.Pr_realdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Pr_realdiv_p2 <- update(m.Pr_realdiv_p, 
                          bf(Pr_per100g ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ realdivLogStd*treatment + (1|block/plot)), 
                          seed = SEED) #3 div

summary(m.Pr_realdiv_p2, prob=0.9)

#remove realdivLogStd*week
m.Pr_realdiv_p31 <- update(m.Pr_realdiv_p, 
                           bf(Pr_per100g ~ realdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ realdivLogStd*treatment + (1|block/plot)), 
                           seed = SEED) #1 div
summary(m.Pr_realdiv_p31, prob=0.9)

#remove treatment*week 
m.Pr_realdiv_p32 <- update(m.Pr_realdiv_p, 
                           bf(Pr_per100g ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot),
                              hu ~ realdivLogStd*treatment + (1|block/plot)), 
                           seed = SEED)
summary(m.Pr_realdiv_p32, prob=0.9)

#remove realdivLogStd*week and treatment*week 
m.Pr_realdiv_p4 <- update(m.Pr_realdiv_p, 
                          bf(Pr_per100g ~ realdivLogStd*treatment + week + (1|block/plot),
                             hu ~ realdivLogStd*treatment + (1|block/plot)), 
                          seed = SEED) #7 div
summary(m.Pr_realdiv_p4, prob=0.9)

#remove week
m.Pr_realdiv_p5 <- update(m.Pr_realdiv_p, 
                          bf(Pr_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd*treatment + (1|block/plot)), 
                          seed = SEED)
summary(m.Pr_realdiv_p5, prob=0.9)


save(m.Pr_realdiv_p, m.Pr_realdiv_p2, m.Pr_realdiv_p31, m.Pr_realdiv_p32, m.Pr_realdiv_p4,
     m.Pr_realdiv_p5,
     file="./statistics/brms/240221_Pr_realdiv_priors.RData")  

rm(m.Pr_realdiv_p, m.Pr_realdiv_p2, m.Pr_realdiv_p31, m.Pr_realdiv_p32, m.Pr_realdiv_p4,
   m.Pr_realdiv_p5, m.Pr_realdiv_p7, m.Pr_realdiv_p8)

#### Om ~ realdiv ####
#looic: p7 is most parsimonous with a elpd that is not significantly worse (more than 2 SE elpd diff), while elpd_diff < 4
#stepwise termremoval at 90% CI: p7, as hu~treatment is not marginally significant, but hu~week is significant (95% CI)
sum(subset(dat, week=="W1")$Om_per100g == 0) #89
sum(subset(dat, week=="W2")$Om_per100g == 0) #42


load("./statistics/brms/231219_Om_realdiv_priors.RData")
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

#for both weeks  
m.Om_realdiv_p <- brm(
  bf(Om_per100g ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd*treatment + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #11 div

summary(m.Om_realdiv_p)
pp_check(m.Om_realdiv_p, ndraws=100)+
  xlim(0,2000)

#remove 3way interaction:
m.Om_realdiv_p2 <- update(m.Om_realdiv_p, 
                          bf(Om_per100g ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~ realdivLogStd*treatment + (1|block/plot)), 
                          seed = SEED) #3 div
summary(m.Om_realdiv_p2)

#remove realdivLogStd*week
m.Om_realdiv_p31 <- update(m.Om_realdiv_p, 
                           bf(Om_per100g ~ realdivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~ realdivLogStd*treatment + (1|block/plot)), 
                           seed = SEED)
summary(m.Om_realdiv_p31)

#remove treatment*week 
m.Om_realdiv_p32 <- update(m.Om_realdiv_p, 
                           bf(Om_per100g ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot),
                              hu ~ realdivLogStd*treatmentk + (1|block/plot)), 
                           seed = SEED)
summary(m.Om_realdiv_p32)

#remove realdivLogStd*week and treatment*week 
m.Om_realdiv_p4 <- update(m.Om_realdiv_p, 
                          bf(Om_per100g ~ realdivLogStd*treatment + week + (1|block/plot),
                             hu ~ realdivLogStd*treatment + (1|block/plot)), 
                          seed = SEED) #7 div
summary(m.Om_realdiv_p4, prob=0.9)

#remove week
m.Om_realdiv_p5 <- update(m.Om_realdiv_p, 
                          bf(Om_per100g ~ realdivLogStd*treatment + (1|block/plot),
                             hu ~ realdivLogStd*treatment + (1|block/plot)), 
                          seed = SEED) #1 div
summary(m.Om_realdiv_p5, prob=0.9)


save(m.Om_realdiv_p, m.Om_realdiv_p2, m.Om_realdiv_p31, m.Om_realdiv_p32, m.Om_realdiv_p4,
     m.Om_realdiv_p5,
     file="./statistics/brms/231219_Om_realdiv_priors.RData")  

rm(m.Om_realdiv_p, m.Om_realdiv_p2, m.Om_realdiv_p31, m.Om_realdiv_p32, m.Om_realdiv_p4,
   m.Om_realdiv_p5)



#### model selection ####
library(brms)
library(ggplot2)

#selection criteria: picking the most parsimonious model which has an elpd_diff > -4 to the model with the best fit
#saving the selected models in a seperate file:

#Ba ~ realdiv
  load("./statistics/brms/240108_Ba_realdiv_priors.RData")
  loo.Ba <- loo(m.Ba_realdiv_p, m.Ba_realdiv_p2, m.Ba_realdiv_p31, m.Ba_realdiv_p32, m.Ba_realdiv_p4,
                m.Ba_realdiv_p5)
  loo.Ba
  
  rm(m.Ba_realdiv_p, m.Ba_realdiv_p2, m.Ba_realdiv_p31, m.Ba_realdiv_p32, m.Ba_realdiv_p4) 

#Fu ~ realdviv
  load("./statistics/brms/231219_Fu_realdiv_priors.RData")  
  loo.Fu <- loo(m.Fu_realdiv_p, m.Fu_realdiv_p2, m.Fu_realdiv_p31, m.Fu_realdiv_p32, m.Fu_realdiv_p4,
                m.Fu_realdiv_p5)
  loo.Fu
  
  rm(m.Fu_realdiv_p, m.Fu_realdiv_p2, m.Fu_realdiv_p31, m.Fu_realdiv_p32, m.Fu_realdiv_p4)

#Pl ~ realdiv
  load("./statistics/brms/231219_Pl_realdiv_priors.RData")  
  loo.Pl <- loo(m.Pl_realdiv_p, m.Pl_realdiv_p2, m.Pl_realdiv_p31, m.Pl_realdiv_p32, m.Pl_realdiv_p4,
                m.Pl_realdiv_p5 )
  loo.Pl
  
  rm(m.Pl_realdiv_p, m.Pl_realdiv_p2, m.Pl_realdiv_p31, m.Pl_realdiv_p32, m.Pl_realdiv_p4)

#Pr ~ realdiv
  load("./statistics/brms/231219_Pr_realdiv_priors.RData")  
  loo.Pr <- loo(m.Pr_realdiv_p, m.Pr_realdiv_p2, m.Pr_realdiv_p31, m.Pr_realdiv_p32, m.Pr_realdiv_p4,
                m.Pr_realdiv_p5)
  loo.Pr
  
  rm(m.Pr_realdiv_p, m.Pr_realdiv_p2, m.Pr_realdiv_p31, m.Pr_realdiv_p32, m.Pr_realdiv_p4,
     m.Pr_realdiv_p5, m.Pr_realdiv_p6, m.Pr_realdiv_p8)

#Om ~ realdiv 
  load("./statistics/brms/231219_Om_realdiv_priors.RData")  
  
  loo.Om <- loo(m.Om_realdiv_p, m.Om_realdiv_p2, m.Om_realdiv_p31, m.Om_realdiv_p32, m.Om_realdiv_p4,
                m.Om_realdiv_p5)
  loo.Om
  
  rm(m.Om_realdiv_p, m.Om_realdiv_p2, m.Om_realdiv_p31, m.Om_realdiv_p32, m.Om_realdiv_p4,
     m.Om_realdiv_p5, m.Om_realdiv_p6, m.Om_realdiv_p8)
  
#save the best fit models:
  save(m.Ba_realdiv_p5, m.Fu_realdiv_p5, m.Pl_realdiv_p5, m.Om_realdiv_p7, m.Pr_realdiv_p7,
    file = "./statistics/brms/240221_TrophDens_realdiv_mselect.RData")
  