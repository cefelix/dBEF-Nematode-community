####
#assessing the enrichment index
library(brms)
library(dplyr)
library(ggplot2)


#### getting our data in shape ####
#excluding 60 species control plots
dat <- subset(dBEF_nem21, sowndiv != 60)


#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog) %>%
  mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
         .after = realdivLog) %>%
  mutate(funcdivStd = funcdiv - mean(funcdiv) / sd(funcdiv),
         .after = funcdiv)

#priors    
beta_coeff_priors <- prior(normal(0,10), class = "b") 



#### EI ~ sowndiv, family = gaussian #####
load("./statistics/brms/240209_EI_sowndiv.RData")

SEED = 22061996

#using a gausian dis  
m.EI.sowndiv_gaus_p <- brm(
  bf(EI ~ sowndivLogStd*treatment*week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = gaussian,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #

pp_check(m.EI.sowndiv_gaus_p, ndraws=100)
pp_check(m.EI.sowndiv_gaus_p, type="stat")
summary(m.EI.sowndiv_gaus_p, prob =0.9)

#remove 3 way interaction:
m.EI.sowndiv_gaus_p2 <- update(m.EI.sowndiv_gaus_p,
                               bf(EI ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + (1|block/plot)),
                               seed = SEED) #7 div
pp_check(m.EI.sowndiv_gaus_p2, ndraws=100)
summary(m.EI.sowndiv_gaus_p2, prob =0.9)

#remove treatment*week:
m.EI.sowndiv_gaus_p31 <- update(m.EI.sowndiv_gaus_p,
                                bf(EI ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot)),
                                seed = SEED) #2 div
pp_check(m.EI.sowndiv_gaus_p31, ndraws=100)
summary(m.EI.sowndiv_gaus_p31, prob =0.9)

#remove sowndivLogStd*week:
m.EI.sowndiv_gaus_p32 <- update(m.EI.sowndiv_gaus_p,
                                bf(EI ~ sowndivLogStd*treatment + treatment*week + (1|block/plot)),
                                seed = SEED) #2 div
pp_check(m.EI.sowndiv_gaus_p32, ndraws=100)
summary(m.EI.sowndiv_gaus_p32, prob =0.9)

#remove sowndivLogStd*week and streatment*week:
m.EI.sowndiv_gaus_p4 <- update(m.EI.sowndiv_gaus_p,
                               bf(EI ~ sowndivLogStd*treatment + week + (1|block/plot)),
                               seed = SEED) #4 div
pp_check(m.EI.sowndiv_gaus_p4, ndraws=100)
summary(m.EI.sowndiv_gaus_p4, prob =0.9) 

#remove week
m.EI.sowndiv_gaus_p5 <- update(m.EI.sowndiv_gaus_p,
                               bf(EI ~ sowndivLogStd*treatment + (1|block/plot)),
                               seed = SEED) #1 div
pp_check(m.EI.sowndiv_gaus_p5, ndraws=100)
summary(m.EI.sowndiv_gaus_p5, prob =0.9)
conditional_effects(m.EI.sowndiv_gaus_p5)


loo_EI_sowndiv <- loo( m.EI.sowndiv_gaus_p,  m.EI.sowndiv_gaus_p2,  m.EI.sowndiv_gaus_p31,  m.EI.sowndiv_gaus_p32,  
                       m.EI.sowndiv_gaus_p4,  m.EI.sowndiv_gaus_p5)
loo_EI_sowndiv #p5 wins

#save the models
save(m.EI.sowndiv_gaus_p, m.EI.sowndiv_gaus_p2, m.EI.sowndiv_gaus_p31, m.EI.sowndiv_gaus_p32,
     m.EI.sowndiv_gaus_p4, m.EI.sowndiv_gaus_p5,
     file = "./statistics/brms/240209_EI_sowndiv.RData")



#EI ~ sowndiv, family = "skew_normal" 
    #using a skew-normal distribution  
    m.EI.sowndiv_skewn_p <- brm(
      bf(EI ~ sowndivLogStd*treatment*week + (1|block/plot)),
      data = dat, 
      prior = beta_coeff_priors,
      family = "skew_normal",
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #
    
    pp_check(m.EI.sowndiv_skewn_p, ndraws=100) #not better than gaussian
    pp_check(m.EI.sowndiv_skewn_p, type="stat")
    summary(m.EI.sowndiv_skewn_p, prob =0.9)


#### EI ~ realdiv, family = gaussian #####
load("./statistics/brms/240209_EI_realdiv.RData")

SEED = 22061996

#using a gausian dis  
m.EI.realdiv_gaus_p <- brm(
  bf(EI ~ realdivLogStd*treatment*week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = gaussian,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #

pp_check(m.EI.realdiv_gaus_p, ndraws=100)
pp_check(m.EI.realdiv_gaus_p, type="stat")
summary(m.EI.realdiv_gaus_p, prob =0.9)

#remove 3 way interaction:
m.EI.realdiv_gaus_p2 <- update(m.EI.realdiv_gaus_p,
                               bf(EI ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + (1|block/plot)),
                               seed = SEED) #7 div
pp_check(m.EI.realdiv_gaus_p2, ndraws=100)
summary(m.EI.realdiv_gaus_p2, prob =0.9)

#remove treatment*week:
m.EI.realdiv_gaus_p31 <- update(m.EI.realdiv_gaus_p,
                                bf(EI ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot)),
                                seed = SEED) #2 div
pp_check(m.EI.realdiv_gaus_p31, ndraws=100)
summary(m.EI.realdiv_gaus_p31, prob =0.9)

#remove realdivLogStd*week:
m.EI.realdiv_gaus_p32 <- update(m.EI.realdiv_gaus_p,
                                bf(EI ~ realdivLogStd*treatment + treatment*week + (1|block/plot)),
                                seed = SEED) #2 div
pp_check(m.EI.realdiv_gaus_p32, ndraws=100)
summary(m.EI.realdiv_gaus_p32, prob =0.9)

#remove realdivLogStd*week and streatment*week:
m.EI.realdiv_gaus_p4 <- update(m.EI.realdiv_gaus_p,
                               bf(EI ~ realdivLogStd*treatment + week + (1|block/plot)),
                               seed = SEED) #4 div
pp_check(m.EI.realdiv_gaus_p4, ndraws=100)
summary(m.EI.realdiv_gaus_p4, prob =0.9) 

#remove week
m.EI.realdiv_gaus_p5 <- update(m.EI.realdiv_gaus_p,
                               bf(EI ~ realdivLogStd*treatment + (1|block/plot)),
                               seed = SEED) #1 div
pp_check(m.EI.realdiv_gaus_p5, ndraws=100)
summary(m.EI.realdiv_gaus_p5, prob =0.9)
conditional_effects(m.EI.realdiv_gaus_p5)

#save the models
save(m.EI.realdiv_gaus_p, m.EI.realdiv_gaus_p2, m.EI.realdiv_gaus_p31, m.EI.realdiv_gaus_p32,
     m.EI.realdiv_gaus_p4, m.EI.realdiv_gaus_p5,
     file = "./statistics/brms/240209_EI_realdiv.RData")

#### save selected models ####
loo.EI.sown <- loo(m.EI.sowndiv_gaus_p, m.EI.sowndiv_gaus_p2, m.EI.sowndiv_gaus_p31, m.EI.sowndiv_gaus_p32,
                   m.EI.sowndiv_gaus_p4, m.EI.sowndiv_gaus_p5)
loo.EI.sown   #p5

loo.EI.real <- loo(m.EI.realdiv_gaus_p, m.EI.realdiv_gaus_p2, m.EI.realdiv_gaus_p31, m.EI.realdiv_gaus_p32,
    m.EI.realdiv_gaus_p4, m.EI.realdiv_gaus_p5) #p5
loo.EI.real   #p5

save(m.EI.sowndiv_gaus_p5, m.EI.realdiv_gaus_p5,
     file = "./statistics/brms/240215_EI_mselect.RData")

