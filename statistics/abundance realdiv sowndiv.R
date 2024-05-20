####loading data and packages ####
# a workflow: https://m-clark.github.io/posts/2021-02-28-practical-bayes-part-i/
# rethinking in brms: https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/interactions.html
library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)

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

#a seed:
SEED = 22061996

#priors    
beta_coeff_priors <- prior(normal(0,10), class = "b")  

#### total_nematodes ~ sowndivLogStd ,fam=negbinomial, -60sp ####
  m.abun_all.sowndiv_p <- brm(total_nematodes ~ sowndivLogStd*treatment*week + offset(log(soilDW)) + (1|block/plot), 
                              data = dat, family = "negbinomial",
                              chains = 3,
                              cores = 3,
                              iter = 4000, warmup = 2000, #(tail ESS too low with iter=2000)
                              seed = SEED,
                              control = list(adapt_delta = 0.99) ) #20 div 

m.abun_all.sowndiv_p %>% pp_check(ndraws=100)
m.abun_all.sowndiv_p %>% summary()

  #remove 3 way interaction:
  m.abun_all.sowndiv_p2 <- update(m.abun_all.sowndiv_p, 
                                  bf(total_nematodes ~ sowndivLogStd*treatment + week*treatment + 
                                       week*sowndivLogStd + offset(log(soilDW)) + (1|block/plot)),
                                  seed = SEED) 
  summary(m.abun_all.sowndiv_p2, prob=0.9)
  
  #remove 2 way interactions with week (1):
  m.abun_all.sowndiv_p31 <- update(m.abun_all.sowndiv_p, 
                                  bf(total_nematodes ~ sowndivLogStd*treatment + 
                                       week*sowndivLogStd + offset(log(soilDW)) + (1|block/plot)),
                                  seed = SEED) #16 div
  
  #remove 2-way interactions with week (2):
  m.abun_all.sowndiv_p32 <- update(m.abun_all.sowndiv_p, 
                                   bf(total_nematodes ~ sowndivLogStd*treatment +  week*treatment + 
                                      offset(log(soilDW)) + (1|block/plot)),
                                   seed = SEED) #37 div
  
  #remove both 2-way interactions with week:
  m.abun_all.sowndiv_p4 <- update(m.abun_all.sowndiv_p, 
                                   bf(total_nematodes ~ sowndivLogStd*treatment + week +
                                        offset(log(soilDW)) + (1|block/plot)),
                                   seed = SEED) #15 div
  
  #remove week:
  m.abun_all.sowndiv_p5 <- update(m.abun_all.sowndiv_p, 
                                  bf(total_nematodes ~ sowndivLogStd*treatment +
                                       offset(log(soilDW)) + (1|block/plot)),
                                  seed = SEED) 
  m.abun_all.sowndiv_p5 %>% pp_check(ndraws=100)
  
  #save models:
  save(m.abun_all.sowndiv_p, m.abun_all.sowndiv_p2, m.abun_all.sowndiv_p31,
       m.abun_all.sowndiv_p32, m.abun_all.sowndiv_p4, m.abun_all.sowndiv_p5,
       file="./statistics/brms/240216_abun_offset_sowndiv.RData")
  
#### total_nematodes ~ realdiv####
  m.abun_all.realdiv_p <- brm(total_nematodes ~ realdivLogStd*treatment*week + offset(log(soilDW)) + (1|block/plot), 
                              data = dat, family = "negbinomial",
                              chains = 3,
                              cores = 3,
                              iter = 4000, warmup = 2000, #(tail ESS too low with iter=2000)
                              seed = SEED,
                              control = list(adapt_delta = 0.99) ) #20 div 
  
  m.abun_all.realdiv_p %>% pp_check(ndraws=100)
  m.abun_all.realdiv_p %>% summary()
  
  #remove 3 way interaction:
  m.abun_all.realdiv_p2 <- update(m.abun_all.realdiv_p, 
                                  bf(total_nematodes ~ realdivLogStd*treatment + week*treatment + 
                                       week*realdivLogStd + offset(log(soilDW)) + (1|block/plot)),
                                  seed = SEED,
                                  iter=6000, warmup=3000) #iter 4000 causes too low bulk/tail ESS 
  summary(m.abun_all.realdiv_p2, prob=0.9)
  
  #remove 2 way interactions with week (1):
  m.abun_all.realdiv_p31 <- update(m.abun_all.realdiv_p, 
                                   bf(total_nematodes ~ realdivLogStd*treatment + 
                                        week*realdivLogStd + offset(log(soilDW)) + (1|block/plot)),
                                   seed = SEED) #10 div
  
  #remove 2-way interactions with week (2):
  m.abun_all.realdiv_p32 <- update(m.abun_all.realdiv_p, 
                                   bf(total_nematodes ~ realdivLogStd*treatment +  week*treatment + 
                                        offset(log(soilDW)) + (1|block/plot)),
                                   seed = SEED) #23 div
  
  #remove both 2-way interactions with week:
  m.abun_all.realdiv_p4 <- update(m.abun_all.realdiv_p, 
                                  bf(total_nematodes ~ realdivLogStd*treatment + week +
                                       offset(log(soilDW)) + (1|block/plot)),
                                  seed = SEED) #35 div
  
  #remove week:
  m.abun_all.realdiv_p5 <- update(m.abun_all.realdiv_p, 
                                  bf(total_nematodes ~ realdivLogStd*treatment +
                                       offset(log(soilDW)) + (1|block/plot)),
                                  seed = SEED) 
  
  #save models:
  save(m.abun_all.realdiv_p, m.abun_all.realdiv_p2, m.abun_all.realdiv_p31,
       m.abun_all.realdiv_p32, m.abun_all.realdiv_p4, m.abun_all.realdiv_p5,
       file="./statistics/brms/240216_abun_offset_realdiv.RData")
  
####dens ~ sowndiv ####
  m.dens_sowndiv_p <- brm(ind_per100g ~ sowndivLogStd*treatment*week + (1|block/plot), 
                              data = dat, family = "lognormal",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000, #(tail ESS too low with iter=2000)
                              seed = SEED,
                              control = list(adapt_delta = 0.99) ) 
  pp_check(m.dens_sowndiv_p)
  
  m.dens_sowndiv_p2 <- update(m.dens_sowndiv_p,
                              bf(ind_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                   treatment*week + (1|block/plot)) ) #12 div
  
  m.dens_sowndiv_p31 <- update(m.dens_sowndiv_p,
                              bf(ind_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot)) ) #9 div
  
  m.dens_sowndiv_p32 <- update(m.dens_sowndiv_p,
                              bf(ind_per100g ~ sowndivLogStd*treatment + treatment*week + (1|block/plot)) ) #10 div
  
  m.dens_sowndiv_p4 <- update(m.dens_sowndiv_p,
                              bf(ind_per100g ~ sowndivLogStd*treatment + week + (1|block/plot)) ) #3 div
  
  m.dens_sowndiv_p5 <- update(m.dens_sowndiv_p,
                              bf(ind_per100g ~ sowndivLogStd*treatment + (1|block/plot)) ) #2 div
 
####dens ~ realdiv ####
  m.dens_realdiv_p <- brm(ind_per100g ~ realdivLogStd*treatment*week + (1|block/plot), 
                          data = dat, family = "lognormal",
                          chains = 3,
                          cores = 3,
                          iter = 2000, warmup = 1000, #(tail ESS too low with iter=2000)
                          seed = SEED,
                          control = list(adapt_delta = 0.99) ) #3 div
  pp_check(m.dens_realdiv_p)
  
  m.dens_realdiv_p2 <- update(m.dens_realdiv_p,
                              bf(ind_per100g ~ realdivLogStd*treatment + realdivLogStd*week + 
                                   treatment*week + (1|block/plot)) ) #28 div
  
  m.dens_realdiv_p31 <- update(m.dens_realdiv_p,
                               bf(ind_per100g ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot)) ) #2 div
  
  m.dens_realdiv_p32 <- update(m.dens_realdiv_p,
                               bf(ind_per100g ~ realdivLogStd*treatment + treatment*week + (1|block/plot)) ) #4 div
  
  m.dens_realdiv_p4 <- update(m.dens_realdiv_p,
                              bf(ind_per100g ~ realdivLogStd*treatment + week + (1|block/plot)) ) #3 div
  
  m.dens_realdiv_p5 <- update(m.dens_realdiv_p,
                              bf(ind_per100g ~ realdivLogStd*treatment + (1|block/plot)) ) #0 div   

#### select best fit based on loo-IC ####
  loo.abun.realdiv <- loo(m.abun_all.realdiv_p, m.abun_all.realdiv_p2, m.abun_all.realdiv_p31,
                          m.abun_all.realdiv_p32, m.abun_all.realdiv_p4, m.abun_all.realdiv_p5)
  loo.abun.realdiv #p5
  summary(m.abun_all.realdiv_p5)
  
  loo.abun.sowndiv <- loo(m.abun_all.sowndiv_p, m.abun_all.sowndiv_p2, m.abun_all.sowndiv_p31,
                          m.abun_all.sowndiv_p32, m.abun_all.sowndiv_p4, m.abun_all.sowndiv_p5)
  loo.abun.sowndiv #p5
  summary(m.abun_all.sowndiv_p5)
  
  loo.dens.sowndiv <- loo( m.dens_sowndiv_p,  m.dens_sowndiv_p2,  m.dens_sowndiv_p32,  
                           m.dens_sowndiv_p31,  m.dens_sowndiv_p4,  m.dens_sowndiv_p5)
  loo.dens.sowndiv #p5
  summary(m.dens_sowndiv_p5)
  
  loo.dens.realdiv <- loo( m.dens_realdiv_p,  m.dens_realdiv_p2,  m.dens_realdiv_p32,  
                           m.dens_realdiv_p31,  m.dens_realdiv_p4,  m.dens_realdiv_p5)
  loo.dens.realdiv #p5
  summary(m.dens_realdiv_p5)
  
  #save best fit models:
  save(m.abun_all.realdiv_p5, m.abun_all.sowndiv_p5, m.dens_sowndiv_p5, m.dens_realdiv_p5,
       file="./statistics/brms/240225_abun_dens_all_mselect.RData")
  

#### OLD (2017 data): negBinom 31 vog: total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), fam=negbinomial, -60sp, year =2017 ####
SEED = 22061996
dat = dBEF_nem17 %>% 
  filter(sowndiv != 60)
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog)


m.abun.nBinomOffS.31vog <- brm(total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), 
                            data = dat, family = "negbinomial",
                            chains = 3,
                            cores = 3,
                            iter = 3000, warmup = 1500,
                            seed = SEED,
                            control = list(adapt_delta = 0.99) ) 
#all good

pp_check(m.abun.nBinomOffS.31vog, ndraws = 100)

#save both nbinom offset models
save(m.abun.nBinomOffS.31,
     m.abun.nBinomOffS.31vog,
     file = "./statistics/brms/231127_abundance_OffSet.RData")
