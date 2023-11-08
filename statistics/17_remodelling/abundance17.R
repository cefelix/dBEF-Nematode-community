library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)
 SEED=22061996



#### modelling abundance as a count with soilDW as an offset variable ####

#https://en.wikipedia.org/wiki/Poisson_regression#:~:text=In%20these%20examples%2C%20exposure%20is,right%20side%20of%20the%20equation.
#https://www.dataquest.io/blog/tutorial-poisson-regression-in-r/
#https://github.com/paul-buerkner/brms/issues/98
#offset() is usable as in other packages

#### Poisson 21: total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), fam=poisson ####
#as the link function for poisson is log(), lets log transform offset variable:
dBEF_nem17$soilDW %>% log() %>% density() %>% plot() #bimodal
#lets check plotwise:
ggplot(data = dBEF_nem17, aes(x = sowndivLog, y = soilDW))+
  geom_point()+
  facet_wrap(~block)

#the response variables density
dBEF_nem17$total_nematodes  %>% density() %>% plot()
rpois(1e3, lambda = 2) %>% density() %>% plot()
SEED = 22061996

m.abun.PoisOffS.21 <- brm(total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), 
                          data = dBEF_nem17, family = "poisson",
                          chains = 3,
                          cores = 3,
                          iter = 3000, warmup = 1500,
                          seed = SEED,
                          control = list(adapt_delta = 0.9) ) 
#74 divergent transitions,
#291 exceeded max_treedepth

m.abun.PoisOffS.22 <- update(m.abun.PoisOffS.21,
                             seed = SEED,
                             control = list(adapt_delta = 0.999,
                                            max_treedepth = 12)) 
#427 transitions exceeded max_treedepth
pp_check(m.abun.PoisOffS.22, ndraws=100) #thats at least a better fit than m.[...].14


m.abun.PoisOffS.23 <- update(m.abun.PoisOffS.22,
                             seed = SEED,
                             control = list(adapt_delta = 0.999,
                                            max_treedepth = 15))
#EES too low

m.abun.PoisOffS.24 <- update(m.abun.PoisOffS.23,
                             seed = SEED,
                             iter=4000, warmup =1500,
                             control = list(adapt_delta = 0.999,
                                            max_treedepth = 15))

pp_check(m.abun.PoisOffS.24, ndraws = 100)
#underpredicting until mean

save(m.abun.PoisOffS.11, m.abun.PoisOffS.21,
     m.abun.PoisOffS.12, m.abun.PoisOffS.22,
     m.abun.PoisOffS.13, m.abun.PoisOffS.23,
     m.abun.PoisOffS.14, m.abun.PoisOffS.24,
     file = "./statistics/17_remodelling/brms/231106_abundance_OffSet.RData")


#### negBinom 11: total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), fam=negbinomial ####
m.abun.nBinomOffS.11 <- brm(total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), 
                            data = dBEF_nem17, family = "negbinomial",
                            chains = 3,
                            cores = 3,
                            iter = 3000, warmup = 1500,
                            seed = SEED,
                            control = list(adapt_delta = 0.9) ) 
#9 divergent transitions

m.abun.nBinomOffS.12 <- update(m.abun.nBinomOffS.11,
                               seed = SEED,
                               control = list(adapt_delta = 0.999,
                                              max_treedepth = 12)) 

m.abun.nBinomOffS.13 <- update(m.abun.nBinomOffS.12,
                               seed = SEED,
                               control = list(adapt_delta = 0.999,
                                              max_treedepth = 15)) 



summary(m.abun.nBinomOffS.13)
pp_check(m.abun.nBinomOffS.13, ndraws = 100)

#### negBinom 21: total_nematodes ~ sowndivLog*treatment + offset(log(soilDW)) + (1|block/plot), fam=negbinomial ####
SEED = 22061996
m.abun17.nBinomOffS.21 <- brm(total_nematodes ~ sowndivLog*treatment + offset(log(soilDW)) + (1|block/plot), 
                            data = dBEF_nem17, family = "negbinomial",
                            chains = 3,
                            cores = 3,
                            iter = 3000, warmup = 1500,
                            seed = SEED,
                            control = list(adapt_delta = 0.9) ) 
#8 divergent transitions

m.abun17.nBinomOffS.22 <- update(m.abun17.nBinomOffS.21,
                               control = list(adapt_delta = 0.999,
                                              max_treedepth = 10)) 

pp_check(m.abun17.nBinomOffS.22, ndraws=100)
summary(m.abun17.nBinomOffS.22)

#### save offset models ####
save(#m.abun.PoisOffS.21, m.abun.PoisOffS.22, 
     #m.abun.PoisOffS.23, m.abun.PoisOffS.24,
     #m.abun.nBinomOffS.11, m.abun.nBinomOffS.12, 
     #m.abun.nBinomOffS.13,
     m.abun17.nBinomOffS.21, m.abun17.nBinomOffS.22,
     file="./statistics/17_remodelling/brms/231108_abundance17_OffSet.RData")



####plot the offset models####
#conditional_effects(m.abun.nBinomOffS.13)

conditional_effects(m.abun17.nBinomOffS.22, prob=0.95)

