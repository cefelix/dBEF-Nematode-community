####data and packages####
library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(hexbin)
library(GGally)

# a seed:
SEED = 22061996

####Hill_q0 11: Hill_q0 ~ sowndivLog*treatment + (1|block/plot), fam=negBinom####
dBEF_nem21 %>% filter(is.na(Hill_q0)==FALSE) %>% pull(Hill_q0) %>% density %>% plot()
#as it is a count, lets try negbinomial:
SEED = 22061996
m.Hill_q0.11 <- brm(Hill_q0 ~ sowndivLog*treatment + (1|block/plot), 
                            data = dBEF_nem21, family = "negbinomial",
                            chains = 3,
                            cores = 3,
                            iter = 2000, warmup = 1000,
                            seed = SEED,
                            control = list(adapt_delta = 0.9) ) 

m.Hill_q0.12 <- update(m.Hill_q0.11,
                       control = list(adapt_delta = 0.99))

pp_check(m.Hill_q0.12, ndraws=100) # thats off

####H1 Shannon ####
dBEF_nem21 %>% filter(is.na(Shannon_H)==FALSE) %>% pull(Shannon_H) %>% density %>% plot()
SEED = 22061996
m.Shannon_H.11 <- brm(Shannon_H ~ sowndivLog*treatment + (1|block/plot), 
                    data = dBEF_nem21, family = "gaussian",
                    chains = 3,
                    cores = 3,
                    iter = 2000, warmup = 1000,
                    seed = SEED,
                    control = list(adapt_delta = 0.9) ) 
  #10 divergent transitions

m.Shannon_H.12 <- update(m.Shannon_H.11,
                         control=list(adapt_delta=0.99))
  #7 divergent transitions
m.Shannon_H.13 <- update(m.Shannon_H.12,
                         control=list(adapt_delta=0.999))

pp_check(m.Shannon_H.13, ndraws=100) # a little overpredictive in quartile below mean

#### hill q1 ####
dBEF_nem21 %>% filter(is.na(Hill_q1)==FALSE) %>% pull(Hill_q1) %>% density %>% plot()
SEED = 22061996
m.Hill_q1.11 <- brm(Hill_q1 ~ sowndivLog*treatment + (1|block/plot), 
                      data = dBEF_nem21, family = "gaussian",
                      chains = 3,
                      cores = 3,
                      iter = 2000, warmup = 1000,
                      seed = SEED,
                      control = list(adapt_delta = 0.9) ) 

m.Hill_q1.12 <- update(m.Hill_q1.11,
                       control = list(adapt_delta = 0.99) )

pp_check(m.Hill_q1.12, ndraws=100)
mcmc_plot(m.Hill_q1.12 )

#### Hill q1 without 60 species ####
dBEF_nem21.no60 <- subset(dBEF_nem21, sowndiv != 60)
  SEED = 22061996
  m.no60.Hill_q1.11 <- brm(Hill_q1 ~ sowndivLog*treatment + (1|block/plot), 
                      data = dBEF_nem21.no60, family = "gaussian",
                      chains = 3,
                      cores = 3,
                      iter = 2000, warmup = 1000,
                      seed = SEED,
                      control = list(adapt_delta = 0.9) ) 
  
  m.no60.Hill_q1.12 <- update(m.no60.Hill_q1.11,
                              control=list(adapt_delta=0.99))
  
  pp_check(m.no60.Hill_q1.12, ndraws=100)



####save ####
save(m.Hill_q0.11, m.Hill_q0.12,
     m.Shannon_H.11, m.Shannon_H.12, m.Shannon_H.13,
     m.Hill_q1.11, m.Hill_q1.12,
     m.no60.Hill_q1.11, m.no60.Hill_q1.12 ,
     file = "./statistics/brms/231108_hillNumbers.RData" )

load(file = "./statistics/brms/231108_hillNumbers.RData" )
conditional_effects(m.Hill_q1.12, prob=0.89)
