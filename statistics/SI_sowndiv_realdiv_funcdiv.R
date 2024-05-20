#assessing the maturity index
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

#an appropriate distribution: 
SI <- dat %>% 
  filter(!is.na(SI)) %>%
  select(SI) 
plot(density(SI$SI))  
rm(SI) #gaus looks best, but also try gamma

ggplot(dat, aes(x=sowndivLog, y=SI, col = treatment))+
  geom_jitter(width = 0.2)


#### SI ~ sowndiv, family = gaussian #####
load("./statistics/brms/240131_SI_sowndiv.RData")

SEED = 22061996

#using a gausian dis  
m.SI.sowndiv_gaus_p <- brm(
  bf(SI ~ sowndivLogStd*treatment*week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = gaussian,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #

pp_check(m.SI.sowndiv_gaus_p, ndraws=100)
summary(m.SI.sowndiv_gaus_p, prob =0.9)

#remove 3 way interaction:
m.SI.sowndiv_gaus_p2 <- update(m.SI.sowndiv_gaus_p,
                               bf(SI ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + (1|block/plot)),
                               seed = SEED) #7 div
pp_check(m.SI.sowndiv_gaus_p2, ndraws=100)
summary(m.SI.sowndiv_gaus_p2, prob =0.9)

#remove treatment*week:
m.SI.sowndiv_gaus_p31 <- update(m.SI.sowndiv_gaus_p,
                                bf(SI ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot)),
                                seed = SEED) #2 div
pp_check(m.SI.sowndiv_gaus_p31, ndraws=100)
summary(m.SI.sowndiv_gaus_p31, prob =0.9)

#remove sowndivLogStd*week:
m.SI.sowndiv_gaus_p32 <- update(m.SI.sowndiv_gaus_p,
                                bf(SI ~ sowndivLogStd*treatment + treatment*week + (1|block/plot)),
                                seed = SEED) #2 div
pp_check(m.SI.sowndiv_gaus_p32, ndraws=100)
summary(m.SI.sowndiv_gaus_p32, prob =0.9)

#remove sowndivLogStd*week and streatment*week:
m.SI.sowndiv_gaus_p4 <- update(m.SI.sowndiv_gaus_p,
                               bf(SI ~ sowndivLogStd*treatment + week + (1|block/plot)),
                               seed = SEED) #4 div
pp_check(m.SI.sowndiv_gaus_p4, ndraws=100)
summary(m.SI.sowndiv_gaus_p4, prob =0.9) 

#remove week
m.SI.sowndiv_gaus_p5 <- update(m.SI.sowndiv_gaus_p,
                               bf(SI ~ sowndivLogStd*treatment + (1|block/plot)),
                               seed = SEED) #1 div
pp_check(m.SI.sowndiv_gaus_p5, ndraws=100)
summary(m.SI.sowndiv_gaus_p5, prob =0.9)
conditional_effects(m.SI.sowndiv_gaus_p5)

loo_SI_sowndiv <- loo( m.SI.sowndiv_gaus_p,  m.SI.sowndiv_gaus_p2,  m.SI.sowndiv_gaus_p31,  m.SI.sowndiv_gaus_p32,  
                       m.SI.sowndiv_gaus_p4,  m.SI.sowndiv_gaus_p5)
loo_SI_sowndiv #p5 wins

#save the models
save(m.SI.sowndiv_gaus_p, m.SI.sowndiv_gaus_p2, m.SI.sowndiv_gaus_p31, m.SI.sowndiv_gaus_p32,
     m.SI.sowndiv_gaus_p4, m.SI.sowndiv_gaus_p5,
     file = "./statistics/brms/240131_SI_sowndiv.RData")


#### SI ~ realdiv, family = gaussian #####
load("./statistics/brms/240131_SI_realdiv.RData")
    SEED = 22061996
    
    #using a gausian dis  
    m.SI.realdiv_gaus_p <- brm(
      bf(SI ~ realdivLogStd*treatment*week + (1|block/plot)),
      data = dat, 
      prior = beta_coeff_priors,
      family = gaussian,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #
    
    pp_check(m.SI.realdiv_gaus_p, ndraws=100)
    summary(m.SI.realdiv_gaus_p, prob =0.9)
    
    #remove 3 way interaction:
    m.SI.realdiv_gaus_p2 <- update(m.SI.realdiv_gaus_p,
                                   bf(SI ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + (1|block/plot)),
                                   seed = SEED) #7 div
    pp_check(m.SI.realdiv_gaus_p2, ndraws=100)
    summary(m.SI.realdiv_gaus_p2, prob =0.9)
    
    #remove treatment*week:
    m.SI.realdiv_gaus_p31 <- update(m.SI.realdiv_gaus_p,
                                    bf(SI ~ realdivLogStd*treatment + realdivLogStd*week + (1|block/plot)),
                                    seed = SEED) #2 div
    pp_check(m.SI.realdiv_gaus_p31, ndraws=100)
    summary(m.SI.realdiv_gaus_p31, prob =0.9)
    
    #remove realdivLogStd*week:
    m.SI.realdiv_gaus_p32 <- update(m.SI.realdiv_gaus_p,
                                    bf(SI ~ realdivLogStd*treatment + treatment*week + (1|block/plot)),
                                    seed = SEED) #2 div
    pp_check(m.SI.realdiv_gaus_p32, ndraws=100)
    summary(m.SI.realdiv_gaus_p32, prob =0.9)
    
    #remove realdivLogStd*week and streatment*week:
    m.SI.realdiv_gaus_p4 <- update(m.SI.realdiv_gaus_p,
                                   bf(SI ~ realdivLogStd*treatment + week + (1|block/plot)),
                                   seed = SEED) #4 div
    pp_check(m.SI.realdiv_gaus_p4, ndraws=100)
    summary(m.SI.realdiv_gaus_p4, prob =0.9) 
    
    #remove week
    m.SI.realdiv_gaus_p5 <- update(m.SI.realdiv_gaus_p,
                                   bf(SI ~ realdivLogStd*treatment + (1|block/plot)),
                                   seed = SEED) #1 div
    pp_check(m.SI.realdiv_gaus_p5, ndraws=100)
    summary(m.SI.realdiv_gaus_p5, prob =0.9)
    conditional_effects(m.SI.realdiv_gaus_p5)
    
    loo_SI_realdiv <- loo( m.SI.realdiv_gaus_p,  m.SI.realdiv_gaus_p2,  m.SI.realdiv_gaus_p31,  m.SI.realdiv_gaus_p32,  
                           m.SI.realdiv_gaus_p4,  m.SI.realdiv_gaus_p5)
    loo_SI_realdiv #p5 wins



save(m.SI.realdiv_gaus_p,  m.SI.realdiv_gaus_p2,  m.SI.realdiv_gaus_p31,  m.SI.realdiv_gaus_p32,  
     m.SI.realdiv_gaus_p4,  m.SI.realdiv_gaus_p5,
     file = "./statistics/brms/240131_SI_realdiv.RData")

#### save the selected models ####
save(m.SI.sowndiv_gaus_p5, m.SI.realdiv_gaus_p5,
     file = "./statistics/brms/240215_SI_mselect.RData")

