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

dat <- dat %>% mutate(Pr4_Om4 = Pr4 + Om4, .after = Om4) %>%
  mutate(Pr5_Om5 = Pr5 + Om5, .after = Pr4_Om4) %>%
  mutate(Pr_Om_per100g = Pr_per100g + Om_per100g, .after = Om_per100g)


datW1 <- subset(dat, week=="W1")
datW2 <- subset(dat, week=="W2")

#priors    
beta_coeff_priors <- prior(normal(0,10), class = "b") 

####Pr4+Om4 ~ realdiv ####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(dat$Pr4_Om4==0) #46 -> hu ~ realdivLogStd*treatment

m.Pr4_Om4_realdiv_p <- brm(
  bf(Pr4_Om4 ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd*treatment  + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 4000, warmup = 2000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #10 div

pp_check(m.Pr4_Om4_realdiv_p, ndraws=100)+
  xlim(0,300)
summary(m.Pr4_Om4_realdiv_p, prob =0.9)

m.Pr4_Om4_realdiv_p2 <- update(m.Pr4_Om4_realdiv_p,
                               bf(Pr4_Om4 ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                                    (1|block/plot),
                                  hu ~ realdivLogStd*treatment  + (1|block/plot)),
                               seed=SEED)
m.Pr4_Om4_realdiv_p2 %>% pp_check(ndraws=100)

m.Pr4_Om4_realdiv_p31 <- update(m.Pr4_Om4_realdiv_p,
                                bf(Pr4_Om4 ~ realdivLogStd*treatment + treatment*week + 
                                     (1|block/plot),
                                   hu ~ realdivLogStd*treatment  + (1|block/plot)),
                                seed=SEED)
m.Pr4_Om4_realdiv_p31 %>% pp_check(ndraws=100)

m.Pr4_Om4_realdiv_p32 <- update(m.Pr4_Om4_realdiv_p,
                                bf(Pr4_Om4 ~ realdivLogStd*treatment + realdivLogStd*week  + 
                                     (1|block/plot),
                                   hu ~ realdivLogStd*treatment  + (1|block/plot)),
                                seed=SEED)
m.Pr4_Om4_realdiv_p32 %>% pp_check(ndraws=100)

m.Pr4_Om4_realdiv_p4 <- update(m.Pr4_Om4_realdiv_p,
                               bf(Pr4_Om4 ~ realdivLogStd*treatment  + week + (1|block/plot),
                                  hu ~ realdivLogStd*treatment  + (1|block/plot)),
                               seed=SEED)
m.Pr4_Om4_realdiv_p4 %>% pp_check(ndraws=100)

m.Pr4_Om4_realdiv_p5 <- update(m.Pr4_Om4_realdiv_p,
                               bf(Pr4_Om4 ~ realdivLogStd*treatment + (1|block/plot),
                                  hu ~ realdivLogStd*treatment  + (1|block/plot)),
                               seed=SEED)
m.Pr4_Om4_realdiv_p5 %>% pp_check(ndraws=100)

loo.Pr4_Om4 <- loo(m.Pr4_Om4_realdiv_p, m.Pr4_Om4_realdiv_p2, m.Pr4_Om4_realdiv_p31, 
                   m.Pr4_Om4_realdiv_p32, m.Pr4_Om4_realdiv_p4, m.Pr4_Om4_realdiv_p5)
loo.Pr4_Om4

save(m.Pr4_Om4_realdiv_p, m.Pr4_Om4_realdiv_p2, m.Pr4_Om4_realdiv_p31, 
     m.Pr4_Om4_realdiv_p32, m.Pr4_Om4_realdiv_p4, m.Pr4_Om4_realdiv_p5,
     file = "./statistics/brms/240205_Pr4_Om4_realdiv.RData")

rm(m.Pr4_Om4_realdiv_p, m.Pr4_Om4_realdiv_p2, m.Pr4_Om4_realdiv_p31, 
   m.Pr4_Om4_realdiv_p32, m.Pr4_Om4_realdiv_p4, m.Pr4_Om4_realdiv_p5)

#### Pr + Om  ~ realdiv ####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(dat$Pr_Om==0) #43 -> hu ~ realdivLogStd*treatment

m.Pr_Om_realdiv_p <- brm(
  bf(Pr_Om_per100g ~ realdivLogStd*treatment*week + (1|block/plot),
     hu ~ realdivLogStd*treatment  + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #10 div

pp_check(m.Pr_Om_realdiv_p, ndraws=100)+
  xlim(0,300)
summary(m.Pr_Om_realdiv_p, prob =0.9)

m.Pr_Om_realdiv_p2 <- update(m.Pr_Om_realdiv_p,
                               bf(Pr_Om_per100g  ~ realdivLogStd*treatment + realdivLogStd*week + treatment*week + 
                                    (1|block/plot),
                                  hu ~ realdivLogStd*treatment  + (1|block/plot)),
                               seed=SEED)
m.Pr_Om_realdiv_p2 %>% pp_check(ndraws=100)

m.Pr_Om_realdiv_p31 <- update(m.Pr_Om_realdiv_p,
                                bf(Pr_Om_per100g ~ realdivLogStd*treatment + treatment*week + 
                                     (1|block/plot),
                                   hu ~ realdivLogStd*treatment  + (1|block/plot)),
                                seed=SEED)
m.Pr_Om_realdiv_p31 %>% pp_check(ndraws=100)

m.Pr_Om_realdiv_p32 <- update(m.Pr_Om_realdiv_p,
                                bf(Pr_Om_per100g ~ realdivLogStd*treatment + realdivLogStd*week  + 
                                     (1|block/plot),
                                   hu ~ realdivLogStd*treatment  + (1|block/plot)),
                                seed=SEED)
m.Pr_Om_realdiv_p32 %>% pp_check(ndraws=100)

m.Pr_Om_realdiv_p4 <- update(m.Pr_Om_realdiv_p,
                               bf(Pr_Om_per100g ~ realdivLogStd*treatment  + week + (1|block/plot),
                                  hu ~ realdivLogStd*treatment  + (1|block/plot)),
                               seed=SEED)
m.Pr_Om_realdiv_p4 %>% pp_check(ndraws=100)

m.Pr_Om_realdiv_p5 <- update(m.Pr_Om_realdiv_p,
                               bf(Pr_Om_per100g ~ realdivLogStd*treatment + (1|block/plot),
                                  hu ~ realdivLogStd*treatment  + (1|block/plot)),
                               seed=SEED)
m.Pr_Om_realdiv_p5 %>% pp_check(ndraws=100)

m.Pr_Om_realdiv_p6 <- update(m.Pr_Om_realdiv_p,
                             bf(Pr_Om_per100g ~ realdivLogStd*treatment + (1|block/plot),
                                hu ~ realdivLogStd + week  + (1|block/plot)),
                             seed=SEED)
m.Pr_Om_realdiv_p6 %>% pp_check(ndraws=100)
summary(m.Pr_Om_realdiv_p6)

m.Pr_Om_realdiv_p7 <- update(m.Pr_Om_realdiv_p,
                             bf(Pr_Om_per100g ~ realdivLogStd*treatment + (1|block/plot),
                                hu ~  week  + (1|block/plot)),
                             seed=SEED)
m.Pr_Om_realdiv_p7 %>% pp_check(ndraws=100)

m.Pr_Om_realdiv_p8 <- update(m.Pr_Om_realdiv_p,
                             bf(Pr_Om_per100g ~ realdivLogStd*treatment + (1|block/plot),
                                hu ~ 1),
                             seed=SEED)
m.Pr_Om_realdiv_p8 %>% pp_check(ndraws=100)


loo.Pr_Om <- loo(m.Pr_Om_realdiv_p, m.Pr_Om_realdiv_p2, m.Pr_Om_realdiv_p31, 
                   m.Pr_Om_realdiv_p32, m.Pr_Om_realdiv_p4, m.Pr_Om_realdiv_p5,
                   m.Pr_Om_realdiv_p6, m.Pr_Om_realdiv_p7, m.Pr_Om_realdiv_p8)
loo.Pr_Om

save(m.Pr_Om_realdiv_p, m.Pr_Om_realdiv_p2, m.Pr_Om_realdiv_p31, 
     m.Pr_Om_realdiv_p32, m.Pr_Om_realdiv_p4, m.Pr_Om_realdiv_p5,
     m.Pr_Om_realdiv_p6, m.Pr_Om_realdiv_p7, m.Pr_Om_realdiv_p8,
     file = "./statistics/brms/240212_Pr_Om_realdiv.RData")

rm(m.Pr_Om_realdiv_p, m.Pr_Om_realdiv_p2, m.Pr_Om_realdiv_p31, 
   m.Pr_Om_realdiv_p32, m.Pr_Om_realdiv_p4, m.Pr_Om_realdiv_p5,
   m.Pr_Om_realdiv_p6, m.Pr_Om_realdiv_p7, m.Pr_Om_realdiv_p8)

#### model selection ####
load("./statistics/brms/240212_Pr_Om_realdiv.RData")

loo.Pr4_Om4 <- loo(m.Pr4_Om4_realdiv_p, m.Pr4_Om4_realdiv_p2, m.Pr4_Om4_realdiv_p31, 
                   m.Pr4_Om4_realdiv_p32, m.Pr4_Om4_realdiv_p4, m.Pr4_Om4_realdiv_p5)
loo.Pr4_Om4 #p5

#posterior predictive check
pp_check(m.Pr4_Om4_realdiv_p5, ndraws = 100)

#conditional effects
conditional_effects(m.Pr4_Om4_realdiv_p5) #pos. for t1, neutral t3, negative t2

#save
save(m.Pr4_Om4_realdiv_p5,
     file = "./statistics/brms/240205_Pr_Om_cp_realdiv_mselect.RData")
