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
  mutate(Pr5_Om5 = Pr5 + Om5, .after = Pr4_Om4)


datW1 <- subset(dat, week=="W1")
datW2 <- subset(dat, week=="W2")

#priors    
beta_coeff_priors <- prior(normal(0,10), class = "b") 

####Pr4+Om4 ~ sowndiv ####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996

sum(dat$Pr4_Om4==0) #46 -> hu ~ sowndivLogStd*treatment

m.Pr4_Om4_sowndiv_p <- brm(
  bf(Pr4_Om4 ~ sowndivLogStd*treatment*week + (1|block/plot),
     hu ~ sowndivLogStd*treatment  + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #10 div

pp_check(m.Pr4_Om4_sowndiv_p, ndraws=100)+
  xlim(0,300)
summary(m.Pr4_Om4_sowndiv_p, prob =0.9)
conditional_effects(m.Pr4_Om4_sowndiv_p)

m.Pr4_Om4_sowndiv_p2 <- update(m.Pr4_Om4_sowndiv_p,
                           bf(Pr4_Om4 ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                (1|block/plot),
                              hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                           seed=SEED)
m.Pr4_Om4_sowndiv_p2 %>% pp_check(ndraws=100)

m.Pr4_Om4_sowndiv_p31 <- update(m.Pr4_Om4_sowndiv_p,
                            bf(Pr4_Om4 ~ sowndivLogStd*treatment + treatment*week + 
                                 (1|block/plot),
                               hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                            seed=SEED)
m.Pr4_Om4_sowndiv_p31 %>% pp_check(ndraws=100)

m.Pr4_Om4_sowndiv_p32 <- update(m.Pr4_Om4_sowndiv_p,
                            bf(Pr4_Om4 ~ sowndivLogStd*treatment + sowndivLogStd*week  + 
                                 (1|block/plot),
                               hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                            seed=SEED)
m.Pr4_Om4_sowndiv_p32 %>% pp_check(ndraws=100)

m.Pr4_Om4_sowndiv_p4 <- update(m.Pr4_Om4_sowndiv_p,
                           bf(Pr4_Om4 ~ sowndivLogStd*treatment  + week + (1|block/plot),
                              hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                           seed=SEED)
m.Pr4_Om4_sowndiv_p4 %>% pp_check(ndraws=100)

m.Pr4_Om4_sowndiv_p5 <- update(m.Pr4_Om4_sowndiv_p,
                           bf(Pr4_Om4 ~ sowndivLogStd*treatment + (1|block/plot),
                              hu ~ sowndivLogStd*treatment  + (1|block/plot)),
                           seed=SEED)
m.Pr4_Om4_sowndiv_p5 %>% pp_check(ndraws=100)

loo.Pr4_Om4 <- loo(m.Pr4_Om4_sowndiv_p, m.Pr4_Om4_sowndiv_p2, m.Pr4_Om4_sowndiv_p31, 
               m.Pr4_Om4_sowndiv_p32, m.Pr4_Om4_sowndiv_p4, m.Pr4_Om4_sowndiv_p5)
loo.Pr4_Om4

save(m.Pr4_Om4_sowndiv_p, m.Pr4_Om4_sowndiv_p2, m.Pr4_Om4_sowndiv_p31, 
     m.Pr4_Om4_sowndiv_p32, m.Pr4_Om4_sowndiv_p4, m.Pr4_Om4_sowndiv_p5,
     file = "./statistics/brms/240205_Pr4_Om4_sowndiv.RData")

rm(m.Pr4_Om4_sowndiv_p, m.Pr4_Om4_sowndiv_p2, m.Pr4_Om4_sowndiv_p31, 
   m.Pr4_Om4_sowndiv_p32, m.Pr4_Om4_sowndiv_p4, m.Pr4_Om4_sowndiv_p5)

