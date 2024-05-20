#assessing the maturity index
library(brms)


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
MI <- dat %>% 
  filter(!is.na(MI)) %>%
  select(MI) 
plot(density(MI$MI))  
    rm(MI) #gaus looks best, but also try gamma


#### MI ~ sowndiv #####
load("./statistics/brms/240130_MI_sowndiv.RData")    
SEED = 22061996

#using a gausian dis  
m.MI.sowndiv_gaus_p <- brm(
  bf(MI ~ sowndivLogStd*treatment*week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = gaussian,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #1 div, 387 exceeded max_treedepth

  pp_check(m.MI.sowndiv_gaus_p, ndraws=100)+
    xlim(0,5)
  summary(m.MI.sowndiv_gaus_p, prob =0.9)
  
  #remove 3 way interaction:
  m.MI.sowndiv_gaus_p2 <- update(m.MI.sowndiv_gaus_p,
                                 bf(MI ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + (1|block/plot)),
                                 seed = SEED) #7 div
  pp_check(m.MI.sowndiv_gaus_p2, ndraws=100)+
    xlim(0,5)
  summary(m.MI.sowndiv_gaus_p2, prob =0.9)
  
  #remove treatment*week:
  m.MI.sowndiv_gaus_p31 <- update(m.MI.sowndiv_gaus_p,
                                  bf(MI ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot)),
                                  seed = SEED) #2 div
      pp_check(m.MI.sowndiv_gaus_p31, ndraws=100)+
        xlim(0,5)
      summary(m.MI.sowndiv_gaus_p31, prob =0.9)
  
  #remove sowndivLogStd*week:
  m.MI.sowndiv_gaus_p32 <- update(m.MI.sowndiv_gaus_p,
                                  bf(MI ~ sowndivLogStd*treatment + treatment*week + (1|block/plot)),
                                  seed = SEED) #2 div
      pp_check(m.MI.sowndiv_gaus_p32, ndraws=100)+
        xlim(0,5)
      summary(m.MI.sowndiv_gaus_p32, prob =0.9)
  
  #remove sowndivLogStd*week and streatment*week:
  m.MI.sowndiv_gaus_p4 <- update(m.MI.sowndiv_gaus_p,
                                  bf(MI ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                  seed = SEED) #4 div
      pp_check(m.MI.sowndiv_gaus_p4, ndraws=100)+
        xlim(0,5)
      summary(m.MI.sowndiv_gaus_p4, prob =0.9) 
  
  #remove week
  m.MI.sowndiv_gaus_p5 <- update(m.MI.sowndiv_gaus_p,
                                 bf(MI ~ sowndivLogStd*treatment + (1|block/plot)),
                                 seed = SEED) #1 div
      pp_check(m.MI.sowndiv_gaus_p5, ndraws=100)+
        xlim(0,5)
      summary(m.MI.sowndiv_gaus_p5, prob =0.9)
      
  loo_MI_sowndiv <- loo( m.MI.sowndiv_gaus_p,  m.MI.sowndiv_gaus_p2,  m.MI.sowndiv_gaus_p31,  m.MI.sowndiv_gaus_p32,  
                         m.MI.sowndiv_gaus_p4,  m.MI.sowndiv_gaus_p5)    
  loo_MI_sowndiv #p5 wins
  pp_check(m.MI.sowndiv_gaus_p5, ndraws=100)
  conditional_effects(m.MI.sowndiv_gaus_p5)
  
  save(m.MI.sowndiv_gaus_p, m.MI.sowndiv_gaus_p2, m.MI.sowndiv_gaus_p31, m.MI.sowndiv_gaus_p32,
       m.MI.sowndiv_gaus_p4, m.MI.sowndiv_gaus_p5,
       file = "./statistics/brms/240130_MI_sowndiv.RData")
  
  
  
