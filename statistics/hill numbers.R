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
SEED = 19111996

#load files: 
#load(file = "./statistics/brms/231127_Hill.RData" )

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
  beta_coeff_priors <- prior(normal(0,20), class = "b")  
  beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
  beta_coeff_priors3 <- prior(normal(0,2), class = "b")  

  
#### Shannon all ~ sowndivLogStd ####
#with default priors:  
  m.all.Shannon.gaus_d <- brm(Hill_q1 ~ sowndivLogStd*treatment + (1|block/plot), 
                        data = dat, family = "gaussian",
                        chains = 3,
                        cores = 3,
                        iter = 2000, warmup = 1000,
                        seed = SEED,
                        control = list(adapt_delta = 0.99) ) #0 div
    pp_check(m.all.Shannon.gaus_d, ndraws=100) #okay
    summary(m.all.Shannon.gaus_d, prob = 0.9)
    
  m.all.Shannon.gamma_d <- brm(Hill_q1 ~ sowndivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.all.Shannon.gamma_d, ndraws=100) #okay
  summary(m.all.Shannon.gamma_d, prob = 0.9)
    

  #with a normal(0,20) prior for beta coefficients
    m.all.Shannon.gaus_p <- brm(Hill_q1 ~ sowndivLogStd*treatment + (1|block/plot), 
                                 data = dat, family = "gaussian",
                                 chains = 3,
                                 cores = 3,
                                 iter = 2000, warmup = 1000,
                                 prior = beta_coeff_priors,
                                 seed = SEED,
                                 control = list(adapt_delta = 0.99) )  
    
    m.all.Shannon.gamma_p <- brm(Hill_q1 ~ sowndivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gamma",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  
    
  #a narrower normal(0,5) prior:
    m.all.Shannon.gaus_p2 <- update( m.all.Shannon.gaus_p,
                                      prior =   beta_coeff_priors2,
                                      seed = SEED)
    
    m.all.Shannon.gamma_p2 <- update( m.all.Shannon.gamma_p,
                                      prior =   beta_coeff_priors2,
                                      seed = SEED)
    
    pp_check(m.all.Shannon.gamma_p2, ndraws=100)
    summary(m.all.Shannon.gamma_p2, prob = 0.9)
    
  #a narrower normal(0,2) prior:
    m.all.Shannon.gaus_p3 <- update( m.all.Shannon.gaus_p2,
                                      prior =   beta_coeff_priors3,
                                      seed = SEED)
    pp_check(m.all.Shannon.gaus_p3, ndraws=100)
    summary(m.all.Shannon.gaus_p3, prob = 0.9)
    
    m.all.Shannon.gamma_p3 <- update( m.all.Shannon.gamma_p2,
                               prior =   beta_coeff_priors3,
                               seed = SEED)
    pp_check(m.all.Shannon.gamma_p3, ndraws=100)
    summary(m.all.Shannon.gamma_p3, prob = 0.9)
  
#compare them: 
  loo(m.all.Shannon.gaus_p3, m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p, m.all.Shannon.gaus_d,
      m.all.Shannon.gamma_p3, m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p, m.all.Shannon.gamma_d)
  #best: all within 2 SE of elpd range, so all 
  
  save(m.all.Shannon.gaus_p3, m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p, m.all.Shannon.gaus_d,
       m.all.Shannon.gamma_p3, m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p, m.all.Shannon.gamma_d,
       file = "./statistics/brms/231214_all_HillQ1sowndiv.RData")
  
  rm(m.all.Shannon.gaus_p3, m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p, m.all.Shannon.gaus_d,
     m.all.Shannon.gamma_p3, m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p, m.all.Shannon.gamma_d)


#### Shannon Ba ~ sowndivLogStd ####
  
  #with default priors:  
  m.Ba.Shannon.gaus_d <- brm(Hill_q1.Ba ~ sowndivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gaussian",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) ) #0 div
  pp_check(m.Ba.Shannon.gaus_d, ndraws=100) #okay
  summary(m.Ba.Shannon.gaus_d, prob = 0.9)
  
  m.Ba.Shannon.gamma_d <- brm(Hill_q1.Ba ~ sowndivLogStd*treatment + (1|block/plot), 
                               data = dat, family = "gamma",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.Ba.Shannon.gamma_d, ndraws=100) #okay
  summary(m.Ba.Shannon.gamma_d, prob = 0.9)
  
  
  #with a normal(0,20) prior for beta coefficients
  m.Ba.Shannon.gaus_p <- brm(Hill_q1.Ba ~ sowndivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gaussian",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              prior = beta_coeff_priors,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) )  
  
  m.Ba.Shannon.gamma_p <- brm(Hill_q1.Ba ~ sowndivLogStd*treatment + (1|block/plot), 
                               data = dat, family = "gamma",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               prior = beta_coeff_priors,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) )  
  
  #a narrower normal(0,5) prior:
  m.Ba.Shannon.gaus_p2 <- update( m.Ba.Shannon.gaus_p,
                                   prior =   beta_coeff_priors2,
                                   seed = SEED)
  
  m.Ba.Shannon.gamma_p2 <- update( m.Ba.Shannon.gamma_p,
                                    prior =   beta_coeff_priors2,
                                    seed = SEED)
  
  pp_check(m.Ba.Shannon.gamma_p2, ndraws=100)
  summary(m.Ba.Shannon.gamma_p2, prob = 0.9)
  
  #a narrower normal(0,2) prior:
  m.Ba.Shannon.gaus_p3 <- update( m.Ba.Shannon.gaus_p2,
                                   prior =   beta_coeff_priors3,
                                   seed = SEED)
  pp_check(m.Ba.Shannon.gaus_p3, ndraws=100)
  summary(m.Ba.Shannon.gaus_p3, prob = 0.9)
  
  m.Ba.Shannon.gamma_p3 <- update( m.Ba.Shannon.gamma_p2,
                                    prior =   beta_coeff_priors3,
                                    seed = SEED)
  pp_check(m.Ba.Shannon.gamma_p3, ndraws=100)
  summary(m.Ba.Shannon.gamma_p3, prob = 0.9)
  
  #compare them: 
  loo(m.Ba.Shannon.gaus_p3, m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p, m.Ba.Shannon.gaus_d,
      m.Ba.Shannon.gamma_p3, m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p, m.Ba.Shannon.gamma_d)
  #best:
  summary(m.Ba.Shannon.gaus_p3, prob=0.9)
  summary(m.Ba.Shannon.gamma_p, prob=0.9)
  
  save(m.Ba.Shannon.gaus_p3, m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p, m.Ba.Shannon.gaus_d,
       m.Ba.Shannon.gamma_p3, m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p, m.Ba.Shannon.gamma_d,
       file = "./statistics/brms/231214_Ba_HillQ1sowndiv_.RData")
  
  rm(m.Ba.Shannon.gaus_p3, m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p, m.Ba.Shannon.gaus_d,
     m.Ba.Shannon.gamma_p3, m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p, m.Ba.Shannon.gamma_d)
  
  #### Shannon Fu ~ sowndivLogStd ####
  
  #with default priors:  
  m.Fu.Shannon.gaus_d <- brm(Hill_q1.Fu ~ sowndivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) ) #0 div
  pp_check(m.Fu.Shannon.gaus_d, ndraws=100) #okay
  summary(m.Fu.Shannon.gaus_d, prob = 0.9)
  
  m.Fu.Shannon.gamma_d <- brm(Hill_q1.Fu ~ sowndivLogStd*treatment + (1|block/plot), 
                               data = dat, family = "gamma",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.Fu.Shannon.gamma_d, ndraws=100) #okay
  summary(m.Fu.Shannon.gamma_d, prob = 0.9)
  
  
  #with a normal(0,20) prior for beta coefficients
  m.Fu.Shannon.gaus_p <- brm(Hill_q1.Fu ~ sowndivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  
  
  m.Fu.Shannon.gamma_p <- brm(Hill_q1.Fu ~ sowndivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              prior = beta_coeff_priors,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) )  
  
  #a narrower normal(0,5) prior:
  m.Fu.Shannon.gaus_p2 <- update( m.Fu.Shannon.gaus_p,
                                  prior =   beta_coeff_priors2,
                                  seed = SEED)
  
  m.Fu.Shannon.gamma_p2 <- update( m.Fu.Shannon.gamma_p,
                                   prior =   beta_coeff_priors2,
                                   seed = SEED)
  
  pp_check(m.Fu.Shannon.gamma_p2, ndraws=100)
  summary(m.Fu.Shannon.gamma_p2, prob = 0.9)
  
  #a narrower normal(0,2) prior:
  m.Fu.Shannon.gaus_p3 <- update( m.Fu.Shannon.gaus_p2,
                                  prior =   beta_coeff_priors3,
                                  seed = SEED)
  pp_check(m.Fu.Shannon.gaus_p3, ndraws=100)
  summary(m.Fu.Shannon.gaus_p3, prob = 0.9)
  
  m.Fu.Shannon.gamma_p3 <- update( m.Fu.Shannon.gamma_p2,
                                   prior =   beta_coeff_priors3,
                                   seed = SEED)
  pp_check(m.Fu.Shannon.gamma_p3, ndraws=100)
  summary(m.Fu.Shannon.gamma_p3, prob = 0.9)
  
  #compare them: 
  loo(m.Fu.Shannon.gaus_p3, m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p, m.Fu.Shannon.gaus_d,
      m.Fu.Shannon.gamma_p3, m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p, m.Fu.Shannon.gamma_d)
  #best:
  summary(m.Fu.Shannon.gaus_p3, prob=0.9)
  summary(m.Fu.Shannon.gamma_p, prob=0.9)
  
  save(m.Fu.Shannon.gaus_p3, m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p, m.Fu.Shannon.gaus_d,
       m.Fu.Shannon.gamma_p3, m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p, m.Fu.Shannon.gamma_d,
       file = "./statistics/brms/231214_Fu_HillQ1sowndiv_.RData")
  
  rm(m.Fu.Shannon.gaus_p3, m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p, m.Fu.Shannon.gaus_d,
     m.Fu.Shannon.gamma_p3, m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p, m.Fu.Shannon.gamma_d)
  
  
  #### Shannon Pl ~ sowndivLogStd ####
  
  #with default priors:  
  m.Pl.Shannon.gaus_d <- brm(Hill_q1.Pl ~ sowndivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) ) #0 div
  pp_check(m.Pl.Shannon.gaus_d, ndraws=100) #okay
  summary(m.Pl.Shannon.gaus_d, prob = 0.9)
  
  m.Pl.Shannon.gamma_d <- brm(Hill_q1.Pl ~ sowndivLogStd*treatment + (1|block/plot), 
                               data = dat, family = "gamma",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.Pl.Shannon.gamma_d, ndraws=100) #okay
  summary(m.Pl.Shannon.gamma_d, prob = 0.9)
  
  
  #with a normal(0,20) prior for beta coefficients
  m.Pl.Shannon.gaus_p <- brm(Hill_q1.Pl ~ sowndivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  
  
  m.Pl.Shannon.gamma_p <- brm(Hill_q1.Pl ~ sowndivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              prior = beta_coeff_priors,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) )  
  
  #a narrower normal(0,5) prior:
  m.Pl.Shannon.gaus_p2 <- update( m.Pl.Shannon.gaus_p,
                                  prior =   beta_coeff_priors2,
                                  seed = SEED)
  
  m.Pl.Shannon.gamma_p2 <- update( m.Pl.Shannon.gamma_p,
                                   prior =   beta_coeff_priors2,
                                   seed = SEED)
  
  pp_check(m.Pl.Shannon.gamma_p2, ndraws=100)
  summary(m.Pl.Shannon.gamma_p2, prob = 0.9)
  
  #a narrower normal(0,2) prior:
  m.Pl.Shannon.gaus_p3 <- update( m.Pl.Shannon.gaus_p2,
                                  prior =   beta_coeff_priors3,
                                  seed = SEED)
  pp_check(m.Pl.Shannon.gaus_p3, ndraws=100)
  summary(m.Pl.Shannon.gaus_p3, prob = 0.9)
  
  m.Pl.Shannon.gamma_p3 <- update( m.Pl.Shannon.gamma_p2,
                                   prior =   beta_coeff_priors3,
                                   seed = SEED)
  pp_check(m.Pl.Shannon.gamma_p3, ndraws=100)
  summary(m.Pl.Shannon.gamma_p3, prob = 0.9)
  
  #compare them: 
  loo(m.Pl.Shannon.gaus_p3, m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p, m.Pl.Shannon.gaus_d,
      m.Pl.Shannon.gamma_p3, m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p, m.Pl.Shannon.gamma_d)
  #best:
  summary(m.Pl.Shannon.gaus_p3, prob=0.9)
  summary(m.Pl.Shannon.gamma_p, prob=0.9)
  
  save(m.Pl.Shannon.gaus_p3, m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p, m.Pl.Shannon.gaus_d,
       m.Pl.Shannon.gamma_p3, m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p, m.Pl.Shannon.gamma_d,
       file = "./statistics/brms/231214_Pl_HillQ1sowndiv_.RData")
  
  rm(m.Pl.Shannon.gaus_p3, m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p, m.Pl.Shannon.gaus_d,
     m.Pl.Shannon.gamma_p3, m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p, m.Pl.Shannon.gamma_d)
  
#### Shannon Pr ~ sowndivLogStd ####
  #with default priors:  
  m.Pr.Shannon.gaus_d <- brm(Hill_q1.Pr ~ sowndivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) ) #0 div
  pp_check(m.Pr.Shannon.gaus_d, ndraws=100) #okay
  summary(m.Pr.Shannon.gaus_d, prob = 0.9)
  
  m.Pr.Shannon.gamma_d <- brm(Hill_q1.Pr ~ sowndivLogStd*treatment + (1|block/plot), 
                               data = dat, family = "gamma",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.Pr.Shannon.gamma_d, ndraws=100) #okay
  summary(m.Pr.Shannon.gamma_d, prob = 0.9)
  
  
  #with a normal(0,20) prior for beta coefficients
  m.Pr.Shannon.gaus_p <- brm(Hill_q1.Pr ~ sowndivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  
  
  m.Pr.Shannon.gamma_p <- brm(Hill_q1.Pr ~ sowndivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              prior = beta_coeff_priors,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) )  
  
  #a narrower normal(0,5) prior:
  m.Pr.Shannon.gaus_p2 <- update( m.Pr.Shannon.gaus_p,
                                  prior =   beta_coeff_priors2,
                                  seed = SEED)
  
  m.Pr.Shannon.gamma_p2 <- update( m.Pr.Shannon.gamma_p,
                                   prior =   beta_coeff_priors2,
                                   seed = SEED)
  
  pp_check(m.Pr.Shannon.gamma_p2, ndraws=100)
  summary(m.Pr.Shannon.gamma_p2, prob = 0.9)
  
  #a narrower normal(0,2) prior:
  m.Pr.Shannon.gaus_p3 <- update( m.Pr.Shannon.gaus_p2,
                                  prior =   beta_coeff_priors3,
                                  seed = SEED)
  pp_check(m.Pr.Shannon.gaus_p3, ndraws=100)
  summary(m.Pr.Shannon.gaus_p3, prob = 0.9)
  
  m.Pr.Shannon.gamma_p3 <- update( m.Pr.Shannon.gamma_p2,
                                   prior =   beta_coeff_priors3,
                                   seed = SEED)
  pp_check(m.Pr.Shannon.gamma_p3, ndraws=100)
  summary(m.Pr.Shannon.gamma_p3, prob = 0.9)
  
  #compare them: 
  loo(m.Pr.Shannon.gaus_p3, m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p, m.Pr.Shannon.gaus_d,
      m.Pr.Shannon.gamma_p3, m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p, m.Pr.Shannon.gamma_d)
  #best:
  summary(m.Pr.Shannon.gaus_p3, prob=0.9)
  summary(m.Pr.Shannon.gamma_p, prob=0.9)
  
  save(m.Pr.Shannon.gaus_p3, m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p, m.Pr.Shannon.gaus_d,
       m.Pr.Shannon.gamma_p3, m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p, m.Pr.Shannon.gamma_d,
       file = "./statistics/brms/231214_Pr_HillQ1sowndiv_.RData")
  
  rm(m.Pr.Shannon.gaus_p3, m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p, m.Pr.Shannon.gaus_d,
     m.Pr.Shannon.gamma_p3, m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p, m.Pr.Shannon.gamma_d)
  
  #### Shannon Om ~ sowndivLogStd ####
  
  #with default priors:  
  m.Om.Shannon.gaus_d <- brm(Hill_q1.Om ~ sowndivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) ) #0 div
  pp_check(m.Om.Shannon.gaus_d, ndraws=100) #okay
  summary(m.Om.Shannon.gaus_d, prob = 0.9)
  
  m.Om.Shannon.gamma_d <- brm(Hill_q1.Om ~ sowndivLogStd*treatment + (1|block/plot), 
                               data = dat, family = "gamma",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.Om.Shannon.gamma_d, ndraws=100) #okay
  summary(m.Om.Shannon.gamma_d, prob = 0.9)
  
  
  #with a normal(0,20) prior for beta coefficients
  m.Om.Shannon.gaus_p <- brm(Hill_q1.Om ~ sowndivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  
  
  m.Om.Shannon.gamma_p <- brm(Hill_q1.Om ~ sowndivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              prior = beta_coeff_priors,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) )  
  
  #a narrower normal(0,5) prior:
  m.Om.Shannon.gaus_p2 <- update( m.Om.Shannon.gaus_p,
                                  prior =   beta_coeff_priors2,
                                  seed = SEED)
  
  m.Om.Shannon.gamma_p2 <- update( m.Om.Shannon.gamma_p,
                                   prior =   beta_coeff_priors2,
                                   seed = SEED)
  
  pp_check(m.Om.Shannon.gamma_p2, ndraws=100)
  summary(m.Om.Shannon.gamma_p2, prob = 0.9)
  
  #a narrower normal(0,2) prior:
  m.Om.Shannon.gaus_p3 <- update( m.Om.Shannon.gaus_p2,
                                  prior =   beta_coeff_priors3,
                                  seed = SEED)
  pp_check(m.Om.Shannon.gaus_p3, ndraws=100)
  summary(m.Om.Shannon.gaus_p3, prob = 0.9)
  
  m.Om.Shannon.gamma_p3 <- update( m.Om.Shannon.gamma_p2,
                                   prior =   beta_coeff_priors3,
                                   seed = SEED)
  pp_check(m.Om.Shannon.gamma_p3, ndraws=100)
  summary(m.Om.Shannon.gamma_p3, prob = 0.9)
  
  #compare them: 
  loo(m.Om.Shannon.gaus_p3, m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p, m.Om.Shannon.gaus_d,
      m.Om.Shannon.gamma_p3, m.Om.Shannon.gamma_p2, m.Om.Shannon.gamma_p, m.Om.Shannon.gamma_d)
  #best:
  summary(m.Om.Shannon.gaus_p3, prob=0.9)
  summary(m.Om.Shannon.gamma_p, prob=0.9)
  
  save(m.Om.Shannon.gaus_p3, m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p, m.Om.Shannon.gaus_d,
       m.Om.Shannon.gamma_p3, m.Om.Shannon.gamma_p2, m.Om.Shannon.gamma_p, m.Om.Shannon.gamma_d,
       file = "./statistics/brms/231214_Om_HillQ1sowndiv_.RData")
  
  rm(m.Om.Shannon.gaus_p3, m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p, m.Om.Shannon.gaus_d,
     m.Om.Shannon.gamma_p3, m.Om.Shannon.gamma_p2, m.Om.Shannon.gamma_p, m.Om.Shannon.gamma_d)
  
  #### Shannon all ~ realdivLogStd ####
  #with default priors:  
  m.all.Shannon.gaus_d <- brm(Hill_q1 ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gaussian",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) ) #0 div
  pp_check(m.all.Shannon.gaus_d, ndraws=100) #okay
  summary(m.all.Shannon.gaus_d, prob = 0.9)
  
  m.all.Shannon.gamma_d <- brm(Hill_q1 ~ realdivLogStd*treatment + (1|block/plot), 
                               data = dat, family = "gamma",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.all.Shannon.gamma_d, ndraws=100) #okay
  summary(m.all.Shannon.gamma_d, prob = 0.9)
  
  
  #with a normal(0,20) prior for beta coefficients
  m.all.Shannon.gaus_p <- brm(Hill_q1 ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gaussian",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              prior = beta_coeff_priors,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) )  
  
  m.all.Shannon.gamma_p <- brm(Hill_q1 ~ realdivLogStd*treatment + (1|block/plot), 
                               data = dat, family = "gamma",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               prior = beta_coeff_priors,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) )  
  
  #a narrower normal(0,5) prior:
  m.all.Shannon.gaus_p2 <- update( m.all.Shannon.gaus_p,
                                   prior =   beta_coeff_priors2,
                                   seed = SEED)
  
  m.all.Shannon.gamma_p2 <- update( m.all.Shannon.gamma_p,
                                    prior =   beta_coeff_priors2,
                                    seed = SEED)
  
  pp_check(m.all.Shannon.gamma_p2, ndraws=100)
  summary(m.all.Shannon.gamma_p2, prob = 0.9)
  
  #a narrower normal(0,2) prior:
  m.all.Shannon.gaus_p3 <- update( m.all.Shannon.gaus_p2,
                                   prior =   beta_coeff_priors3,
                                   seed = SEED)
  pp_check(m.all.Shannon.gaus_p3, ndraws=100)
  summary(m.all.Shannon.gaus_p3, prob = 0.9)
  
  m.all.Shannon.gamma_p3 <- update( m.all.Shannon.gamma_p2,
                                    prior =   beta_coeff_priors3,
                                    seed = SEED)
  pp_check(m.all.Shannon.gamma_p3, ndraws=100)
  summary(m.all.Shannon.gamma_p3, prob = 0.9)
  
  #compare them: 
  loo(m.all.Shannon.gaus_p3, m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p, m.all.Shannon.gaus_d,
      m.all.Shannon.gamma_p3, m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p, m.all.Shannon.gamma_d)
  #best: all within 2 SE of elpd range, so all 
  
  save(m.all.Shannon.gaus_p3, m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p, m.all.Shannon.gaus_d,
       m.all.Shannon.gamma_p3, m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p, m.all.Shannon.gamma_d,
       file = "./statistics/brms/231214_all_HillQ1realdiv.RData")
  
  
  
  #### Shannon Ba ~ realdivLogStd ####
  
  #with default priors:  
  m.Ba.Shannon.gaus_d <- brm(Hill_q1.Ba ~ realdivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) ) #0 div
  pp_check(m.Ba.Shannon.gaus_d, ndraws=100) #okay
  summary(m.Ba.Shannon.gaus_d, prob = 0.9)
  
  m.Ba.Shannon.gamma_d <- brm(Hill_q1.Ba ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.Ba.Shannon.gamma_d, ndraws=100) #okay
  summary(m.Ba.Shannon.gamma_d, prob = 0.9)
  
  
  #with a normal(0,20) prior for beta coefficients
  m.Ba.Shannon.gaus_p <- brm(Hill_q1.Ba ~ realdivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  
  
  m.Ba.Shannon.gamma_p <- brm(Hill_q1.Ba ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              prior = beta_coeff_priors,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) )  
  
  #a narrower normal(0,5) prior:
  m.Ba.Shannon.gaus_p2 <- update( m.Ba.Shannon.gaus_p,
                                  prior =   beta_coeff_priors2,
                                  seed = SEED)
  
  m.Ba.Shannon.gamma_p2 <- update( m.Ba.Shannon.gamma_p,
                                   prior =   beta_coeff_priors2,
                                   seed = SEED)
  
  pp_check(m.Ba.Shannon.gamma_p2, ndraws=100)
  summary(m.Ba.Shannon.gamma_p2, prob = 0.9)
  
  #a narrower normal(0,2) prior:
  m.Ba.Shannon.gaus_p3 <- update( m.Ba.Shannon.gaus_p2,
                                  prior =   beta_coeff_priors3,
                                  seed = SEED)
  pp_check(m.Ba.Shannon.gaus_p3, ndraws=100)
  summary(m.Ba.Shannon.gaus_p3, prob = 0.9)
  
  m.Ba.Shannon.gamma_p3 <- update( m.Ba.Shannon.gamma_p2,
                                   prior =   beta_coeff_priors3,
                                   seed = SEED)
  pp_check(m.Ba.Shannon.gamma_p3, ndraws=100)
  summary(m.Ba.Shannon.gamma_p3, prob = 0.9)
  
  #compare them: 
  loo(m.Ba.Shannon.gaus_p3, m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p, m.Ba.Shannon.gaus_d,
      m.Ba.Shannon.gamma_p3, m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p, m.Ba.Shannon.gamma_d)
  #best:
  summary(m.Ba.Shannon.gaus_p3, prob=0.9)
  summary(m.Ba.Shannon.gamma_p, prob=0.9)
  
  save(m.Ba.Shannon.gaus_p3, m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p, m.Ba.Shannon.gaus_d,
       m.Ba.Shannon.gamma_p3, m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p, m.Ba.Shannon.gamma_d,
       file = "./statistics/brms/231214_Ba_HillQ1realdiv_.RData")
  
  rm(m.Ba.Shannon.gaus_p3, m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p, m.Ba.Shannon.gaus_d,
     m.Ba.Shannon.gamma_p3, m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p, m.Ba.Shannon.gamma_d)
  
  #### Shannon Fu ~ realdivLogStd ####
  
  #with default priors:  
  m.Fu.Shannon.gaus_d <- brm(Hill_q1.Fu ~ realdivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) ) #0 div
  pp_check(m.Fu.Shannon.gaus_d, ndraws=100) #okay
  summary(m.Fu.Shannon.gaus_d, prob = 0.9)
  
  m.Fu.Shannon.gamma_d <- brm(Hill_q1.Fu ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.Fu.Shannon.gamma_d, ndraws=100) #okay
  summary(m.Fu.Shannon.gamma_d, prob = 0.9)
  
  
  #with a normal(0,20) prior for beta coefficients
  m.Fu.Shannon.gaus_p <- brm(Hill_q1.Fu ~ realdivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  
  
  m.Fu.Shannon.gamma_p <- brm(Hill_q1.Fu ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              prior = beta_coeff_priors,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) )  
  
  #a narrower normal(0,5) prior:
  m.Fu.Shannon.gaus_p2 <- update( m.Fu.Shannon.gaus_p,
                                  prior =   beta_coeff_priors2,
                                  seed = SEED)
  
  m.Fu.Shannon.gamma_p2 <- update( m.Fu.Shannon.gamma_p,
                                   prior =   beta_coeff_priors2,
                                   seed = SEED)
  
  pp_check(m.Fu.Shannon.gamma_p2, ndraws=100)
  summary(m.Fu.Shannon.gamma_p2, prob = 0.9)
  
  #a narrower normal(0,2) prior:
  m.Fu.Shannon.gaus_p3 <- update( m.Fu.Shannon.gaus_p2,
                                  prior =   beta_coeff_priors3,
                                  seed = SEED)
  pp_check(m.Fu.Shannon.gaus_p3, ndraws=100)
  summary(m.Fu.Shannon.gaus_p3, prob = 0.9)
  
  m.Fu.Shannon.gamma_p3 <- update( m.Fu.Shannon.gamma_p2,
                                   prior =   beta_coeff_priors3,
                                   seed = SEED)
  pp_check(m.Fu.Shannon.gamma_p3, ndraws=100)
  summary(m.Fu.Shannon.gamma_p3, prob = 0.9)
  
  #compare them: 
  loo(m.Fu.Shannon.gaus_p3, m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p, m.Fu.Shannon.gaus_d,
      m.Fu.Shannon.gamma_p3, m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p, m.Fu.Shannon.gamma_d)
  #best:
  summary(m.Fu.Shannon.gaus_p3, prob=0.9)
  summary(m.Fu.Shannon.gamma_p, prob=0.9)
  
  save(m.Fu.Shannon.gaus_p3, m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p, m.Fu.Shannon.gaus_d,
       m.Fu.Shannon.gamma_p3, m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p, m.Fu.Shannon.gamma_d,
       file = "./statistics/brms/231214_Fu_HillQ1realdiv_.RData")
  
  rm(m.Fu.Shannon.gaus_p3, m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p, m.Fu.Shannon.gaus_d,
     m.Fu.Shannon.gamma_p3, m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p, m.Fu.Shannon.gamma_d)
  
  
  #### Shannon Pl ~ realdivLogStd ####
  
  #with default priors:  
  m.Pl.Shannon.gaus_d <- brm(Hill_q1.Pl ~ realdivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) ) #0 div
  pp_check(m.Pl.Shannon.gaus_d, ndraws=100) #okay
  summary(m.Pl.Shannon.gaus_d, prob = 0.9)
  
  m.Pl.Shannon.gamma_d <- brm(Hill_q1.Pl ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.Pl.Shannon.gamma_d, ndraws=100) #okay
  summary(m.Pl.Shannon.gamma_d, prob = 0.9)
  
  
  #with a normal(0,20) prior for beta coefficients
  m.Pl.Shannon.gaus_p <- brm(Hill_q1.Pl ~ realdivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  
  
  m.Pl.Shannon.gamma_p <- brm(Hill_q1.Pl ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              prior = beta_coeff_priors,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) )  
  
  #a narrower normal(0,5) prior:
  m.Pl.Shannon.gaus_p2 <- update( m.Pl.Shannon.gaus_p,
                                  prior =   beta_coeff_priors2,
                                  seed = SEED)
  
  m.Pl.Shannon.gamma_p2 <- update( m.Pl.Shannon.gamma_p,
                                   prior =   beta_coeff_priors2,
                                   seed = SEED)
  
  pp_check(m.Pl.Shannon.gamma_p2, ndraws=100)
  summary(m.Pl.Shannon.gamma_p2, prob = 0.9)
  
  #a narrower normal(0,2) prior:
  m.Pl.Shannon.gaus_p3 <- update( m.Pl.Shannon.gaus_p2,
                                  prior =   beta_coeff_priors3,
                                  seed = SEED)
  pp_check(m.Pl.Shannon.gaus_p3, ndraws=100)
  summary(m.Pl.Shannon.gaus_p3, prob = 0.9)
  
  m.Pl.Shannon.gamma_p3 <- update( m.Pl.Shannon.gamma_p2,
                                   prior =   beta_coeff_priors3,
                                   seed = SEED)
  pp_check(m.Pl.Shannon.gamma_p3, ndraws=100)
  summary(m.Pl.Shannon.gamma_p3, prob = 0.9)
  
  #compare them: 
  loo(m.Pl.Shannon.gaus_p3, m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p, m.Pl.Shannon.gaus_d,
      m.Pl.Shannon.gamma_p3, m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p, m.Pl.Shannon.gamma_d)
  #best:
  summary(m.Pl.Shannon.gaus_p3, prob=0.9)
  summary(m.Pl.Shannon.gamma_p, prob=0.9)
  
  save(m.Pl.Shannon.gaus_p3, m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p, m.Pl.Shannon.gaus_d,
       m.Pl.Shannon.gamma_p3, m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p, m.Pl.Shannon.gamma_d,
       file = "./statistics/brms/231214_Pl_HillQ1_.RData")
  
  rm(m.Pl.Shannon.gaus_p3, m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p, m.Pl.Shannon.gaus_d,
     m.Pl.Shannon.gamma_p3, m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p, m.Pl.Shannon.gamma_d)
  
  #### Shannon Pr ~ realdivLogStd ####
  #with default priors:  
  m.Pr.Shannon.gaus_d <- brm(Hill_q1.Pr ~ realdivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) ) #0 div
  pp_check(m.Pr.Shannon.gaus_d, ndraws=100) #okay
  summary(m.Pr.Shannon.gaus_d, prob = 0.9)
  
  m.Pr.Shannon.gamma_d <- brm(Hill_q1.Pr ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.Pr.Shannon.gamma_d, ndraws=100) #okay
  summary(m.Pr.Shannon.gamma_d, prob = 0.9)
  
  
  #with a normal(0,20) prior for beta coefficients
  m.Pr.Shannon.gaus_p <- brm(Hill_q1.Pr ~ realdivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  
  
  m.Pr.Shannon.gamma_p <- brm(Hill_q1.Pr ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              prior = beta_coeff_priors,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) )  
  
  #a narrower normal(0,5) prior:
  m.Pr.Shannon.gaus_p2 <- update( m.Pr.Shannon.gaus_p,
                                  prior =   beta_coeff_priors2,
                                  seed = SEED)
  
  m.Pr.Shannon.gamma_p2 <- update( m.Pr.Shannon.gamma_p,
                                   prior =   beta_coeff_priors2,
                                   seed = SEED)
  
  pp_check(m.Pr.Shannon.gamma_p2, ndraws=100)
  summary(m.Pr.Shannon.gamma_p2, prob = 0.9)
  
  #a narrower normal(0,2) prior:
  m.Pr.Shannon.gaus_p3 <- update( m.Pr.Shannon.gaus_p2,
                                  prior =   beta_coeff_priors3,
                                  seed = SEED)
  pp_check(m.Pr.Shannon.gaus_p3, ndraws=100)
  summary(m.Pr.Shannon.gaus_p3, prob = 0.9)
  
  m.Pr.Shannon.gamma_p3 <- update( m.Pr.Shannon.gamma_p2,
                                   prior =   beta_coeff_priors3,
                                   seed = SEED)
  pp_check(m.Pr.Shannon.gamma_p3, ndraws=100)
  summary(m.Pr.Shannon.gamma_p3, prob = 0.9)
  
  #compare them: 
  loo(m.Pr.Shannon.gaus_p3, m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p, m.Pr.Shannon.gaus_d,
      m.Pr.Shannon.gamma_p3, m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p, m.Pr.Shannon.gamma_d)
  #best:
  summary(m.Pr.Shannon.gaus_p3, prob=0.9)
  summary(m.Pr.Shannon.gamma_p, prob=0.9)
  
  save(m.Pr.Shannon.gaus_p3, m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p, m.Pr.Shannon.gaus_d,
       m.Pr.Shannon.gamma_p3, m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p, m.Pr.Shannon.gamma_d,
       file = "./statistics/brms/231214_Pr_HillQ1realdiv_.RData")
  
  rm(m.Pr.Shannon.gaus_p3, m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p, m.Pr.Shannon.gaus_d,
     m.Pr.Shannon.gamma_p3, m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p, m.Pr.Shannon.gamma_d)
  
  #### Shannon Om ~ realdivLogStd ####
  
  #with default priors:  
  m.Om.Shannon.gaus_d <- brm(Hill_q1.Om ~ realdivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) ) #0 div
  pp_check(m.Om.Shannon.gaus_d, ndraws=100) #okay
  summary(m.Om.Shannon.gaus_d, prob = 0.9)
  
  m.Om.Shannon.gamma_d <- brm(Hill_q1.Om ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) ) #0 div  
  pp_check(m.Om.Shannon.gamma_d, ndraws=100) #okay
  summary(m.Om.Shannon.gamma_d, prob = 0.9)
  
  
  #with a normal(0,20) prior for beta coefficients
  m.Om.Shannon.gaus_p <- brm(Hill_q1.Om ~ realdivLogStd*treatment + (1|block/plot), 
                             data = dat, family = "gaussian",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  
  
  m.Om.Shannon.gamma_p <- brm(Hill_q1.Om ~ realdivLogStd*treatment + (1|block/plot), 
                              data = dat, family = "gamma",
                              chains = 3,
                              cores = 3,
                              iter = 2000, warmup = 1000,
                              prior = beta_coeff_priors,
                              seed = SEED,
                              control = list(adapt_delta = 0.99) )  
  
  #a narrower normal(0,5) prior:
  m.Om.Shannon.gaus_p2 <- update( m.Om.Shannon.gaus_p,
                                  prior =   beta_coeff_priors2,
                                  seed = SEED)
  
  m.Om.Shannon.gamma_p2 <- update( m.Om.Shannon.gamma_p,
                                   prior =   beta_coeff_priors2,
                                   seed = SEED)
  
  pp_check(m.Om.Shannon.gamma_p2, ndraws=100)
  summary(m.Om.Shannon.gamma_p2, prob = 0.9)
  
  #a narrower normal(0,2) prior:
  m.Om.Shannon.gaus_p3 <- update( m.Om.Shannon.gaus_p2,
                                  prior =   beta_coeff_priors3,
                                  seed = SEED)
  pp_check(m.Om.Shannon.gaus_p3, ndraws=100)
  summary(m.Om.Shannon.gaus_p3, prob = 0.9)
  
  m.Om.Shannon.gamma_p3 <- update( m.Om.Shannon.gamma_p2,
                                   prior =   beta_coeff_priors3,
                                   seed = SEED)
  pp_check(m.Om.Shannon.gamma_p3, ndraws=100)
  summary(m.Om.Shannon.gamma_p3, prob = 0.9)
  
  #compare them: 
  loo(m.Om.Shannon.gaus_p3, m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p, m.Om.Shannon.gaus_d,
      m.Om.Shannon.gamma_p3, m.Om.Shannon.gamma_p2, m.Om.Shannon.gamma_p, m.Om.Shannon.gamma_d)
  #best:
  summary(m.Om.Shannon.gaus_p3, prob=0.9)
  summary(m.Om.Shannon.gamma_p, prob=0.9)
  
  save(m.Om.Shannon.gaus_p3, m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p, m.Om.Shannon.gaus_d,
       m.Om.Shannon.gamma_p3, m.Om.Shannon.gamma_p2, m.Om.Shannon.gamma_p, m.Om.Shannon.gamma_d,
       file = "./statistics/brms/231214_Om_HillQ1realdiv_.RData")
  
  rm(m.Om.Shannon.gaus_p3, m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p, m.Om.Shannon.gaus_d,
     m.Om.Shannon.gamma_p3, m.Om.Shannon.gamma_p2, m.Om.Shannon.gamma_p, m.Om.Shannon.gamma_d)



