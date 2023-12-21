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
  beta_coeff_priors <- prior(normal(0,10), class = "b")  
  
#### Shannon all ~ sowndivLogStd: gaus_p4 (looic: gauss_p5) ####
  #all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
  # stepwise elimiantion: gaus_p4, as week is significant (95% CI)
  #looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
    #differ significantly (less than 2 SE) from the best fit model
  
  #using a gaussian distribution with a normal(0,10) prior for beta coefficients
    m.all.Shannon.gaus_p <- brm(Hill_q1 ~ sowndivLogStd*treatment*week + (1|block/plot), 
                                 data = dat, family = "gaussian",
                                 chains = 3,
                                 cores = 3,
                                 iter = 2000, warmup = 1000,
                                 prior = beta_coeff_priors,
                                 seed = SEED,
                                 control = list(adapt_delta = 0.99) )  #4 div
    summary(m.all.Shannon.gaus_p, prob=0.9 )
    
    #remove 3 way interaction:
    m.all.Shannon.gaus_p2 <- update(m.all.Shannon.gaus_p, 
                                    bf(Hill_q1 ~ sowndivLogStd*treatment + treatment*week + 
                                         sowndivLogStd*week + (1|block/plot)),
                                    seed = SEED) # 8 div
    summary(m.all.Shannon.gaus_p2, prob=0.9)
    
    #remove treatment*week
    m.all.Shannon.gaus_p31 <- update(m.all.Shannon.gaus_p, 
                                    bf(Hill_q1 ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                         (1|block/plot)),
                                    seed = SEED) #3 div
    summary(m.all.Shannon.gaus_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.all.Shannon.gaus_p32 <- update(m.all.Shannon.gaus_p, 
                                     bf(Hill_q1 ~ sowndivLogStd*treatment + treatment*week + 
                                          (1|block/plot)),
                                     seed = SEED) # 4 div
    summary(m.all.Shannon.gaus_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.all.Shannon.gaus_p4 <- update(m.all.Shannon.gaus_p, 
                                     bf(Hill_q1 ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                     seed = SEED)
    summary(m.all.Shannon.gaus_p4, prob=0.95 )
      #stepwise simplification: choose p4, as week is significant (95% CI)
    
    #remove week:
    m.all.Shannon.gaus_p5 <- update(m.all.Shannon.gaus_p, 
                                    bf(Hill_q1 ~ sowndivLogStd*treatment + (1|block/plot)),
                                    seed = SEED)
    summary(m.all.Shannon.gaus_p5, prob=0.9 )
    
    loo(m.all.Shannon.gaus_p5,m.all.Shannon.gaus_p4, m.all.Shannon.gaus_p32, m.all.Shannon.gaus_p31,
        m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p)
    
  #using a gamma distribution with a normal(0,10) prior for beta coefficients
    m.all.Shannon.gamma_p <- brm(Hill_q1 ~ sowndivLogStd*treatment*week + (1|block/plot), 
                             data = dat, family = "gamma",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  # 8 div
    summary(m.all.Shannon.gamma_p, prob=0.9)
    
    #remove 3 way interaction:
    m.all.Shannon.gamma_p2 <- update(m.all.Shannon.gamma_p, 
                                    bf(Hill_q1 ~ sowndivLogStd*treatment + treatment*week + 
                                         sowndivLogStd*week + (1|block/plot)),
                                    seed = SEED) # 5 div
    summary(m.all.Shannon.gamma_p2, prob=0.9)
    
    #remove treatment*week
    m.all.Shannon.gamma_p31 <- update(m.all.Shannon.gamma_p, 
                                     bf(Hill_q1 ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                          (1|block/plot)),
                                     seed = SEED) #12 div
    summary(m.all.Shannon.gamma_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.all.Shannon.gamma_p32 <- update(m.all.Shannon.gamma_p, 
                                     bf(Hill_q1 ~ sowndivLogStd*treatment + treatment*week + 
                                          (1|block/plot)),
                                     seed = SEED) #4 div
    summary(m.all.Shannon.gamma_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.all.Shannon.gamma_p4 <- update(m.all.Shannon.gamma_p, 
                                    bf(Hill_q1 ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                    seed = SEED)
    summary(m.all.Shannon.gamma_p4, prob=0.9 )
    #stepwise simplification keep week, as it is marginally significant (CI=90%)
    pp_check(m.all.Shannon.gamma_p4, ndraws=100)
    
    
    #remove week:
    m.all.Shannon.gamma_p5 <- update(m.all.Shannon.gamma_p, 
                                    bf(Hill_q1 ~ sowndivLogStd*treatment + (1|block/plot)),
                                    seed = SEED)
    summary(m.all.Shannon.gamma_p5, prob=0.9 )
    
    #compare
    all.loo <- loo(m.all.Shannon.gamma_p5, m.all.Shannon.gamma_p4, m.all.Shannon.gamma_p32,m.all.Shannon.gamma_p31,
                   m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p,
                   m.all.Shannon.gaus_p5,m.all.Shannon.gaus_p4, m.all.Shannon.gaus_p32, m.all.Shannon.gaus_p31,
                   m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p)
    all.loo
    #using a gaussian distribution shows better fit
    
    save(m.all.Shannon.gamma_p5, m.all.Shannon.gamma_p4, m.all.Shannon.gamma_p32,m.all.Shannon.gamma_p31,
        m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p,
        m.all.Shannon.gaus_p5,m.all.Shannon.gaus_p4, m.all.Shannon.gaus_p32, m.all.Shannon.gaus_p31,
        m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p,
        file="./statistics/brms/231221_hill_all_sowndiv.RData")
    
    rm(m.all.Shannon.gamma_p5, m.all.Shannon.gamma_p4, m.all.Shannon.gamma_p32,m.all.Shannon.gamma_p31,
      m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p,
      m.all.Shannon.gaus_p5,m.all.Shannon.gaus_p4, m.all.Shannon.gaus_p32, m.all.Shannon.gaus_p31,
      m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p)

#### Shannon Ba ~ sowndivLogStd:  ####
    #all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
    # stepwise elimiantion: gaus_p4, as week is significant (95% CI)
    #looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
    #differ significantly (less than 2 SE) from the best fit model
    
    #using a gaussian distribution with a normal(0,10) prior for beta coefficients
    m.Ba.Shannon.gaus_p <- brm(Hill_q1.Ba ~ sowndivLogStd*treatment*week + (1|block/plot), 
                                data = dat, family = "gaussian",
                                chains = 3,
                                cores = 3,
                                iter = 2000, warmup = 1000,
                                prior = beta_coeff_priors,
                                seed = SEED,
                                control = list(adapt_delta = 0.99) )  #4 div
    summary(m.Ba.Shannon.gaus_p, prob=0.9 )
    
    #remove 3 way interaction:
    m.Ba.Shannon.gaus_p2 <- update(m.Ba.Shannon.gaus_p, 
                                    bf(Hill_q1.Ba ~ sowndivLogStd*treatment + treatment*week + 
                                         sowndivLogStd*week + (1|block/plot)),
                                    seed = SEED) # 8 div
    summary(m.Ba.Shannon.gaus_p2, prob=0.9)
    
    #remove treatment*week
    m.Ba.Shannon.gaus_p31 <- update(m.Ba.Shannon.gaus_p, 
                                     bf(Hill_q1.Ba ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                          (1|block/plot)),
                                     seed = SEED) #3 div
    summary(m.Ba.Shannon.gaus_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.Ba.Shannon.gaus_p32 <- update(m.Ba.Shannon.gaus_p, 
                                     bf(Hill_q1.Ba ~ sowndivLogStd*treatment + treatment*week + 
                                          (1|block/plot)),
                                     seed = SEED) # 4 div
    summary(m.Ba.Shannon.gaus_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.Ba.Shannon.gaus_p4 <- update(m.Ba.Shannon.gaus_p, 
                                    bf(Hill_q1.Ba ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                    seed = SEED)
    summary(m.Ba.Shannon.gaus_p4, prob=0.95 )
    #stepwise simplification: choose p4, as week is significant (95% CI)
    
    #remove week:
    m.Ba.Shannon.gaus_p5 <- update(m.Ba.Shannon.gaus_p, 
                                    bf(Hill_q1.Ba ~ sowndivLogStd*treatment + (1|block/plot)),
                                    seed = SEED)
    summary(m.Ba.Shannon.gaus_p5, prob=0.9 )
    
    loo(m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
        m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p)
    
    #using a gamma distribution with a normal(0,10) prior for beta coefficients
    m.Ba.Shannon.gamma_p <- brm(Hill_q1.Ba ~ sowndivLogStd*treatment*week + (1|block/plot), 
                                 data = dat, family = "gamma",
                                 chains = 3,
                                 cores = 3,
                                 iter = 2000, warmup = 1000,
                                 prior = beta_coeff_priors,
                                 seed = SEED,
                                 control = list(adapt_delta = 0.99) )  # 8 div
    summary(m.Ba.Shannon.gamma_p, prob=0.9)
    
    #remove 3 way interaction:
    m.Ba.Shannon.gamma_p2 <- update(m.Ba.Shannon.gamma_p, 
                                     bf(Hill_q1.Ba ~ sowndivLogStd*treatment + treatment*week + 
                                          sowndivLogStd*week + (1|block/plot)),
                                     seed = SEED) # 5 div
    summary(m.Ba.Shannon.gamma_p2, prob=0.9)
    
    #remove treatment*week
    m.Ba.Shannon.gamma_p31 <- update(m.Ba.Shannon.gamma_p, 
                                      bf(Hill_q1.Ba ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                           (1|block/plot)),
                                      seed = SEED) #12 div
    summary(m.Ba.Shannon.gamma_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.Ba.Shannon.gamma_p32 <- update(m.Ba.Shannon.gamma_p, 
                                      bf(Hill_q1.Ba ~ sowndivLogStd*treatment + treatment*week + 
                                           (1|block/plot)),
                                      seed = SEED) #4 div
    summary(m.Ba.Shannon.gamma_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.Ba.Shannon.gamma_p4 <- update(m.Ba.Shannon.gamma_p, 
                                     bf(Hill_q1.Ba ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                     seed = SEED)
    summary(m.Ba.Shannon.gamma_p4, prob=0.9 )
    #stepwise simplification keep week, as it is marginally significant (CI=90%)
    pp_check(m.Ba.Shannon.gamma_p4, ndraws=100)
    
    
    #remove week:
    m.Ba.Shannon.gamma_p5 <- update(m.Ba.Shannon.gamma_p, 
                                     bf(Hill_q1.Ba ~ sowndivLogStd*treatment + (1|block/plot)),
                                     seed = SEED)
    summary(m.Ba.Shannon.gamma_p5, prob=0.9 )
    
    #compare
    Ba.loo <- loo(m.Ba.Shannon.gamma_p5, m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
                   m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
                   m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
                   m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p)
    Ba.loo
    #using a gaussian distribution shows better fit
    
    save(m.Ba.Shannon.gamma_p5, m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
         m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
         m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
         m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p,
         file="./statistics/brms/231221_hill_Ba_sowndiv.RData")
    
    rm(m.Ba.Shannon.gamma_p5, m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
       m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
       m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
       m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p)
    
#### Shannon Fu ~ sowndivLogStd:  ####
    #all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
    # stepwise elimiantion: gaus_p4, as week is significant (95% CI)
    #looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
    #differ significantly (less than 2 SE) from the best fit model
    
    #using a gaussian distribution with a normal(0,10) prior for beta coefficients
    m.Fu.Shannon.gaus_p <- brm(Hill_q1.Fu ~ sowndivLogStd*treatment*week + (1|block/plot), 
                               data = dat, family = "gaussian",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               prior = beta_coeff_priors,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) )  #4 div
    summary(m.Fu.Shannon.gaus_p, prob=0.9 )
    
    #remove 3 way interaction:
    m.Fu.Shannon.gaus_p2 <- update(m.Fu.Shannon.gaus_p, 
                                   bf(Hill_q1.Fu ~ sowndivLogStd*treatment + treatment*week + 
                                        sowndivLogStd*week + (1|block/plot)),
                                   seed = SEED) # 8 div
    summary(m.Fu.Shannon.gaus_p2, prob=0.9)
    
    #remove treatment*week
    m.Fu.Shannon.gaus_p31 <- update(m.Fu.Shannon.gaus_p, 
                                    bf(Hill_q1.Fu ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                         (1|block/plot)),
                                    seed = SEED) #3 div
    summary(m.Fu.Shannon.gaus_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.Fu.Shannon.gaus_p32 <- update(m.Fu.Shannon.gaus_p, 
                                    bf(Hill_q1.Fu ~ sowndivLogStd*treatment + treatment*week + 
                                         (1|block/plot)),
                                    seed = SEED) # 4 div
    summary(m.Fu.Shannon.gaus_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.Fu.Shannon.gaus_p4 <- update(m.Fu.Shannon.gaus_p, 
                                   bf(Hill_q1.Fu ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                   seed = SEED)
    summary(m.Fu.Shannon.gaus_p4, prob=0.95 )
    #stepwise simplification: choose p4, as week is significant (95% CI)
    
    #remove week:
    m.Fu.Shannon.gaus_p5 <- update(m.Fu.Shannon.gaus_p, 
                                   bf(Hill_q1.Fu ~ sowndivLogStd*treatment + (1|block/plot)),
                                   seed = SEED)
    summary(m.Fu.Shannon.gaus_p5, prob=0.9 )
    
    loo(m.Fu.Shannon.gaus_p5,m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
        m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p)
    
    #using a gamma distribution with a normal(0,10) prior for beta coefficients
    m.Fu.Shannon.gamma_p <- brm(Hill_q1.Fu ~ sowndivLogStd*treatment*week + (1|block/plot), 
                                data = dat, family = "gamma",
                                chains = 3,
                                cores = 3,
                                iter = 2000, warmup = 1000,
                                prior = beta_coeff_priors,
                                seed = SEED,
                                control = list(adapt_delta = 0.99) )  # 8 div
    summary(m.Fu.Shannon.gamma_p, prob=0.9)
    
    #remove 3 way interaction:
    m.Fu.Shannon.gamma_p2 <- update(m.Fu.Shannon.gamma_p, 
                                    bf(Hill_q1.Fu ~ sowndivLogStd*treatment + treatment*week + 
                                         sowndivLogStd*week + (1|block/plot)),
                                    seed = SEED) # 5 div
    summary(m.Fu.Shannon.gamma_p2, prob=0.9)
    
    #remove treatment*week
    m.Fu.Shannon.gamma_p31 <- update(m.Fu.Shannon.gamma_p, 
                                     bf(Hill_q1.Fu ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                          (1|block/plot)),
                                     seed = SEED) #12 div
    summary(m.Fu.Shannon.gamma_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.Fu.Shannon.gamma_p32 <- update(m.Fu.Shannon.gamma_p, 
                                     bf(Hill_q1.Fu ~ sowndivLogStd*treatment + treatment*week + 
                                          (1|block/plot)),
                                     seed = SEED) #4 div
    summary(m.Fu.Shannon.gamma_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.Fu.Shannon.gamma_p4 <- update(m.Fu.Shannon.gamma_p, 
                                    bf(Hill_q1.Fu ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                    seed = SEED)
    summary(m.Fu.Shannon.gamma_p4, prob=0.9 )
    #stepwise simplification keep week, as it is marginally significant (CI=90%)
    pp_check(m.Fu.Shannon.gamma_p4, ndraws=100)
    
    
    #remove week:
    m.Fu.Shannon.gamma_p5 <- update(m.Fu.Shannon.gamma_p, 
                                    bf(Hill_q1.Fu ~ sowndivLogStd*treatment + (1|block/plot)),
                                    seed = SEED)
    summary(m.Fu.Shannon.gamma_p5, prob=0.9 )
    
    #compare
    Fu.loo <- loo(m.Fu.Shannon.gamma_p5, m.Fu.Shannon.gamma_p4, m.Fu.Shannon.gamma_p32,m.Fu.Shannon.gamma_p31,
                  m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p,
                  m.Fu.Shannon.gaus_p5,m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
                  m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p)
    Fu.loo
    #using a gaussian distribution shows better fit
    
    save(m.Fu.Shannon.gamma_p5, m.Fu.Shannon.gamma_p4, m.Fu.Shannon.gamma_p32,m.Fu.Shannon.gamma_p31,
         m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p,
         m.Fu.Shannon.gaus_p5,m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
         m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p,
         file="./statistics/brms/231221_hill_Fu_sowndiv.RData")
    
    rm(m.Fu.Shannon.gamma_p5, m.Fu.Shannon.gamma_p4, m.Fu.Shannon.gamma_p32,m.Fu.Shannon.gamma_p31,
       m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p,
       m.Fu.Shannon.gaus_p5,m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
       m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p)

####Shannon Pl ~ sowndivLogStd:  ####
    #all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
    # stepwise elimiantion: gaus_p4, as week is significant (95% CI)
    #looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
    #differ significantly (less than 2 SE) from the best fit model
    
    #using a gaussian distribution with a normal(0,10) prior for beta coefficients
    m.Pl.Shannon.gaus_p <- brm(Hill_q1.Pl ~ sowndivLogStd*treatment*week + (1|block/plot), 
                               data = dat, family = "gaussian",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               prior = beta_coeff_priors,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) )  #4 div
    summary(m.Pl.Shannon.gaus_p, prob=0.9 )
    
    #remove 3 way interaction:
    m.Pl.Shannon.gaus_p2 <- update(m.Pl.Shannon.gaus_p, 
                                   bf(Hill_q1.Pl ~ sowndivLogStd*treatment + treatment*week + 
                                        sowndivLogStd*week + (1|block/plot)),
                                   seed = SEED) # 8 div
    summary(m.Pl.Shannon.gaus_p2, prob=0.9)
    
    #remove treatment*week
    m.Pl.Shannon.gaus_p31 <- update(m.Pl.Shannon.gaus_p, 
                                    bf(Hill_q1.Pl ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                         (1|block/plot)),
                                    seed = SEED) #3 div
    summary(m.Pl.Shannon.gaus_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.Pl.Shannon.gaus_p32 <- update(m.Pl.Shannon.gaus_p, 
                                    bf(Hill_q1.Pl ~ sowndivLogStd*treatment + treatment*week + 
                                         (1|block/plot)),
                                    seed = SEED) # 4 div
    summary(m.Pl.Shannon.gaus_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.Pl.Shannon.gaus_p4 <- update(m.Pl.Shannon.gaus_p, 
                                   bf(Hill_q1.Pl ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                   seed = SEED)
    summary(m.Pl.Shannon.gaus_p4, prob=0.95 )
    #stepwise simplification: choose p4, as week is significant (95% CI)
    
    #remove week:
    m.Pl.Shannon.gaus_p5 <- update(m.Pl.Shannon.gaus_p, 
                                   bf(Hill_q1.Pl ~ sowndivLogStd*treatment + (1|block/plot)),
                                   seed = SEED)
    summary(m.Pl.Shannon.gaus_p5, prob=0.9 )
    
    loo(m.Pl.Shannon.gaus_p5,m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
        m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p)
    
    #using a gamma distribution with a normal(0,10) prior for beta coefficients
    m.Pl.Shannon.gamma_p <- brm(Hill_q1.Pl ~ sowndivLogStd*treatment*week + (1|block/plot), 
                                data = dat, family = "gamma",
                                chains = 3,
                                cores = 3,
                                iter = 2000, warmup = 1000,
                                prior = beta_coeff_priors,
                                seed = SEED,
                                control = list(adapt_delta = 0.99) )  # 8 div
    summary(m.Pl.Shannon.gamma_p, prob=0.9)
    
    #remove 3 way interaction:
    m.Pl.Shannon.gamma_p2 <- update(m.Pl.Shannon.gamma_p, 
                                    bf(Hill_q1.Pl ~ sowndivLogStd*treatment + treatment*week + 
                                         sowndivLogStd*week + (1|block/plot)),
                                    seed = SEED) # 5 div
    summary(m.Pl.Shannon.gamma_p2, prob=0.9)
    
    #remove treatment*week
    m.Pl.Shannon.gamma_p31 <- update(m.Pl.Shannon.gamma_p, 
                                     bf(Hill_q1.Pl ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                          (1|block/plot)),
                                     seed = SEED) #12 div
    summary(m.Pl.Shannon.gamma_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.Pl.Shannon.gamma_p32 <- update(m.Pl.Shannon.gamma_p, 
                                     bf(Hill_q1.Pl ~ sowndivLogStd*treatment + treatment*week + 
                                          (1|block/plot)),
                                     seed = SEED) #4 div
    summary(m.Pl.Shannon.gamma_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.Pl.Shannon.gamma_p4 <- update(m.Pl.Shannon.gamma_p, 
                                    bf(Hill_q1.Pl ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                    seed = SEED)
    summary(m.Pl.Shannon.gamma_p4, prob=0.9 )
    #stepwise simplification keep week, as it is marginally significant (CI=90%)
    pp_check(m.Pl.Shannon.gamma_p4, ndraws=100)
    
    
    #remove week:
    m.Pl.Shannon.gamma_p5 <- update(m.Pl.Shannon.gamma_p, 
                                    bf(Hill_q1.Pl ~ sowndivLogStd*treatment + (1|block/plot)),
                                    seed = SEED)
    summary(m.Pl.Shannon.gamma_p5, prob=0.9 )
    
    #compare
    Pl.loo <- loo(m.Pl.Shannon.gamma_p5, m.Pl.Shannon.gamma_p4, m.Pl.Shannon.gamma_p32,m.Pl.Shannon.gamma_p31,
                  m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p,
                  m.Pl.Shannon.gaus_p5,m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
                  m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p)
    Pl.loo
    #using a gaussian distribution shows better fit
    
    save(m.Pl.Shannon.gamma_p5, m.Pl.Shannon.gamma_p4, m.Pl.Shannon.gamma_p32,m.Pl.Shannon.gamma_p31,
         m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p,
         m.Pl.Shannon.gaus_p5,m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
         m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p,
         file="./statistics/brms/231221_hill_Pl_sowndiv.RData")
    
    rm(m.Pl.Shannon.gamma_p5, m.Pl.Shannon.gamma_p4, m.Pl.Shannon.gamma_p32,m.Pl.Shannon.gamma_p31,
       m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p,
       m.Pl.Shannon.gaus_p5,m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
       m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p)
    
#### Shannon Pr ~ sowndivLogStd:  ####
    #all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
    # stepwise elimiantion: gaus_p4, as week is significant (95% CI)
    #looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
    #differ significantly (less than 2 SE) from the best fit model
    
    #using a gaussian distribution with a normal(0,10) prior for beta coefficients
    m.Pr.Shannon.gaus_p <- brm(Hill_q1.Pr ~ sowndivLogStd*treatment*week + (1|block/plot), 
                               data = dat, family = "gaussian",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               prior = beta_coeff_priors,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) )  #4 div
    summary(m.Pr.Shannon.gaus_p, prob=0.9 )
    
    #remove 3 way interaction:
    m.Pr.Shannon.gaus_p2 <- update(m.Pr.Shannon.gaus_p, 
                                   bf(Hill_q1.Pr ~ sowndivLogStd*treatment + treatment*week + 
                                        sowndivLogStd*week + (1|block/plot)),
                                   seed = SEED) # 8 div
    summary(m.Pr.Shannon.gaus_p2, prob=0.9)
    
    #remove treatment*week
    m.Pr.Shannon.gaus_p31 <- update(m.Pr.Shannon.gaus_p, 
                                    bf(Hill_q1.Pr ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                         (1|block/plot)),
                                    seed = SEED) #3 div
    summary(m.Pr.Shannon.gaus_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.Pr.Shannon.gaus_p32 <- update(m.Pr.Shannon.gaus_p, 
                                    bf(Hill_q1.Pr ~ sowndivLogStd*treatment + treatment*week + 
                                         (1|block/plot)),
                                    seed = SEED) # 4 div
    summary(m.Pr.Shannon.gaus_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.Pr.Shannon.gaus_p4 <- update(m.Pr.Shannon.gaus_p, 
                                   bf(Hill_q1.Pr ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                   seed = SEED)
    summary(m.Pr.Shannon.gaus_p4, prob=0.95 )
    #stepwise simplification: choose p4, as week is significant (95% CI)
    
    #remove week:
    m.Pr.Shannon.gaus_p5 <- update(m.Pr.Shannon.gaus_p, 
                                   bf(Hill_q1.Pr ~ sowndivLogStd*treatment + (1|block/plot)),
                                   seed = SEED)
    summary(m.Pr.Shannon.gaus_p5, prob=0.9 )
    
    loo(m.Pr.Shannon.gaus_p5,m.Pr.Shannon.gaus_p4, m.Pr.Shannon.gaus_p32, m.Pr.Shannon.gaus_p31,
        m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p)
    
    #using a gamma distribution with a normal(0,10) prior for beta coefficients
    m.Pr.Shannon.gamma_p <- brm(Hill_q1.Pr ~ sowndivLogStd*treatment*week + (1|block/plot), 
                                data = dat, family = "gamma",
                                chains = 3,
                                cores = 3,
                                iter = 2000, warmup = 1000,
                                prior = beta_coeff_priors,
                                seed = SEED,
                                control = list(adapt_delta = 0.99) )  # 8 div
    summary(m.Pr.Shannon.gamma_p, prob=0.9)
    
    #remove 3 way interaction:
    m.Pr.Shannon.gamma_p2 <- update(m.Pr.Shannon.gamma_p, 
                                    bf(Hill_q1.Pr ~ sowndivLogStd*treatment + treatment*week + 
                                         sowndivLogStd*week + (1|block/plot)),
                                    seed = SEED) # 5 div
    summary(m.Pr.Shannon.gamma_p2, prob=0.9)
    
    #remove treatment*week
    m.Pr.Shannon.gamma_p31 <- update(m.Pr.Shannon.gamma_p, 
                                     bf(Hill_q1.Pr ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                          (1|block/plot)),
                                     seed = SEED) #12 div
    summary(m.Pr.Shannon.gamma_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.Pr.Shannon.gamma_p32 <- update(m.Pr.Shannon.gamma_p, 
                                     bf(Hill_q1.Pr ~ sowndivLogStd*treatment + treatment*week + 
                                          (1|block/plot)),
                                     seed = SEED) #4 div
    summary(m.Pr.Shannon.gamma_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.Pr.Shannon.gamma_p4 <- update(m.Pr.Shannon.gamma_p, 
                                    bf(Hill_q1.Pr ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                    seed = SEED)
    summary(m.Pr.Shannon.gamma_p4, prob=0.9 )
    #stepwise simplification keep week, as it is marginally significant (CI=90%)
    pp_check(m.Pr.Shannon.gamma_p4, ndraws=100)
    
    
    #remove week:
    m.Pr.Shannon.gamma_p5 <- update(m.Pr.Shannon.gamma_p, 
                                    bf(Hill_q1.Pr ~ sowndivLogStd*treatment + (1|block/plot)),
                                    seed = SEED)
    summary(m.Pr.Shannon.gamma_p5, prob=0.9 )
    
    #compare
    Pr.loo <- loo(m.Pr.Shannon.gamma_p5, m.Pr.Shannon.gamma_p4, m.Pr.Shannon.gamma_p32,m.Pr.Shannon.gamma_p31,
                  m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p,
                  m.Pr.Shannon.gaus_p5,m.Pr.Shannon.gaus_p4, m.Pr.Shannon.gaus_p32, m.Pr.Shannon.gaus_p31,
                  m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p)
    Pr.loo
    #using a gaussian distribution shows better fit
    
    save(m.Pr.Shannon.gamma_p5, m.Pr.Shannon.gamma_p4, m.Pr.Shannon.gamma_p32,m.Pr.Shannon.gamma_p31,
         m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p,
         m.Pr.Shannon.gaus_p5,m.Pr.Shannon.gaus_p4, m.Pr.Shannon.gaus_p32, m.Pr.Shannon.gaus_p31,
         m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p,
         file="./statistics/brms/231221_hill_Pr_sowndiv.RData")
    
    rm(m.Pr.Shannon.gamma_p5, m.Pr.Shannon.gamma_p4, m.Pr.Shannon.gamma_p32,m.Pr.Shannon.gamma_p31,
       m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p,
       m.Pr.Shannon.gaus_p5,m.Pr.Shannon.gaus_p4, m.Pr.Shannon.gaus_p32, m.Pr.Shannon.gaus_p31,
       m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p)
    
#### Shannon Om ~ sowndivLogStd:  ####
    #all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
    # stepwise elimiantion: gaus_p4, as week is significant (95% CI)
    #looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
    #differ significantly (less than 2 SE) from the best fit model
    
    #using a gaussian distribution with a normal(0,10) prior for beta coefficients
    m.Om.Shannon.gaus_p <- brm(Hill_q1.Om ~ sowndivLogStd*treatment*week + (1|block/plot), 
                               data = dat, family = "gaussian",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               prior = beta_coeff_priors,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) )  #4 div
    summary(m.Om.Shannon.gaus_p, prob=0.9 )
    
    #remove 3 way interaction:
    m.Om.Shannon.gaus_p2 <- update(m.Om.Shannon.gaus_p, 
                                   bf(Hill_q1.Om ~ sowndivLogStd*treatment + treatment*week + 
                                        sowndivLogStd*week + (1|block/plot)),
                                   seed = SEED) # 8 div
    summary(m.Om.Shannon.gaus_p2, prob=0.9)
    
    #remove treatment*week
    m.Om.Shannon.gaus_p31 <- update(m.Om.Shannon.gaus_p, 
                                    bf(Hill_q1.Om ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                         (1|block/plot)),
                                    seed = SEED) #3 div
    summary(m.Om.Shannon.gaus_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.Om.Shannon.gaus_p32 <- update(m.Om.Shannon.gaus_p, 
                                    bf(Hill_q1.Om ~ sowndivLogStd*treatment + treatment*week + 
                                         (1|block/plot)),
                                    seed = SEED) # 4 div
    summary(m.Om.Shannon.gaus_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.Om.Shannon.gaus_p4 <- update(m.Om.Shannon.gaus_p, 
                                   bf(Hill_q1.Om ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                   seed = SEED)
    summary(m.Om.Shannon.gaus_p4, prob=0.95 )
    #stepwise simplification: choose p4, as week is significant (95% CI)
    
    #remove week:
    m.Om.Shannon.gaus_p5 <- update(m.Om.Shannon.gaus_p, 
                                   bf(Hill_q1.Om ~ sowndivLogStd*treatment + (1|block/plot)),
                                   seed = SEED)
    summary(m.Om.Shannon.gaus_p5, prob=0.9 )
    
    loo(m.Om.Shannon.gaus_p5,m.Om.Shannon.gaus_p4, m.Om.Shannon.gaus_p32, m.Om.Shannon.gaus_p31,
        m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p)
    
    #using a gamma distribution with a normal(0,10) prior for beta coefficients
    m.Om.Shannon.gamma_p <- brm(Hill_q1.Om ~ sowndivLogStd*treatment*week + (1|block/plot), 
                                data = dat, family = "gamma",
                                chains = 3,
                                cores = 3,
                                iter = 2000, warmup = 1000,
                                prior = beta_coeff_priors,
                                seed = SEED,
                                control = list(adapt_delta = 0.99) )  # 8 div
    summary(m.Om.Shannon.gamma_p, prob=0.9)
    
    #remove 3 way interaction:
    m.Om.Shannon.gamma_p2 <- update(m.Om.Shannon.gamma_p, 
                                    bf(Hill_q1.Om ~ sowndivLogStd*treatment + treatment*week + 
                                         sowndivLogStd*week + (1|block/plot)),
                                    seed = SEED) # 5 div
    summary(m.Om.Shannon.gamma_p2, prob=0.9)
    
    #remove treatment*week
    m.Om.Shannon.gamma_p31 <- update(m.Om.Shannon.gamma_p, 
                                     bf(Hill_q1.Om ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                          (1|block/plot)),
                                     seed = SEED) #12 div
    summary(m.Om.Shannon.gamma_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.Om.Shannon.gamma_p32 <- update(m.Om.Shannon.gamma_p, 
                                     bf(Hill_q1.Om ~ sowndivLogStd*treatment + treatment*week + 
                                          (1|block/plot)),
                                     seed = SEED) #4 div
    summary(m.Om.Shannon.gamma_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.Om.Shannon.gamma_p4 <- update(m.Om.Shannon.gamma_p, 
                                    bf(Hill_q1.Om ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                    seed = SEED)
    summary(m.Om.Shannon.gamma_p4, prob=0.9 )
    #stepwise simplification keep week, as it is marginally significant (CI=90%)
    pp_check(m.Om.Shannon.gamma_p4, ndraws=100)
    
    
    #remove week:
    m.Om.Shannon.gamma_p5 <- update(m.Om.Shannon.gamma_p, 
                                    bf(Hill_q1.Om ~ sowndivLogStd*treatment + (1|block/plot)),
                                    seed = SEED)
    summary(m.Om.Shannon.gamma_p5, prob=0.9 )
    
    #compare
    Om.loo <- loo(m.Om.Shannon.gamma_p5, m.Om.Shannon.gamma_p4, m.Om.Shannon.gamma_p32,m.Om.Shannon.gamma_p31,
                  m.Om.Shannon.gamma_p2, m.Om.Shannon.gamma_p,
                  m.Om.Shannon.gaus_p5,m.Om.Shannon.gaus_p4, m.Om.Shannon.gaus_p32, m.Om.Shannon.gaus_p31,
                  m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p)
    Om.loo
    #using a gaussian distribution shows better fit
    
    save(m.Om.Shannon.gamma_p5, m.Om.Shannon.gamma_p4, m.Om.Shannon.gamma_p32,m.Om.Shannon.gamma_p31,
         m.Om.Shannon.gamma_p2, m.Om.Shannon.gamma_p,
         m.Om.Shannon.gaus_p5,m.Om.Shannon.gaus_p4, m.Om.Shannon.gaus_p32, m.Om.Shannon.gaus_p31,
         m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p,
         file="./statistics/brms/231221_hill_Om_sowndiv.RData")
    
    rm(m.Om.Shannon.gamma_p5, m.Om.Shannon.gamma_p4, m.Om.Shannon.gamma_p32,m.Om.Shannon.gamma_p31,
       m.Om.Shannon.gamma_p2, m.Om.Shannon.gamma_p,
       m.Om.Shannon.gaus_p5,m.Om.Shannon.gaus_p4, m.Om.Shannon.gaus_p32, m.Om.Shannon.gaus_p31,
       m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p)
    
    

  
 