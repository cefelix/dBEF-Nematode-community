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
  
  load("./statistics/brms/240214_hill_all_sowndiv.RData")
  
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
    summary(m.all.Shannon.gaus_p4, prob=0.9 )
      #stepwise simplification: choose p4, as week is marginally significant (90% CI)
    
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
        file="./statistics/brms/240214_hill_all_sowndiv.RData")
    
    rm(m.all.Shannon.gamma_p5, m.all.Shannon.gamma_p4, m.all.Shannon.gamma_p32,m.all.Shannon.gamma_p31,
      m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p,
      m.all.Shannon.gaus_p5,m.all.Shannon.gaus_p4, m.all.Shannon.gaus_p32, m.all.Shannon.gaus_p31,
      m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p)

#### Shannon Ba ~ sowndivLogStd: gamma_p4 (looic: gamma_p4)  ####
    #all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
    # stepwise elimiantion: gaus_p4, as week is significant (95% CI)
    #looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
    #differ significantly (less than 2 SE) from the best fit model
    load("./statistics/brms/240214_hill_Ba_sowndiv.RData")
    
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
                                     seed = SEED) #10 div
    summary(m.Ba.Shannon.gaus_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.Ba.Shannon.gaus_p32 <- update(m.Ba.Shannon.gaus_p, 
                                     bf(Hill_q1.Ba ~ sowndivLogStd*treatment + treatment*week + 
                                          (1|block/plot)),
                                     seed = SEED) # 9 div
    summary(m.Ba.Shannon.gaus_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.Ba.Shannon.gaus_p4 <- update(m.Ba.Shannon.gaus_p, 
                                    bf(Hill_q1.Ba ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                    seed = SEED)
    summary(m.Ba.Shannon.gaus_p4, prob=0.95 ) #77 divergent transitions
    #stepwise simplification: choose p4, as week is significant (95% CI)
    
    #remove week:
    m.Ba.Shannon.gaus_p5 <- update(m.Ba.Shannon.gaus_p, 
                                    bf(Hill_q1.Ba ~ sowndivLogStd*treatment + (1|block/plot)),
                                    seed = SEED)
    summary(m.Ba.Shannon.gaus_p5, prob=0.9 ) #0 div trans
    
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
                                 control = list(adapt_delta = 0.99) )  #14 div
    summary(m.Ba.Shannon.gamma_p, prob=0.9)
    
    #remove 3 way interaction:
    m.Ba.Shannon.gamma_p2 <- update(m.Ba.Shannon.gamma_p, 
                                     bf(Hill_q1.Ba ~ sowndivLogStd*treatment + treatment*week + 
                                          sowndivLogStd*week + (1|block/plot)),
                                     seed = SEED) # 22 div
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
                                     seed = SEED) #14 div
    summary(m.Ba.Shannon.gamma_p4, prob=0.95 )
    #stepwise simplification: use gamma_p4 keep week, as it is significant (CI=95%)
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
    #use gamma p4, as it is most parsimonous model with an elpd_diff < 2 SE 
    
    save(m.Ba.Shannon.gamma_p5, m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
         m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
         m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
         m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p,
         file="./statistics/brms/240214_hill_Ba_sowndiv.RData")
    
    rm(m.Ba.Shannon.gamma_p5, m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
       m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
       m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
       m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p)
    
#### Shannon Fu ~ sowndivLogStd: gamma_p4 (looic: gamma_p4) ####
    #all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
    # stepwise elimiantion: gaus_p4, as week is significant (95% CI)
    #looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
    #differ significantly (less than 2 SE) from the best fit model
    load(file="./statistics/brms/240214_hill_Fu_sowndiv.RData")
      
    #using a gaussian distribution with a normal(0,10) prior for beta coefficients
    m.Fu.Shannon.gaus_p <- brm(Hill_q1.Fu ~ sowndivLogStd*treatment*week + (1|block/plot), 
                               data = dat, family = "gaussian",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               prior = beta_coeff_priors,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) )  #1 div
    summary(m.Fu.Shannon.gaus_p, prob=0.9 )
    
    #remove 3 way interaction:
    m.Fu.Shannon.gaus_p2 <- update(m.Fu.Shannon.gaus_p, 
                                   bf(Hill_q1.Fu ~ sowndivLogStd*treatment + treatment*week + 
                                        sowndivLogStd*week + (1|block/plot)),
                                   seed = SEED) # 12 div
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
    summary(m.Fu.Shannon.gaus_p4, prob=0.9 ) 
    #keep week as it is marginally significant (CI 90%)
  
    #remove week:
    m.Fu.Shannon.gaus_p5 <- update(m.Fu.Shannon.gaus_p, 
                                   bf(Hill_q1.Fu ~ sowndivLogStd*treatment + (1|block/plot)),
                                   seed = SEED)
    summary(m.Fu.Shannon.gaus_p5, prob=0.9 )
    #stepwise simplification
  
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
    summary(m.Fu.Shannon.gamma_p4, prob=0.95 )
    #stepwise simplification keep week, as it is significant (CI=95%)
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
    #using a gamma distribution shows better fit
    
    save(m.Fu.Shannon.gamma_p5, m.Fu.Shannon.gamma_p4, m.Fu.Shannon.gamma_p32,m.Fu.Shannon.gamma_p31,
         m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p,
         m.Fu.Shannon.gaus_p5,m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
         m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p,
         file="./statistics/brms/240214_hill_Fu_sowndiv.RData")
    
    rm(m.Fu.Shannon.gamma_p5, m.Fu.Shannon.gamma_p4, m.Fu.Shannon.gamma_p32,m.Fu.Shannon.gamma_p31,
       m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p,
       m.Fu.Shannon.gaus_p5,m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
       m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p)

####Shannon Pl ~ sowndivLogStd: gaus_p5 (looic: gaus_p5)  ####
    #all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
    # stepwise elimiantion: gaus_p4, as week is significant (95% CI)
    #looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
    #differ significantly (less than 2 SE) from the best fit model
    load(file="./statistics/brms/240214_hill_Pl_sowndiv.RData")
    
    #using a gaussian distribution with a normal(0,10) prior for beta coefficients
    m.Pl.Shannon.gaus_p <- brm(Hill_q1.Pl ~ sowndivLogStd*treatment*week + (1|block/plot), 
                               data = dat, family = "gaussian",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               prior = beta_coeff_priors,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) )  #3 div
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
   
    
    #remove week:
    m.Pl.Shannon.gaus_p5 <- update(m.Pl.Shannon.gaus_p, 
                                   bf(Hill_q1.Pl ~ sowndivLogStd*treatment + (1|block/plot)),
                                   seed = SEED)
    summary(m.Pl.Shannon.gaus_p5, prob=0.9 )
    #stepwise simplification: choose p5, as none of additional terms is significant
    pp_check(m.Pl.Shannon.gaus_p5, ndraws = 100)
    
    
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
         file="./statistics/brms/240214_hill_Pl_sowndiv.RData")
    
    rm(m.Pl.Shannon.gamma_p5, m.Pl.Shannon.gamma_p4, m.Pl.Shannon.gamma_p32,m.Pl.Shannon.gamma_p31,
       m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p,
       m.Pl.Shannon.gaus_p5,m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
       m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p)
    
#### Shannon PrOm ~ sowndivLogStd: gamma_p (looic: gamma_p5)  ####
    #all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
    # stepwise elimiantion: gaus_p4, as week is significant (95% CI)
    #looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
    #differ significantly (less than 2 SE) from the best fit model
    load(file="./statistics/brms/240214_hill_PrOm_sowndiv.RData")
    
    
    #using a gaussian distribution with a normal(0,10) prior for beta coefficients
    m.PrOm.Shannon.gaus_p <- brm(Hill_q1.PrOm ~ sowndivLogStd*treatment*week + (1|block/plot), 
                               data = dat, family = "gaussian",
                               chains = 3,
                               cores = 3,
                               iter = 2000, warmup = 1000,
                               prior = beta_coeff_priors,
                               seed = SEED,
                               control = list(adapt_delta = 0.99) )  #4 div
    summary(m.PrOm.Shannon.gaus_p, prob=0.95 )
    #keep 3-fold interaction as it is significant (CI = 0.95)
    
    #remove 3 way interaction:
    m.PrOm.Shannon.gaus_p2 <- update(m.PrOm.Shannon.gaus_p, 
                                   bf(Hill_q1.PrOm ~ sowndivLogStd*treatment + treatment*week + 
                                        sowndivLogStd*week + (1|block/plot)),
                                   seed = SEED) # 8 div
    summary(m.PrOm.Shannon.gaus_p2, prob=0.9)
    
    #remove treatment*week
    m.PrOm.Shannon.gaus_p31 <- update(m.PrOm.Shannon.gaus_p, 
                                    bf(Hill_q1.PrOm ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                         (1|block/plot)),
                                    seed = SEED) #3 div
    summary(m.PrOm.Shannon.gaus_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.PrOm.Shannon.gaus_p32 <- update(m.PrOm.Shannon.gaus_p, 
                                    bf(Hill_q1.PrOm ~ sowndivLogStd*treatment + treatment*week + 
                                         (1|block/plot)),
                                    seed = SEED) # 4 div
    summary(m.PrOm.Shannon.gaus_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.PrOm.Shannon.gaus_p4 <- update(m.PrOm.Shannon.gaus_p, 
                                   bf(Hill_q1.PrOm ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                   seed = SEED)
    summary(m.PrOm.Shannon.gaus_p4, prob=0.9 )
    #stepwise simplification (when not considering 3-fold interactions): choose p4, as week is marginally significant (90% CI)
    
    #remove week:
    m.PrOm.Shannon.gaus_p5 <- update(m.PrOm.Shannon.gaus_p, 
                                   bf(Hill_q1.PrOm ~ sowndivLogStd*treatment + (1|block/plot)),
                                   seed = SEED)
    summary(m.PrOm.Shannon.gaus_p5, prob=0.9 )
    pp_check(m.PrOm.Shannon.gaus_p, ndraws=100)
    
    loo(m.PrOm.Shannon.gaus_p5,m.PrOm.Shannon.gaus_p4, m.PrOm.Shannon.gaus_p32, m.PrOm.Shannon.gaus_p31,
        m.PrOm.Shannon.gaus_p2, m.PrOm.Shannon.gaus_p)
    
    #using a gamma distribution with a normal(0,10) prior for beta coefficients
    m.PrOm.Shannon.gamma_p <- brm(Hill_q1.PrOm ~ sowndivLogStd*treatment*week + (1|block/plot), 
                                data = dat, family = "gamma",
                                chains = 3,
                                cores = 3,
                                iter = 2000, warmup = 1000,
                                prior = beta_coeff_priors,
                                seed = SEED,
                                control = list(adapt_delta = 0.99) )  # 1 div
    summary(m.PrOm.Shannon.gamma_p, prob=0.9)
    pp_check(m.PrOm.Shannon.gamma_p, ndraws = 100) #visually much better fit than gaus_p
    
    #keep 3-way interaction
    
    #remove 3 way interaction:
    m.PrOm.Shannon.gamma_p2 <- update(m.PrOm.Shannon.gamma_p, 
                                    bf(Hill_q1.PrOm ~ sowndivLogStd*treatment + treatment*week + 
                                         sowndivLogStd*week + (1|block/plot)),
                                    seed = SEED) # 5 div
    summary(m.PrOm.Shannon.gamma_p2, prob=0.9)
    
    #remove treatment*week
    m.PrOm.Shannon.gamma_p31 <- update(m.PrOm.Shannon.gamma_p, 
                                     bf(Hill_q1.PrOm ~ sowndivLogStd*treatment + sowndivLogStd*week + 
                                          (1|block/plot)),
                                     seed = SEED) #12 div
    summary(m.PrOm.Shannon.gamma_p31, prob=0.9 )
    
    #remove sowndiv*week
    m.PrOm.Shannon.gamma_p32 <- update(m.PrOm.Shannon.gamma_p, 
                                     bf(Hill_q1.PrOm ~ sowndivLogStd*treatment + treatment*week + 
                                          (1|block/plot)),
                                     seed = SEED) #4 div
    summary(m.PrOm.Shannon.gamma_p32, prob=0.9 )
    
    #remove both treatment*week and sowndiv*week:
    m.PrOm.Shannon.gamma_p4 <- update(m.PrOm.Shannon.gamma_p, 
                                    bf(Hill_q1.PrOm ~ sowndivLogStd*treatment + week + (1|block/plot)),
                                    seed = SEED)
    summary(m.PrOm.Shannon.gamma_p4, prob=0.9 )
    #stepwise simplification without 3-way interaction: keep week, as it is marginally significant (CI=90%)
    
    pp_check(m.PrOm.Shannon.gamma_p, ndraws=100)
    
    
    #remove week:
    m.PrOm.Shannon.gamma_p5 <- update(m.PrOm.Shannon.gamma_p, 
                                    bf(Hill_q1.PrOm ~ sowndivLogStd*treatment + (1|block/plot)),
                                    seed = SEED)
    summary(m.PrOm.Shannon.gamma_p5, prob=0.9 )
    
    #compare
    PrOm.loo <- loo(m.PrOm.Shannon.gamma_p5, m.PrOm.Shannon.gamma_p4, m.PrOm.Shannon.gamma_p32,m.PrOm.Shannon.gamma_p31,
                  m.PrOm.Shannon.gamma_p2, m.PrOm.Shannon.gamma_p,
                  m.PrOm.Shannon.gaus_p5,m.PrOm.Shannon.gaus_p4, m.PrOm.Shannon.gaus_p32, m.PrOm.Shannon.gaus_p31,
                  m.PrOm.Shannon.gaus_p2, m.PrOm.Shannon.gaus_p)
    PrOm.loo
    #using a gaussian distribution shows better fit
    
    save(m.PrOm.Shannon.gamma_p5, m.PrOm.Shannon.gamma_p4, m.PrOm.Shannon.gamma_p32,m.PrOm.Shannon.gamma_p31,
         m.PrOm.Shannon.gamma_p2, m.PrOm.Shannon.gamma_p,
         m.PrOm.Shannon.gaus_p5,m.PrOm.Shannon.gaus_p4, m.PrOm.Shannon.gaus_p32, m.PrOm.Shannon.gaus_p31,
         m.PrOm.Shannon.gaus_p2, m.PrOm.Shannon.gaus_p,
         file="./statistics/brms/240214_hill_PrOm_sowndiv.RData")
    
    rm(m.PrOm.Shannon.gamma_p5, m.PrOm.Shannon.gamma_p4, m.PrOm.Shannon.gamma_p32,m.PrOm.Shannon.gamma_p31,
       m.PrOm.Shannon.gamma_p2, m.PrOm.Shannon.gamma_p,
       m.PrOm.Shannon.gaus_p5,m.PrOm.Shannon.gaus_p4, m.PrOm.Shannon.gaus_p32, m.PrOm.Shannon.gaus_p31,
       m.PrOm.Shannon.gaus_p2, m.PrOm.Shannon.gaus_p)
    

#### model selection ####
  library(brms)
  library(ggplot2)
    
  #selection criteria: picking the most parsimonious model which has an elpd_diff > -4 to the model with the best fit
  #saving the selected models in a seperate file:
    
#all trophic guilds: gaus_p5
  load(file="./statistics/brms/240214_hill_all_sowndiv.RData")
    
    all.loo <- loo(m.all.Shannon.gamma_p5, m.all.Shannon.gamma_p4, m.all.Shannon.gamma_p32,m.all.Shannon.gamma_p31,
                   m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p,
                   m.all.Shannon.gaus_p5,m.all.Shannon.gaus_p4, m.all.Shannon.gaus_p32, m.all.Shannon.gaus_p31,
                   m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p)
    all.loo #gaus p5
    pp_check(m.all.Shannon.gaus_p5, ndraws = 100)
    conditional_effects(m.all.Shannon.gaus_p5)
    
    rm(m.all.Shannon.gamma_p5, m.all.Shannon.gamma_p4, m.all.Shannon.gamma_p32,m.all.Shannon.gamma_p31,
       m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p,
       m.all.Shannon.gaus_p5,
       m.all.Shannon.gaus_p4, m.all.Shannon.gaus_p32, m.all.Shannon.gaus_p31,
       m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p)
    
#Ba: gamma_p5
    load(file="./statistics/brms/240214_hill_Ba_sowndiv.RData")
    Ba.loo <- loo(m.Ba.Shannon.gamma_p5, m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
                  m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
                  m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
                  m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p)
    Ba.loo #gamma_p5
    pp_check(m.Ba.Shannon.gamma_p5, ndraws = 100)
    conditional_effects(m.Ba.Shannon.gamma_p5)
    
    rm(m.Ba.Shannon.gamma_p5, 
       m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
       m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
       #m.Ba.Shannon.gaus_p5,
       m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
       m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p)
    
    
    
    
#Fu ~ sowndiv: gaus_p5
    load("./statistics/brms/240214_hill_Fu_sowndiv.RData")
    Fu.loo <- loo(m.Fu.Shannon.gamma_p5, m.Fu.Shannon.gamma_p4, m.Fu.Shannon.gamma_p32,m.Fu.Shannon.gamma_p31,
                  m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p,
                  m.Fu.Shannon.gaus_p5,m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
                  m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p)
    Fu.loo
    pp_check(m.Fu.Shannon.gaus_p5, ndraws=100)
    conditional_effects(m.Fu.Shannon.gaus_p5)
    
    rm(m.Fu.Shannon.gamma_p5, 
      m.Fu.Shannon.gamma_p4, m.Fu.Shannon.gamma_p32,m.Fu.Shannon.gamma_p31,
      m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p,
      m.Fu.Shannon.gaus_p5,
      m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
      m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p)
    
#Pl ~ sowndiv: gamma_p5
    load("./statistics/brms/240214_hill_Pl_sowndiv.RData")
    Pl.loo <- loo(m.Pl.Shannon.gamma_p5, m.Pl.Shannon.gamma_p4, m.Pl.Shannon.gamma_p32,m.Pl.Shannon.gamma_p31,
                  m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p,
                  m.Pl.Shannon.gaus_p5,m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
                  m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p)
    Pl.loo
    pp_check(m.Pl.Shannon.gamma_p5, ndraws = 100)
    conditional_effects(m.Pl.Shannon.gamma_p5)
    
    rm(m.Pl.Shannon.gamma_p5, 
       m.Pl.Shannon.gamma_p4, m.Pl.Shannon.gamma_p32,m.Pl.Shannon.gamma_p31,
       m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p,
       #m.Pl.Shannon.gaus_p5,
       m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
       m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p)
    
#PrOm ~ sowndiv: gamma_p5
    load("./statistics/brms/240214_hill_PrOm_sowndiv.RData")
    PrOm.loo <- loo(m.PrOm.Shannon.gamma_p5, m.PrOm.Shannon.gamma_p4, m.PrOm.Shannon.gamma_p32,m.PrOm.Shannon.gamma_p31,
                  m.PrOm.Shannon.gamma_p2, m.PrOm.Shannon.gamma_p,
                  m.PrOm.Shannon.gaus_p5,m.PrOm.Shannon.gaus_p4, m.PrOm.Shannon.gaus_p32, m.PrOm.Shannon.gaus_p31,
                  m.PrOm.Shannon.gaus_p2, m.PrOm.Shannon.gaus_p)
    PrOm.loo
    pp_check(m.PrOm.Shannon.gamma_p5, ndraws = 100)
    conditional_effects(m.PrOm.Shannon.gamma_p5)
    
    rm(#m.PrOm.Shannon.gamma_p5, 
      m.PrOm.Shannon.gamma_p4, m.PrOm.Shannon.gamma_p32,m.PrOm.Shannon.gamma_p31,
      m.PrOm.Shannon.gamma_p2, m.PrOm.Shannon.gamma_p,
      m.PrOm.Shannon.gaus_p5,
      m.PrOm.Shannon.gaus_p4, m.PrOm.Shannon.gaus_p32, m.PrOm.Shannon.gaus_p31,
      m.PrOm.Shannon.gaus_p2, m.PrOm.Shannon.gaus_p)
    
#
#save selected models in .RData file:
    save(m.all.Shannon.gaus_p5, 
         m.Ba.Shannon.gamma_p5,
         m.Fu.Shannon.gaus_p5,
         m.Pl.Shannon.gamma_p5,
         m.PrOm.Shannon.gamma_p5,
         file = "./statistics/brms/240214_Hill_sowndiv_mselect.RData")    
    
    
    
    
    

  
 