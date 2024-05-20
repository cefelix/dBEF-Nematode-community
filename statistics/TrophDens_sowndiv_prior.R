library(brms)
library(rstan)
library(ggplot2)
library(emmeans)
library(bayestestR)
library(dplyr)

# fitting trophic group densities ~ sowndiv

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
   
    
    #the narrowest prior is still basically flat (at realistic values): 
    ggplot(data.frame(density(rnorm(1e5, 0, 10)) ), 
           aes(x=x, y=y))+
      geom_point()+
      xlim(-1.5,1.5)
    
    exp(4.5) 
    #this would mean that increasing plant diversity by 1 SD would lead to 348% more individuals per dry weight of soil
    


####Ba ~ sowndiv ####
  #stepwise elimination: p8, as none of the additional terms/interactions was at least marginally significant (90%CI)
  #looic: p7, as it the most parsimonous with a elpd_diff of less than 4 and elpd_diff less than 2 se_diff  
  
  #see here https://users.aalto.fi/~ave/CV-FAQ.html#12_What_is_the_interpretation_of_ELPD__elpd_loo__elpd_diff
  load("./statistics/brms/240108_Ba_sowndiv_priors.RData")
  
  beta_coeff_priors <- prior(normal(0,10), class = "b")  
  SEED = 22061996
  sum(subset(dat, week=="W1")$Ba_per100g == 0) #9
  sum(subset(dat, week=="W2")$Ba_per100g == 0) #2
  
  #for both weeks  
  m.Ba_sowndiv_p <- brm(
    bf(Ba_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
       hu ~ 1 ),
    data = dat, 
    prior = beta_coeff_priors,
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #1 div
  
  pp_check(m.Ba_sowndiv_p, ndraws=100)+
    xlim(0,2000)
  summary(m.Ba_sowndiv_p, prob =0.9)
  
  #as different orientation of sowndivLogStd:treatment2:weekW2 and sowndivLogStd:treatment3:weekW2 -->
  #check pairwise with emmeans:
  emtrends(m.Ba_sowndiv_p, specs = c("treatment", "week"), var="sowndivLogStd") %>%
    summary() #CIs overlap
  
  #remove 3way interaction:
  m.Ba_sowndiv_p2 <- update(m.Ba_sowndiv_p, 
                            bf(Ba_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                 (1|block/plot),
                               hu ~ 1), 
                            seed = SEED) #2 div
  summary(m.Ba_sowndiv_p2, prob =0.9)
  
  emtrends(m.Ba_sowndiv_p2, specs = c("treatment", "week"), var="sowndivLogStd") %>%
    summary() #CIs overlap
  
  #remove sowndivLogStd*week
  m.Ba_sowndiv_p31 <- update(m.Ba_sowndiv_p, 
                             bf(Ba_per100g ~ sowndivLogStd*treatment + treatment*week + (1|block/plot),
                                hu ~ 1 + (1|block/plot)), 
                             seed = SEED) #4 div
  
  summary(m.Ba_sowndiv_p31, prob=0.9)
  emtrends(m.Ba_sowndiv_p31, specs = c("treatment", "week"), var="sowndivLogStd") %>%
    summary()
  
  #remove treatment*week 
  m.Ba_sowndiv_p32 <- update(m.Ba_sowndiv_p, 
                             bf(Ba_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot),
                                hu ~ 1), 
                             seed = SEED)
  summary(m.Ba_sowndiv_p32, prob=0.9)
  
  emtrends(m.Ba_sowndiv_p32, specs = c("treatment", "week"), var="sowndivLogStd") %>%
    summary()
  
  #remove sowndivLogStd:week
  m.Ba_sowndiv_p4 <- update(m.Ba_sowndiv_p, 
                            bf(Ba_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                               hu ~ 1), 
                            seed = SEED)
  summary(m.Ba_sowndiv_p4, prob=0.9)
  pp_check(m.Ba_sowndiv_p4, ndraws=100)+xlim(0,500)
  
  #remove week
  m.Ba_sowndiv_p5 <- update(m.Ba_sowndiv_p, 
                            bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                               hu ~ 1), 
                            seed = SEED)
  summary(m.Ba_sowndiv_p5, prob=0.9)
  pp_check(m.Ba_sowndiv_p5, ndraws=100)+xlim(0,500)
  #compare models
  
  loo.Ba <- loo(m.Ba_sowndiv_p, m.Ba_sowndiv_p2, m.Ba_sowndiv_p31, m.Ba_sowndiv_p32, m.Ba_sowndiv_p4,
                m.Ba_sowndiv_p5)
  loo.Ba
  
  save(m.Ba_sowndiv_p, m.Ba_sowndiv_p2, m.Ba_sowndiv_p31, m.Ba_sowndiv_p32, m.Ba_sowndiv_p4,
       m.Ba_sowndiv_p5,
       file="./statistics/brms/240108_Ba_sowndiv_priors.RData")
  
  rm(m.Ba_sowndiv_p, m.Ba_sowndiv_p2, m.Ba_sowndiv_p31, m.Ba_sowndiv_p32, m.Ba_sowndiv_p4,
     m.Ba_sowndiv_p5)  
  

####Fu ~ sowndiv ####
    #stepwise elimination of non-significant terms: p4, as week is marginally significant
    #looic: p5, as it is the most parsimonious model and all elpd_diff lay in range of 2 SE_diff
    
    load("./statistics/brms/231219_Fu_sowndiv_priors.RData")
    beta_coeff_priors <- prior(normal(0,10), class = "b")  
    SEED = 22061996
    
    sum(subset(dat, week=="W1")$Fu_per100g == 0) #3
    sum(subset(dat, week=="W2")$Fu_per100g == 0) #1
      #use hu~1, as too little zeros to estimate anything
    
    #for both weeks  
    m.Fu_sowndiv_p <- brm(
      bf(Fu_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
         hu ~ 1),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #all good 
    
    summary(m.Fu_sowndiv_p, prob =0.9)
    pp_check(m.Fu_sowndiv_p, ndraws=100)+
      xlim(0,2000)
    
    #remove 3way interaction:
    m.Fu_sowndiv_p2 <- update(m.Fu_sowndiv_p, 
                              bf(Fu_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                   (1|block/plot),
                                 hu ~ 1), 
                              seed = SEED) #4 div
    summary(m.Fu_sowndiv_p2, prob =0.9)
    
    #remove sowndivLogStd*week
    m.Fu_sowndiv_p31 <- update(m.Fu_sowndiv_p, 
                               bf(Fu_per100g ~ sowndivLogStd*treatment + treatment*week + (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED) # 1 div
    summary(m.Fu_sowndiv_p31, prob =0.9)
    
    #remove treatment*week 
    m.Fu_sowndiv_p32 <- update(m.Fu_sowndiv_p, 
                               bf(Fu_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED)
    summary(m.Fu_sowndiv_p32, prob =0.9)
    
    
    #remove sowndivLogStd*week and treatment*week 
    m.Fu_sowndiv_p4 <- update(m.Fu_sowndiv_p, 
                              bf(Fu_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                                 hu ~ 1), 
                              seed = SEED)
    summary(m.Fu_sowndiv_p4, prob=0.9) #stepwise: keep week, as it is marginally siginficant at 90% CI
    
    #remove week
    m.Fu_sowndiv_p5 <- update(m.Fu_sowndiv_p, 
                              bf(Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                 hu ~ 1), 
                              seed = SEED)
    summary(m.Fu_sowndiv_p5, prob=0.9)
    
    
    loo.Fu <- loo(m.Fu_sowndiv_p, m.Fu_sowndiv_p2, m.Fu_sowndiv_p31, m.Fu_sowndiv_p32, m.Fu_sowndiv_p4,
                  m.Fu_sowndiv_p5)
    
    loo.Fu
    
    save(m.Fu_sowndiv_p, m.Fu_sowndiv_p2, m.Fu_sowndiv_p31, m.Fu_sowndiv_p32, m.Fu_sowndiv_p4,
         m.Fu_sowndiv_p5,
         file="./statistics/brms/231219_Fu_sowndiv_priors.RData")  
    
    rm(m.Fu_sowndiv_p, m.Fu_sowndiv_p2, m.Fu_sowndiv_p31, m.Fu_sowndiv_p32,
       m.Fu_sowndiv_p5)
    

####Pl ~ sowndiv ####
    #stepwise elimination of non-significant terms: p4, as week is significant
    #looic: p4, as it is the most parsimonous model with a elpd_se_diff of less than 2 SE
    
    load("./statistics/brms/231219_Pl_sowndiv_priors.RData")
    beta_coeff_priors <- prior(normal(0,10), class = "b")  
    SEED = 22061996
    
    #for both weeks  
    m.Pl_sowndiv_p <- brm(
      bf(Pl_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
         hu ~ 1),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #5 div
    summary(m.Pl_sowndiv_p, prob =0.9)
    
    pp_check(m.Pl_sowndiv_p, ndraws=100)+
      xlim(0,2000)
    
    #remove 3way interaction:
    m.Pl_sowndiv_p2 <- update(m.Pl_sowndiv_p, 
                              bf(Pl_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                   (1|block/plot),
                                 hu ~ 1), 
                              seed = SEED)
    summary(m.Pl_sowndiv_p2, prob =0.9)
    
    
    #remove sowndivLogStd*week
    m.Pl_sowndiv_p31 <- update(m.Pl_sowndiv_p, 
                               bf(Pl_per100g ~ sowndivLogStd*treatment + treatment*week + (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED)
    summary(m.Pl_sowndiv_p31, prob =0.9)
    
    #remove treatment*week 
    m.Pl_sowndiv_p32 <- update(m.Pl_sowndiv_p, 
                               bf(Pl_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED)
    summary(m.Pl_sowndiv_p32, prob =0.9)
    
    #remove sowndivLogStd*week and treatment*week 
    m.Pl_sowndiv_p4 <- update(m.Pl_sowndiv_p, 
                              bf(Pl_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                                 hu ~ 1), 
                              seed = SEED)
    summary(m.Pl_sowndiv_p4, prob=0.95)
    #stepwise simplification; keep week as it is significant when considering a 95% CI
    
    #remove week
    m.Pl_sowndiv_p5 <- update(m.Pl_sowndiv_p, 
                              bf(Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                 hu ~ 1), 
                              seed = SEED)
    summary(m.Pl_sowndiv_p5, prob=0.9)
    
    
    loo.Pl <- loo(m.Pl_sowndiv_p, m.Pl_sowndiv_p2, m.Pl_sowndiv_p31, m.Pl_sowndiv_p32, m.Pl_sowndiv_p4,
                  m.Pl_sowndiv_p5 )
    loo.Pl
    
    save(m.Pl_sowndiv_p, m.Pl_sowndiv_p2, m.Pl_sowndiv_p31, m.Pl_sowndiv_p32, m.Pl_sowndiv_p4,
         m.Pl_sowndiv_p5,
         file="./statistics/brms/231219_Pl_sowndiv_priors.RData")  
    
    rm(m.Pl_sowndiv_p, m.Pl_sowndiv_p2, m.Pl_sowndiv_p31, m.Pl_sowndiv_p32, 
        m.Pl_sowndiv_p5)
        

####Pr ~ sowndiv ####
    #
    
    load("./statistics/brms/231219_Pr_sowndiv_priors.RData")
    beta_coeff_priors <- prior(normal(0,10), class = "b")  
    SEED = 22061996
    
    #for both weeks  
    m.Pr_sowndiv_p <- brm(
      bf(Pr_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #all good 
    
    summary(m.Pr_sowndiv_p, prob=0.9)
    pp_check(m.Pr_sowndiv_p, ndraws=100)+
      xlim(0,2000)
    
    #remove 3way interaction:
    m.Pr_sowndiv_p2 <- update(m.Pr_sowndiv_p, 
                              bf(Pr_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                   (1|block/plot),
                                 hu ~ sowndivLogStd*treatment + (1|block/plot)), 
                              seed = SEED)
    
    summary(m.Pr_sowndiv_p2, prob=0.9)
    
    #remove sowndivLogStd*week
    m.Pr_sowndiv_p31 <- update(m.Pr_sowndiv_p, 
                               bf(Pr_per100g ~ sowndivLogStd*treatment + treatment*week + (1|block/plot),
                                  hu ~ sowndivLogStd*treatment + (1|block/plot)), 
                               seed = SEED)
    
    summary(m.Pr_sowndiv_p31, prob=0.9)
    
    #remove treatment*week 
    m.Pr_sowndiv_p32 <- update(m.Pr_sowndiv_p, 
                               bf(Pr_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot),
                                  hu ~ sowndivLogStd*treatment + (1|block/plot)), 
                               seed = SEED)
    summary(m.Pr_sowndiv_p32, prob=0.9)
    
    #remove sowndivLogStd*week and treatment*week 
    m.Pr_sowndiv_p4 <- update(m.Pr_sowndiv_p, 
                              bf(Pr_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                                 hu ~ sowndivLogStd*treatment + (1|block/plot)), 
                              seed = SEED)
    summary(m.Pr_sowndiv_p4, prob=0.9)
    
    #remove week
    m.Pr_sowndiv_p5 <- update(m.Pr_sowndiv_p, 
                              bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                 hu ~ sowndivLogStd*treatment + (1|block/plot)), 
                              seed = SEED)
    summary(m.Pr_sowndiv_p5, prob=0.9)
    
    loo.Pr <- loo(m.Pr_sowndiv_p, m.Pr_sowndiv_p2, m.Pr_sowndiv_p31, m.Pr_sowndiv_p32, m.Pr_sowndiv_p4,
                  m.Pr_sowndiv_p5)
    
    loo.Pr
    
    save(m.Pr_sowndiv_p, m.Pr_sowndiv_p2, m.Pr_sowndiv_p31, m.Pr_sowndiv_p32, m.Pr_sowndiv_p4,
         m.Pr_sowndiv_p5,
         file="./statistics/brms/240221_Pr_sowndiv_priors.RData")  
    
    rm(m.Pr_sowndiv_p, m.Pr_sowndiv_p2, m.Pr_sowndiv_p31, m.Pr_sowndiv_p32, m.Pr_sowndiv_p4,
       m.Pr_sowndiv_p5)
    

####Om ~ sowndiv ####
   #looic: p7 is most parsimonous with a elpd that is not significantly worse (more than 2 SE elpd diff), while elpd_diff < 4
   #stepwise term removal at 90% CI: p6, as hu~treatment is marginally significant
   
   load("./statistics/brms/231219_Om_sowndiv_priors.RData")
   beta_coeff_priors <- prior(normal(0,10), class = "b")  
   SEED = 22061996
   
   #for both weeks  
   m.Om_sowndiv_p <- brm(
     bf(Om_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
        hu ~ sowndivLogStd*treatment + (1|block/plot)),
     data = dat, 
     prior = beta_coeff_priors,
     family = hurdle_lognormal,
     chains = 3,
     cores = 3,
     iter = 2000, warmup = 1000,
     seed = SEED,
     control = list(adapt_delta=0.99)) #all good 
   
   summary(m.Om_sowndiv_p)
   pp_check(m.Om_sowndiv_p, ndraws=100)+
     xlim(0,100)
   
   #remove 3way interaction:
   m.Om_sowndiv_p2 <- update(m.Om_sowndiv_p, 
                             bf(Om_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                  (1|block/plot),
                                hu ~ sowndivLogStd*treatment + (1|block/plot)), 
                             seed = SEED)
   summary(m.Om_sowndiv_p2)
   
   #remove sowndivLogStd*week
   m.Om_sowndiv_p31 <- update(m.Om_sowndiv_p, 
                              bf(Om_per100g ~ sowndivLogStd*treatment + treatment*week + (1|block/plot),
                                 hu ~ sowndivLogStd*treatment + (1|block/plot)), 
                              seed = SEED)
   summary(m.Om_sowndiv_p31)
   
   #remove treatment*week 
   m.Om_sowndiv_p32 <- update(m.Om_sowndiv_p, 
                              bf(Om_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot),
                                 hu ~ sowndivLogStd*treatment + (1|block/plot)), 
                              seed = SEED)
   summary(m.Om_sowndiv_p32)
   
   #remove sowndivLogStd*week and treatment*week 
   m.Om_sowndiv_p4 <- update(m.Om_sowndiv_p, 
                             bf(Om_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                                hu ~ sowndivLogStd*treatment + (1|block/plot)), 
                             seed = SEED)
   summary(m.Om_sowndiv_p4, prob=0.9)
   
   #remove week
   m.Om_sowndiv_p5 <- update(m.Om_sowndiv_p, 
                             bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                hu ~ sowndivLogStd*treatment + (1|block/plot)), 
                             seed = SEED)
   summary(m.Om_sowndiv_p5, prob=0.9)
   
  
   summary(m.Om_sowndiv_p8, prob=0.9)
   
   
   loo.Om <- loo(m.Om_sowndiv_p, m.Om_sowndiv_p2, m.Om_sowndiv_p31, m.Om_sowndiv_p32, m.Om_sowndiv_p4,
                 m.Om_sowndiv_p5)
   pp_check(m.Om_sowndiv_p5, ndraws = 100)
   
   loo.Om
   
   save(m.Om_sowndiv_p, m.Om_sowndiv_p2, m.Om_sowndiv_p31, m.Om_sowndiv_p32, m.Om_sowndiv_p4,
        m.Om_sowndiv_p5,
        file="./statistics/brms/240221_Om_sowndiv_priors.RData")  
   
   rm(m.Om_sowndiv_p, m.Om_sowndiv_p2, m.Om_sowndiv_p31, m.Om_sowndiv_p32, m.Om_sowndiv_p4,
      m.Om_sowndiv_p5, m.Om_sowndiv_p7, m.Om_sowndiv_p8)
  
   
#### model selection ####
   library(brms)
   library(ggplot2)
   
   #selection criteria: picking the most parsimonious model which has an elpd_diff > -4 to the model with the best fit
   #saving the selected models in a seperate file:
   
   #Ba ~ sowndiv
   load("./statistics/brms/240108_Ba_sowndiv_priors.RData")
   loo.Ba <- loo(m.Ba_sowndiv_p, m.Ba_sowndiv_p2, m.Ba_sowndiv_p31, m.Ba_sowndiv_p32, m.Ba_sowndiv_p4,
                 m.Ba_sowndiv_p5)
   loo.Ba
   
   rm(m.Ba_sowndiv_p, m.Ba_sowndiv_p2, m.Ba_sowndiv_p31, m.Ba_sowndiv_p32, m.Ba_sowndiv_p4) 
   
   #Fu ~ realdviv
   load("./statistics/brms/231219_Fu_sowndiv_priors.RData")  
   loo.Fu <- loo(m.Fu_sowndiv_p, m.Fu_sowndiv_p2, m.Fu_sowndiv_p31, m.Fu_sowndiv_p32, m.Fu_sowndiv_p4,
                 m.Fu_sowndiv_p5)
   loo.Fu
   
   rm(m.Fu_sowndiv_p, m.Fu_sowndiv_p2, m.Fu_sowndiv_p31, m.Fu_sowndiv_p32, m.Fu_sowndiv_p4)
   
   #Pl ~ sowndiv
   load("./statistics/brms/231219_Pl_sowndiv_priors.RData")  
   loo.Pl <- loo(m.Pl_sowndiv_p, m.Pl_sowndiv_p2, m.Pl_sowndiv_p31, m.Pl_sowndiv_p32, m.Pl_sowndiv_p4,
                 m.Pl_sowndiv_p5 )
   loo.Pl
   
   rm(m.Pl_sowndiv_p, m.Pl_sowndiv_p2, m.Pl_sowndiv_p31, m.Pl_sowndiv_p32, m.Pl_sowndiv_p4)
   
   #Pr ~ sowndiv
   load("./statistics/brms/240221_Pr_sowndiv_priors.RData")  
   loo.Pr <- loo(m.Pr_sowndiv_p, m.Pr_sowndiv_p2, m.Pr_sowndiv_p31, m.Pr_sowndiv_p32, m.Pr_sowndiv_p)
   loo.Pr
   
   rm(m.Pr_sowndiv_p, m.Pr_sowndiv_p2, m.Pr_sowndiv_p31, m.Pr_sowndiv_p32, m.Pr_sowndiv_p4,
      m.Pr_sowndiv_p5, m.Pr_sowndiv_p6, m.Pr_sowndiv_p8)
   
   #Om ~ sowndiv 
   load("./statistics/brms/240221_Om_sowndiv_priors.RData")  
   
   loo.Om <- loo(m.Om_sowndiv_p, m.Om_sowndiv_p2, m.Om_sowndiv_p31, m.Om_sowndiv_p32, m.Om_sowndiv_p4,
                 m.Om_sowndiv_p5)
   loo.Om
   
   rm(m.Om_sowndiv_p, m.Om_sowndiv_p2, m.Om_sowndiv_p31, m.Om_sowndiv_p32, m.Om_sowndiv_p4)
   
   #save the best fit models:
   save(m.Ba_sowndiv_p5, m.Fu_sowndiv_p5, m.Pl_sowndiv_p5, m.Om_sowndiv_p5, m.Pr_sowndiv_p5,
        file = "./statistics/brms/240221_TrophDens_sowndiv_mselect.RData")
   
####emmeans####
  
   emt = emtrends( m.Pr_sowndiv_p, specs = c("treatment", "week"), var="sowndivLogStd")
   summary(emt, point.est=mean, level = .9) 
   emt.pairs <- pairs(emt)
   summary(emt.pairs, point.est=mean, level = .9)
   bayestestR::p_direction(emt.pairs)
   
   Pr.emm <- emmeans(m.Pr_sowndiv_p, ~ sowndivLogStd * treatment * week)
   mvcontrast(Pr.emm, "pairwise", mult.name = c("treatment", "week"))
   