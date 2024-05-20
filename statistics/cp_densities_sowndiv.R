library(brms)
library(rstan)
library(ggplot2)
library(emmeans)
library(bayestestR)
library(dplyr)


####cp1 ~ sowndiv####
beta_coeff_priors <- prior(normal(0,10), class = "b")  
SEED = 22061996
sum(subset(dat, week=="W1")$cp1_per100g == 0) #61
sum(subset(dat, week=="W2")$cp1_per100g == 0) #48

#an appropriate distribution: 
cp1 <- dat %>% 
  filter(!is.na(cp1_per100g)) %>%
  select(cp1_per100g) 
plot(density(cp1$cp1_per100g))  
rm(SI)


#for both weeks  
m.cp1_sowndiv_p <- brm(
  bf(cp1_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
     hu ~  sowndivLogStd + treatment + week + (1|block/plot)),
  data = dat, 
  prior = beta_coeff_priors,
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #2div

pp_check(m.cp1_sowndiv_p, ndraws=100)+
  xlim(0,200)
summary(m.cp1_sowndiv_p, prob =0.9)

#check pairwise with emmeans:
emtrends(m.cp1_sowndiv_p, specs = c("treatment", "week"), var="sowndivLogStd") %>%
  summary() 

#remove 3way interaction:
m.cp1_sowndiv_p2 <- update(m.cp1_sowndiv_p, 
                          bf(cp1_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                               (1|block/plot),
                             hu ~  sowndivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED) #2 div
summary(m.cp1_sowndiv_p2, prob =0.9)

emtrends(m.cp1_sowndiv_p2, specs = c("treatment", "week"), var="sowndivLogStd") %>%
  summary() #CIs overlap

#remove sowndivLogStd*week
m.cp1_sowndiv_p31 <- update(m.cp1_sowndiv_p, 
                           bf(cp1_per100g ~ sowndivLogStd*treatment + treatment*week + (1|block/plot),
                              hu ~  sowndivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED) #4 div

summary(m.cp1_sowndiv_p31, prob=0.9)
emtrends(m.cp1_sowndiv_p31, specs = c("treatment", "week"), var="sowndivLogStd") %>%
  summary()

#remove treatment*week 
m.cp1_sowndiv_p32 <- update(m.cp1_sowndiv_p, 
                           bf(cp1_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot),
                              hu ~  sowndivLogStd + treatment + week + (1|block/plot)), 
                           seed = SEED)
summary(m.cp1_sowndiv_p32, prob=0.9)

emtrends(m.cp1_sowndiv_p32, specs = c("treatment", "week"), var="sowndivLogStd") %>%
  summary()

#remove sowndivLogStd:week
m.cp1_sowndiv_p4 <- update(m.cp1_sowndiv_p, 
                          bf(cp1_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                             hu ~  sowndivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED) #9 div
summary(m.cp1_sowndiv_p4, prob=0.9)
pp_check(m.cp1_sowndiv_p4, ndraws=100)+xlim(0,500)

#remove week
m.cp1_sowndiv_p5 <- update(m.cp1_sowndiv_p, 
                          bf(cp1_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                             hu ~  sowndivLogStd + treatment + week + (1|block/plot)), 
                          seed = SEED) #0 div
summary(m.cp1_sowndiv_p5, prob=0.9)
pp_check(m.cp1_sowndiv_p5, ndraws=100)+xlim(0,500)

#remove hu ~ sowndivLogStd
m.cp1_sowndiv_p6 <- update(m.cp1_sowndiv_p, 
                           bf(cp1_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                              hu ~ treatment + week + (1|block/plot)), 
                           seed = SEED)
summary(m.cp1_sowndiv_p6, prob=0.9)
pp_check(m.cp1_sowndiv_p6, ndraws=100)+xlim(0,500)

#remove hu ~ week
m.cp1_sowndiv_p7 <- update(m.cp1_sowndiv_p, 
                           bf(cp1_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                              hu ~ treatment + (1|block/plot)), 
                           seed = SEED)
summary(m.cp1_sowndiv_p7, prob=0.9) 
pp_check(m.cp1_sowndiv_p7, ndraws=100)+xlim(0,500)

#remove hu ~ treatment
m.cp1_sowndiv_p8 <- update(m.cp1_sowndiv_p, 
                           bf(cp1_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                              hu ~ 1), 
                           seed = SEED) #1 div
summary(m.cp1_sowndiv_p8, prob=0.9)
pp_check(m.cp1_sowndiv_p8, ndraws=100)+xlim(0,500)
emtrends(m.cp1_sowndiv_p8, specs = c("treatment"), var="sowndivLogStd") %>%
  summary() 


#compare models

loo.cp1 <- loo(m.cp1_sowndiv_p, m.cp1_sowndiv_p2, m.cp1_sowndiv_p31, m.cp1_sowndiv_p32, m.cp1_sowndiv_p4,
              m.cp1_sowndiv_p5, m.cp1_sowndiv_p6, m.cp1_sowndiv_p7, m.cp1_sowndiv_p8)
loo.cp1
#choose p8

save(m.cp1_sowndiv_p, m.cp1_sowndiv_p2, m.cp1_sowndiv_p31, m.cp1_sowndiv_p32, m.cp1_sowndiv_p4,
     m.cp1_sowndiv_p5, m.cp1_sowndiv_p6, m.cp1_sowndiv_p7, m.cp1_sowndiv_p8,
     file="./statistics/brms/240131_cp1_sowndiv_priors.RData")

#### cp2 ~ sowndiv ####
    beta_coeff_priors <- prior(normal(0,10), class = "b")  
    SEED = 22061996
    sum(subset(dat, week=="W1")$cp2_per100g == 0) #0
    sum(subset(dat, week=="W2")$cp2_per100g == 0) #1
    
    #an appropriate distribution: 
    cp2 <- dat %>% 
      filter(!is.na(cp2_per100g)) %>%
      select(cp2_per100g) 
    plot(density(cp2$cp2_per100g))  
    rm(cp2)
    
    
    #for both weeks  
    m.cp2_sowndiv_p <- brm(
      bf(cp2_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
         hu ~ 1 ),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #13 div, bulk/tail ESS too low
    
    pp_check(m.cp2_sowndiv_p, ndraws=100)+
      xlim(0,200)
    summary(m.cp2_sowndiv_p, prob =0.9)
    
    #check pairwise with emmeans:
    emtrends(m.cp2_sowndiv_p, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary() 
    
    #remove 3way interaction:
    m.cp2_sowndiv_p2 <- update(m.cp2_sowndiv_p, 
                               bf(cp2_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                    (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED) #2 div
    summary(m.cp2_sowndiv_p2, prob =0.9)
    
    emtrends(m.cp2_sowndiv_p2, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary() #CIs overlap
    
    #remove sowndivLogStd*week
    m.cp2_sowndiv_p31 <- update(m.cp2_sowndiv_p, 
                                bf(cp2_per100g ~ sowndivLogStd*treatment + treatment*week + (1|block/plot),
                                   hu ~ 1 + (1|block/plot)), 
                                seed = SEED) #4 div
    
    summary(m.cp2_sowndiv_p31, prob=0.9)
    emtrends(m.cp2_sowndiv_p31, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary()
    
    #remove treatment*week 
    m.cp2_sowndiv_p32 <- update(m.cp2_sowndiv_p, 
                                bf(cp2_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot),
                                   hu ~ 1), 
                                seed = SEED)
    summary(m.cp2_sowndiv_p32, prob=0.9)
    
    emtrends(m.cp2_sowndiv_p32, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary()
    
    #remove sowndivLogStd:week
    m.cp2_sowndiv_p4 <- update(m.cp2_sowndiv_p, 
                               bf(cp2_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED)
    summary(m.cp2_sowndiv_p4, prob=0.9)
    pp_check(m.cp2_sowndiv_p4, ndraws=100)+xlim(0,500)
    
    #remove week
    m.cp2_sowndiv_p5 <- update(m.cp2_sowndiv_p, 
                               bf(cp2_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED)
    summary(m.cp2_sowndiv_p5, prob=0.9)
    pp_check(m.cp2_sowndiv_p5, ndraws=100)+xlim(0,500)
    
    emtrends(m.cp2_sowndiv_p5, specs = c("treatment"), var="sowndivLogStd") %>%
      summary()
    
    #compare models
    
    loo.cp2 <- loo(m.cp2_sowndiv_p, m.cp2_sowndiv_p2, m.cp2_sowndiv_p31, m.cp2_sowndiv_p32, m.cp2_sowndiv_p4,
                   m.cp2_sowndiv_p5)
    loo.cp2
    
    save(m.cp2_sowndiv_p, m.cp2_sowndiv_p2, m.cp2_sowndiv_p31, m.cp2_sowndiv_p32, m.cp2_sowndiv_p4,
         m.cp2_sowndiv_p5,
         file="./statistics/brms/240131_cp2_sowndiv_priors.RData")

#### cp3 ~ sowndiv ####
    beta_coeff_priors <- prior(normal(0,10), class = "b")  
    SEED = 22061996
    sum(subset(dat, week=="W1")$cp3_per100g == 0) #4
    sum(subset(dat, week=="W2")$cp3_per100g == 0) #0
    
    #an appropriate distribution: 
    cp3 <- dat %>% 
      filter(!is.na(cp3_per100g)) %>%
      select(cp3_per100g) 
    plot(density(cp3$cp3_per100g))  
    rm(cp3)
    
    
    #for both weeks  
    m.cp3_sowndiv_p <- brm(
      bf(cp3_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
         hu ~ 1 ),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #13 div, bulk/tail ESS too low
    
    pp_check(m.cp3_sowndiv_p, ndraws=100)+
      xlim(0,200)
    summary(m.cp3_sowndiv_p, prob =0.9)
    
    #as different orientation of sowndivLogStd:treatment2:weekW2 and sowndivLogStd:treatment3:weekW2 -->
    #check pairwise with emmeans:
    emtrends(m.cp3_sowndiv_p, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary() #CIs overlap
    
    #remove 3way interaction:
    m.cp3_sowndiv_p2 <- update(m.cp3_sowndiv_p, 
                               bf(cp3_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                    (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED) #2 div
    summary(m.cp3_sowndiv_p2, prob =0.9)
    
    emtrends(m.cp3_sowndiv_p2, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary() #CIs overlap
    
    #remove sowndivLogStd*week
    m.cp3_sowndiv_p31 <- update(m.cp3_sowndiv_p, 
                                bf(cp3_per100g ~ sowndivLogStd*treatment + treatment*week + (1|block/plot),
                                   hu ~ 1 + (1|block/plot)), 
                                seed = SEED) #4 div
    
    summary(m.cp3_sowndiv_p31, prob=0.9)
    emtrends(m.cp3_sowndiv_p31, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary()
    
    #remove treatment*week 
    m.cp3_sowndiv_p32 <- update(m.cp3_sowndiv_p, 
                                bf(cp3_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot),
                                   hu ~ 1), 
                                seed = SEED)
    summary(m.cp3_sowndiv_p32, prob=0.9)
    
    emtrends(m.cp3_sowndiv_p32, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary()
    
    #remove sowndivLogStd:week
    m.cp3_sowndiv_p4 <- update(m.cp3_sowndiv_p, 
                               bf(cp3_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED)
    summary(m.cp3_sowndiv_p4, prob=0.9)
    pp_check(m.cp3_sowndiv_p4, ndraws=100)+xlim(0,500)
    
    #remove week
    m.cp3_sowndiv_p5 <- update(m.cp3_sowndiv_p, 
                               bf(cp3_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED)
    summary(m.cp3_sowndiv_p5, prob=0.9)
    pp_check(m.cp3_sowndiv_p5, ndraws=100)+xlim(0,500)
    #compare models
    
    loo.cp3 <- loo(m.cp3_sowndiv_p, m.cp3_sowndiv_p2, m.cp3_sowndiv_p31, m.cp3_sowndiv_p32, m.cp3_sowndiv_p4,
                   m.cp3_sowndiv_p5)
    loo.cp3
    
    save(m.cp3_sowndiv_p, m.cp3_sowndiv_p2, m.cp3_sowndiv_p31, m.cp3_sowndiv_p32, m.cp3_sowndiv_p4,
         m.cp3_sowndiv_p5,
         file="./statistics/brms/240131_cp3_sowndiv_priors.RData")

    
    
#### cp4 ~ sowndiv ####
    beta_coeff_priors <- prior(normal(0,10), class = "b")  
    SEED = 22061996
    sum(subset(dat, week=="W1")$cp4_per100g == 0) #13
    sum(subset(dat, week=="W2")$cp4_per100g == 0) #2
    
    #an appropriate distribution: 
    cp4 <- dat %>% 
      filter(!is.na(cp4_per100g)) %>%
      select(cp4_per100g) 
    plot(density(cp4$cp4_per100g))  
    rm(cp4)
    
    
    #for both weeks  
    m.cp4_sowndiv_p <- brm(
      bf(cp4_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
         hu ~ 1 ),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #13 div, bulk/tail ESS too low
    
    pp_check(m.cp4_sowndiv_p, ndraws=100)+
      xlim(0,200)
    summary(m.cp4_sowndiv_p, prob =0.9)
    
    #as different orientation of sowndivLogStd:treatment2:weekW2 and sowndivLogStd:treatment3:weekW2 -->
    #check pairwise with emmeans:
    emtrends(m.cp4_sowndiv_p, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary() #CIs overlap
    
    #remove 3way interaction:
    m.cp4_sowndiv_p2 <- update(m.cp4_sowndiv_p, 
                               bf(cp4_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + treatment*week + 
                                    (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED) #2 div
    summary(m.cp4_sowndiv_p2, prob =0.9)
    
    emtrends(m.cp4_sowndiv_p2, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary() #CIs overlap
    
    #remove sowndivLogStd*week
    m.cp4_sowndiv_p31 <- update(m.cp4_sowndiv_p, 
                                bf(cp4_per100g ~ sowndivLogStd*treatment + treatment*week + (1|block/plot),
                                   hu ~ 1 + (1|block/plot)), 
                                seed = SEED) #4 div
    
    summary(m.cp4_sowndiv_p31, prob=0.9)
    emtrends(m.cp4_sowndiv_p31, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary()
    
    #remove treatment*week 
    m.cp4_sowndiv_p32 <- update(m.cp4_sowndiv_p, 
                                bf(cp4_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot),
                                   hu ~ 1), 
                                seed = SEED)
    summary(m.cp4_sowndiv_p32, prob=0.9)
    
    emtrends(m.cp4_sowndiv_p32, specs = c("treatment", "week"), var="sowndivLogStd") %>%
      summary()
    
    #remove sowndivLogStd:week
    m.cp4_sowndiv_p4 <- update(m.cp4_sowndiv_p, 
                               bf(cp4_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED)
    summary(m.cp4_sowndiv_p4, prob=0.9)
    pp_check(m.cp4_sowndiv_p4, ndraws=100)+xlim(0,500)
    
    #remove week
    m.cp4_sowndiv_p5 <- update(m.cp4_sowndiv_p, 
                               bf(cp4_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                  hu ~ 1), 
                               seed = SEED)
    summary(m.cp4_sowndiv_p5, prob=0.9)
    pp_check(m.cp4_sowndiv_p5, ndraws=100)+xlim(0,1000)
    emtrends(m.cp4_sowndiv_p5, specs = c("treatment"), var="sowndivLogStd") %>%
      summary()
    
    #compare models
    
    loo.cp4 <- loo(m.cp4_sowndiv_p, m.cp4_sowndiv_p2, m.cp4_sowndiv_p31, m.cp4_sowndiv_p32, m.cp4_sowndiv_p4,
                   m.cp4_sowndiv_p5)
    loo.cp4
    
    save(m.cp4_sowndiv_p, m.cp4_sowndiv_p2, m.cp4_sowndiv_p31, m.cp4_sowndiv_p32, m.cp4_sowndiv_p4,
         m.cp4_sowndiv_p5,
         file="./statistics/brms/240131_cp4_sowndiv_priors.RData")
    

#### cp5 ~ sowndiv ####
    beta_coeff_priors <- prior(normal(0,10), class = "b")  
    SEED = 22061996
    sum(subset(dat, week=="W1")$cp5_per100g == 0) #105
    sum(subset(dat, week=="W2")$cp5_per100g == 0) #96
    
    #an appropriate distribution: 
    cp5 <- dat %>% 
      filter(!is.na(cp5_per100g)) %>%
      select(cp5_per100g) 
    plot(density(cp5$cp5_per100g))  
    
    #this cannot be fit
    
#### select best fit models ####
    load("./statistics/brms/240131_cp4_sowndiv_priors.RData")
    loo.cp4 <- loo(m.cp4_sowndiv_p, m.cp4_sowndiv_p2, m.cp4_sowndiv_p31, m.cp4_sowndiv_p32, m.cp4_sowndiv_p4,
                   m.cp4_sowndiv_p5)
    loo.cp4 #p5
    
    load("./statistics/brms/240131_cp3_sowndiv_priors.RData")
    loo.cp3 <- loo(m.cp3_sowndiv_p, m.cp3_sowndiv_p2, m.cp3_sowndiv_p31, m.cp3_sowndiv_p32, m.cp3_sowndiv_p4,
                   m.cp3_sowndiv_p5)
    loo.cp3 #p5
     
    load("./statistics/brms/240131_cp2_sowndiv_priors.RData")
    loo.cp2 <- loo(m.cp2_sowndiv_p, m.cp2_sowndiv_p2, m.cp2_sowndiv_p31, m.cp2_sowndiv_p32, m.cp2_sowndiv_p4,
                   m.cp2_sowndiv_p5)
    loo.cp2 #p5
    
    load("./statistics/brms/240131_cp1_sowndiv_priors.RData")
    loo.cp1 <- loo(m.cp1_sowndiv_p, m.cp1_sowndiv_p2, m.cp1_sowndiv_p31, m.cp1_sowndiv_p32, m.cp1_sowndiv_p4,
                   m.cp1_sowndiv_p5)
    loo.cp1 #p5
    
    
    save(m.cp1_sowndiv_p5, m.cp2_sowndiv_p5, m.cp3_sowndiv_p5, m.cp4_sowndiv_p5,
         file="./statistics/brms/240205_cp_sowndiv_mselect.RData")

   