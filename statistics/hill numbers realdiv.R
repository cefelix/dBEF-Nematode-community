####data and packages####
library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(hexbin)
library(GGally)
library(modelbased)
#rule of thumb for p_dir: possibly existent at p_dir > 95%

#looic model selection: 
#if elpd_diff < 4 models have very similar predictive performance and SE underestimation doesnt matter
#--> choose the most parsimonious model with elpd_diff > 4 when compared to best fit model
#background FAQ 16 of https://users.aalto.fi/~ave/CV-FAQ.html#16_How_to_interpret_in_Standard_error_(SE)_of_elpd_difference_(elpd_diff)


# a seed:
SEED = 19111996

#load files: 
#load(file = "./statistics/brms/231127_Hill.RData" )

#exclude 60 sp.:
dat <- subset(dBEF_nem21, realdiv != 60) 
#standardize:  
dat <- dat %>% mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
                      .after = realdivLog) %>%
  mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
         .after = realdivLog)

datW1 <- subset(dat, week=="W1")
  datW1.test <- subset(dat, week=="W1" & treatment == "1" & sowndiv == 16)
  nrow(datW1.test) #7/8 reps
  rm(datW1.test)
datW2 <- subset(dat, week=="W2")

#priors    
beta_coeff_priors <- prior(normal(0,10), class = "b")  

#### Shannon all ~ realdivLogStd: gaus_p4 (looic: gauss_p5) ####
#all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
# stepwise elimiantion: gaus_p4, as for t2w1-t2w2 : p_dir = 93.73 % < 95 %
#looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
#differ significantly (less than 2 SE) from the best fit model
load("./statistics/brms/231221_hill_all_realdiv.RData")

#using a gaussian distribution with a normal(0,10) prior for beta coefficients
m.all.Shannon.gaus_p <- brm(Hill_q1 ~ realdivLogStd*treatment*week + (1|block/plot), 
                            data = dat, family = "gaussian",
                            chains = 3,
                            cores = 3,
                            iter = 2000, warmup = 1000,
                            prior = beta_coeff_priors,
                            seed = SEED,
                            control = list(adapt_delta = 0.99) )  #4 div
summary(m.all.Shannon.gaus_p, prob=0.9 )

#check with emtrends:
    emt = emtrends( m.all.Shannon.gaus_p, specs = c("treatment", "week"), var="realdivLogStd")
    summary(emt, point.est=mean, level = .9) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .9)
    bayestestR::p_direction(emt.pairs) #t2w1-t2w2 : p_dir = 93.73 % > 90 %
    #Thus --> simplify further

#remove 3 way interaction:
m.all.Shannon.gaus_p2 <- update(m.all.Shannon.gaus_p, 
                                bf(Hill_q1 ~ realdivLogStd*treatment + treatment*week + 
                                     realdivLogStd*week + (1|block/plot)),
                                seed = SEED) # 8 div
summary(m.all.Shannon.gaus_p2, prob=0.9)

#remove treatment*week
m.all.Shannon.gaus_p31 <- update(m.all.Shannon.gaus_p, 
                                 bf(Hill_q1 ~ realdivLogStd*treatment + realdivLogStd*week + 
                                      (1|block/plot)),
                                 seed = SEED) #3 div
summary(m.all.Shannon.gaus_p31, prob=0.9 )

#remove realdiv*week
m.all.Shannon.gaus_p32 <- update(m.all.Shannon.gaus_p, 
                                 bf(Hill_q1 ~ realdivLogStd*treatment + treatment*week + 
                                      (1|block/plot)),
                                 seed = SEED) # 4 div
summary(m.all.Shannon.gaus_p32, prob=0.9 )

#remove both treatment*week and realdiv*week:
m.all.Shannon.gaus_p4 <- update(m.all.Shannon.gaus_p, 
                                bf(Hill_q1 ~ realdivLogStd*treatment + week + (1|block/plot)),
                                seed = SEED)
summary(m.all.Shannon.gaus_p4, prob=0.95 )
#stepwise simplification: choose p4, as week is significant (95% CI)
pp_check(m.all.Shannon.gaus_p4, ndraws=100)

#remove week:
m.all.Shannon.gaus_p5 <- update(m.all.Shannon.gaus_p, 
                                bf(Hill_q1 ~ realdivLogStd*treatment + (1|block/plot)),
                                seed = SEED)
summary(m.all.Shannon.gaus_p5, prob=0.9 )

loo(m.all.Shannon.gaus_p5,m.all.Shannon.gaus_p4, m.all.Shannon.gaus_p32, m.all.Shannon.gaus_p31,
    m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p)

#using a gamma distribution with a normal(0,10) prior for beta coefficients
m.all.Shannon.gamma_p <- brm(Hill_q1 ~ realdivLogStd*treatment*week + (1|block/plot), 
                             data = dat, family = "gamma",
                             chains = 3,
                             cores = 3,
                             iter = 2000, warmup = 1000,
                             prior = beta_coeff_priors,
                             seed = SEED,
                             control = list(adapt_delta = 0.99) )  # 7 div
summary(m.all.Shannon.gamma_p, prob=0.9)



#remove 3 way interaction:
m.all.Shannon.gamma_p2 <- update(m.all.Shannon.gamma_p, 
                                 bf(Hill_q1 ~ realdivLogStd*treatment + treatment*week + 
                                      realdivLogStd*week + (1|block/plot)),
                                 seed = SEED) # 7 div
summary(m.all.Shannon.gamma_p2, prob=0.9)

#remove treatment*week
m.all.Shannon.gamma_p31 <- update(m.all.Shannon.gamma_p, 
                                  bf(Hill_q1 ~ realdivLogStd*treatment + realdivLogStd*week + 
                                       (1|block/plot)),
                                  seed = SEED) #8 div
summary(m.all.Shannon.gamma_p31, prob=0.9 )

#remove realdiv*week
m.all.Shannon.gamma_p32 <- update(m.all.Shannon.gamma_p, 
                                  bf(Hill_q1 ~ realdivLogStd*treatment + treatment*week + 
                                       (1|block/plot)),
                                  seed = SEED) #8 div
summary(m.all.Shannon.gamma_p32, prob=0.9 )

#remove both treatment*week and realdiv*week:
m.all.Shannon.gamma_p4 <- update(m.all.Shannon.gamma_p, 
                                 bf(Hill_q1 ~ realdivLogStd*treatment + week + (1|block/plot)),
                                 seed = SEED) #4 div
summary(m.all.Shannon.gamma_p4, prob=0.9 )
#stepwise simplification keep week, as it is marginally significant (CI=90%)
pp_check(m.all.Shannon.gamma_p4, ndraws=100)


#remove week:
m.all.Shannon.gamma_p5 <- update(m.all.Shannon.gamma_p, 
                                 bf(Hill_q1 ~ realdivLogStd*treatment + (1|block/plot)),
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
     file="./statistics/brms/231221_hill_all_realdiv.RData")

rm(m.all.Shannon.gamma_p5, m.all.Shannon.gamma_p4, m.all.Shannon.gamma_p32,m.all.Shannon.gamma_p31,
   m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p,
   m.all.Shannon.gaus_p5,m.all.Shannon.gaus_p4, m.all.Shannon.gaus_p32, m.all.Shannon.gaus_p31,
   m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p)

#### Shannon Ba ~ realdivLogStd: gaus_p (looic: gaus_p5, gamma_p5 is almost similar) ####
#all gaussian models show approximately same fit than the gamma models (elpd_diff > 2 SE_diff)
#stepwise elimiantion: gaus_p, as p_dir of t2w1-t2w2 is 99.87 % (> 95 %) 
#looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
#differ significantly (less than 2 SE) from the best fit model
load("./statistics/brms/231221_hill_Ba_realdiv.RData")

#using a gaussian distribution with a normal(0,10) prior for beta coefficients
m.Ba.Shannon.gaus_p <- brm(Hill_q1.Ba ~ realdivLogStd*treatment*week + (1|block/plot), 
                           data = dat, family = "gaussian",
                           chains = 3,
                           cores = 3,
                           iter = 2000, warmup = 1000,
                           prior = beta_coeff_priors,
                           seed = SEED,
                           control = list(adapt_delta = 0.99) )  #3 div
summary(m.Ba.Shannon.gaus_p, prob=0.9 )
  #3-way interaction is significant, so lets compare all week-treatment combinations:
  emt = emtrends( m.Ba.Shannon.gaus_p, specs = c("treatment", "week"), var="realdivLogStd")
  summary(emt, point.est=mean, level = .95) #t2 in w1 differs significantly from t2 in w2 (95% CI)
  emt.pairs <- pairs(emt)
  summary(emt.pairs, point.est=mean, level = .9)
  bayestestR::p_direction(emt.pairs) #probability of direction 99.87 % between t2w1 and t2w2
  
  #comparing week-treatment combinations using modelbased library (doesnt assess slopes but means)
  #estimate_contrasts(m.Ba.Shannon.gaus_p, contrast = c("treatment", "week"), ci=0.9)
  
  pp_check(m.Ba.Shannon.gaus_p, ndraws=100)


#keep three-way interaction, as it is 

#remove 3 way interaction:
m.Ba.Shannon.gaus_p2 <- update(m.Ba.Shannon.gaus_p, 
                               bf(Hill_q1.Ba ~ realdivLogStd*treatment + treatment*week + 
                                    realdivLogStd*week + (1|block/plot)),
                               seed = SEED) # 14 div
summary(m.Ba.Shannon.gaus_p2, prob=0.9)

#remove treatment*week
m.Ba.Shannon.gaus_p31 <- update(m.Ba.Shannon.gaus_p, 
                                bf(Hill_q1.Ba ~ realdivLogStd*treatment + realdivLogStd*week + 
                                     (1|block/plot)),
                                seed = SEED) #12 div
summary(m.Ba.Shannon.gaus_p31, prob=0.95 ) 
#keep realdiv:week interaction, as it is significant (CI=95%)
pp_check(m.Ba.Shannon.gaus_p31, ndraws=100)

#not to consider when following stepwise simplification:
    #remove realdiv*week
    m.Ba.Shannon.gaus_p32 <- update(m.Ba.Shannon.gaus_p, 
                                    bf(Hill_q1.Ba ~ realdivLogStd*treatment + treatment*week + 
                                         (1|block/plot)),
                                    seed = SEED) # 4 div
    summary(m.Ba.Shannon.gaus_p32, prob=0.9 )
    
    #remove both treatment*week and realdiv*week:
    m.Ba.Shannon.gaus_p4 <- update(m.Ba.Shannon.gaus_p, 
                                   bf(Hill_q1.Ba ~ realdivLogStd*treatment + week + (1|block/plot)),
                                   seed = SEED)
    summary(m.Ba.Shannon.gaus_p4, prob=0.95 )
    #stepwise simplification: choose p4, as week is significant (95% CI)
    
    #remove week:
    m.Ba.Shannon.gaus_p5 <- update(m.Ba.Shannon.gaus_p, 
                                   bf(Hill_q1.Ba ~ realdivLogStd*treatment + (1|block/plot)),
                                   seed = SEED)
    summary(m.Ba.Shannon.gaus_p5, prob=0.9 )
    
    loo(m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
        m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p)

#using a gamma distribution with a normal(0,10) prior for beta coefficients
m.Ba.Shannon.gamma_p <- brm(Hill_q1.Ba ~ realdivLogStd*treatment*week + (1|block/plot), 
                            data = dat, family = "gamma",
                            chains = 3,
                            cores = 3,
                            iter = 2000, warmup = 1000,
                            prior = beta_coeff_priors,
                            seed = SEED,
                            control = list(adapt_delta = 0.99) )  # 18 div
summary(m.Ba.Shannon.gamma_p, prob=0.9)
pp_check(m.Ba.Shannon.gamma_p, ndraws=100)
pp_check(m.Ba.Shannon.gaus_p, ndraws=100)


#as 3-way interaction is significant, compare all week-treatment combinations:
    emt = emtrends( m.Ba.Shannon.gamma_p, specs = c("treatment", "week"), var="realdivLogStd")
    summary(emt, point.est=mean, level = .95) #t2 in w1 differs significantly from t2 in w2 (95% CI)
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .9)
    bayestestR::p_direction(emt.pairs) #probability of direction 99.83 % between t2w1 and t2w2


#remove 3 way interaction:
m.Ba.Shannon.gamma_p2 <- update(m.Ba.Shannon.gamma_p, 
                                bf(Hill_q1.Ba ~ realdivLogStd*treatment + treatment*week + 
                                     realdivLogStd*week + (1|block/plot)),
                                seed = SEED) # 7 div
summary(m.Ba.Shannon.gamma_p2, prob=0.9)

#remove treatment*week
m.Ba.Shannon.gamma_p31 <- update(m.Ba.Shannon.gamma_p, 
                                 bf(Hill_q1.Ba ~ realdivLogStd*treatment + realdivLogStd*week + 
                                      (1|block/plot)),
                                 seed = SEED) #12 div
summary(m.Ba.Shannon.gamma_p31, prob=0.9 )
#if not considering 3-way interactions, keep p31 as realdiv:week is marginally significant!

#remove realdiv*week
m.Ba.Shannon.gamma_p32 <- update(m.Ba.Shannon.gamma_p, 
                                 bf(Hill_q1.Ba ~ realdivLogStd*treatment + treatment*week + 
                                      (1|block/plot)),
                                 seed = SEED) #4 div
summary(m.Ba.Shannon.gamma_p32, prob=0.9 )

#remove both treatment*week and realdiv*week:
m.Ba.Shannon.gamma_p4 <- update(m.Ba.Shannon.gamma_p, 
                                bf(Hill_q1.Ba ~ realdivLogStd*treatment + week + (1|block/plot)),
                                seed = SEED)
summary(m.Ba.Shannon.gamma_p4, prob=0.9 )
#stepwise simplification keep week, as it is marginally significant (CI=90%)
pp_check(m.Ba.Shannon.gamma_p4, ndraws=100)


#remove week:
m.Ba.Shannon.gamma_p5 <- update(m.Ba.Shannon.gamma_p, 
                                bf(Hill_q1.Ba ~ realdivLogStd*treatment + (1|block/plot)),
                                seed = SEED)
summary(m.Ba.Shannon.gamma_p5, prob=0.9 )
pp_check(m.Ba.Shannon.gamma_p5, ndraws=100)

#compare
Ba.loo <- loo(m.Ba.Shannon.gamma_p5, m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
              m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
              m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
              m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p)
Ba.loo
#it does not matter whether you use gaussian or gamma family:
  #at lower values for Shannon's H, gamma predicts better
  #at higher values, gaus predicts better

save(m.Ba.Shannon.gamma_p5, m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
     m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
     m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
     m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p,
     file="./statistics/brms/231221_hill_Ba_realdiv.RData")

rm(m.Ba.Shannon.gamma_p5, m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
   m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
   m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
   m.Ba.Shannon.gaus_p2#, m.Ba.Shannon.gaus_p
   )

#### Shannon Fu ~ realdivLogStd: gamma_p4 (looic: gamma_p5)  ####
#all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
# stepwise elimiantion: gamma_p4
#looic: gamma_p5, as it is the most parsimonous model to not
#differ significantly (less than 2 SE) from the best fit model
load(     file="./statistics/brms/231221_hill_Fu_realdiv.RData")

#using a gaussian distribution with a normal(0,10) prior for beta coefficients
m.Fu.Shannon.gaus_p <- brm(Hill_q1.Fu ~ realdivLogStd*treatment*week + (1|block/plot), 
                           data = dat, family = "gaussian",
                           chains = 3,
                           cores = 3,
                           iter = 2000, warmup = 1000,
                           prior = beta_coeff_priors,
                           seed = SEED,
                           control = list(adapt_delta = 0.99) )  #5 div
summary(m.Fu.Shannon.gaus_p, prob=0.9 )

#remove 3 way interaction:
m.Fu.Shannon.gaus_p2 <- update(m.Fu.Shannon.gaus_p, 
                               bf(Hill_q1.Fu ~ realdivLogStd*treatment + treatment*week + 
                                    realdivLogStd*week + (1|block/plot)),
                               seed = SEED) # 9 div
summary(m.Fu.Shannon.gaus_p2, prob=0.9)

#remove treatment*week
m.Fu.Shannon.gaus_p31 <- update(m.Fu.Shannon.gaus_p, 
                                bf(Hill_q1.Fu ~ realdivLogStd*treatment + realdivLogStd*week + 
                                     (1|block/plot)),
                                seed = SEED) #2 div
summary(m.Fu.Shannon.gaus_p31, prob=0.9 )

#remove realdiv*week
m.Fu.Shannon.gaus_p32 <- update(m.Fu.Shannon.gaus_p, 
                                bf(Hill_q1.Fu ~ realdivLogStd*treatment + treatment*week + 
                                     (1|block/plot)),
                                seed = SEED) # 14 div
summary(m.Fu.Shannon.gaus_p32, prob=0.9 )

#remove both treatment*week and realdiv*week:
m.Fu.Shannon.gaus_p4 <- update(m.Fu.Shannon.gaus_p, 
                               bf(Hill_q1.Fu ~ realdivLogStd*treatment + week + (1|block/plot)),
                               seed = SEED) #29 div
summary(m.Fu.Shannon.gaus_p4, prob=0.9 )
#stepwise simplification: choose p4, as week is marginally significant (90% CI)
  pp_check(m.Fu.Shannon.gaus_p4, ndraws=100) 
    #okayish, there are two many samples only one Fu-species

#remove week:
m.Fu.Shannon.gaus_p5 <- update(m.Fu.Shannon.gaus_p, 
                               bf(Hill_q1.Fu ~ realdivLogStd*treatment + (1|block/plot)),
                               seed = SEED) #0 div
summary(m.Fu.Shannon.gaus_p5, prob=0.9 )

loo(m.Fu.Shannon.gaus_p5,m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
    m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p)

#using a gamma distribution with a normal(0,10) prior for beta coefficients
m.Fu.Shannon.gamma_p <- brm(Hill_q1.Fu ~ realdivLogStd*treatment*week + (1|block/plot), 
                            data = dat, family = "gamma",
                            chains = 3,
                            cores = 3,
                            iter = 2000, warmup = 1000,
                            prior = beta_coeff_priors,
                            seed = SEED,
                            control = list(adapt_delta = 0.99) )  # 2 div
summary(m.Fu.Shannon.gamma_p, prob=0.9)

#remove 3 way interaction:
m.Fu.Shannon.gamma_p2 <- update(m.Fu.Shannon.gamma_p, 
                                bf(Hill_q1.Fu ~ realdivLogStd*treatment + treatment*week + 
                                     realdivLogStd*week + (1|block/plot)),
                                seed = SEED) # 5 div
summary(m.Fu.Shannon.gamma_p2, prob=0.9)

#remove treatment*week
m.Fu.Shannon.gamma_p31 <- update(m.Fu.Shannon.gamma_p, 
                                 bf(Hill_q1.Fu ~ realdivLogStd*treatment + realdivLogStd*week + 
                                      (1|block/plot)),
                                 seed = SEED) #12 div
summary(m.Fu.Shannon.gamma_p31, prob=0.9 )

#remove realdiv*week
m.Fu.Shannon.gamma_p32 <- update(m.Fu.Shannon.gamma_p, 
                                 bf(Hill_q1.Fu ~ realdivLogStd*treatment + treatment*week + 
                                      (1|block/plot)),
                                 seed = SEED) #4 div
summary(m.Fu.Shannon.gamma_p32, prob=0.9 )

#remove both treatment*week and realdiv*week:
m.Fu.Shannon.gamma_p4 <- update(m.Fu.Shannon.gamma_p, 
                                bf(Hill_q1.Fu ~ realdivLogStd*treatment + week + (1|block/plot)),
                                seed = SEED) #19 div
summary(m.Fu.Shannon.gamma_p4, prob=0.9 )
#stepwise simplification keep week, as it is marginally significant (CI=90%)
pp_check(m.Fu.Shannon.gamma_p4, ndraws=100) #looks better than gaussian fit


#remove week:
m.Fu.Shannon.gamma_p5 <- update(m.Fu.Shannon.gamma_p, 
                                bf(Hill_q1.Fu ~ realdivLogStd*treatment + (1|block/plot)),
                                seed = SEED) #2 div
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
     file="./statistics/brms/231221_hill_Fu_realdiv.RData")

rm(m.Fu.Shannon.gamma_p5, #m.Fu.Shannon.gamma_p4, 
   m.Fu.Shannon.gamma_p32,m.Fu.Shannon.gamma_p31,
   m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p,
   m.Fu.Shannon.gaus_p5,m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
   m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p)

####Shannon Pl ~ realdivLogStd: gaus_p (looic: gaus_p5) ####
#all gaussian models show slightly better fit than the gamma models (elpd_diff < 2 SE_diff)
# stepwise elimiantion: gaus_p, as for t2w1-t2w2 : p_dir=96.43% (>90%)
#looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
#differ significantly (less than 2 SE) from the best fit model
load(     file="./statistics/brms/231221_hill_Pl_realdiv.RData")

#using a gaussian distribution with a normal(0,10) prior for beta coefficients
m.Pl.Shannon.gaus_p <- brm(Hill_q1.Pl ~ realdivLogStd*treatment*week + (1|block/plot), 
                           data = dat, family = "gaussian",
                           chains = 3,
                           cores = 3,
                           iter = 2000, warmup = 1000,
                           prior = beta_coeff_priors,
                           seed = SEED,
                           control = list(adapt_delta = 0.99) )  #5 div
summary(m.Pl.Shannon.gaus_p, prob=0.95 ) 
#3-way interaction is significant (CI=95%)
pp_check(m.Pl.Shannon.gaus_p, ndraws = 100)

#check with emmeans:
  emt = emtrends( m.Pl.Shannon.gaus_p, specs = c("treatment", "week"), var="realdivLogStd")
  summary(emt, point.est=mean, level = .9) 
  emt.pairs <- pairs(emt)
  summary(emt.pairs, point.est=mean, level = .9) #t2w1 vs t2w2 doesnt cut 0 at 90%HPD
  bayestestR::p_direction(emt.pairs) #t2 w1 vs t2w2 : p_dir=96.43% (>95%)

#remove 3 way interaction:
m.Pl.Shannon.gaus_p2 <- update(m.Pl.Shannon.gaus_p, 
                               bf(Hill_q1.Pl ~ realdivLogStd*treatment + treatment*week + 
                                    realdivLogStd*week + (1|block/plot)),
                               seed = SEED) # 5 div
summary(m.Pl.Shannon.gaus_p2, prob=0.9)

#remove treatment*week
m.Pl.Shannon.gaus_p31 <- update(m.Pl.Shannon.gaus_p, 
                                bf(Hill_q1.Pl ~ realdivLogStd*treatment + realdivLogStd*week + 
                                     (1|block/plot)),
                                seed = SEED) #3 div
summary(m.Pl.Shannon.gaus_p31, prob=0.9 )

#remove realdiv*week
m.Pl.Shannon.gaus_p32 <- update(m.Pl.Shannon.gaus_p, 
                                bf(Hill_q1.Pl ~ realdivLogStd*treatment + treatment*week + 
                                     (1|block/plot)),
                                seed = SEED) # 4 div
summary(m.Pl.Shannon.gaus_p32, prob=0.9 )

#remove both treatment*week and realdiv*week:
m.Pl.Shannon.gaus_p4 <- update(m.Pl.Shannon.gaus_p, 
                               bf(Hill_q1.Pl ~ realdivLogStd*treatment + week + (1|block/plot)),
                               seed = SEED) #11 div
summary(m.Pl.Shannon.gaus_p4, prob=0.9 )
#stepwise simplification: choose p4, as week is marginally significant (90% CI)

#remove week:
m.Pl.Shannon.gaus_p5 <- update(m.Pl.Shannon.gaus_p, 
                               bf(Hill_q1.Pl ~ realdivLogStd*treatment + (1|block/plot)),
                               seed = SEED)
summary(m.Pl.Shannon.gaus_p5, prob=0.9 )

loo(m.Pl.Shannon.gaus_p5,m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
    m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p)

#using a gamma distribution with a normal(0,10) prior for beta coefficients
m.Pl.Shannon.gamma_p <- brm(Hill_q1.Pl ~ realdivLogStd*treatment*week + (1|block/plot), 
                            data = dat, family = "gamma",
                            chains = 3,
                            cores = 3,
                            iter = 2000, warmup = 1000,
                            prior = beta_coeff_priors,
                            seed = SEED,
                            control = list(adapt_delta = 0.99) )  # 14 div
summary(m.Pl.Shannon.gamma_p, prob=0.9)
#3 way interaction marginally significant, check emtrends:
  emt = emtrends( m.Pl.Shannon.gamma_p, specs = c("treatment", "week"), var="realdivLogStd")
  summary(emt, point.est=mean, level = .9) 
  emt.pairs <- pairs(emt)
  summary(emt.pairs, point.est=mean, level = .9)
  bayestestR::p_direction(emt.pairs) #t2w1 vs t2w2: p_dir = 95.4

#remove 3 way interaction:
m.Pl.Shannon.gamma_p2 <- update(m.Pl.Shannon.gamma_p, 
                                bf(Hill_q1.Pl ~ realdivLogStd*treatment + treatment*week + 
                                     realdivLogStd*week + (1|block/plot)),
                                seed = SEED) # 15 div
summary(m.Pl.Shannon.gamma_p2, prob=0.9)

#remove treatment*week
m.Pl.Shannon.gamma_p31 <- update(m.Pl.Shannon.gamma_p, 
                                 bf(Hill_q1.Pl ~ realdivLogStd*treatment + realdivLogStd*week + 
                                      (1|block/plot)),
                                 seed = SEED) #15 div
summary(m.Pl.Shannon.gamma_p31, prob=0.9 )

#remove realdiv*week
m.Pl.Shannon.gamma_p32 <- update(m.Pl.Shannon.gamma_p, 
                                 bf(Hill_q1.Pl ~ realdivLogStd*treatment + treatment*week + 
                                      (1|block/plot)),
                                 seed = SEED) #4 div
summary(m.Pl.Shannon.gamma_p32, prob=0.9 )

#remove both treatment*week and realdiv*week:
m.Pl.Shannon.gamma_p4 <- update(m.Pl.Shannon.gamma_p, 
                                bf(Hill_q1.Pl ~ realdivLogStd*treatment + week + (1|block/plot)),
                                seed = SEED)
summary(m.Pl.Shannon.gamma_p4, prob=0.9 )

#remove week:
m.Pl.Shannon.gamma_p5 <- update(m.Pl.Shannon.gamma_p, 
                                bf(Hill_q1.Pl ~ realdivLogStd*treatment + (1|block/plot)),
                                seed = SEED)
summary(m.Pl.Shannon.gamma_p5, prob=0.9 )

#compare
Pl.loo <- loo(m.Pl.Shannon.gamma_p5, m.Pl.Shannon.gamma_p4, m.Pl.Shannon.gamma_p32,m.Pl.Shannon.gamma_p31,
              m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p,
              m.Pl.Shannon.gaus_p5,m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
              m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p)
Pl.loo
#using a gaussian distribution shows slightly better fit

save(m.Pl.Shannon.gamma_p5, m.Pl.Shannon.gamma_p4, m.Pl.Shannon.gamma_p32,m.Pl.Shannon.gamma_p31,
     m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p,
     m.Pl.Shannon.gaus_p5,m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
     m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p,
     file="./statistics/brms/231221_hill_Pl_realdiv.RData")

rm(m.Pl.Shannon.gamma_p5, m.Pl.Shannon.gamma_p4, m.Pl.Shannon.gamma_p32,m.Pl.Shannon.gamma_p31,
   m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p,
   m.Pl.Shannon.gaus_p5,m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
   m.Pl.Shannon.gaus_p2#, m.Pl.Shannon.gaus_p
   )

#### Shannon Pr ~ realdivLogStd: gamma_p (looic: gamma_p5) ####
#all gaussian models show significantly better fit than the gamma models (elpd_diff > 2 SE_diff)
# stepwise elimiantion: gamma_p, as prob_dir between t3w1 and t3w2 is 91.33% (>90%)
#looic: gauss_p5, as it is the most parsimonous model whose looic doesnt 
#differ significantly (less than 2 SE) from the best fit model
load(file="./statistics/brms/231221_hill_Pr_realdiv.RData")

#using a gaussian distribution with a normal(0,10) prior for beta coefficients
m.Pr.Shannon.gaus_p <- brm(Hill_q1.Pr ~ realdivLogStd*treatment*week + (1|block/plot), 
                           data = dat, family = "gaussian",
                           chains = 3,
                           cores = 3,
                           iter = 2000, warmup = 1000,
                           prior = beta_coeff_priors,
                           seed = SEED,
                           control = list(adapt_delta = 0.99) )  #9 div
summary(m.Pr.Shannon.gaus_p, prob=0.9 )
  
  #do not do this:
  #almost marginally significant, look at emtrends:
    emt = emtrends( m.Pr.Shannon.gaus_p, specs = c("treatment", "week"), var="realdivLogStd")
    summary(emt, point.est=mean, level = .9) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .9)
    bayestestR::p_direction(emt.pairs) #t1w1-t1w2: 93.57% (>90%) 
  
#remove 3 way interaction:
m.Pr.Shannon.gaus_p2 <- update(m.Pr.Shannon.gaus_p, 
                               bf(Hill_q1.Pr ~ realdivLogStd*treatment + treatment*week + 
                                    realdivLogStd*week + (1|block/plot)),
                               seed = SEED) # 10 div
summary(m.Pr.Shannon.gaus_p2, prob=0.9)

#remove treatment*week
m.Pr.Shannon.gaus_p31 <- update(m.Pr.Shannon.gaus_p, 
                                bf(Hill_q1.Pr ~ realdivLogStd*treatment + realdivLogStd*week + 
                                     (1|block/plot)),
                                seed = SEED) #16 div
summary(m.Pr.Shannon.gaus_p31, prob=0.9 )

#remove realdiv*week
m.Pr.Shannon.gaus_p32 <- update(m.Pr.Shannon.gaus_p, 
                                bf(Hill_q1.Pr ~ realdivLogStd*treatment + treatment*week + 
                                     (1|block/plot)),
                                seed = SEED) # 4 div
summary(m.Pr.Shannon.gaus_p32, prob=0.9 )

#remove both treatment*week and realdiv*week:
m.Pr.Shannon.gaus_p4 <- update(m.Pr.Shannon.gaus_p, 
                               bf(Hill_q1.Pr ~ realdivLogStd*treatment + week + (1|block/plot)),
                               seed = SEED) #9 div
summary(m.Pr.Shannon.gaus_p4, prob=0.9 )
#stepwise simplification: choose p4, as week is marginally significant (95% CI)

    #remove week:
    m.Pr.Shannon.gaus_p5 <- update(m.Pr.Shannon.gaus_p, 
                                   bf(Hill_q1.Pr ~ realdivLogStd*treatment + (1|block/plot)),
                                   seed = SEED)
    summary(m.Pr.Shannon.gaus_p5, prob=0.9 )

loo(m.Pr.Shannon.gaus_p5,m.Pr.Shannon.gaus_p4, m.Pr.Shannon.gaus_p32, m.Pr.Shannon.gaus_p31,
    m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p)

#using a gamma distribution with a normal(0,10) prior for beta coefficients
m.Pr.Shannon.gamma_p <- brm(Hill_q1.Pr ~ realdivLogStd*treatment*week + (1|block/plot), 
                            data = dat, family = "gamma",
                            chains = 3,
                            cores = 3,
                            iter = 2000, warmup = 1000,
                            prior = beta_coeff_priors,
                            seed = SEED,
                            control = list(adapt_delta = 0.99) )  # 8 div
summary(m.Pr.Shannon.gamma_p, prob=0.9) #
#keep p, as 3-way interaction is marginally significant
pp_check(m.Pr.Shannon.gamma_p, ndraws=100)
#look at emtrends:
    emt = emtrends(m.Pr.Shannon.gamma_p, specs = c("treatment", "week"), var="realdivLogStd")
    summary(emt, point.est=mean, level = .9) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .9)
    bayestestR::p_direction(emt.pairs) #t3w1 vs t3w2: p_dir=91.33%

#remove 3 way interaction:
m.Pr.Shannon.gamma_p2 <- update(m.Pr.Shannon.gamma_p, 
                                bf(Hill_q1.Pr ~ realdivLogStd*treatment + treatment*week + 
                                     realdivLogStd*week + (1|block/plot)),
                                seed = SEED) # 15 div
summary(m.Pr.Shannon.gamma_p2, prob=0.9)

#remove treatment*week
m.Pr.Shannon.gamma_p31 <- update(m.Pr.Shannon.gamma_p, 
                                 bf(Hill_q1.Pr ~ realdivLogStd*treatment + realdivLogStd*week + 
                                      (1|block/plot)),
                                 seed = SEED) #9 div
summary(m.Pr.Shannon.gamma_p31, prob=0.9 )

#remove realdiv*week
m.Pr.Shannon.gamma_p32 <- update(m.Pr.Shannon.gamma_p, 
                                 bf(Hill_q1.Pr ~ realdivLogStd*treatment + treatment*week + 
                                      (1|block/plot)),
                                 seed = SEED) #4 div
summary(m.Pr.Shannon.gamma_p32, prob=0.9 )

#remove both treatment*week and realdiv*week:
m.Pr.Shannon.gamma_p4 <- update(m.Pr.Shannon.gamma_p, 
                                bf(Hill_q1.Pr ~ realdivLogStd*treatment + week + (1|block/plot)),
                                seed = SEED) # 4 div
summary(m.Pr.Shannon.gamma_p4, prob=0.9 )

#stepwise simplification keep week, as it is marginally significant (CI=90%)
pp_check(m.Pr.Shannon.gamma_p4, ndraws=100)


#remove week:
m.Pr.Shannon.gamma_p5 <- update(m.Pr.Shannon.gamma_p, 
                                bf(Hill_q1.Pr ~ realdivLogStd*treatment + (1|block/plot)),
                                seed = SEED)
summary(m.Pr.Shannon.gamma_p5, prob=0.9 )

#compare
Pr.loo <- loo(m.Pr.Shannon.gamma_p5, m.Pr.Shannon.gamma_p4, m.Pr.Shannon.gamma_p32,m.Pr.Shannon.gamma_p31,
              m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p,
              m.Pr.Shannon.gaus_p5,m.Pr.Shannon.gaus_p4, m.Pr.Shannon.gaus_p32, m.Pr.Shannon.gaus_p31,
              m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p)
Pr.loo
#using a gamma distribution shows much better fit

save(m.Pr.Shannon.gamma_p5, m.Pr.Shannon.gamma_p4, m.Pr.Shannon.gamma_p32,m.Pr.Shannon.gamma_p31,
     m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p,
     m.Pr.Shannon.gaus_p5,m.Pr.Shannon.gaus_p4, m.Pr.Shannon.gaus_p32, m.Pr.Shannon.gaus_p31,
     m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p,
     file="./statistics/brms/231221_hill_Pr_realdiv.RData")

rm(m.Pr.Shannon.gamma_p5, m.Pr.Shannon.gamma_p4, m.Pr.Shannon.gamma_p32,m.Pr.Shannon.gamma_p31,
   m.Pr.Shannon.gamma_p2, #m.Pr.Shannon.gamma_p,
   m.Pr.Shannon.gaus_p5,m.Pr.Shannon.gaus_p4, m.Pr.Shannon.gaus_p32, m.Pr.Shannon.gaus_p31,
   m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p)

#### Shannon Om ~ realdivLogStd: no fitting possible  ####
ggplot(data = dat, aes(y=Sample, x=Hill_q1.Om))+
  geom_point()

sum(dat$Hill_q1.Om == 1 | dat$Hill_q1.Om == 2) #196 of 240 
  #-> this wont fit neither a gaussian nor a gamma distribution


load( file="./statistics/brms/231221_hill_Om_realdiv.RData")

#using a gaussian distribution with a normal(0,10) prior for beta coefficients
m.Om.Shannon.gaus_p <- brm(Hill_q1.Om ~ realdivLogStd*treatment*week + (1|block/plot), 
                           data = dat, family = "gaussian",
                           chains = 3,
                           cores = 3,
                           iter = 2000, warmup = 1000,
                           prior = beta_coeff_priors,
                           seed = SEED,
                           control = list(adapt_delta = 0.99) )  #6 div
summary(m.Om.Shannon.gaus_p, prob=0.9 )

#remove 3 way interaction:
m.Om.Shannon.gaus_p2 <- update(m.Om.Shannon.gaus_p, 
                               bf(Hill_q1.Om ~ realdivLogStd*treatment + treatment*week + 
                                    realdivLogStd*week + (1|block/plot)),
                               seed = SEED) # 8 div
summary(m.Om.Shannon.gaus_p2, prob=0.9)

#remove treatment*week
m.Om.Shannon.gaus_p31 <- update(m.Om.Shannon.gaus_p, 
                                bf(Hill_q1.Om ~ realdivLogStd*treatment + realdivLogStd*week + 
                                     (1|block/plot)),
                                seed = SEED) #3 div
summary(m.Om.Shannon.gaus_p31, prob=0.9 )

#remove realdiv*week
m.Om.Shannon.gaus_p32 <- update(m.Om.Shannon.gaus_p, 
                                bf(Hill_q1.Om ~ realdivLogStd*treatment + treatment*week + 
                                     (1|block/plot)),
                                seed = SEED) # 4 div
summary(m.Om.Shannon.gaus_p32, prob=0.9 )

#remove both treatment*week and realdiv*week:
m.Om.Shannon.gaus_p4 <- update(m.Om.Shannon.gaus_p, 
                               bf(Hill_q1.Om ~ realdivLogStd*treatment + week + (1|block/plot)),
                               seed = SEED)
summary(m.Om.Shannon.gaus_p4, prob=0.9 )
#stepwise simplification: choose p4, as week is significant (95% CI)

#remove week:
m.Om.Shannon.gaus_p5 <- update(m.Om.Shannon.gaus_p, 
                               bf(Hill_q1.Om ~ realdivLogStd*treatment + (1|block/plot)),
                               seed = SEED)
summary(m.Om.Shannon.gaus_p5, prob=0.9 )
#stepwise simplification: choose p5


loo(m.Om.Shannon.gaus_p5,m.Om.Shannon.gaus_p4, m.Om.Shannon.gaus_p32, m.Om.Shannon.gaus_p31,
    m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p)

#using a gamma distribution with a normal(0,10) prior for beta coefficients
m.Om.Shannon.gamma_p <- brm(Hill_q1.Om ~ realdivLogStd*treatment*week + (1|block/plot), 
                            data = dat, family = "gamma",
                            chains = 3,
                            cores = 3,
                            iter = 2000, warmup = 1000,
                            prior = beta_coeff_priors,
                            seed = SEED,
                            control = list(adapt_delta = 0.99) )  # 4 div
summary(m.Om.Shannon.gamma_p, prob=0.9)
pp_check(m.Om.Shannon.gamma_p, ndraws = 100)

#remove 3 way interaction:
m.Om.Shannon.gamma_p2 <- update(m.Om.Shannon.gamma_p, 
                                bf(Hill_q1.Om ~ realdivLogStd*treatment + treatment*week + 
                                     realdivLogStd*week + (1|block/plot)),
                                seed = SEED) # 5 div
summary(m.Om.Shannon.gamma_p2, prob=0.9)

#remove treatment*week
m.Om.Shannon.gamma_p31 <- update(m.Om.Shannon.gamma_p, 
                                 bf(Hill_q1.Om ~ realdivLogStd*treatment + realdivLogStd*week + 
                                      (1|block/plot)),
                                 seed = SEED) #12 div
summary(m.Om.Shannon.gamma_p31, prob=0.9 )

#remove realdiv*week
m.Om.Shannon.gamma_p32 <- update(m.Om.Shannon.gamma_p, 
                                 bf(Hill_q1.Om ~ realdivLogStd*treatment + treatment*week + 
                                      (1|block/plot)),
                                 seed = SEED) #4 div
summary(m.Om.Shannon.gamma_p32, prob=0.9 )

#remove both treatment*week and realdiv*week:
m.Om.Shannon.gamma_p4 <- update(m.Om.Shannon.gamma_p, 
                                bf(Hill_q1.Om ~ realdivLogStd*treatment + week + (1|block/plot)),
                                seed = SEED)
summary(m.Om.Shannon.gamma_p4, prob=0.9 )
#stepwise simplification keep week, as it is marginally significant (CI=90%)
pp_check(m.Om.Shannon.gamma_p4, ndraws=100)


#remove week:
m.Om.Shannon.gamma_p5 <- update(m.Om.Shannon.gamma_p, 
                                bf(Hill_q1.Om ~ realdivLogStd*treatment + (1|block/plot)),
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
     file="./statistics/brms/231221_hill_Om_realdiv.RData")

rm(m.Om.Shannon.gamma_p5, m.Om.Shannon.gamma_p4, m.Om.Shannon.gamma_p32,m.Om.Shannon.gamma_p31,
   m.Om.Shannon.gamma_p2, m.Om.Shannon.gamma_p,
   m.Om.Shannon.gaus_p5,m.Om.Shannon.gaus_p4, m.Om.Shannon.gaus_p32, m.Om.Shannon.gaus_p31,
   m.Om.Shannon.gaus_p2, m.Om.Shannon.gaus_p)

#### model selection ####
library(brms)
library(ggplot2)

#selection criteria: picking the most parsimonious model which has an elpd_diff > -4 to the model with the best fit
#saving the selected models in a seperate file:

#all ~ realdiv: gaus_p5
load("./statistics/brms/231221_hill_all_realdiv.RData")
  all.loo <- loo(m.all.Shannon.gamma_p5, m.all.Shannon.gamma_p4, m.all.Shannon.gamma_p32,m.all.Shannon.gamma_p31,
                m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p,
                m.all.Shannon.gaus_p5,m.all.Shannon.gaus_p4, m.all.Shannon.gaus_p32, m.all.Shannon.gaus_p31,
                m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p)
  all.loo
  rm(m.all.Shannon.gamma_p5, m.all.Shannon.gamma_p4, m.all.Shannon.gamma_p32,m.all.Shannon.gamma_p31,
     m.all.Shannon.gamma_p2, m.all.Shannon.gamma_p,
     #m.all.Shannon.gaus_p5,
     m.all.Shannon.gaus_p4, m.all.Shannon.gaus_p32, m.all.Shannon.gaus_p31,
     m.all.Shannon.gaus_p2, m.all.Shannon.gaus_p)


#Ba ~ realdiv: gaus_p5
load("./statistics/brms/231221_hill_Ba_realdiv.RData")
Ba.loo <- loo(m.Ba.Shannon.gamma_p5, m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
                        m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
                        m.Ba.Shannon.gaus_p5,m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
                        m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p)
Ba.loo
rm(m.Ba.Shannon.gamma_p5, m.Ba.Shannon.gamma_p4, m.Ba.Shannon.gamma_p32,m.Ba.Shannon.gamma_p31,
    m.Ba.Shannon.gamma_p2, m.Ba.Shannon.gamma_p,
    #m.Ba.Shannon.gaus_p5,
   m.Ba.Shannon.gaus_p4, m.Ba.Shannon.gaus_p32, m.Ba.Shannon.gaus_p31,
    m.Ba.Shannon.gaus_p2, m.Ba.Shannon.gaus_p)

#Fu ~ realdiv: gaus_p5
  load("./statistics/brms/231221_hill_Fu_realdiv.RData")
  Fu.loo <- loo(m.Fu.Shannon.gamma_p5, m.Fu.Shannon.gamma_p4, m.Fu.Shannon.gamma_p32,m.Fu.Shannon.gamma_p31,
                m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p,
                m.Fu.Shannon.gaus_p5,m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
                m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p)
  Fu.loo
  rm(#m.Fu.Shannon.gamma_p5, 
     m.Fu.Shannon.gamma_p4, m.Fu.Shannon.gamma_p32,m.Fu.Shannon.gamma_p31,
     m.Fu.Shannon.gamma_p2, m.Fu.Shannon.gamma_p,
     m.Fu.Shannon.gaus_p5,
     m.Fu.Shannon.gaus_p4, m.Fu.Shannon.gaus_p32, m.Fu.Shannon.gaus_p31,
     m.Fu.Shannon.gaus_p2, m.Fu.Shannon.gaus_p)

#Pl ~ realdiv: gaus_p5
  load("./statistics/brms/231221_hill_Pl_realdiv.RData")
  Pl.loo <- loo(m.Pl.Shannon.gamma_p5, m.Pl.Shannon.gamma_p4, m.Pl.Shannon.gamma_p32,m.Pl.Shannon.gamma_p31,
                m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p,
                m.Pl.Shannon.gaus_p5,m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
                m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p)
  Pl.loo
  rm(m.Pl.Shannon.gamma_p5, 
    m.Pl.Shannon.gamma_p4, m.Pl.Shannon.gamma_p32,m.Pl.Shannon.gamma_p31,
    m.Pl.Shannon.gamma_p2, m.Pl.Shannon.gamma_p,
    #m.Pl.Shannon.gaus_p5,
    m.Pl.Shannon.gaus_p4, m.Pl.Shannon.gaus_p32, m.Pl.Shannon.gaus_p31,
    m.Pl.Shannon.gaus_p2, m.Pl.Shannon.gaus_p)

#Pr ~ realdiv: gamma_p5
  load("./statistics/brms/231221_hill_Pr_realdiv.RData")
  Pr.loo <- loo(m.Pr.Shannon.gamma_p5, m.Pr.Shannon.gamma_p4, m.Pr.Shannon.gamma_p32,m.Pr.Shannon.gamma_p31,
                m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p,
                m.Pr.Shannon.gaus_p5,m.Pr.Shannon.gaus_p4, m.Pr.Shannon.gaus_p32, m.Pr.Shannon.gaus_p31,
                m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p)
  Pr.loo
  rm(#m.Pr.Shannon.gamma_p5, 
    m.Pr.Shannon.gamma_p4, m.Pr.Shannon.gamma_p32,m.Pr.Shannon.gamma_p31,
    m.Pr.Shannon.gamma_p2, m.Pr.Shannon.gamma_p,
    m.Pr.Shannon.gaus_p5,
    m.Pr.Shannon.gaus_p4, m.Pr.Shannon.gaus_p32, m.Pr.Shannon.gaus_p31,
    m.Pr.Shannon.gaus_p2, m.Pr.Shannon.gaus_p)

#Om ~ realdiv: no model can predict accurately enough, as most samples have only one or two identified omnivorous individuals in them
  
#save selected models in .RData file:
  save(m.all.Shannon.gaus_p5, 
       m.Ba.Shannon.gaus_p5,
       m.Fu.Shannon.gamma_p5,
       m.Pl.Shannon.gaus_p5,
       m.Pr.Shannon.gamma_p5,
       file = "./statistics/brms/240116_Hill_realdiv_mselect.RData")
  
  
  
  
  
  
  

