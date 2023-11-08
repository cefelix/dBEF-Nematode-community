####loading data and packages ####
# a workflow: https://m-clark.github.io/posts/2021-02-28-practical-bayes-part-i/
# rethinking in brms: https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/interactions.html
library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
dBEF_nem %>%
  str()
load(file = "./statistics/brms/231004_Lm.RData")

dBEF_nem$treatment <- dBEF_nem$treatment %>%
  as.factor()

#### 2.1 combined analysis for data from 2017 and 2021 ####

#examining the individual abundances of 2017 and 2021 in one model, with year as a continuos variable
#this should not be used due to different sampling depths in 2017 (25g FW, 5cm) and 2021 (25g, 10cm)

#distribution: lognormal
par(mfrow=c(1,2))
hist(dBEF_nem$ind_per100g, breaks = seq(min(dBEF_nem$ind_per100g), max(dBEF_nem$ind_per100g), length.out=30))
rlnorm(1000, meanlog = log(10), sdlog = log(2.25)) %>%
  hist()
par(mfrow=c(1,1))

#so lets add a column with log-transformed individual abundances:
dBEF_nem <- dBEF_nem %>%
  mutate(Log_ind_per100g = log(ind_per100g), .after = ind_per100g)
hist(dBEF_nem$Log_ind_per100g)

#lets see the data
ggplot(dBEF_nem, aes(x= log(sowndiv), y=ind_per100g, col=treatment))+
  geom_point()

#m1: abun ~ sowndiv+(1|block) fam=lognormal
#m2: log(abun) ~ sowndiv+(1|block), fam=gaussian
#m3: log(abun) ~ sowndiv+treatment+(1|block), fam=gaussian
#m4: log(abun) ~ sowndiv+SH+PH+(1|block), fam= gaussian, SH+PH numeric
#m5: log(abun) ~ sowndiv+SH+(1|block), fam= gaussian, SH factor
#m6: log(abun) ~ sowndiv*SH+(1|block), fam= gaussian, SH numeric

m1.1 <- brm(ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "lognormal",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000)
summary(m1.1) #5 divergent transitions

m1.2 <- brm(ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "lognormal",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta=0.99))
summary(m1.2) #still 1 divergent transition

m1.3 <- brm(ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "lognormal",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta=0.9999))
summary(m1.3) # no divergent transitions
#
pp_check(m1.3)


m2.1 <- brm(Log_ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000)
summary(m2.1) # 3 divergent transitions

m2.2 <- brm(Log_ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.99))
summary(m2.2) # 1 divergent transition

m2.3 <- brm(Log_ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.9999))

summary(m2.3) #184 transitions after warmup exceeding max. treedepth h -> increase h above 10

m2.4 <- brm(Log_ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.9999,
                           max_treedepth = 11))
summary(m2.4)

pp_check(m2.4)

m3.1 <- brm(Log_ind_per100g ~ sowndiv + treatment + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000)
summary(m3.1) #6 divergent transitions

m3.2 <- brm(Log_ind_per100g ~ sowndiv + treatment + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.99))
summary(m3.2)

pp_check(m3.2)


#include SH and PH as numeric variables
m4.1 <- brm(Log_ind_per100g ~ sowndiv + SH + PH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000)
summary(m4.1)
pp_check(m4.1)

m4.2 <- brm(Log_ind_per100g ~ sowndiv + SH + PH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.99))
summary(m4.2)

m4.3 <- brm(Log_ind_per100g ~ sowndiv + SH + PH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.9999))
summary(m4.3)

pp_check(m4.2)



#include SH and PH as a factors
dBEF_nem$SH <- dBEF_nem$SH %>% as.factor()
dBEF_nem$PH <- dBEF_nem$PH %>% as.factor()

m5.1 <- brm(Log_ind_per100g ~ sowndiv + SH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 10))
summary(m5.1)

#but actually, we want to see whether the effect of sowndiv depends on SH
#thus: interaction
dBEF_nem$SH <- as.numeric(dBEF_nem$SH)

m6.1 <- brm(Log_ind_per100g ~ sowndiv * SH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000)
summary(m6.1) # 4 divergent transitions

m6.2 <- brm(Log_ind_per100g ~ sowndiv * SH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.99))
summary(m6.2) #1 divergent transition

m6.3 <- brm(Log_ind_per100g ~ sowndiv * SH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.9999))
summary(m6.3)

pp_check(m6.3)

#### 2.2 saving brm outputs ####
save(m1.1, m1.2,m1.3, 
     m2.1, m2.2, m2.3, m2.4,
     m3.1, m3.2,
     m4.1, m4.2, m4.3,
     m5.1,
     m6.1, m6.2, m6.3,
     file = "./statistics/brms/231004_Lm.RData")

#### 3.1 analysing 2017 data separately ####
dBEF_nem17 <- subset(dBEF_nem, year == 2017 )

#distribution: lognormal
par(mfrow=c(1,2))
hist(dBEF_nem17$ind_per100g, breaks = seq(min(dBEF_nem17$ind_per100g), max(dBEF_nem17$ind_per100g), length.out=30))
rlnorm(1000, meanlog = log(10), sdlog = log(1.5)) %>%
  hist()
par(mfrow=c(1,1))

#plot the data:
ggplot(dBEF_nem17, aes(x=log(sowndiv), y=ind_per100g, col=treatment) )+
  geom_jitter()
  #treatment 3 (+SH, +PH) show positiv relationship
  #treatment 2 (+Sh, -PH) shows positive trend over 1, 4, 16 species, but very low abundances at 60 sp
  #treatment 1 (-SH, -PH) same as treatment 2

m17.11





#### 4.1 analysing 2021 data separately ####
dBEF_nem21 <- subset(dBEF_nem, year == 2021)
dBEF_nem21$SH <- as.factor(dBEF_nem21$SH)
dBEF_nem21$PH <- as.factor(dBEF_nem21$PH)

#distribution: lognormal
par(mfrow=c(1,2))
hist(dBEF_nem21$ind_per100g, breaks = seq(min(dBEF_nem21$ind_per100g), max(dBEF_nem21$ind_per100g), length.out=30))
rlnorm(1000, meanlog = log(10), sdlog = log(2)) %>%
  hist()
par(mfrow=c(1,1))

#plot the data:
ggplot(dBEF_nem21, aes(x=log(sowndiv), y=ind_per100g, col=treatment) )+
  geom_jitter()

#model terms:
#m21.1: abun ~ sowndiv*treatment + (1|block), fam=lognormal
#m21.2: abun ~ sowndiv*SH*PH + (1|block), fam=lognormal 
  #m1 and m2 should be the same, as treatment is an interaction of SH and PH (?)
# { m21.3: log(abun) ~ sowndiv*treatment + (1|block), fam=gaussian
#   m21.4: abun ~ sowndiv*SH + (1|block), fam=lognormal
#   m21.5: abun ~ sowndiv*PH + (1|block), fam=lognormal } not done yet

#loading the fit models from previous session:
load("./statistics/brms/231017_abun.RData")

#m21.1
  
  #m21.10 arbitrarily leaves out the random factor block
  m21.10 <- brm(ind_per100g ~ sowndiv * treatment, 
                data = dBEF_nem21, family = "lognormal",
                chains = 3,
                cores = 3,
                iter = 2000, warmup = 1000, 
                control = list(adapt_delta =0.99))
  summary(m21.10)
  

  m21.11 <- brm(ind_per100g ~ sowndiv * treatment + (1|block), 
              data = dBEF_nem21, family = "lognormal",
              chains = 3,
              cores = 3,
              iter = 2000, warmup = 1000)
  summary(m21.11) #8 divergent transitions
  
  m21.12 <- update(m21.11, control = list(adapt_delta=0.99))
    summary(m21.12)
  mcmc_plot(m21.12, type = "combo")  
  mcmc_plot(m21.12, type = "violin")
  mcmc_plot(m21.12, type = "rhat")
  
  #posterior predictive check:
  pp_check(m21.12, ndraws = 100)
    #check if model predicts minimum/median/max values well:
    pp_check(m21.12, ndraws=100, 
             type = "stat", stat="min") #min is captured well
    pp_check(m21.12, ndraws=100, 
             type = "stat", stat="median") #model underestimates median
    pp_check(m21.12, ndraws=100, 
             type = "stat", stat="max") #model overestimates max
    
  #bayesian R2 
  bayes_R2(m21.12) #with random effects
  bayes_R2(m21.12, re_formula = NA) #without random effects
  
  
  
#m21.2
  m21.21 <- brm(ind_per100g ~ sowndiv * SH * PH + (1|block), 
                 data = dBEF_nem21, family = "lognormal",
                 chains = 3,
                 cores = 3,
                 iter = 2000, warmup = 1000)
  summary(m21.21) #apparantly this threefold interaction is a bad idea
  
  m21.22 <- update(m21.21, control = list(adapt_delta=0.99))
    summary(m21.22)
  
  m21.23 <- update(m21.21, control = list(adapt_delta=0.99))
    summary(m21.23)
    
  pp_check(m21.23, ndraws=100)
  
  
  
####4.2 saving 2021s models####
  save(m21.10, m21.11,m21.12, 
       m21.21, m21.22, m21.23,
       file = "./statistics/brms/231017_abun.RData")
  
  
####5 modelling only counts####
  dBEF_nem$total_nematodes %>% density() %>% plot()
m.counts.11 <- brm(total_nematodes ~ sowndiv*treatment + (1|block), 
                   data = dBEF_nem21, family = "poisson",
                   chains = 3,
                   cores = 3,
                   iter = 2000, warmup = 1000 )

m.counts.13 <- update(m.counts.11,
                      control = list(adapt_delta = 0.95,
                                     max_treedepth=12))
pp_check(m.counts.11)
pp_check(m.counts.13, ndraws=100) #bimodal, not fitting at all
#same for 17 data
m.counts.11a <- brm(total_nematodes ~ sowndiv*treatment + offset(log(soilDW)) + (1|block), 
                   data = dBEF_nem17, family = "poisson",
                   chains = 3,
                   cores = 3,
                   iter = 3000, warmup = 1500,
                   control = list(adapt_delta = 0.9))
####somewhat working model with count data ####
m.counts.11c <- update(m.counts.11b, 
                       iter=4000, warmup=1500,
                       control = list(adapt_delta = 0.99,
                                      max_treedepth = 15))
pp_check(m.counts.11c) #3 divergent transitions 
  #underpredicting low counts, too high probability mass in the centre

m.counts.11d <- update(m.counts.11c,
                       control = list(adapt_delta=0.995))
pairs(m.counts.11d) %>% plot()

save(m.counts.11,
     m.counts.11a, m.counts.11b, m.counts.11c, m.counts.11d,
     m.counts.12, 
     m.counts.13,
     file = "./statistics/brms/231102_TotAbun_counts.RData")
#exclude 

dBEF_nem21$total_nematodes %>% density() %>% plot()
dBEF_nem21$ind_per100g %>% density() %>% plot()


#### modelling abundance as a count with soilDW as an offset variable ####

#https://en.wikipedia.org/wiki/Poisson_regression#:~:text=In%20these%20examples%2C%20exposure%20is,right%20side%20of%20the%20equation.
#https://www.dataquest.io/blog/tutorial-poisson-regression-in-r/
#https://github.com/paul-buerkner/brms/issues/98
  #offset() is usable as in other packages

#### Poisson 21: total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), fam=poisson ####
#as the link function for poisson is log(), lets log transform offset variable:
dBEF_nem21$soilDW %>% log() %>% density() %>% plot() #bimodal
  #lets check plotwise:
  ggplot(data = dBEF_nem21, aes(x = sowndivLog, y = soilDW))+
    geom_point()+
    facet_wrap(~block)
  
#the response variables density
dBEF_nem21$total_nematodes  %>% density() %>% plot()
rpois(1e3, lambda = 2) %>% density() %>% plot()
SEED = 22061996

m.abun.PoisOffS.21 <- brm(total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), 
                          data = dBEF_nem21, family = "poisson",
                          chains = 3,
                          cores = 3,
                          iter = 3000, warmup = 1500,
                          seed = SEED,
                          control = list(adapt_delta = 0.9) ) 
                          #74 divergent transitions,
                          #291 exceeded max_treedepth

m.abun.PoisOffS.22 <- update(m.abun.PoisOffS.21,
                             seed = SEED,
                             control = list(adapt_delta = 0.999,
                                            max_treedepth = 12)) 
                              #427 transitions exceeded max_treedepth
pp_check(m.abun.PoisOffS.22, ndraws=100) #thats at least a better fit than m.[...].14


m.abun.PoisOffS.23 <- update(m.abun.PoisOffS.22,
                             seed = SEED,
                             control = list(adapt_delta = 0.999,
                                            max_treedepth = 15))
                              #EES too low

m.abun.PoisOffS.24 <- update(m.abun.PoisOffS.23,
                             seed = SEED,
                             iter=4000, warmup =1500,
                             control = list(adapt_delta = 0.999,
                                            max_treedepth = 15))

pp_check(m.abun.PoisOffS.24, ndraws = 100)
#underpredicting until mean

save(m.abun.PoisOffS.11, m.abun.PoisOffS.21,
     m.abun.PoisOffS.12, m.abun.PoisOffS.22,
     m.abun.PoisOffS.13, m.abun.PoisOffS.23,
     m.abun.PoisOffS.14, m.abun.PoisOffS.24,
     file = "./statistics/brms/231106_abundance_OffSet.RData")


#### negBinom 11: total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), fam=negbinomial ####
m.abun.nBinomOffS.11 <- brm(total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), 
                          data = dBEF_nem21, family = "negbinomial",
                          chains = 3,
                          cores = 3,
                          iter = 3000, warmup = 1500,
                          seed = SEED,
                          control = list(adapt_delta = 0.9) ) 
                          #9 divergent transitions

m.abun.nBinomOffS.12 <- update(m.abun.nBinomOffS.11,
                             seed = SEED,
                             control = list(adapt_delta = 0.999,
                                            max_treedepth = 12)) 

m.abun.nBinomOffS.13 <- update(m.abun.nBinomOffS.12,
                               seed = SEED,
                               control = list(adapt_delta = 0.999,
                                              max_treedepth = 15)) 



summary(m.abun.nBinomOffS.13)
pp_check(m.abun.nBinomOffS.13, ndraws = 100)

#### negBinom 21: total_nematodes ~ sowndivLog*treatment + offset(log(soilDW)) + (1|block/plot), fam=negbinomial ####
SEED = 22061996
m.abun.nBinomOffS.21 <- brm(total_nematodes ~ sowndivLog*treatment + offset(log(soilDW)) + (1|block/plot), 
                            data = dBEF_nem21, family = "negbinomial",
                            chains = 3,
                            cores = 3,
                            iter = 3000, warmup = 1500,
                            seed = SEED,
                            control = list(adapt_delta = 0.9) ) 
  #14 divergent transitions

m.abun.nBinomOffS.22 <- update(m.abun.nBinomOffS.21,
                               seed = SEED,
                               control = list(adapt_delta = 0.999,
                                              max_treedepth = 12)) 

pp_check(m.abun.nBinomOffS.22, ndraws=100)
summary(m.abun.nBinomOffS.22)

#### save offset models ####
save(m.abun.PoisOffS.21, m.abun.PoisOffS.22, 
     m.abun.PoisOffS.23, m.abun.PoisOffS.24,
     m.abun.nBinomOffS.11, m.abun.nBinomOffS.12, 
     m.abun.nBinomOffS.13,
     m.abun.nBinomOffS.21, m.abun.nBinomOffS.22,
     file="./statistics/brms/231108_abundance_OffSet.RData")

load(file="./statistics/brms/231108_abundance_OffSet.RData")


####plot the offset models####
conditional_effects(m.abun.nBinomOffS.13)

conditional_effects(m.abun.nBinomOffS.22, prob=0.85)

