####loading data and packages ####
# a workflow: https://m-clark.github.io/posts/2021-02-28-practical-bayes-part-i/
# rethinking in brms: https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/interactions.html
library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)

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
  xlim(-2.5,1.5)

  exp(1.5) 
#this would mean that increasing plant diversity by 1 SD would lead to 348% more individuals 


#### Poisson 21: total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), fam=poisson ####

#the response variables density
dBEF_nem21$total_nematodes  %>% density() %>% plot()
rpois(1e3, lambda = 2) %>% density() %>% plot()
SEED = 22061996

abun.Pois.21 <- brm(total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), 
                          data = dat, family = "poisson",
                          chains = 3,
                          cores = 3,
                          iter = 3000, warmup = 1500,
                          seed = SEED,
                          control = list(adapt_delta = 0.9) ) 
                          #74 divergent transitions,
                          #291 exceeded max_treedepth
pp_check(abun.Pois.21, ndraws=100)

abun.Pois.21b <- update(abun.Pois.21,
                             seed = SEED,
                             control = list(adapt_delta = 0.999,
                                            max_treedepth = 12)) 
                              #427 transitions exceeded max_treedepth
pp_check(abun.Pois.21b, ndraws=100) 

# the poisson distribtution doesnt really fit the data (overdispersion),
# lets use the more general negative binomial distribution instead



#### nBinom 21: total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), fam=negbinomial, -60sp ####
m.abun.nBinomOffS.31 <- brm(total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), 
                          data = dat, family = "negbinomial",
                          chains = 3,
                          cores = 3,
                          iter = 3000, warmup = 1500,
                          seed = SEED,
                          control = list(adapt_delta = 0.99) ) 
                          #all good



pp_check(m.abun.nBinomOffS.31, ndraws = 100)

#### negBinom 31 vog: total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), fam=negbinomial, -60sp, year =2017 ####
SEED = 22061996
dat = dBEF_nem17 %>% 
  filter(sowndiv != 60)
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog)


m.abun.nBinomOffS.31vog <- brm(total_nematodes ~ sowndivLogStd*treatment + offset(log(soilDW)) + (1|block/plot), 
                            data = dat, family = "negbinomial",
                            chains = 3,
                            cores = 3,
                            iter = 3000, warmup = 1500,
                            seed = SEED,
                            control = list(adapt_delta = 0.99) ) 
#all good

pp_check(m.abun.nBinomOffS.31vog, ndraws = 100)

#save both nbinom offset models
save(m.abun.nBinomOffS.31,
     m.abun.nBinomOffS.31vog,
     file = "./statistics/brms/231127_abundance_OffSet.RData")

#### save offset models ####
save(m.abun.PoisOffS.21, m.abun.PoisOffS.22, 
     m.abun.PoisOffS.23, m.abun.PoisOffS.24,
     m.abun.nBinomOffS.11, m.abun.nBinomOffS.12, 
     m.abun.nBinomOffS.13,
     m.abun.nBinomOffS.21, m.abun.nBinomOffS.22,
     file="./statistics/brms/231108_abundance_OffSet.RData")

load(file="./statistics/brms/231108_abundance_OffSet.RData")



