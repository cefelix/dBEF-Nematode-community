####copy####
library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(hexbin)
library(GGally)
library(emmeans)

#week 1: sampled B1 and B2, dry
#week 2: sampled B3 and B4, wet

#exclude 60 sp.:
dat <- subset(dBEF_nem21, sowndiv != 60) 
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog) %>%
  mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
         .after = realdivLog)






####Ba ~ realdiv, stepwise simplification ####

SEED = 22061996  

#most complex, realdiv: 
m.Ba.3way_a <- brm(
  bf(Ba_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + treatment:week + 
       realdivLogStd:treatment:week + (1|block/plot),
     hu ~ week + realdivLogStd + week:realdivLogStd + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #3 divergent transitions


m.Ba.3way_a <- update(m.Ba.3way_a ,
                      control=list(adapt_delta=0.999, 
                                   max_treedepth=12))  #all good
summary(m.Ba.3way_a, prob=0.9 )
emt = emtrends(m.Ba.3way_a, specs = c("treatment", "week"), var="realdivLogStd")
summary(emt, point.est=mean, level = .95) 

#remove 3 way interaction:
m.Ba.2way_a <- update(m.Ba.3way_a, .~. -realdivLogStd:treatment:week,
                      control = list(adapt_delta=0.99,
                                     max_treedepth=10)) #10 div trans
summary(m.Ba.2way_a, prob=0.9)

#remove week:treatment
m.Ba.2way_b <- update(m.Ba.2way_a, .~. -treatment:week,
                      control = list(adapt_delta=0.99,
                                     max_treedepth=10))
summary(m.Ba.2way_b, prob=0.9) 
pp_check(m.Ba.2way_b, ndraws = 100)

#remove hu~week:realdiv  
m.Ba.2way_b2 <- brm(
  bf(Ba_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + (1|block/plot),
     hu ~ week + realdivLogStd + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

summary(m.Ba.2way_b2, prob=0.9) #best
pp_check(m.Ba.2way_b2, ndraws = 100)

#remove hu~week
m.Ba.2way_b3 <- brm(
  bf(Ba_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + (1|block/plot),
     hu ~ realdivLogStd + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #6 diverg transitions

summary(m.Ba.2way_b3, prob=0.9)

#remove hu~realdivLogStd
m.Ba.2way_b4 <- brm(
  bf(Ba_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + (1|block/plot),
     hu ~ week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

summary(m.Ba.2way_b4, prob=0.9)

#remove hu ~ week )
m.Ba.2way_b5 <- brm(
  bf(Ba_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #2 div transitions

####best fit Ba~realdiv model####
m.Ba.2way_b6 <- update(m.Ba.2way_b5, 
                       control = list(adapt_delta=0.999))
summary(m.Ba.2way_b6, prob=0.9)


#remove both week:treatment and week:realdivLogStd  
m.Ba.2way_c <- update(m.Ba.2way_b, .~. -week:realdivLogStd,
                      control = list(adapt_delta=0.99,
                                     max_treedepth=10)) #4 div trans
summary(m.Ba.2way_c)

#remove week  
m.Ba.2way_d <- update(m.Ba.2way_c, .~. -week,
                      control = list(adapt_delta=0.99,
                                     max_treedepth=10)) #all good
summary(m.Ba.2way_d)

#remove hu~week:realdivLogStd
m.Ba.2way_d2 <- brm(
  bf(Ba_per100g ~ realdivLogStd + treatment + 
       realdivLogStd:treatment + (1|block/plot),
     hu ~ week + realdivLogStd + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

summary(m.Ba.2way_d2)

#remove hu~week
m.Ba.2way_d3 <- brm(
  bf(Ba_per100g ~ realdivLogStd + treatment + 
       realdivLogStd:treatment + (1|block/plot),
     hu ~  realdivLogStd + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
summary(m.Ba.2way_d3, prob=0.9)
pp_check(m.Ba.2way_d3, ndraws=100)

#remove hu~realdivLogStd
m.Ba.2way_d4 <- brm(
  bf(Ba_per100g ~ realdivLogStd + treatment + 
       realdivLogStd:treatment + (1|block/plot),
     hu ~  week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
summary(m.Ba.2way_d4)

#remove hu~week
m.Ba.2way_e <- brm(
  bf(Ba_per100g ~ realdivLogStd + treatment + 
       realdivLogStd:treatment + (1|block/plot),
     hu ~  1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #1 div

summary(m.Ba.2way_e)
m.Ba.2way_e <- update(m.Ba.2way_e, control = list(adapt_delta=0.999))


#get slopes:
emt = emtrends(m.Ba.3way_a, specs = c("treatment", "week"), var="realdivLogStd")
summary(emt, point.est=mean, level = .95) 
emt.pairs <- pairs(emt)
summary(emt.pairs, point.est=mean, level = .95)
bayestestR::p_direction(emt.pairs) #probability of direction
#t1w1 - t1w2 78.60%
#t2w1 - t2w2 70.33%
#t3w1 - t3w2 88.27%

#none "significant" -> exclude 3 way interaction:


#### Fu ~ realdiv, hu~1, stepwise simplification####

#hu~1, as there ar only 4 samples with zero fungivores, which is too little to assess prob(Fu==0)
#!running only one chain to speed up fitting!

#3way interaction  
m.Fu.3way_a <- brm(
  bf(Fu_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + treatment:week +
       realdivLogStd:treatment:week + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 1,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

summary(m2.Fu.3way_a, prob=0.9)


#remove 3 way interaction:
m.Fu.2way_a <- update(m.Fu.3way_a, .~. -realdivLogStd:treatment:week) #5 divergs
summary(m2.Fu.2way_a, prob=0.9) #no effect of realdivLogStd:week


####best fit Fu~realdiv model####
#remove realdivLogStd:week
m.Fu.2way_b <- update(m.Fu.2way_a, .~. -realdivLogStd:week) #all good
summary(m.Fu.2way_b, prob=0.9)
pp_check(m.Fu.2way_b, ndraws=100)

#add hu~week*treatment 
m.Fu.2way_b2 <- brm(
  bf(Fu_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + treatment:week + (1|block/plot),
     hu ~ week + treatment + week:treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 1,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #8 divs, 45 exceed max_treedepth 

#remove hu~week:treatment 
m.Fu.2way_b3 <- brm(
  bf(Fu_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + treatment:week + (1|block/plot),
     hu ~ week + treatment + (1|block/plot) ),
  data = dat, 
  family = hurdle_lognormal,
  chains = 1,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #3 divs
summary(m2.Fu.2way_b3)

#remove hu~week 
m.Fu.2way_b4 <- brm(
  bf(Fu_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + treatment:week + (1|block/plot),
     hu ~ week + (1|block/plot) ),
  data = dat, 
  family = hurdle_lognormal,
  chains = 1,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #bulk ESS too low
summary(m2.Fu.2way_b4)

#remove treatment:week
m.Fu.2way_c <- update(m.Fu.2way_a, .~. -treatment:week) #125 divergent
summary(m.Fu.2way_c, prob=0.9)

#remove both   treatment:week AND realdivLogStd:week
m.Fu.2way_d <- update(m.Fu.2way_b, .~. -treatment:week) #3 diverg
summary(m.Fu.2way_d, prob=0.9)



summary(m2.Fu.3way_a, prob=0.9)

#### Pl ~ realdiv, stepwise simplification ####    
sum(dat$Pl_per100g == 0) #2 -> hurdle ~ 1

#1 chain for speed
m.Pl.3way_a <- brm(
  bf(Pl_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + treatment:week + 
       realdivLogStd:treatment:week + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 1,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #2 div trans
summary(m.Pl.3way_a, prob=.9)

#remove 3-way interaction:
m.Pl.2way_a <- update(m.Pl.3way_a, .~. -realdivLogStd:treatment:week ) #1 div trans
summary(m.Pl.2way_a , prob=0.9)

#remove treatment:week
m.Pl.2way_b <- update(m.Pl.2way_a, .~. -treatment:week, 
                      chains=3) #2 divs
summary(m.Pl.2way_b , prob=0.9)

#remove realdivLogStd:week
m.Pl.2way_b2 <- update(m.Pl.2way_a, .~. -realdivLogStd:week) #2 divergent, worse than b
summary(m.Pl.2way_b2 , prob=0.9)

#remove both treatment:week AND realdivLogStd:week
m.Pl.2way_c <- update(m.Pl.2way_b, .~. -realdivLogStd:week,
                      iter = 3000, warmup = 1500,) #9 divergents

m.Pl.2way_c2 <- update( m.Pl.2way_c, control=list(adapt_delta=0.999)) #11
summary(m.Pl.2way_c , prob=0.9)  

####best fit Pl~realdiv model####
#remove week
m.Pl.2way_d <- update(m.Pl.2way_c, .~. -week,
                      chains=3) #all good
summary(m.Pl.2way_d )
pp_check(m.Pl.2way_d, ndraws=100)

#### Pr ~ realdiv, stepwise simplification ####    
sum(dat$Pr_per100g == 0) #57 -> hurdle ~ term

b <- prior(normal(0,20), class = "b")

#1 chain for speed
m.Pr.3way_a <- brm(
  bf(Pr_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + treatment:week + 
       realdivLogStd:treatment:week + (1|block/plot),
     hu ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + treatment:week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  #prior = b,
  chains = 1,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #2 div transitions
summary(m.Pr.3way_a, prob=.9)  

pp_check(m.Pr.3way_a, ndraws=100 )

#remove 3 way interaction
m.Pr.2way_a <- update(m.Pr.3way_a, .~. -realdivLogStd:treatment:week,
                      chains=3) #1 div trans, tail ESS too low

summary(m.Pr.2way_a, prob=0.9)  

#remove hu~treatment:week
m.Pr.2way_a2 <- brm(
  bf(Pr_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + treatment:week + (1|block/plot),
     hu ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #8 div, tail ESS too low

summary(m.Pr.2way_a2, prob=0.9)  

#remove treatment:week
m.Pr.2way_b <- update(m.Pr.2way_a2, .~. -treatment:week) #17 div
summary(m.Pr.2way_b, prob=0.9)

#remove hu~realdivLogStd:treatment
m.Pr.2way_b2 <- brm(
  bf(Pr_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week +  (1|block/plot),
     hu ~ realdivLogStd + treatment + week +  realdivLogStd:week + 
       (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #1 div, tail ESS too low  
summary(m.Pr.2way_b2, prob=0.9)

#remove hu~realdivLogStd:week
m.Pr.2way_b3 <- brm(
  bf(Pr_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week +  (1|block/plot),
     hu ~ realdivLogStd + treatment + week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #6 div transitions
summary(m.Pr.2way_b3, prob=0.9)

#remove hu~realdivLogStd
m.Pr.2way_b4 <- brm(
  bf(Pr_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week +  (1|block/plot),
     hu ~ treatment + week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #1 div
summary(m.Pr.2way_b4, prob=0.9)

#remove hu~week
m.Pr.2way_b5 <- brm(
  bf(Pr_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week +  (1|block/plot),
     hu ~ treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #35, Tail ESS too low
summary(m.Pr.2way_b5, prob=0.9)

#remove hu~treatment
m.Pr.2way_b6 <- brm(
  bf(Pr_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week +  (1|block/plot),
     hu ~ week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #13 div
summary(m.Pr.2way_b6, prob=0.9)

#remove treatment:week from m.Pr.2way_b4
m.Pr.2way_c4 <- update(m.Pr.2way_b4, .~. - treatment:week) #13 div, tail ESS too low
summary(m.Pr.2way_c4, prob=0.9)

#remove week
m.Pr.2way_c42 <- update(m.Pr.2way_c4, .~. - week) 
summary(m.Pr.2way_c42, prob=0.9)



#remove treatment:week from m.Pr.2way_b5
m.Pr.2way_c5 <- update(m.Pr.2way_b5, .~. - treatment:week) #6 div
summary(m.Pr.2way_c5, prob=0.9)


#remove treatment:week from m.Pr.2way_b6
m.Pr.2way_c6 <- update(m.Pr.2way_b6, .~. - treatment:week) #233
summary(m.Pr.2way_c6, prob=0.9)

#remove week:realdivLogStd
m.Pr.2way_d <- update(m.Pr.2way_c6, .~. - realdivLogStd:week) #27 div
summary(m.Pr.2way_d, prob=0.9)

####best Pr~realdiv model####
#remove week 
m.Pr.2way_d2 <- update(m.Pr.2way_d, .~. -week, 
                       control=list(adapt_delta=0.999)) #0 div trans
summary(m.Pr.2way_d2, prob=0.9)


#remove hu~week
m.Pr.2way_e <- brm(
  bf(Pr_per100g ~ realdivLogStd + treatment  + 
       realdivLogStd:treatment +  (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #1 div transition

m.Pr.2way_e2 <-update(m.Pr.2way_e, 
                      control=list(adapt_delta=0.999))

summary(m.Pr.2way_e, prob=0.9)
loo(m.Pr.2way_e, m.Pr.2way_d2)

####Om~realdiv, stepwise simplification####  
sum(dat$Om_per100g ==0) #131
#1 chain to be faster

m.Om.3way_a <- brm(
  bf(Om_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + treatment:week + 
       realdivLogStd:treatment:week + (1|block/plot),
     hu ~ realdivLogStd + treatment + week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 1,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #2 div transitions
summary(m.Om.3way_a, prob=.9) 

#remove 3way interaction
m.Om.2way_a <- update(m.Om.3way_a, .~. -realdivLogStd:treatment:week) #3 div
summary(m.Om.2way_a, prob=.9) 

#remove hu~realdivLogStd
m.Om.2way_a2 <- brm(
  bf(Om_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + treatment:week + (1|block/plot),
     hu ~ treatment + week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 1,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #4 div transitions
summary(m.Om.2way_a2, prob=.9) 

#remove treatment:week
m.Om.2way_b <- update(m.Om.2way_a2, .~. -treatment:week) #7div
summary(m.Om.2way_b, prob=.9) 

#remove realdivLogStd:week
m.Om.2way_c <- update(m.Om.2way_b, .~. -realdivLogStd:week) #all good
summary(m.Om.2way_c, prob=.9) 


####Om best fit ####
#remove week
m.Om.2way_d <- update(m.Om.2way_c, .~. -week) #all good
summary(m.Om.2way_d, prob=.9) 
pp_check(m.Om.2way_d, ndraw=100)+
  xlim(0,180)


####save models####

#best realdiv models
save(m.Ba.2way_b6,
     
     file="./statistics/brms/231205_densBEST_realdiv.RData")

#best realdiv models
save(file="./statistics/brms/231205_densBEST_realdiv.RData")

#all Ba~realdiv models        
save(m.Ba.2way_a, m.Ba.2way_b,
     m.Ba.2way_b2, m.Ba.2way_b3, m.Ba.2way_b4, m.Ba.2way_b5, m.Ba.2way_b6,
     m.Ba.2way_c,
     m.Ba.2way_d, m.Ba.2way_d2, m.Ba.2way_d3,m.Ba.2way_d4,
     m.Ba.2way_e,
     m.Ba.3way_a,
     file="./statistics/brms/231205_densBa_realdiv.RData")

#all Fu~realdiv models   
save(m.Fu.2way_a, m.Fu.2way_b,
     m.Fu.2way_b2, m.Fu.2way_b3, m.Fu.2way_b4,
     m.Fu.2way_c, m.Fu.2way_d,
     m.Fu.3way_a,
     file="./statistics/brms/231205_densFu_realdiv.RData")

#all Pl~realdiv models
save(m.Pl.2way_a, m.Pl.2way_b,
     m.Pl.2way_b2,
     m.Pl.2way_c, m.Pl.2way_d,
     m.Pl.3way_a,
     file="./statistics/brms/231205_densPl_realdiv.RData")

#all Pr~realdiv models
save(m.Pr.2way_a, m.Pr.2way_a2,
     m.Pr.2way_b, m.Pr.2way_b2, m.Pr.2way_b3, m.Pr.2way_b4, m.Pr.2way_b5, m.Pr.2way_b6,
     m.Pr.2way_c4, m.Pr.2way_c42, m.Pr.2way_c5,
     m.Pr.2way_c6, 
     m.Pr.2way_d, m.Pr.2way_d2, m.Pr.2way_e, m.Pr.2way_e2,
     m.Pr.3way_a,
     
     file="./statistics/brms/231205_densPr_realdiv.RData")

#all Om~realdiv
save(m.Om.2way_a, m.Om.2way_a2,
     m.Om.2way_b, 
     m.Om.2way_c,
     m.Om.2way_d,
     m.Om.3way_a,
     file="./statistics/brms/231205_densOm_realdiv.RData")





summary(m.Fu.3way_a)
#get slopes:
emt = emtrends(m.Fu.3way_a, specs = c("treatment", "week"), var="realdivLogStd")
summary(emt, point.est=mean, level = .95) 
emt.pairs <- pairs(emt)
summary(emt.pairs, point.est=mean, level = .9)
bayestestR::p_direction(emt.pairs) #probability of direction  
#t1w1 - t1w2: 75.03%
#t2w1 - t2w2: 52.43%
#t3w1 - t3w2: 68.17%
#no "significant differences" -> get rid of 3way interaction

m.Fu.3way_realdiv <- brm(
  bf(Fu_per100g ~ realdivLogStd + treatment + week + 
       realdivLogStd:treatment + realdivLogStd:week + treatment:week +
       realdivLogStd:treatment:week + (1|block/plot),
     hu ~ week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #3 divergent transitions

pp_check(m.Fu.3way_realdiv, ndraws=100)+
  xlim(0,2000)


