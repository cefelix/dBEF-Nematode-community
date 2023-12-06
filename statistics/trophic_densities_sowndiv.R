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


#### exploratrory compare the weeks: ####

dat2 <- dat %>% 
  mutate(sowndivLogStd = as.factor(sowndivLogStd))


#Ba:
  ggplot(dat2, aes(x=sowndivLogStd, y=Ba_per100g, col=treatment, shape=week) )+
    #geom_boxplot(alpha=0.4, outlier.shape = NA)+
    geom_point(position = position_dodge(width=0.7), alpha=0.9)+
    scale_x_discrete(name = "sown plant diversity",
                       labels = c("1", "2", "4", "8", "16"))
  
  #prob of zero in week 1: 7.89 %
    sum(subset(dat, week == "W1")$Ba_per100g == 0) /
      nrow(subset(dat, week == "W1"))             
  # "           in week 2: 1.75%
    sum(subset(dat, week == "W2")$Ba_per100g == 0) /
      nrow(subset(dat, week == "W2"))  
  #ratio 3.46
  
#Fu:
  ggplot(dat2, aes(x=sowndivLogStd, y=Fu_per100g, col=treatment, shape=week) )+
    geom_point(position = position_dodge(width=0.7), alpha=0.9)+
    #geom_boxplot(alpha=0.4, outlier.shape = NA)+
    scale_x_discrete(name = "sown plant diversity",
                     labels = c("1", "2", "4", "8", "16"))
  
  #prob of zero in week 1: 2.63 %
    sum(subset(dat, week == "W1")$Fu_per100g == 0) / 
      nrow(subset(dat, week == "W1")) 
  # "           in week 2: 0.88 %
    sum(subset(dat, week == "W2")$Fu_per100g == 0) /
      nrow(subset(dat, week == "W2")) 
  #ratio 2.30
    
#Pl:
  ggplot(dat2, aes(x=sowndivLogStd, y=Pl_per100g, col=treatment, shape=week) )+
    geom_point(position = position_dodge(width=0.7), alpha=0.9)+
    #geom_boxplot(alpha=0.4, outlier.shape = NA)+
    scale_x_discrete(name = "sown plant diversity",
                     labels = c("1", "2", "4", "8", "16"))  
  
  #prob of zero in week 1: 1.75 %
    sum(subset(dat, week == "W1")$Pl_per100g == 0) / 
      nrow(subset(dat, week == "W1")) 
  # "           in week 2: 0 %
    sum(subset(dat, week == "W2")$Pl_per100g == 0) /
      nrow(subset(dat, week == "W2")) 
  #ratio: inf
    
#Pr:
  ggplot(dat2, aes(x=sowndivLogStd, y=Pr_per100g, col=treatment, shape=week) )+
    geom_point(position = position_dodge(width=0.7), alpha=0.9)+
    #geom_boxplot(alpha=0.4, outlier.shape = NA)+
    scale_x_discrete(name = "sown plant diversity",
                     labels = c("1", "2", "4", "8", "16"))  
  #prob of zero in week 1: 38.60 %
    sum(subset(dat, week == "W1")$Pr_per100g == 0) / 
      nrow(subset(dat, week == "W1"))
  # "           in week 2: 11.40 %
    sum(subset(dat, week == "W2")$Pr_per100g == 0) /
      nrow(subset(dat, week == "W2")) 
  #ratio 2.60  
  
#Om:
  ggplot(dat2, aes(x=sowndivLogStd, y=Om_per100g, col=treatment, shape=week) )+
    geom_point(position = position_dodge(width=0.7), alpha=0.9)+
    #geom_boxplot(alpha=0.4, outlier.shape = NA)+
    stat_summary( fun.y = "mean" )+
    scale_x_discrete(name = "sown plant diversity",
                     labels = c("1", "2", "4", "8", "16"))   
  
  #prob of zero in week 1: 78.07 %
    sum(subset(dat, week == "W1")$Om_per100g == 0) / 
      nrow(subset(dat, week == "W1")) 
  # "           in week 2: 36.84 %
    sum(subset(dat, week == "W2")$Om_per100g == 0) /
      nrow(subset(dat, week == "W2")) 
  #ratio 1.63
    
####Ba ~ sowndiv, stepwise simplification ####

SEED = 22061996  

#most complex, sowndiv: 
  m.Ba.3way_a <- brm(
    bf(Ba_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
         sowndivLogStd:treatment:week + (1|block/plot),
       hu ~ week + sowndivLogStd + week:sowndivLogStd + (1|block/plot)),
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
  #### emmeans term check 3fold ####
  #3fold interaction
  summary(m.Ba.3way_a, prob=0.9 )
  emt = emtrends(m.Ba.3way_a, specs = c("treatment", "week"), var="sowndivLogStd")
  summary(emt, point.est=mean, level = .9) 
  emt.pairs <- pairs(emt)
  summary(emt.pairs, point.est=mean, level = .90)
  bayestestR::p_direction(emt.pairs)
  # t1 week1 - t1 week 2: 79.10%
  # t2 week1 - t2 week 2: 72.07%
  # t3 week1 - t3 week 2: 88.20%
  contrast(emt)
  
  #combinations of treatment and week
  emm1 <- emmeans(m.Ba.3way_a, "week")
  emm2 <- emmeans(m.Ba.3way_a, ~treatment*week)
  emm3 <- emm2[1] + emm1
  summary(emm3)
  contrast(emm3)
  
  
#remove 3 way interaction:
  m.Ba.2way_a <- update(m.Ba.3way_a, .~. -sowndivLogStd:treatment:week,
                        control = list(adapt_delta=0.99,
                                       max_treedepth=10)) #10 div trans
  summary(m.Ba.2way_a, prob=0.9)
  
  ####emmeans term check 2way interactions####
  
#remove week:treatment
  m.Ba.2way_b <- update(m.Ba.2way_a, .~. -treatment:week,
                        control = list(adapt_delta=0.99,
                                       max_treedepth=10))
  summary(m.Ba.2way_b, prob=0.9) 
  pp_check(m.Ba.2way_b, ndraws = 100)
 
    #remove hu~week:sowndiv  
      m.Ba.2way_b2 <- brm(
        bf(Ba_per100g ~ sowndivLogStd + treatment + week + 
             sowndivLogStd:treatment + sowndivLogStd:week + (1|block/plot),
           hu ~ week + sowndivLogStd + (1|block/plot)),
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
        bf(Ba_per100g ~ sowndivLogStd + treatment + week + 
             sowndivLogStd:treatment + sowndivLogStd:week + (1|block/plot),
           hu ~ sowndivLogStd + (1|block/plot)),
        data = dat, 
        family = hurdle_lognormal,
        chains = 3,
        cores = 3,
        iter = 2000, warmup = 1000,
        seed = SEED,
        control = list(adapt_delta=0.99)) #2 div trans
      
      summary(m.Ba.2way_b3, prob=0.9)
      
      #remove hu~sowndivLogStd
      m.Ba.2way_b4 <- brm(
        bf(Ba_per100g ~ sowndivLogStd + treatment + week + 
             sowndivLogStd:treatment + sowndivLogStd:week + (1|block/plot),
           hu ~ week + (1|block/plot)),
        data = dat, 
        family = hurdle_lognormal,
        chains = 3,
        cores = 3,
        iter = 2000, warmup = 1000,
        seed = SEED,
        control = list(adapt_delta=0.99)) #1 div
      
      summary(m.Ba.2way_b4, prob=0.9)
      
      #remove hu ~ sowndivLogStd + week + sowndivLogStd:week + (1|block/plot)
      m.Ba.2way_b5 <- brm(
        bf(Ba_per100g ~ sowndivLogStd + treatment + week + 
             sowndivLogStd:treatment + sowndivLogStd:week + (1|block/plot),
           hu ~ 1),
        data = dat, 
        family = hurdle_lognormal,
        chains = 3,
        cores = 3,
        iter = 2000, warmup = 1000,
        seed = SEED,
        control = list(adapt_delta=0.99)) #4 div trans
      
      summary(m.Ba.2way_b5, prob=0.9)
      
      #assess effects of week:sowndivLogStd with emmeans 
        emt = emtrends(m.Ba.2way_b5, specs = c("treatment", "week"), var="sowndivLogStd")
        summary(emt, point.est=mean, level = .90) 
        emt.pairs <- pairs(emt)
        summary(emt.pairs, point.est=mean, level = .90)
        bayestestR::p_direction(emt.pairs)
        #treat1 w1- treat1 w2: 91.5%
        #treat2 w1- treat2 w2: 91.5%
        #treat2 w1- treat2 w2: 91.5%

           
  
#remove both week:sowndivLogStd  
  m.Ba.2way_c <- update(m.Ba.2way_b5, .~. -week:sowndivLogStd,
                         control = list(adapt_delta=0.99,
                                        max_treedepth=10))
  summary(m.Ba.2way_c)
  
  

####best fit Ba ~ sowndiv model ####    
#remove week  
  m.Ba.2way_c2 <- update(m.Ba.2way_c, .~. -week,
                        control = list(adapt_delta=0.99,
                                       max_treedepth=10))
  summary(m.Ba.2way_c2, prob=0.9)
  
#### Fu ~ sowndiv, hu~1, stepwise simplification####
    
  #hu~1, as there ar only 4 samples with zero fungivores, which is too little to assess prob(Fu==0)
  #!running only one chain to speed up fitting!
  
  #3way interaction  
  m.Fu.3way_a <- brm(
    bf(Fu_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + treatment:week +
         sowndivLogStd:treatment:week + (1|block/plot),
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
  m.Fu.2way_a <- update(m2.Fu.3way_a, .~. -sowndivLogStd:treatment:week) #5 divergs
  summary(m2.Fu.2way_a, prob=0.9) #no effect of sowndivLogStd:week
  
  
  ####best fit Fu~sowndiv model####
  #remove sowndivLogStd:week
  m.Fu.2way_b <- update(m2.Fu.2way_a, .~. -sowndivLogStd:week) #all good
  summary(m2.Fu.2way_b, prob=0.9)
  pp_check(m2.Fu.2way_b, ndraws=100)
  
      #add hu~week*treatment 
      m2.Fu.2way_b2 <- brm(
        bf(Fu_per100g ~ sowndivLogStd + treatment + week + 
             sowndivLogStd:treatment + treatment:week + (1|block/plot),
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
        bf(Fu_per100g ~ sowndivLogStd + treatment + week + 
             sowndivLogStd:treatment + treatment:week + (1|block/plot),
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
        bf(Fu_per100g ~ sowndivLogStd + treatment + week + 
             sowndivLogStd:treatment + treatment:week + (1|block/plot),
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
  
  #remove both   treatment:week AND sowndivLogStd:week
  m.Fu.2way_d <- update(m.Fu.2way_b, .~. -treatment:week) #3 diverg
  summary(m.Fu.2way_d, prob=0.9)
  
  
  
  summary(m2.Fu.3way_a, prob=0.9)
  
#### Pl ~ sowndiv, stepwise simplification ####    
  sum(dat$Pl_per100g == 0) #2 -> hurdle ~ 1
  
  #1 chain for speed
  m.Pl.3way_a <- brm(
    bf(Pl_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
         sowndivLogStd:treatment:week + (1|block/plot),
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
  m.Pl.2way_a <- update(m.Pl.3way_a, .~. -sowndivLogStd:treatment:week ) #1 div trans
  summary(m.Pl.2way_a , prob=0.9)
  
  #remove treatment:week
  m.Pl.2way_b <- update(m.Pl.2way_a, .~. -treatment:week, 
                        chains=3) #2 divs
  summary(m.Pl.2way_b , prob=0.9)
  
  #remove sowndivLogStd:week
  m.Pl.2way_b2 <- update(m.Pl.2way_a, .~. -sowndivLogStd:week) #2 divergent, worse than b
  summary(m.Pl.2way_b2 , prob=0.9)
  
  #remove both treatment:week AND sowndivLogStd:week
  m.Pl.2way_c <- update(m.Pl.2way_b, .~. -sowndivLogStd:week,
                        iter = 3000, warmup = 1500,) #26 divergents
  
  m.Pl.2way_c2 <- update( m.Pl.2way_c, control=list(adapt_delta=0.999)) #11
  summary(m.Pl.2way_c2 , prob=0.9)  
 
  ####best fit Pl~sowndiv model####
  #remove week
  m.Pl.2way_d <- update(m.Pl.2way_c, .~. -week,
                        chains=3) #all good
  summary(m.Pl.2way_d )
  pp_check(m.Pl.2way_d, ndraws=100)
  
#### Pr ~ sowndiv, stepwise simplification ####    
  sum(dat$Pr_per100g == 0) #57 -> hurdle ~ term
  
  beta_coeff_priors <- prior(normal(0,20), class = "b")
  get_prior(bf(Pr_per100g ~ sowndivLogStd + treatment + week + 
                 sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
                 sowndivLogStd:treatment:week + (1|block/plot),
               hu ~ 1),
            data = dat, 
            family = hurdle_lognormal)
  
  #1 chain for speed
  m.Pr.3way_a <- brm(
    bf(Pr_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
         sowndivLogStd:treatment:week + (1|block/plot),
       hu ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + (1|block/plot)),
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
  m.Pr.2way_a <- update(m.Pr.3way_a, .~. -sowndivLogStd:treatment:week,
                        chains=3) #1 div trans, tail ESS too low
    
  summary(m.Pr.2way_a, prob=0.9)  
  
  #remove hu~treatment:week
  m.Pr.2way_a2 <- brm(
    bf(Pr_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + (1|block/plot),
       hu ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + (1|block/plot)),
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
  
  #remove hu~sowndivLogStd:treatment
  m.Pr.2way_b2 <- brm(
    bf(Pr_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week +  (1|block/plot),
       hu ~ sowndivLogStd + treatment + week +  sowndivLogStd:week + 
         (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #1 div, tail ESS too low  
  summary(m.Pr.2way_b2, prob=0.9)
  
  #remove hu~sowndivLogStd:week
  m.Pr.2way_b3 <- brm(
    bf(Pr_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week +  (1|block/plot),
       hu ~ sowndivLogStd + treatment + week + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #6 div transitions
  summary(m.Pr.2way_b3, prob=0.9)
  
  #remove hu~sowndivLogStd
  m.Pr.2way_b4 <- brm(
    bf(Pr_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week +  (1|block/plot),
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
    bf(Pr_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week +  (1|block/plot),
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
    bf(Pr_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week +  (1|block/plot),
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
  m.Pr.2way_c4 <- m.Pr.2way_c
  rm(m.Pr.2way_c)
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
  
  #remove week:sowndivLogStd
  m.Pr.2way_d <- update(m.Pr.2way_c6, .~. - sowndivLogStd:week) #27 div
  summary(m.Pr.2way_d, prob=0.9)
  
  ####best Pr~sowndiv model####
  #remove week 
  m.Pr.2way_d2 <- update(m.Pr.2way_d, .~. -week, 
                         control=list(adapt_delta=0.999)) #0 div trans
  summary(m.Pr.2way_d2, prob=0.9)
  
  
  #remove hu~week
  m.Pr.2way_e <- brm(
    bf(Pr_per100g ~ sowndivLogStd + treatment  + 
         sowndivLogStd:treatment +  (1|block/plot),
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
  
####Om~sowndiv, stepwise simplification####  
  sum(dat$Om_per100g ==0) #131
  #1 chain to be faster
  
  m.Om.3way_a <- brm(
    bf(Om_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
         sowndivLogStd:treatment:week + (1|block/plot),
       hu ~ sowndivLogStd + treatment + week + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 1,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #2 div transitions
  summary(m.Om.3way_a, prob=.9) 
  
  #remove 3way interaction
  m.Om.2way_a <- update(m.Om.3way_a, .~. -sowndivLogStd:treatment:week) #3 div
  summary(m.Om.2way_a, prob=.9) 
  
  #remove hu~sowndivLogStd
  m.Om.2way_a2 <- brm(
    bf(Om_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + (1|block/plot),
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
  
  #remove sowndivLogStd:week
  m.Om.2way_c <- update(m.Om.2way_b, .~. -sowndivLogStd:week) #all good
  summary(m.Om.2way_c, prob=.9) 
  
  
  ####Om best fit ####
  #remove week
  m.Om.2way_d <- update(m.Om.2way_c, .~. -week) #all good
  summary(m.Om.2way_d, prob=.9) 
  pp_check(m.Om.2way_d, ndraw=100)+
    xlim(0,180)

  
####save models####

#best sowndiv models
    save(m.Ba.2way_c2,
         m.Fu.2way_b ,
         m.Pl.2way_d ,
         m.Pr.2way_d2,
         m.Om.2way_d ,
         
    file="./statistics/brms/231205_densBEST_sowndiv.RData")
  
  load(file="./statistics/brms/231205_densBEST_sowndiv.RData")
  
#all Ba~sowndiv models        
 save(m.Ba.2way_a, m.Ba.2way_b,
      m.Ba.2way_b2, m.Ba.2way_b3, m.Ba.2way_b4, m.Ba.2way_b5,
      m.Ba.2way_c, m.Ba.2way_c2,
      m.Ba.3way_a,
      file="./statistics/brms/231205_densBa_sowndiv.RData")
    
 #all Fu~sowndiv models   
  save(m.Fu.2way_a, m.Fu.2way_b,
       m.Fu.2way_b2, m.Fu.2way_b3, m.Fu.2way_b4,
       m.Fu.2way_c, m.Fu.2way_d,
       m.Fu.3way_a,
       file="./statistics/brms/231205_densFu_sowndiv.RData")
  
  #all Pl~sowndiv models
    save(m.Pl.2way_a, m.Pl.2way_b,
         m.Pl.2way_b2,
         m.Pl.2way_c, m.Pl.2way_d,
         m.Pl.3way_a,
         file="./statistics/brms/231205_densPl_sowndiv.RData")
    
  #all Pr~sowndiv models
    save(m.Pr.2way_a, m.Pr.2way_a2,
         m.Pr.2way_b, m.Pr.2way_b2, m.Pr.2way_b3, m.Pr.2way_b4, m.Pr.2way_b5, m.Pr.2way_b6,
         m.Pr.2way_c4, m.Pr.2way_c42, m.Pr.2way_c5, m.Pr.2way_c52,
         m.Pr.2way_c6, m.Pr.2way_c62,
         m.Pr.2way_d, m.Pr.2way_d2, m.Pr.2way_e, m.Pr.2way_e2,
         m.Pr.3way_a,
         
         file="./statistics/brms/231205_densPr_sowndiv.RData")
    
  #all Om~sowndiv
    save(m.Om.2way_a, m.Om.2way_a2,
         m.Om.2way_b, 
         m.Om.2way_c,
         m.Om.2way_d,
         m.Om.3way_a,
         file="./statistics/brms/231205_densOm_sowndiv.RData")

    
        
    
    
    summary(m.Fu.3way_a)
    #get slopes:
    emt = emtrends(m.Fu.3way_a, specs = c("treatment", "week"), var="sowndivLogStd")
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
    

    
#### plot predictions from 2way-models ####    
    load("./statistics/brms/231201_FuBa_hu_week.RData")  
    
    #get predictions:
    predictions.Ba.sown <- conditional_effects(m.Ba.2way_b)[[4]]
    predictions.Ba.real <- conditional_effects(m.Ba.2way_b_realdiv)[[4]]
    predictions.Fu.sown <- conditional_effects(m.Fu.2way_b)[[4]]
    predictions.Fu.real <- conditional_effects(m.Fu.2way_b_realdiv)[[4]]   
    
    #our data as in the models
    dat <- subset(dBEF_nem21, sowndiv != 60) #60 sp. plots
    dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                          .after = sowndivLog) %>%
      mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
             .after = realdivLog)
    
    #Ba sowndiv
    predictions = predictions.Ba.sown 
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    BREAKS = unique(dat$sowndivLogStd)
    ggplot(dat, aes(x=sowndivLogStd, y=Ba_per100g) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.2, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=predictions, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "Ba per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")
    
    #Ba realdiv
    predictions = predictions.Ba.real 
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    BREAKS = unique(dat$sowndivLogStd)
    ggplot(dat, aes(x=realdivLogStd, y=Ba_per100g) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.075, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=predictions, aes(x= realdivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "realized plant diversity", breaks = BREAKS,
                         labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "Ba per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")
    
    #Fu sowndiv
    predictions = predictions.Fu.sown 
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    BREAKS = unique(dat$sowndivLogStd)
    ggplot(dat, aes(x=sowndivLogStd, y=Fu_per100g) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.2, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=predictions, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "Fu per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")
    
    #Fu realdiv
    predictions = predictions.Fu.real
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    BREAKS = round(unique(dat$realdivLogStd), digits = 1)
    ggplot(dat, aes(x=realdivLogStd, y=Fu_per100g) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.05, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=predictions, aes(x= realdivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "realized plant diversity", breaks = BREAKS)+#,
                         #labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "Fu per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")
    
####plot predictions from 3-way models####
    
    conditions <- make_conditions(m.Fu.3way_a, "week")
    predictions.Fu3way.sown <- conditional_effects(m.Fu.3way_a, "sowndivLogStd:treatment", conditions = conditions)[[1]]
    
    
    #Fu sowndiv
    predictions = predictions.Fu3way.sown
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    BREAKS = -1*unique(dat$sowndivLogStd)
    ggplot(dat, aes(x=sowndivLogStd, y=Fu_per100g) )+
      #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
      #                                  fill=treatment), 
      #            alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.05, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=predictions, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype=week, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "realized plant diversity", breaks = BREAKS,
        labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "Fu per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")
    
    summary(m.Fu.3way_a)
    
    
    
    
####check priors####
  get_prior( bf(Fu_per100g ~ realdivLogStd + treatment + week + 
           realdivLogStd:treatment + (1|block/plot),
         hu ~ week + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal)  
    
    