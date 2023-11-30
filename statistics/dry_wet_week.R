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
                      .after = sowndivLog)


#### exploratrory compare the weeks: ####

dat2 <- dat %>% 
  mutate(sowndivLogStd = as.factor(sowndivLogStd))


#Ba:
  ggplot(dat2, aes(x=sowndivLogStd, y=Ba_per100g, col=treatment, shape=week) )+
    #geom_boxplot(alpha=0.4, outlier.shape = NA)+
    geom_point(position = position_dodge(width=0.7), alpha=0.9)+
    scale_x_discrete(name = "sown plant diversity",
                       labels = c("1", "2", "4", "8", "16"))
  
#Fu:
  ggplot(dat2, aes(x=sowndivLogStd, y=Fu_per100g, col=treatment, shape=week) )+
    geom_point(position = position_dodge(width=0.7), alpha=0.9)+
    #geom_boxplot(alpha=0.4, outlier.shape = NA)+
    scale_x_discrete(name = "sown plant diversity",
                     labels = c("1", "2", "4", "8", "16"))
  
#Pl:
  ggplot(dat2, aes(x=sowndivLogStd, y=Pl_per100g, col=treatment, shape=week) )+
    geom_point(position = position_dodge(width=0.7), alpha=0.9)+
    #geom_boxplot(alpha=0.4, outlier.shape = NA)+
    scale_x_discrete(name = "sown plant diversity",
                     labels = c("1", "2", "4", "8", "16"))  
  
#Pr:
  ggplot(dat2, aes(x=sowndivLogStd, y=Pr_per100g, col=treatment, shape=week) )+
    geom_point(position = position_dodge(width=0.7), alpha=0.9)+
    #geom_boxplot(alpha=0.4, outlier.shape = NA)+
    scale_x_discrete(name = "sown plant diversity",
                     labels = c("1", "2", "4", "8", "16"))  
  
#Om:
  ggplot(dat2, aes(x=sowndivLogStd, y=Om_per100g, col=treatment, shape=week) )+
    geom_point(position = position_dodge(width=0.7), alpha=0.9)+
    #geom_boxplot(alpha=0.4, outlier.shape = NA)+
    scale_x_discrete(name = "sown plant diversity",
                     labels = c("1", "2", "4", "8", "16"))    
  
####bacterivores 3way ####
  
#most complex: 
  m.Ba.3way_a <- brm(
    bf(Ba_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
         sowndivLogStd:treatment:week + (1|block/plot),
       hu ~ week + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #1 divergent transition
  
  m.Ba.3way_a <- update(m.Ba.3way_a ,
                         control=list(adapt_delta=0.999))  #all good
  
  #get slopes:
  emt = emtrends(m.Ba.3way_a2, specs = c("treatment", "week"), var="sowndivLogStd")
  summary(emt, point.est=mean, level = .95) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .95)
    bayestestR::p_direction(emt.pairs) #probability of direction
      #t1w1 - t1w2 78.60%
      #t2w1 - t2w2 70.33%
      #t3w1 - t3w2 88.27%
    
    #none "significant" -> exclude 3 way interaction:

####bacterivores 2way 1 ####        
    m.Ba.2way_a <- brm(
      bf(Ba_per100g ~ sowndivLogStd + treatment + week + 
           sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + (1|block/plot),
         hu ~ week + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #1 diverg
    
    m.Ba.2way_a2 <- update(m.Ba.2way_a,
                             control=list(adapt_delta=0.999)) #986 exceeded max_treedepth 
    
    m.Ba.2way_a2 <- update(m.Ba.2way_a,
                             control=list(adapt_delta=0.999,
                                          max_treedepth=12)) #4divergent transitions
    
    m.Ba.2way_a2 <- update(m.Ba.2way_a,
                          control=list(adapt_delta=0.9999,
                                       max_treedepth=12)) #all good
    
    pp_check(m.Ba.2way_a2, ndraws=100)+ 
      xlim(0,1000)
    
    #get slopes:
    emt = emtrends(m.Ba.2way_a2, specs = c("treatment", "week"), var="sowndivLogStd")
    summary(emt, point.est=mean, level = .95) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .95)
    bayestestR::p_direction(emt.pairs) #probability of direction
  
####bacterivores 2way 2 ####        
    m.Ba.2way_b <- brm(
      bf(Ba_per100g ~ sowndivLogStd + treatment + week + 
           sowndivLogStd:treatment + (1|block/plot),
         hu ~ week + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) 

    
    pp_check(m.Ba.2way_b, ndraws=100)+ 
      xlim(0,1000)
    
    #get slopes:
    emt = emtrends(m.Ba.2way_b, specs = c("treatment", "week"), var="sowndivLogStd")
    summary(emt, point.est=mean, level = .95) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .95)
    bayestestR::p_direction(emt.pairs) #probability of direction  
    
    
    
    
#compare and decide    
    loo(m.Ba.2way_b, m.Ba.2way_a2, m.Ba.3way_a2)
  
  
  
#### fungivores 3 way ####
  
  m.Fu.3way_a <- brm(
    bf(Fu_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + treatment:week +
         sowndivLogStd:treatment:week + (1|block/plot),
       hu ~ week + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #all good
    
    #get slopes:
    emt = emtrends(m.Fu.3way_a, specs = c("treatment", "week"), var="sowndivLogStd")
    summary(emt, point.est=mean, level = .95) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .95)
    bayestestR::p_direction(emt.pairs) #probability of direction  
    #t1w1 - t1w2: 75.03%
    #t2w1 - t2w2: 52.43%
    #t3w1 - t3w2: 68.17%
    #no "significant differences" -> get rid of 3way interaction
    
#### fungivores 2 way 1.1 ####
    
    m.Fu.2way_a <- brm(
      bf(Fu_per100g ~ sowndivLogStd + treatment + week + 
           sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + (1|block/plot),
         hu ~ week + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #all good
    
    #get slopes:
    emt = emtrends(m.Fu.2way_a, specs = c("treatment", "week"), var="sowndivLogStd")
    summary(emt, point.est=mean, level = .95) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .95)
    bayestestR::p_direction(emt.pairs) #probability of direction      
    
#### fungivores 2 way 1.2 ####
    
    m.Fu.2way_a2 <- brm(
      bf(Fu_per100g ~ sowndivLogStd + treatment + week + 
           sowndivLogStd:treatment + sowndivLogStd:week + (1|block/plot),
         hu ~ week + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #6 divergent transitions
    
    m.Fu.2way_a2 <- update(m.Fu.2way_a2, 
                           control = list(adapt_delta=0.999)) #all good
    
    #get slopes:
    emt = emtrends(m.Fu.2way_a2, specs = c("treatment", "week"), var="sowndivLogStd")
    summary(emt, point.est=mean, level = .95) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .95)
    bayestestR::p_direction(emt.pairs) #probability of direction         
    
#### fungivores 2 way 1.3 ####
    
    m.Fu.2way_a3 <- brm(
      bf(Fu_per100g ~ sowndivLogStd + treatment + week + 
           sowndivLogStd:treatment + sowndivLogStd:week + (1|block/plot),
         hu ~ week + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #a6 divergent transitions
    
      m.Fu.2way_a3 <- update(m.Fu.2way_a3, 
                             control = list(adapt_delta=0.999)) #all good
      
      pp_check(m.Fu.2way_a3, ndraws=100)+
        xlim(0,1000)
    
    #get slopes:
    emt = emtrends(m.Fu.2way_a3, specs = c("treatment", "week"), var="sowndivLogStd")
    summary(emt, point.est=mean, level = .95) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .95)
    bayestestR::p_direction(emt.pairs) #probability of direction         
    
#### fungivores 2 way 2 ####
    
    m.Fu.2way_b <- brm(
      bf(Fu_per100g ~ sowndivLogStd + treatment + week + 
           sowndivLogStd:treatment + (1|block/plot),
         hu ~ week + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #2 divergent transitions
    
    m.Fu.2way_b <- update(m.Fu.2way_b,
                          control = list(adapt_delta=0.999)) #all good
    
    summary(m.Fu.2way_b)
    pp_check(m.Fu.2way_b2, ndraws=100)+
      xlim(0,1000) #but rather bad fit
    
    #get slopes:
    emt = emtrends(m.Fu.2way_b, specs = c("treatment"), var="sowndivLogStd")
    summary(emt, point.est=mean, level = .95) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .95)
    bayestestR::p_direction(emt.pairs) #probability of direction 
  
    
    ####save fu and ba models####
    save(m.Ba.3way_a,
         m.Ba.2way_a,
         m.Ba.2way_b,
      
         m.Fu.3way_a,
         m.Fu.2way_a,
         m.Fu.2way_a2,
         m.Fu.2way_a3,
         m.Fu.2way_b, 
         file="./statistics/brms/231130_FuBa_hu_week.RData")
    
    
    