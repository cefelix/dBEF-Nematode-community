#cp densities

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


#cp1:
ggplot(dat2, aes(x=sowndivLogStd, y=cp1_per100g, col=treatment, shape=week) )+
  #geom_boxplot(alpha=0.4, outlier.shape = NA)+
  geom_point(position = position_dodge(width=0.7), alpha=0.9)+
  scale_x_discrete(name = "sown plant diversity",
                   labels = c("1", "2", "4", "8", "16"))

#prob of zero in week 1: 53.51 %
sum(subset(dat, week == "W1")$cp1_per100g == 0) / #61
  nrow(subset(dat, week == "W1"))             
# "           in week 2: 42.11%
sum(subset(dat, week == "W2")$cp1_per100g == 0) / #48
  nrow(subset(dat, week == "W2"))  


#cp2:
ggplot(dat2, aes(x=sowndivLogStd, y=cp2_per100g, col=treatment, shape=week) )+
  geom_point(position = position_dodge(width=0.7), alpha=0.9)+
  #geom_boxplot(alpha=0.4, outlier.shape = NA)+
  scale_x_discrete(name = "sown plant diversity",
                   labels = c("1", "2", "4", "8", "16"))

#prob of zero in week 1: 0 %
sum(subset(dat, week == "W1")$cp2_per100g == 0) / 
  nrow(subset(dat, week == "W1")) 
# "           in week 2: 0.88 %
sum(subset(dat, week == "W2")$cp2_per100g == 0) /
  nrow(subset(dat, week == "W2")) 


#cp3:
ggplot(dat2, aes(x=sowndivLogStd, y=cp3_per100g, col=treatment, shape=week) )+
  geom_point(position = position_dodge(width=0.7), alpha=0.9)+
  #geom_boxplot(alpha=0.4, outlier.shape = NA)+
  scale_x_discrete(name = "sown plant diversity",
                   labels = c("1", "2", "4", "8", "16"))  

#prob of zero in week 1: 3.51 %
sum(subset(dat, week == "W1")$cp3_per100g == 0) / 
  nrow(subset(dat, week == "W1")) 
# "           in week 2: 0 %
sum(subset(dat, week == "W2")$cp3_per100g == 0) /
  nrow(subset(dat, week == "W2")) 
#ratio: inf

#cp4:
ggplot(dat2, aes(x=sowndivLogStd, y=cp4_per100g, col=treatment, shape=week) )+
  geom_point(position = position_dodge(width=0.7), alpha=0.9)+
  #geom_boxplot(alpha=0.4, outlier.shape = NA)+
  scale_x_discrete(name = "sown plant diversity",
                   labels = c("1", "2", "4", "8", "16"))  
#prob of zero in week 1: 11.40 %
sum(subset(dat, week == "W1")$cp4_per100g == 0) / 
  nrow(subset(dat, week == "W1"))
# "           in week 2: 1.74 %
sum(subset(dat, week == "W2")$cp4_per100g == 0) /
  nrow(subset(dat, week == "W2")) 
#ratio 2.60  

#cp5:
ggplot(dat2, aes(x=sowndivLogStd, y=cp5_per100g, col=treatment, shape=week) )+
  geom_point(position = position_dodge(width=0.7), alpha=0.9)+
  #geom_boxplot(alpha=0.4, outlier.shape = NA)+
  scale_x_discrete(name = "sown plant diversity",
                   labels = c("1", "2", "4", "8", "16"))   

#prob of zero in week 1: 92.11 %
sum(subset(dat, week == "W1")$cp5_per100g == 0) / 
  nrow(subset(dat, week == "W1")) 
# "           in week 2: 84.21 %
sum(subset(dat, week == "W2")$cp5_per100g == 0) /
  nrow(subset(dat, week == "W2")) 
#ratio 1.63


####cp1 ~ sowndiv, stepwise simplification ####

SEED = 22061996  

#most complex, sowndiv: 
m.cp1.3way_a <- brm(
  bf(cp1_per100g ~ sowndivLogStd + treatment + week + 
       sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
       sowndivLogStd:treatment:week + (1|block/plot),
     hu ~ week + sowndivLogStd + week:sowndivLogStd + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #18 div trans

summary(m.cp1.3way_a, prob=0.9)

  #remove 3 way
  m.cp1.2way_a <- update(m.cp1.3way_a, .~. -sowndivLogStd:treatment:week) #1 div
  summary(m.cp1.2way_a, prob=0.9)

  #remove hu~ week:sowndivLogStd
    m.cp1.2way_b <- brm(
      bf(cp1_per100g ~ sowndivLogStd + treatment + week + 
           sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + (1|block/plot),
         hu ~ week + sowndivLogStd + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #3 div trans
    summary(m.cp1.2way_b, prob=0.9)
    
  #remove sowndivLogStd:week 
  m.cp1.2way_b2 <- update( m.cp1.2way_b, .~. - sowndivLogStd:week ) #9 div
  summary(m.cp1.2way_b2, prob=0.9)
  
  #remove hu~week
  m.cp1.2way_c <- brm(
    bf(cp1_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + treatment:week + (1|block/plot),
       hu ~ sowndivLogStd + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #
  summary(m.cp1.2way_c, prob=0.9) #8 div trans
  
  #remove hu~sowndiv
  m.cp1.2way_d <- brm(
    bf(cp1_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + treatment:week + (1|block/plot),
       hu ~ 1),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #2 div
  
  ####best fit cp1####
  m.cp1.2way_d <- update(m.cp1.2way_d,
                         control = list(adapt_delta=0.999))
  summary(m.cp1.2way_d, prob=0.9)
  
  
  
####cp2 ~ sowndiv, stepwise simplification ####

#most complex, sowndiv: 
m.cp2.3way_a <- brm(
  bf(cp2_per100g ~ sowndivLogStd + treatment + week + 
       sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
       sowndivLogStd:treatment:week + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good
summary(m.cp2.3way_a, prob=0.9)

#remove 3 way
m.cp2.2way_a <- update(m.cp2.3way_a, .~.-sowndivLogStd:treatment:week)
summary(m.cp2.2way_a, prob=0.9)

#remove treatment:week
m.cp2.2way_a2 <- update(m.cp2.2way_a, .~.-treatment:week) #19 div
summary(m.cp2.2way_a2, prob=0.9)

#remove sowndivLogStd:week
m.cp2.2way_a3 <- update(m.cp2.2way_a2, .~.-sowndivLogStd:week) #1 div
summary(m.cp2.2way_a3, prob=0.9) 

#remove week
m.cp2.2way_a4 <- update(m.cp2.2way_a3, .~.-week) #1 div trans


####best fit cp2 model####
m.cp2.2way_a4 <- update(m.cp2.2way_a4, 
                        control=list(adapt_delta=0.999))
summary(m.cp2.2way_a4, prob=0.9) 
 

####cp3 ~ sowndiv, stepwise simplification ####

  #most complex, sowndiv: 
  m.cp3.3way_a <- brm(
    bf(cp3_per100g ~ sowndivLogStd + treatment + week + 
         sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
         sowndivLogStd:treatment:week + (1|block/plot),
       hu ~ sowndivLogStd + treatment + treatment:sowndivLogStd + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #1 div
  summary(m.cp3.3way_a, prob=0.9)
  
  #remove 3way interaction
   m.cp3.2way_a <- update(m.cp3.3way_a, .~. -sowndivLogStd:treatment:week ) #all god
   summary(m.cp3.2way_a, prob=0.9)
   
  #remove hu~ treatment:sowndivLogStd 
   m.cp3.2way_a2 <- brm(
     bf(cp3_per100g ~ sowndivLogStd + treatment + week + 
          sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + (1|block/plot),
        hu ~ sowndivLogStd + treatment + (1|block/plot)),
     data = dat, 
     family = hurdle_lognormal,
     chains = 3,
     cores = 3,
     iter = 2000, warmup = 1000,
     seed = SEED,
     control = list(adapt_delta=0.99))  #1 div 
   summary(m.cp3.2way_a2, prob=0.9)
   
   #remove hu~sowndivLogStd 
   m.cp3.2way_a3 <- brm(
     bf(cp3_per100g ~ sowndivLogStd + treatment + week + 
          sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + (1|block/plot),
        hu ~  treatment + (1|block/plot)),
     data = dat, 
     family = hurdle_lognormal,
     chains = 3,
     cores = 3,
     iter = 2000, warmup = 1000,
     seed = SEED,
     control = list(adapt_delta=0.99))
   summary(m.cp3.2way_a3, prob=0.9)
   
   #remove hu~treatment
   m.cp3.2way_a4 <- brm(
     bf(cp3_per100g ~ sowndivLogStd + treatment + week + 
          sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + (1|block/plot),
        hu ~  1),
     data = dat, 
     family = hurdle_lognormal,
     chains = 3,
     cores = 3,
     iter = 2000, warmup = 1000,
     seed = SEED,
     control = list(adapt_delta=0.99))
   summary(m.cp3.2way_a4, prob=0.9) #12 div trans, tail ESS too low
   
   #remove treatment:week 
   m.cp3.2way_b <- update(m.cp3.2way_a4, .~. -treatment:week) #5 div
   summary(m.cp3.2way_b, prob=0.9) 
   
   #remove sowndivLogStd:week
   m.cp3.2way_c <- update(m.cp3.2way_b, .~. -sowndivLogStd:week) #8 div
   summary(m.cp3.2way_c, prob=0.9) 

   #remove treatment:week and sowndivLogStd
   m.cp3.2way_c2 <- update(m.cp3.2way_c, .~. -treatment:week) #4 div
   
   ####best fit cp3 model####
   m.cp3.2way_c2 <- update(m.cp3.2way_c2,  control = list(adapt_delta=0.999)) #all good
   summary(m.cp3.2way_c2, prob=0.9) 
   
   
####cp4 ~ sowndiv, stepwise simplification ####

#most complex, sowndiv: 
m.cp4.3way_a <- brm(
  bf(cp4_per100g ~ sowndivLogStd + treatment + week + 
       sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
       sowndivLogStd:treatment:week + (1|block/plot),
     hu ~ sowndivLogStd + treatment + week + 
       sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
       sowndivLogStd:treatment:week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
summary(m.cp4.3way_a, prob=0.9)


####cp5 ~ sowndiv, stepwise simplification ####

#most complex, sowndiv: 
m.cp5.3way_a <- brm(
  bf(cp5_per100g ~ sowndivLogStd + treatment + week + 
       sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
       sowndivLogStd:treatment:week + (1|block/plot),
     hu ~ sowndivLogStd + treatment + week + 
       sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + 
       sowndivLogStd:treatment:week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #15 div
summary(m.cp5.3way_a, prob=0.9)



#### save models####
#best fit:
save(m.cp1.2way_d, 
     m.cp2.2way_a4,
     m.cp3.2way_c2, 
     m.cp4.3way_a,
     m.cp5.3way_a,
     file ="./statistics/brms/231206_densBEST_CP_sowndiv.RData")


#cp1 

save(m.cp1.2way_a, m.cp1.2way_b, m.cp1.2way_b2, 
     m.cp1.2way_c, m.cp1.2way_d, m.cp1.3way_a, 
     file ="./statistics/brms/231206_densCP1_sowndiv.RData")

#cp2
save(m.cp2.2way_a, m.cp2.2way_a2, m.cp2.2way_a3, m.cp2.2way_a4,
file ="./statistics/brms/231206_densCP2_sowndiv.RData")

#cp3 
save(m.cp3.2way_a, m.cp3.2way_a2, m.cp3.2way_a3, m.cp3.2way_a4,
     m.cp3.2way_b,
     m.cp3.2way_c, m.cp3.2way_c2,
     m.cp3.3way_a, m.cp3.3way_b, 
     
  file ="./statistics/brms/231206_densCP3_sowndiv.RData")

