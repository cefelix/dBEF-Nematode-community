#modelling sowndiv*treatment effects on density of each trophic guild
  #problem: we have quite a few samples where we have zero nematodes of a certain group
  #thus, we need to use distributions which allow to contain zeros: 
  #https://www.andrewheiss.com/blog/2022/05/09/hurdle-lognormal-gaussian-brms/

####data and packages####
library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(hexbin)
library(GGally)
library(magrittr)

#set a seed:
SEED = 22061996

#to run: 11a, 22b


#### exploration #### 
  dBEF_nemSH1 <- subset(dBEF_nem, SH == 1)
  dBEF_nemSH5 <- subset(dBEF_nem, SH == 5)
  dBEF_nemSH15 <- subset(dBEF_nem, SH == 15)
  dBEF_nemSH19 <- subset(dBEF_nem, SH == 19)  
  
p.1 <- ggplot(dBEF_nemSH1, aes(y = Fu_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  labs(title = "SH 1")+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")

p.5 <- ggplot(dBEF_nemSH5, aes(y = Fu_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")+
  labs(title = "SH 5")

p.15 <- ggplot(dBEF_nemSH15, aes(y = Fu_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 15")+
  geom_smooth(method="lm")

p.19 <- ggplot(dBEF_nemSH19, aes(y = Fu_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 19")+
  geom_smooth(method="lm")

grid.arrange(p.1, p.5, p.15, p.19)
          rm(p.1, p.5, p.15, p.19)




#### exploration - 2021's data####
  dBEF_nem21_t1 <- subset(dBEF_nem21, treatment == 1) 
  dBEF_nem21_t2 <- subset(dBEF_nem21, treatment == 2)
  dBEF_nem21_t3 <- subset(dBEF_nem21, treatment == 1)

  
  sum(dBEF_nem21$Fu_per100g == 0) #4 of 240 samples have zero fungivores

#a jitter plot with an OLS regression line
p.all <- ggplot(dBEF_nem21, aes(x = log(sowndiv), y = Fu_per100g))+
          geom_jitter(aes(col=block))+
          labs(title = "all treatments")

p.t1 <- ggplot(dBEF_nem21_t1, aes(x = log(sowndiv), y = Fu_per100g))+
          geom_jitter(aes(col=block))+
          geom_smooth(method="lm")+
          labs(title = "-SH -PH")

p.t2 <- ggplot(dBEF_nem21_t2, aes(x = log(sowndiv), y = Fu_per100g))+
          geom_jitter(aes(col=block))+
          geom_smooth(method="lm")+
          labs(title = "+SH -PH")

p.t3 <- ggplot(dBEF_nem21_t2, aes(x = log(sowndiv), y = Fu_per100g))+
          geom_jitter(aes(col=block))+
          geom_smooth(method="lm")+
          labs(title = "+SH +PH")

grid.arrange(p.all, p.t1, p.t2, p.t3)
  rm(dBEF_nem21_t1, dBEF_nem21_t2, dBEF_nem21_t3,
     p.all, p.t1, p.t2, p.t3)



#### 11a hurdle: Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~1, family = hurdle_lognormal ####
m.Fu.hurdle11a <- brm(
  bf(Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = "hurdle_lognormal",
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Fu.hurdle11, ndraws = 100)+
  xlim(0, 3000) #this concentrates a lot of probability mass on the mode

m.Fu.hurdle12a <- update(m.Fu.hurdle11,
                        control = list(adapt_delta=0.999))





#### 21a  hurdle: Fu_per100g ~ sowndivLog*treatment + (1|block/plot), hu~1, family = hurdle_lognormal ####
m.Fu.hurdle21a <- brm(
  bf(Fu_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = "hurdle_lognormal",
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Fu.hurdle21a, ndraws = 100)+
  xlim(0, 3000)


#### 21b  hurdle: Fu_per100g ~ sowndivLog*treatment + (1|block/plot), hu~term, family = hurdle_lognormal ####

m.Fu.hurdle21b <- brm(
  bf(Fu_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ sowndivLog*treatment + (1|block/plot)),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
  #6 divergent transitions
  #1946 exceeded max_treedepth
  #bulk/tail EES too low

m.Fu.hurdle22b <- update(m.Fu.hurdle21b,
                         iter = 3000, warmup = 1500,
                         control = list(adapt_delta = 0.999,
                                        max_treedepth=12))  

loo(m.Fu.hurdle21a, m.Fu.hurdle21b)

####31a hurdle: Fu_per100g ~ sowndivLog*treatment + (1|block/plot), hu~1, family = hurdle_lognormal, no 60 sp plots ####
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60)

m.Fu.hurdle31a <- brm(
  bf(Fu_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Fu.hurdle31a, ndraws=100)+
  xlim(0,3000)
summary(m.Fu.hurdle31a)


####31b hurdle: Fu_per100g ~ sowndivLog*treatment + (1|block/plot), hu~term, family = hurdle_lognormal, no 60 sp plots ####
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60)

m.Fu.hurdle31b <- brm(
  bf(Fu_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ sowndivLog*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Fu.hurdle31b, ndraws=100)
summary(m.Fu.hurdle31b)

conditional_effects(m.Fu.hurdle31b)

#### 31c hurdle: Fu_per100g ~ sowndivLog*treatment + (1|block/plot), hu~treatment, family = hurdle_lognormal, no 60 sp plots ####
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60)

m.Fu.hurdle31c <- brm(
  bf(Fu_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Fu.hurdle31c, ndraws=100)

#### 31d hurdle: Fu_per100g ~ sowndivLog*treatment + (1|block/plot), hu~sowndiv, family = hurdle_lognormal, no 60 sp plots ####
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60)

m.Fu.hurdle31d <- brm(
  bf(Fu_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ sowndivLog + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Fu.hurdle31d, ndraws=100)

#### compare the different hurdle model predictions: ####

#1 to 60 sp
predictions21a <- conditional_effects(m.Fu.hurdle21a)[[3]]
predictions21b <- conditional_effects(m.Fu.hurdle21b)[[3]]

ggplot(dBEF_nem21, aes(x=sowndivLog, y=Fu_per100g, col=treatment) )+
  geom_jitter(width=0.3, alpha=0.5)+
  geom_smooth(data=predictions21a, aes(x= sowndivLog, y=estimate__, col=treatment), 
              linetype="dashed")+
  geom_smooth(data=predictions21b, aes(x= sowndivLog, y=estimate__, col=treatment), 
              linetype="solid", alpha=0.7)

# 1 to 16 sp
predictions31a <- conditional_effects(m.Fu.hurdle31a)[[3]]
  predictions31b <- conditional_effects(m.Fu.hurdle31b)[[3]]
  predictions31c <- conditional_effects(m.Fu.hurdle31c)[[3]]
  predictions31d <- conditional_effects(m.Fu.hurdle31d)[[3]]

ggplot(dBEF_nem21, aes(x=sowndivLog, y=Fu_per100g, col=treatment) )+
  geom_jitter(width=0.3, alpha=0.5)+
  geom_smooth(data=predictions31c, aes(x= sowndivLog, y=estimate__, col=treatment), 
              linetype="dashed")+ #hu~treatment
  geom_smooth(data=predictions31d, aes(x= sowndivLog, y=estimate__, col=treatment), 
              linetype="solid", alpha=0.7)#+ #hu~diversity

  geom_smooth(data=predictions31b, aes(x= sowndivLog, y=estimate__, col=treatment), 
              linetype="dotted") #hu~term


#60 vs 16 sp max
ggplot(dBEF_nem21, aes(x=sowndivLog, y=Fu_per100g, col=treatment) )+
  geom_jitter(width=0.3, alpha=0.5)+
  geom_smooth(data=predictions31a, aes(x= sowndivLog, y=estimate__, col=treatment), 
              linetype="solid", alpha=0.7)+
  geom_smooth(data=predictions21a, aes(x= sowndivLog, y=estimate__, col=treatment), 
              linetype="dashed", alpha=0.7)

#hurdle~term gives very weird results (like negative densities at low diversity)
#maybe a justification to use hu~1 is, that we have actually very little data with zeros,
#but those few datapoints have a huge influence on the regression line
loo(m.Fu.hurdle31a, m.Fu.hurdle31b,
    m.Fu.hurdle31c, m.Fu.hurdle31d) 
      #31a:2883.3, 31b:2871.7
      #31c:2885.3, 31d:2887.7



####saving models####
save(#m.Fu.hurdle11a, m.Fu.hurdle12a,
     m.Fu.hurdle21a,
     m.Fu.hurdle21b, #m.Fu.hurdle22b,
     m.Fu.hurdle31a,
     m.Fu.hurdle31b,
     m.Fu.hurdle31c,
     m.Fu.hurdle31d,
     file="./statistics/brms/231117_Fu.RData")

load(file="./statistics/brms/231107_Fu_hurdle.RData")
conditional_effects(m.Fu.hurdle51)






