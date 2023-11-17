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

#set a seed:
SEED = 22061996

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
  dBEF_nemSH1 <- subset(dBEF_nem, SH == 1)
  dBEF_nemSH5 <- subset(dBEF_nem, SH == 5)
  dBEF_nemSH15 <- subset(dBEF_nem, SH == 15)
  dBEF_nemSH19 <- subset(dBEF_nem, SH == 19)
  
  dBEF_nem21 <- subset(dBEF_nem, year == 2021)



#### exploration #### 
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

#### 1 - both years ####
  #load the fit models:
  load("./statistics/brms/231101_Fu_allData.RData")

hist(dBEF_nem$Fu_per100gLog)
m.0.Fu11 <- brm(Fu_per100gLog ~ sowndiv*treatment + year + (1|block),
                     data = dBEF_nem, family = "gaussian",
                     chains = 3,
                     cores = 3,
                     iter = 2000, warmup = 1000,
                     control = list(adapt_delta=0.99)) 
pp_check(m.0.Fu11, ndraws=100)

#lets make a yearblock variable, and use a random slope:
#https://www.r-bloggers.com/2019/09/bayesian-linear-mixed-models-random-intercepts-slopes-and-missing-data/

dBEF_nem <- dBEF_nem %>%
  mutate(yearblock = as.factor(paste(dBEF_nem$year, dBEF_nem$block, sep="")),
         .after = block)

m.0.Fu12 <- brm(Fu_per100gLog ~ sowndiv*SH + (Fu_per100gLog|yearblock),
                data = dBEF_nem, family = "gaussian",
                chains = 3,
                cores = 3,
                iter = 2000, warmup = 1000,
                control = list(adapt_delta=0.99)) 
summary(m.0.Fu12)
pp_check(m.0.Fu12, ndraws=100)
 m.0.Fu13 <- update()




#### 2 - 2021's data####
dBEF_nem21 <- subset(dBEF_nem, year==2021)
  #one subset for each treatment (to plot) 
  dBEF_nem21_t1 <- subset(dBEF_nem21, treatment == 1) 
  dBEF_nem21_t2 <- subset(dBEF_nem21, treatment == 2)
  dBEF_nem21_t3 <- subset(dBEF_nem21, treatment == 1)

  
#exploration:
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



#### hurdle11: Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_lognormal ####
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


####saving models####
save(m.Fu.hurdle11a, m.Fu.hurdle12a,
     file="./statistics/brms/231117_Fu.RData")

load(file="./statistics/brms/231107_Fu_hurdle.RData")
conditional_effects(m.Fu.hurdle51)






