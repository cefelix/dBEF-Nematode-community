####data and packages####
#problem: we have quite a few samples where we have zero nematodes of a certain group
#thus, we need to use distributions which allow to contain zeros: 
#https://www.andrewheiss.com/blog/2022/05/09/hurdle-lognormal-gaussian-brms/


library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
dBEF_nem %>%
  str()

####2017's data####
dBEF_nem17 <- subset(dBEF_nem, year==2017)

#m17.11


#### 2 - 2021's data####
dBEF_nem21 <- subset(dBEF_nem, year==2021)


#### 2.1 - cp1 ####
#examine the data
par(mfrow = c(1,1))
hist(dBEF_nem$cp1_per100g, breaks = seq(min(dBEF_nem$cp1_per100g), max(dBEF_nem$cp1_per100g), length.out=30))
#how many samples have zero cp1 nematodes: 
sum(dBEF_nem21$cp1_per100g == 0) #116 of 240

ggplot(dBEF_nem21, aes(x = log(sowndiv), y = cp1_per100g))+
  geom_point()

#m21.11
m21.cp1.0 <- brm(cp1_per100g ~ sowndiv * treatment, 
                 data = dBEF_nem21, family = "lognormal",
                 chains = 3,
                 cores = 3,
                 iter = 2000, warmup = 1000, 
                 control = list(adapt_delta =0.99))

#### 2.2 - cp2 ####
hist(dBEF_nem21$cp2_per100g, breaks = seq(min(dBEF_nem21$cp2_per100g), max(dBEF_nem21$cp2_per100g), length.out=30))
sum(dBEF_nem$cp2_per100g == 0) # 1 sample with zero cp2 nematodes


#### 2.3 - cp3 #### 
hist(dBEF_nem21$cp3_per100g, breaks = seq(min(dBEF_nem21$cp3_per100g), max(dBEF_nem21$cp3_per100g), length.out=30))
sum(dBEF_nem21$cp3_per100g == 0) # 5 sample2 with zero cp3 nematodes

#### 2.4 - cp4 ####
hist(dBEF_nem21$cp4_per100g, breaks = seq(min(dBEF_nem21$cp4_per100g), max(dBEF_nem21$cp4_per100g), length.out=30))
sum(dBEF_nem$cp4_per100g == 0) # 19 samples with zero cp4 nematodes

#### 2.5 - cp5 ####
hist(dBEF_nem21$cp5_per100g, breaks = seq(min(dBEF_nem21$cp5_per100g), max(dBEF_nem21$cp5_per100g), length.out=30))
sum(dBEF_nem21$cp5_per100g == 0) # 213 sample with zero cp2 nematodes



