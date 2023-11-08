####data and packages####
  library(brms)
  library(dplyr)
  library(tidyr)
  library(bayesplot)
  library(ggplot2)
  library(gridExtra)
  library(hexbin)
  library(GGally)

# a seed:
SEED = 22061996

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
  dBEF_nemSH1 <- subset(dBEF_nem, SH == 1)
  dBEF_nemSH5 <- subset(dBEF_nem, SH == 5)
  dBEF_nemSH15 <- subset(dBEF_nem, SH == 15)
  dBEF_nemSH19 <- subset(dBEF_nem, SH == 19)
  
  dBEF_nem21 <- subset(dBEF_nem, year==2021)


#### exploration ####
p.1 <- ggplot(dBEF_nemSH1, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  labs(title = "SH 1")+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")

p.5 <- ggplot(dBEF_nemSH5, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")+
  labs(title = "SH 5")

p.15 <- ggplot(dBEF_nemSH15, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 15")+
  geom_smooth(method="lm")

p.19 <- ggplot(dBEF_nemSH19, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 19")+
  geom_smooth(method="lm")

grid.arrange(p.1, p.5, p.15, p.19)
rm(p.1, p.5, p.15, p.19)


#### OUTDATED: Ba_per100gLog ~ sowndiv*treatment + (1|block) ####
m.Ba.21 <- brm(Ba_per100gLog ~ sowndiv*treatment + (1|block),
               data = dBEF_nem21, family = "gaussian",
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               control = list(adapt_delta=0.9))
pp_check(m.Ba.21, ndraws=21) #9 divergent transitions, ESS too low
m.Ba.21b <- update(m.Ba.21, 
                   control = list(adapt_delta=0.99))
pp_check(m.Ba.21b, ndraws=100) #bad fit, lets standardize


#### 11 hurdle: Ba_per100gLog ~ sowndivLog*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Ba.hurdle11 <- brm(
  bf(Ba_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Ba.hurdle11, ndraws=100)

#### 21 hurdle: Ba_per100g ~ sowndivLog*treatment + (1|block/plot), family = hurdle_lognormal ####
m.Ba.hurdle21 <- brm(
  bf(Ba_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Ba.hurdle21, ndraws=100)+
  xlim(0,1000)


#### 31 hurdle: Ba_per100gLog ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Ba.hurdle31 <- brm(
  bf(Ba_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Ba.hurdle31, ndraws=100)

#### 41 hurdle: Ba_per100gLog ~ sowndiv*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Ba.hurdle41 <- brm(
  bf(Ba_per100gLog.hurdle ~ sowndiv*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Ba.hurdle41, ndraws=100)
#### save hurdle models ####
save(m.Ba.hurdle11, #log(Ba) ~ sowndivLog, fam=hurdle_gaussian
     m.Ba.hurdle21, #Ba  ~ sowndivLog, fam=hurdle_lognormal
     m.Ba.hurdle31, #log(Ba) ~ sowndivLogStd, fam=hurdle_gaussian
     m.Ba.hurdle41, #log(Ba) ~ sowndiv, fam=hurdle_gaussian
     file="./statistics/brms/231108_Ba_hurdle.RData")

load(file="./statistics/brms/231107_Ba_hurdle.RData")
conditional_effects(m.Ba.hurdle11)
