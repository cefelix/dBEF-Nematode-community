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


#### exploration ####
  dBEF_nemSH1 <- subset(dBEF_nem, SH == 1)
  dBEF_nemSH5 <- subset(dBEF_nem, SH == 5)
  dBEF_nemSH15 <- subset(dBEF_nem, SH == 15)
  dBEF_nemSH19 <- subset(dBEF_nem, SH == 19)
  
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
rm(p.1, p.5, p.15, p.19,
   dBEF_nemSH19, dBEF_nemSH15, dBEF_nemSH5, dBEF_nemSH1)



#### 11a hurdle: Ba_per100g ~ sowndivLog*treatment + (1|block/plot), family = hurdle_lognormal ####
m.Ba.hurdle11a <- brm(
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

pp_check(m.Ba.hurdle11a, ndraws=100)+
  xlim(0,1000)


#### save hurdle models ####
save(m.Ba.hurdle21, #Ba  ~ sowndivLog, fam=hurdle_lognormal
     file="./statistics/brms/231117_Ba.RData")

load(file="./statistics/brms/231107_Ba_hurdle.RData")
conditional_effects(m.Ba.hurdle11)
