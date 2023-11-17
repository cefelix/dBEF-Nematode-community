
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
SEED = 19111996

#### exploration ####
  dBEF_nemSH1 <- subset(dBEF_nem, SH == 1)
  dBEF_nemSH5 <- subset(dBEF_nem, SH == 5)
  dBEF_nemSH15 <- subset(dBEF_nem, SH == 15)
  dBEF_nemSH19 <- subset(dBEF_nem, SH == 19)
  
p.1 <- ggplot(dBEF_nemSH1, aes(y = Pl_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  labs(title = "SH 1")+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")

p.5 <- ggplot(dBEF_nemSH5, aes(y = Pl_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")+
  labs(title = "SH 5")

p.15 <- ggplot(dBEF_nemSH15, aes(y = Pl_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 15")+
  geom_smooth(method="lm")

p.19 <- ggplot(dBEF_nemSH19, aes(y = Pl_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 19")+
  geom_smooth(method="lm")

grid.arrange(p.1, p.5, p.15, p.19)
rm(p.1, p.5, p.15, p.19)


#### 11 Pl_per100gLog ~ sowndivLogStd * treatment + (1|block/plot), fam=gaussian ####

m.Pl.11 <- brm(Pl_per100gLog ~ sowndivLogStd*treatment + (1|block/plot),
                   data = dBEF_nem21, family = "gaussian",
                   seed = SEED,
                   chains = 3,
                   cores = 3,
                   iter = 2000, warmup = 1000,
                   control = list(adapt_delta=0.9)) #7 divergent transitions

m.Pl.12 <- update(m.Pl.11,
                  control=list(adapt_delta=0.99))

pp_check(m.Pl.12, ndraws = 100)


#### 21 Pl_per100gZeroC ~ sowndiv * treatment + (1|block), fam=lognormal ####

m.Pl.21 <- brm(Pl_per100gZeroC ~ sowndiv*treatment + (1|block),
               data = dBEF_nem21, family = "lognormal",
               seed = SEED,
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               control = list(adapt_delta=0.9)) #2 divergent transitions 

m.Pl.22 <- update(m.Pl.21,
                  control = list(adapt_delta=0.99))

pp_check(m.Pl.22, ndraw=100)


#### 11a hurdle: Pl_per100g ~ sowndivLogStd * treatment + (1|block/plot), fam=hurdle_lognormal ####
m.Pl.hurdle11a <- brm(
  bf(Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pl.hurdle31, ndraws=100)+
  xlim(0,3000)

summary(m.Pl.hurdle11a)


#### saving models ####
save(m.Pl.11, m.Pl.12,
     m.Pl.21, m.Pl.22,
     m.Pl.hurdle11a, #Pl ~ sowndivLogStd, fam=hurdle_lognormal
     file="./statistics/brms/231117_Pl.RData")

load(file="./statistics/brms/231107_Pl_hurdle.RData")

conditional_effects(m.Pl.hurdle21)

