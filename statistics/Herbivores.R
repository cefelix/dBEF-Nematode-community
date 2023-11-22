
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





#### 11a hurdle: Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~1, family = hurdle_lognormal ####
SEED = 22061996
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

pp_check(m.Pl.hurdle11a, ndraws=100)+
  xlim(0,1000)

#### 11b hurdle: Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~term, family = hurdle_lognormal ####
SEED = 22061996
m.Pl.hurdle11b <- brm(
  bf(Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dBEF_nem21, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
  #1 divergent trasnitions, 1000+ exceeded max_treedepth
  #rhat=1.24, bulk EES too low

m.Pl.hurdle12b <- update(m.Pl.hurdle11b, 
                         iter = 3000, warmup = 1500,
                         control=list(adapt_delta=0.999,
                                      max_treedepth=12))
#54 divergent transitions, bulk EES too low


#### 31a hurdle: Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~1, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 22061996
dat <- subset(dBEF_nem21, sowndiv != 60)

m.Pl.hurdle31a <- brm(
  bf(Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pl.hurdle11a, ndraws=100)+
  xlim(0,1000)


####31a_Pl 2-5 hurdle: Pl2-5 ~ sowndivLogStd*treatment + (1|block/plot), hu~1, family = hurdle_lognormal, no 60 sp plots ####
SEED = 22061996
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60)

sum(dat$Pl2!=0) #215 of 228
sum(dat$Pl3!=0) #218
sum(dat$Pl4!=0) #118
sum(dat$Pl5!=0) #4

m.Pl2.hurdle31a <- brm(
  bf(Pl2 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good


m.Pl3.hurdle31a <- brm(
  bf(Pl3 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good


m.Pl4.hurdle31a <- brm(
  bf(Pl4 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #1 divergent transition

m.Pl4.hurdle32a <- update(m.Pl4.hurdle31a,
                          control=list(adapt_delta=0.999)) #all good

m.Pl5.hurdle31a <- brm(
  bf(Pl5 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) 
  #3000 transitions exceeded max_treedepth, bulk ESS too high, max Rhat=3.51
  #no point with 4 samples != 0

#### 31b hurdle: Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~term, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 22061996
dat <- subset(dBEF_nem21, sowndiv != 60)

m.Pl.hurdle31b <- brm(
  bf(Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
  #22 divergent transitions, bulk EES too low

m.Pl.hurdle32b <- update(m.Pl.hurdle31b,
                         iter = 3000, warmup = 1500,
                         control = list(adapt_delta=0.999))
#18 divergent transitions, 3800+ exceeded max_treedepth


####31b_Pl 2-5 hurdle: Pl2-5 ~ sowndivLogStd*treatment + (1|block/plot), hu~term, family = hurdle_lognormal, no 60 sp plots ####
SEED = 22061996
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60)

sum(dat$Pl2!=0) #215 of 228
sum(dat$Pl3!=0) #218
sum(dat$Pl4!=0) #118
sum(dat$Pl5!=0) #4

m.Pl2.hurdle31b <- brm(
  bf(Pl2 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good


m.Pl3.hurdle31b <- brm(
  bf(Pl3 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good


m.Pl4.hurdle31b <- brm(
  bf(Pl4 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good

m.Pl5.hurdle31b <- brm(
  bf(Pl5 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))
  #3000 transitions exceeded max_treedepth, max Rhat=3.69, bulk/tail ESS too low


#### saving models ####
save(m.Pl.hurdle11a, 
     m.Pl.hurdle11b, m.Pl.hurdle12b,
     m.Pl.hurdle31a,
       m.Pl2.hurdle31a,
       m.Pl3.hurdle31a,
       m.Pl4.hurdle31a, m.Pl4.hurdle32a,
       m.Pl5.hurdle31a, #pointless as only 4 samples !=0
     
     m.Pl.hurdle31b, m.Pl.hurdle32b, 
       m.Pl2.hurdle31b,
       m.Pl3.hurdle31b,
       m.Pl4.hurdle31b,
       m.Pl5.hurdle31b, #pointless as only 4 samples !=0
     file="./statistics/brms/231122_Pl.RData")

load(file="./statistics/brms/231122_Pl.RData")

conditional_effects(m.Pl.hurdle21)

