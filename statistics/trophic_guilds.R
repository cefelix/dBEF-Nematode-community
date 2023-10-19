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

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
dBEF_nem %>%
  str()
dBEF_nem$treatment <- as.factor(dBEF_nem$treatment)
dBEF_nem <- dBEF_nem %>%
  mutate(col.sowndiv = as.factor(sowndiv), .after=sowndiv)

#### 0 - all data #### 
#inspect data
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


####1 - 2017's data####
dBEF_nem17 <- subset(dBEF_nem, year==2017)

#m17.11


#### 2 - 2021's data####
dBEF_nem21 <- subset(dBEF_nem, year==2021)
  #one subset for each treatment (to plot) 
  dBEF_nem21_t1 <- subset(dBEF_nem21, treatment == 1) 
  dBEF_nem21_t2 <- subset(dBEF_nem21, treatment == 2)
  dBEF_nem21_t3 <- subset(dBEF_nem21, treatment == 1)
  
#### 2.01 - are 5 trophic group density variables independent from each other? ####  
corr21 <- dBEF_nem21 %>% 
  select(c(Fu_per100g, Ba_per100g, Pl_per100g, Pr_per100g, Om_per100g))
ggpairs(corr21) #Fu and Ba seem to be correlated
rm(corr21)

  
#### 2.11 - Fungivores exploration ####   
  #examine the data 
  hist(dBEF_nem21$Fu_per100g, breaks = seq(min(dBEF_nem21$Fu_per100g), 
                                           max(dBEF_nem21$Fu_per100g), 
                                           length.out=30))
  ggplot(dBEF_nem21, aes(x=Fu_per100g))+
    geom_density()

#how many samples have zero bacterivores:
sum(dBEF_nem21$Fu_per100g == 0) #4 of 240

#a jitter plot with an OLS regression line
p.all <- ggplot(dBEF_nem21, aes(x = log(sowndiv), y = Fu_per100g))+
          geom_jitter()+
          labs(title = "all treatments")

p.t1 <- ggplot(dBEF_nem21_t1, aes(x = log(sowndiv), y = Fu_per100g))+
          geom_jitter()+
          geom_smooth(method="lm")+
          labs(title = "-SH -PH")

p.t2 <- ggplot(dBEF_nem21_t2, aes(x = log(sowndiv), y = Fu_per100g))+
          geom_jitter()+
          geom_smooth(method="lm")+
          labs(title = "+SH -PH")

p.t3 <- ggplot(dBEF_nem21_t2, aes(x = log(sowndiv), y = Fu_per100g))+
          geom_jitter()+
          geom_smooth(method="lm")+
          labs(title = "+SH +PH")

grid.arrange(p.all, p.t1, p.t2, p.t3)
  rm(dBEF_nem21_t1, dBEF_nem21_t2, dBEF_nem21_t3,
     p.all, p.t1, p.t2, p.t3)

#### 2.12 - Fungivores modelling ####
  # m.Fu.11:  density ~ sowndiv * treatment + (1|block), fam=lognormal, data = noZeros  
  # m.Fu.11b: log(density) ~ sowndiv * treatment + (1|block), fam=gaussian , data = noZeros  
  # m.Fu.12:  log(density) ~ sowndiv * treatment + (1|block), fam=gaussian , data = all
  # m.Fu.13:  density ~ sowndiv * treatment, fam= hurdle_lognormal, data = 
  
#log transforming Fu densities, if bigger than zero:  
dBEF_nem21 <- dBEF_nem21 %>% 
  mutate(Fu_per100gLog = ifelse(Fu_per100g == 0, 0, log(Fu_per100g) ) , 
         .after = Fu_per100g)
  #a df excluding zeros:
  dBEF_nem21_noZero <- subset(dBEF_nem21, Fu_per100g != 0)
  #load the existing models:
  load("./statistics/brms/231019_TrophicGuilds.RData")
  

# first, a log transformed model excluding all zero values:
m.Fu.11 <- brm(Fu_per100g ~ sowndiv*treatment + (1|block),
               data = dBEF_nem21_noZero, family = "lognormal",
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               control = list(adapt_delta=0.99))
summary(m.Fu.11)
  pp_check(m.Fu.11, draws=100)      #model underpredicts
  mcmc_plot(m.Fu.11, type="combo")  #
  mcmc_plot(m.Fu.11, type="violin")
  
  mcmc_plot(m.Fu.11, type="neff") #at least >10% of posterior samples, better >50%
  mcmc_plot(m.Fu.11, type="trace") #should bounce around a value randomly
  
#and the same, but using Log-transformed fu density and family = gaussian  
m.Fu.11b <- update(m.Fu.11, newdata = dBEF_nem21_noZero,
                   Fu_per100gLog ~ sowndiv*treatment + (1|block),
                   family="gaussian")
summary(m.Fu.11b)
pp_check(m.Fu.11b) #model underestimates, except for very high / very low values
mcmc_plot(m.Fu.11b, type="pairs",
          off_diag_fun="hex",
          diag_fun="dens")
  

# second, a log transformed model with all values: 
m.Fu.12a <- brm(Fu_per100gLog ~ sowndiv*treatment + (1|block),
               data = dBEF_nem21, family = "gaussian",
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               control = list(adapt_delta=0.9))
summary(m.Fu.12)
  plot(m.Fu.12)
  pp_check(m.Fu.12) #underestimating the amount of zeros
                    #shifting the peak to the left
                    #overestimating right tail frequencies
  mcmc_plot(m.Fu.12, type="neff") #bigger than 0.1
  mcmc_plot(m.Fu.12a, type = "pairs",
            diag_fun = "dens",
            off_diag_fun = "hex",
            fixed = TRUE)
  mcmc_intervals(m.Fu.12)
  
#third, a log transformed hurdle model (without random factor)
m.Fu.13 <- brm(Fu_per100g ~ sowndiv*treatment,
             data = dBEF_nem21, family = "hurdle_lognormal",
             chains = 3,
             cores = 3,
             iter = 2000, warmup = 1000,
             control = list(adapt_delta=0.99))
summary(m.Fu.13)
pp_check(m.Fu.13)
mcmc_plot(m.Fu.13, type = "pairs",
          off_diag_fun="hex",
          diag_fun="dens")

#and to compare with Random factors: 
m.Fu.13b <- update(m.Fu.13, Fu_per100g ~ sowndiv*treatment + (1|block), newdata = dBEF_nem21)
#pp_check(m.Fu.13b)

#save models

save(m.Fu.11, m.Fu.11b,
     m.Fu.12, 
     m.Fu.13, m.Fu.13b,
     file = "./statistics/brms/231019_TrophicGuilds.RData")

####2.21 - Bacterivores 

