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

#add log-transformed densities for each trophic guild:
dBEF_nem <- dBEF_nem %>%
  #log transform Ba densities
  mutate(Ba_per100gLog = ifelse(Ba_per100g == 0, 0.001, Ba_per100g), 
         .after=Ba_per100g) %>%
  mutate(Ba_per100gLog = log(Ba_per100gLog)) %>%
  #log Fu 
  mutate(Fu_per100gLog = ifelse(Fu_per100g == 0, 0.001, Fu_per100g), 
         .after=Fu_per100g) %>%
  mutate(Fu_per100gLog = log(Fu_per100gLog)) %>%
  #log Pr
  mutate(Pr_per100gLog = ifelse(Pr_per100g == 0, 0.001, Pr_per100g), 
         .after=Pr_per100g) %>%
  mutate(Pr_per100gLog = log(Pr_per100gLog)) %>% 
  #log Pl
  mutate(Pl_per100gLog = ifelse(Pl_per100g == 0, 0.001, Pl_per100g), 
       .after=Pl_per100g) %>%
  mutate(Pl_per100gLog = log(Pl_per100gLog)) %>% 
  #log Om
  mutate(Om_per100gLog = ifelse(Om_per100g == 0, 0.001, Om_per100g), 
         .after=Pl_per100g) %>%
  mutate(Om_per100gLog = log(Om_per100gLog))


  
  





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

#### 0.11 - all data, log(Fu) modelling####
hist(dBEF_nem$Fu_per100gLog)
m.0.Fu11 <- brm(Fu_per100gLog ~ sowndiv*treatment + year + (1|block),
                     data = dBEF_nem, family = "gaussian",
                     chains = 3,
                     cores = 3,
                     iter = 2000, warmup = 1000,
                     control = list(adapt_delta=0.99)) 
pp_check(m.0.Fu11, ndraws=100)

#lets make a yearblock variable:
dBEF_nem <- dBEF_nem %>%
  mutate(yearblock = as.factor(paste(dBEF_nem$year, dBEF_nem$block, sep="")),
         .after = block)

m.0.Fu12 <- brm(Fu_per100gLog ~ sowndiv*SH + (1|yearblock),
                data = dBEF_nem, family = "gaussian",
                chains = 3,
                cores = 3,
                iter = 2000, warmup = 1000,
                control = list(adapt_delta=0.99)) 
pp_check(m.0.Fu12, ndraws=100)

m.0.Fu13 <- brm(Fu_per100gLog ~ sowndiv*SH*PH + (1|yearblock),
                data = dBEF_nem, family = "gaussian",
                chains = 3,
                cores = 3,
                iter = 2000, warmup = 1000,
                control = list(adapt_delta=0.99)) 
pp_check(m.0.Fu13, ndraws=100) # 1396 transitions that exceeded max_treedepth

m.0.Fu14 <- update(m.0.Fu13,
                   control=list(max_treedepth = 12))
pp_check(m.0.Fu14, ndraw=100)

m.0.Fu15 <- brm(Fu_per100gLog ~ sowndiv*SH*PH + SWC_gravimetric + (1|yearblock),
                   data = dBEF_nem, family = "gaussian",
                   chains = 3,
                   cores = 3,
                   iter = 2000, warmup = 1000,
                   control = list(adapt_delta=0.99)) 
pp_check(m.0.Fu15, ndraws=100) #1691 transitions exceeding max_treedepth

m.0.Fu16 <- update(m.0.Fu15,
                   control=list(max_treedepth=100))
pp_check(m.0.Fu16, ndraw=100)


#save the models above:
save(m.0.Fu11,
     m.0.Fu12,
     m.0.Fu13,
     m.0.Fu14,
     m.0.Fu15,
     m.0.Fu16,
     file = "./statistics/brms/231101_Fu_allData.RData")


#### 0.21 - all data, standardized Fu modelling ####

Fu_scaled <- scale(dBEF_nem$Fu_per100g)[,1]
dBEF_nem <- dBEF_nem %>%
  mutate(Fu_per100gStd = Fu_scaled, 
         .after =  Fu_per100g)

a <- dBEF_nem$Fu_per100g %>%
  scale() 
a[,1] %>% str()

dBEF_nem %>%
  pull(Fu_per100gStd) %>%
  density() %>%
  plot()

#%>%
  unlist %>%
  as.numeric() %>%

dBEF_nem$Fu_per100g








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
dBEF_nem21_noZero <- subset(dBEF_nem21, Fu_per100g != 0)

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
  
#### 2.11b - Fungivores random factors #### 
  #soil dry weight:
  ggplot(dBEF_nem21, aes(x=soilDW, y = Fu_per100g, col = block))+
    geom_point()
  #gravimetric water content
  ggplot(dBEF_nem21, aes(x=SWC_gravimetric, y = Fu_per100g, col = block))+
    geom_point() #the B4 sample to the very left is B4A13D2
  #gravimetric water content against DW
  ggplot(dBEF_nem21, aes(x=SWC_gravimetric, y = soilDW, col = block))+
    geom_point()
  
  

#### 2.12 - Fungivores modelling ####
  # m.Fu.11:  density ~ sowndiv * treatment + (1|block), fam=lognormal, data = noZeros  
  # m.Fu.11b: log(density) ~ sowndiv * treatment + (1|block), fam=gaussian , data = noZeros  
  # m.Fu.12a: log(density) ~ sowndiv * treatment + (1|block), fam=gaussian , data = all
  # m.Fu.12b: log(density) ~ sowndiv * treatment + soilWC + (1|block), fam=gaussian, data = all
  # m.Fu.12c: log(density) ~ sowndiv * treatment, fam=gaussian , data = all
  # m.Fu.13:  density ~ sowndiv * treatment, fam= hurdle_lognormal, data = 
  
  #the already fit models:
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
  

#### 2.12a - Fu_per100gLog ~ sowndiv*treatment + (1|block): ####
m.Fu.12 <- brm(Fu_per100gLog ~ sowndiv*treatment + (1|block),
               data = dBEF_nem21, family = "gaussian",
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               control = list(adapt_delta=0.9))
summary(m.Fu.12)
#model diagnostics:  
  plot(m.Fu.12)
  pp_check(m.Fu.12, ndraws=100) #overpredicting at very low values 
                    #not accounting for zeros properly
                    #
  mcmc_plot(m.Fu.12, type="neff") #bigger than 0.1
  mcmc_plot(m.Fu.12, type = "pairs",
            diag_fun = "dens",
            off_diag_fun = "hex",
            fixed = TRUE)
  mcmc_intervals(m.Fu.12)
  
  #model plotting:
  #new data to create regression curve:
  #thats only block1
  nd <- tibble(sowndiv = seq(from = 1, 
                             to = 60, 
                             length.out = 30) %>% rep(., times = 3),
               treatment = rep(1:3, each = 30),
               block = rep("B1", 90))
  #fitted values
  f <-
    fitted(m.Fu.12, newdata = nd) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(treatment = as.factor(treatment)) # treatment has to a factor to be plotted
  
  #plotting fitted stuff
  ggplot(dBEF_nem21, aes(x=log(sowndiv)))+
    geom_jitter(aes(y = Fu_per100g, color = treatment), width = 0.125, alpha=0.5)+
    geom_smooth(data = f,
                aes(y = exp(Estimate),  #exp estimate to not plot logged data
                    color = treatment),
                stat = "identity")+
    scale_x_continuous()+
    ylab("Fungivores per 100g DW")+
    theme_classic()
  
  
#### 2.12a2 - Fu_per100gLog ~ sowndiv*treatment + SWC_gravimetric ####  
  m.Fu.12a2 <- brm(Fu_per100gLog ~ sowndiv*treatment + SWC_gravimetric,
                 data = dBEF_nem21, family = "gaussian",
                 chains = 3,
                 cores = 3,
                 iter = 2000, warmup = 1000,
                 control = list(adapt_delta=0.9))  
  pp_check(m.Fu.12a2, ndraws=100) #include block:
  m.Fu.12a4 <- update(m.Fu.12a3,
                      control = list(adapt_delta = 0.99))
  
#### 2.12a5 Fu_per100gLog ~ SWC_gravimetric ####  
  m.Fu.12a6 <- brm(Fu_per100gLog ~ SWC_gravimetric + (1|block),
                   data = dBEF_nem21, family = "gaussian",
                   chains = 3,
                   cores = 3,
                   iter = 2000, warmup = 1000,
                   control = list(adapt_delta=0.9))
  m.Fu.12a7 <- update(m.Fu.12a6, 
                      control = list(adapt_delta=0.99))
  
  pp_check(m.Fu.12a7, ndraws=100)
  
#### 2.12b - Fu_per100gLog ~ sowndiv*treatment + (Fu_per100gLog|block): ####
  m.Fu.12b1 <- brm(Fu_per100gLog ~ sowndiv*treatment + (Fu_per100gLog|block),
                 data = dBEF_nem21, family = "gaussian",
                 chains = 3,
                 cores = 3,
                 iter = 4000, warmup = 1000,
                 control = list(adapt_delta=0.99))  
  m.Fu.12b2 <- update(m.Fu.12b1, control=list(max_treedepth=12))
  #too little data?
  pp_check(m.Fu.12b1, ndraws=100)
  summary(m.Fu.12b1)

#### 2.12b - Fu_per100gLog ~ sowndiv*treatment+SWC+(1|block) ####
  m.Fu.12b <- brm(Fu_per100gLog ~ sowndiv * treatment + SWC_gravimetric + (1|block),
                  data = dBEF_nem21, family = "gaussian",
                  chains = 3,
                  cores = 3,
                  iter = 2000, warmup = 1000,
                  control = list(adapt_delta=0.9))
  #6 divergent transitions, increase delta:
  m.Fu.12b. <- update(m.Fu.12b, control = list(adapt_delta = 0.99))
  
  summary(m.Fu.12b.)
  pp_check(m.Fu.12b.) #thats a bad fit, lets remove block:
  
  m.Fu.12b2 <- brm(Fu_per100gLog ~ sowndiv * treatment + SWC_gravimetric,
                  data = dBEF_nem21, family = "gaussian",
                  chains = 3,
                  cores = 3,
                  iter = 2000, warmup = 1000,
                  control = list(adapt_delta=0.9))
  summary(m.Fu.12b2)
  pp_check(m.Fu.12b2)
  
  
  
  

#### 2.12c - log transformed model with all values, but without RE block: #### 
#plotting it without block:
  m.Fu.12c <- brm(Fu_per100gLog ~ sowndiv*treatment,
                 data = dBEF_nem21, family = "gaussian",
                 chains = 3,
                 cores = 3,
                 iter = 2000, warmup = 1000,
                 control = list(adapt_delta=0.9))
  
  summary(m.Fu.12c)
  pp_check(m.Fu.12c, ndraws=100) 
  #not accounting for block makes the 
  #model concentrate probability mass in the tails tails
  
#model plotting:
  #new data to create regression curve:
  nd <- tibble(sowndiv = seq(from = 1, 
                             to = 60, 
                             length.out = 30) %>% rep(., times = 3),
               treatment = rep(1:3, each = 30),
               )
  #fitted values
  f <-
    fitted(m.Fu.12c, newdata = nd) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(treatment = as.factor(treatment)) # treatment has to a factor to be plotted
  
  #plotting fitted stuff
  ggplot(dBEF_nem21, aes(x=log(sowndiv)))+
    geom_jitter(aes(y = Fu_per100g, color = treatment), width = 0.125, alpha = 0.6)+
    geom_smooth(data = f,
                aes(y = exp(Estimate),  #exp estimate to not plot logged data
                    color = treatment),
                stat = "identity")+
    scale_x_continuous()+
    ylab("Fungivores per 100g DW")+
    theme_classic()
  
 
  f$sowndiv <- f$sowndiv %>% as.numeric()
  f$treatment %>% str()
  
  
  
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

#### 2.13 model prediction accuracy comparison using loo ####
m.Fu.11 <- add_criterion(m.Fu.11, c("loo", "waic"))
m.Fu.11b <- add_criterion(m.Fu.11b, c("loo", "waic"))
m.Fu.12a <- add_criterion(m.Fu.12a, c("loo", "waic"), moment_match = TRUE)
m.Fu.13 <- add_criterion(m.Fu.13, c("loo", "waic"))
m.Fu.13b <- add_criterion(m.Fu.13b, c("loo", "waic"))

waic(m.Fu.11)
waic(m.Fu.11b)
waic(m.Fu.12)
waic(m.Fu.13)
waic(m.Fu.13b)

loo(m.Fu.11)
loo(m.Fu.11b)
loo(m.Fu.12a)

#all models:
model_weights(m.Fu.11, m.Fu.11b, m.Fu.12, m.Fu.13, m.Fu.13b,
              weights = "loo") %>%
  round(digits = 3)

l <- loo_compare(m.Fu.11, m.Fu.11b, m.Fu.12, m.Fu.13, m.Fu.13b, 
            criterion = "loo") 
print(l, simplify = F)

#no zero data:
model_weights(m.Fu.11, m.Fu.11b,
              weights = "loo") %>%
  round(digits = 3)


#save models

save(m.Fu.11, m.Fu.11b,
     m.Fu.12, 
     m.Fu.12a2, m.Fu.12a3, m.Fu.12a4, m.Fu.12a5, m.Fu.12a6, m.Fu.12a7,
     m.Fu.12b, m.Fu.12c, m.Fu.12b., m.Fu.12b1, m.Fu.12b2,
     m.Fu.13, m.Fu.13b,
     file = "./statistics/brms/231101_TrophicGuilds.RData")

####2.21 - Bacterivores 

