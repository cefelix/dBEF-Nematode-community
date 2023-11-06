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


#### OLD: Fu_per100gLog ~ sowndiv*treatment + (1|block) AND ~sowndiv*treatment + (1|block/plot) ####
m.Fu.12 <- brm(Fu_per100gLog ~ sowndiv*treatment + (1|block),
               data = dBEF_nem21, family = "gaussian",
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               control = list(adapt_delta=0.9))
  
pp_check(m.Fu.12, ndraws=100) #overpredicting at low values 
                    #not accounting for zeros properly
                    #has a good fit, but ignores plot nested in block

# ...
#lets NEST plot in block:
# ...
m.Fu.12_nestRE <- brm(Fu_per100gLog ~ sowndiv*treatment + (1|block/plot),
                         data = dBEF_nem21, family = "gaussian",
                         chains = 3,
                         cores = 3,
                         iter = 2000, warmup = 1000,
                         control = list(adapt_delta=0.9))
#2 divergent transitions:
m.Fu.12_nestRE_b <- update(m.Fu.12_nestRE,
                           control = list(adapt_delta=0.99))

pp_check(m.Fu.12_nestRE_b, ndraws = 100) #fit is worse than without nesting plot in block

save(m.Fu.12,
     m.Fu.12_nestRE, m.Fu.12_nestRE_b,
     file = "./statistics/brms/231103_Fu_outdated.RData")

  
#### 2.12b - Fu_per100gLog ~ sowndiv*treatment + (Fu_per100gLog|block): ####
  m.Fu.12b1 <- brm(Fu_per100gLog ~ sowndiv*treatment + (Fu_per100gLog|block),
                 data = dBEF_nem21, family = "gaussian",
                 chains = 3,
                 cores = 3,
                 iter = 4000, warmup = 1000,
                 control = list(adapt_delta=0.99))  
  m.Fu.12b2 <- update(m.Fu.12b1, control=list(max_treedepth=12))
  #too little data?
  pp_check(m.Fu.12b2, ndraws=100)
  summary(m.Fu.12b2)

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
  pp_check(m.Fu.12b., ndraws=100) #thats a bad fit
  
  
  


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

#### hurdle11: Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_lognormal ####
m.Fu.hurdle11 <- brm(
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

m.Fu.hurdle12 <- update(m.Fu.hurdle11,
                        control = list(adapt_delta=0.999))





#### hurdle21: Fu_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_lognormal####
#WARNING: this is a WRONG family for our data, as we already log-transformed the positive 
  #Fu densities prior to the model fitting!


dBEF_nem21$Fu_per100gLog.hurdle %>% density() %>% plot()


m.Fu.hurdle21 <- brm(
  bf(Fu_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = "hurdle_lognormal",
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Fu.hurdle21, ndraws = 100) 
  #okayish, quite a lot of probability mass in the right tail though!

#### hurdle31: Fu_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_gaussian####

m.Fu.hurdle31 <- brm(
  bf(Fu_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Fu.hurdle31, ndraws=100)


####saving Fu hurdle models####
save(m.Fu.hurdle11,
     m.Fu.hurdle21,
     m.Fu.hurdle31,
     file="./statistics/brms/231103_Fu_hurdle.RData")





