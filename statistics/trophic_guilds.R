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

####1 - 2017's data####
dBEF_nem17 <- subset(dBEF_nem, year==2017)

#m17.11


#### 2 - 2021's data####
dBEF_nem21 <- subset(dBEF_nem, year==2021)
  #one subset for each treatment (to plot) 
  dBEF_nem21_t1 <- subset(dBEF_nem21, treatment == 1) 
  dBEF_nem21_t2 <- subset(dBEF_nem21, treatment == 2)
  dBEF_nem21_t3 <- subset(dBEF_nem21, treatment == 1)
  
  
#### 2.01 - are trophic guild density variables independent from each other? ####  
corr21 <- dBEF_nem21 %>% 
  select(c(Fu_per100g, Ba_per100g, Pl_per100g, Pr_per100g, Om_per100g))
ggpairs(corr21) #Fu and Ba seem to be correlated
rm(corr21)

  
#### 2.11 - Fungivores exploration ####   
  #examine the data 
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
  
#exploration of random factors
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
  load("./statistics/brms/231101_TrophicGuilds.RData")
  
  

#### 2.12a - Fu_per100gLog ~ sowndiv*treatment + (1|block): ####
m.Fu.12 <- brm(Fu_per100gLog ~ sowndiv*treatment + (1|block),
               data = dBEF_nem21, family = "gaussian",
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               control = list(adapt_delta=0.9))
  
pp_check(m.Fu.12, ndraws=100) #overpredicting at low values 
                    #not accounting for zeros properly
  
#model plotting:
  #new data to create regression curve:
  #thats only block1
  nd <- tibble(sowndiv = seq(from = 1, 
                             to = 60, 
                             length.out = 30) %>% rep(., times = 3),
               treatment = rep(1:3, each = 30),
               block = rep("B4", 90))
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
  pp_check(m.Fu.12b., ndraws=100) #thats a bad fit, lets remove block:
  
  
  


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




####2.21 - Bacterivores ####
  #explore:
  hist(dBEF_nem$Ba_per100gLog)

#### Ba_per100gLog ~ sowndiv*treatment + (1|block) ####
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


#### 2.3 - Herbivores ####
#### 2.31 - Pl_per100gLog ~ sowndiv*treatment + (1|block) ####
m.Pl.31 <- brm(Pl_per100gLog ~ sowndiv*treatment + (1|block),
                  data = dBEF_nem21, family = "gaussian",
                  chains = 3,
                  cores = 3,
                  iter = 2000, warmup = 1000,
                  control = list(adapt_delta=0.9))
pp_check(m.Pl.31, ndraws = 100)


pr_uniform = prior(uniform(-100, 100), lb=-100, ub=100, class = "b")
m.Pl.31.unifPrior <- brm(Pl_per100gLog ~ sowndiv*treatment + (1|block),
                        data = dBEF_nem21, 
                        iter = 1000,
                        sample_prior = "only",
                        prior=pr_uniform)
pp_check(m.Pl.31.unifPrior, ndraws = 100)+
  xlim(-10, 20)

pr_gaussian = prior(normal(0, 5), class="b")
m.Pl.31.unifPrior <- brm(Pl_per100gLog ~ sowndiv*treatment + (1|block),
                         data = dBEF_nem21, 
                         iter = 1000,
                         sample_prior = "only",
                         prior=pr_gaussian)
pp_check(m.Pl.31.unifPrior, ndraws = 100)+
  xlim(-10, 20)


#lets set the default priors manually:
get_prior(Pl_per100gLog ~ sowndiv*treatment + (1|block),
          data = dBEF_nem21, family = "gaussian")
          #student_t(3, 0, 2.5) means a student t distribution with:
          #3 degrees of freedom, center of 0, range of -2.5 to +2.5
priors <- c(
  prior(uniform(-100, 100), class="b"),
  prior("student_t(3, 0, 2.5)", class = "sd"),
  prior("student_t(3, 0, 2.5)", class = "sigma"),
  prior("student_t(3, 5.6, 2.5)", class = "intercept")
) 

curve(dnorm(x, 0, 1), from = -5, to=5)
curve(dt(x, df=3), from=-5, to=5, col="red", add = TRUE)
curve(dt(x, df=10), from=-5, to=5, col="green", add = TRUE)
curve(dnorm(x, 0, 0.5), from = -5, to=5, col="pink", add=TRUE)

#lets check 
prior

priors <- 

#check out priors:
library(pubh)
rnorm(1e4, mean=0, sd=4) %>% density() %>% plot()
rnorm(1e4, mean=0, sd=4) %>% inv_logit() %>% density() %>% plot()


#### 2.32 - Pl_per100g ~ sowndiv * treatment + (1|block), fam=lognormal ####
#add small constant to zeros:
dBEF_nem21 <- dBEF_nem21 %>% 
  mutate(Pl_per100g  = ifelse(Pl_per100g == 0, 0.01, Pl_per100g ))

m.Pl.32 <- brm(Pl_per100g ~ sowndiv*treatment + (1|block),
               data = dBEF_nem21, family = "lognormal",
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               control = list(adapt_delta=0.9)) #2 divergent transitions -->
#increase delta:
m.Pl.33 <- update(m.Pl.32,
                  control = list(adapt_delta=0.99))
pp_check(m.Pl.33, ndraw=100)





####outdated####
#excluding zeros:
dBEF_nem21_noZero <- subset(dBEF_nem21, Fu_per100g > 0)
# a log transformed model excluding all zero values:
m.Fu.11 <- brm(Fu_per100g ~ sowndiv*treatment + (1|block),
               data = dBEF_nem21_noZero, family = "lognormal",
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               control = list(adapt_delta=0.99))
summary(m.Fu.11)
pp_check(m.Fu.11, ndraws=100)+
  xlim(0, 3000) #model assigns higher probability to right tail
pp_check(m.Fu.11, ndraws=100)+
  xlim(0, 1000) #and much higher probability at the peak

#and the same, but using Log-transformed fu density and family = gaussian  
m.Fu.11b <- update(m.Fu.11, newdata = dBEF_nem21_noZero,
                   Fu_per100gLog ~ sowndiv*treatment + (1|block),
                   family="gaussian")
summary(m.Fu.11b)
pp_check(m.Fu.11b) #model underestimates, except for very high / very low values
mcmc_plot(m.Fu.11b, type="pairs",
          off_diag_fun="hex",
          diag_fun="dens")

##---

#2.12a2 - Fu_per100gLog ~ sowndiv*treatment + SWC_gravimetric #
m.Fu.12a2 <- brm(Fu_per100gLog ~ sowndiv*treatment + SWC_gravimetric,
                 data = dBEF_nem21, family = "gaussian",
                 chains = 3,
                 cores = 3,
                 iter = 2000, warmup = 1000,
                 control = list(adapt_delta=0.9))
summary(m.Fu.12a2)
pp_check(m.Fu.12a2, ndraws=100) #include block:
m.Fu.12a4 <- update(m.Fu.12a3,
                    control = list(adapt_delta = 0.99))

#---

# 2.12a5 Fu_per100gLog ~ SWC_gravimetric #
m.Fu.12a6 <- brm(Fu_per100gLog ~ SWC_gravimetric + (1|block),
                 data = dBEF_nem21, family = "gaussian",
                 chains = 3,
                 cores = 3,
                 iter = 2000, warmup = 1000,
                 control = list(adapt_delta=0.9))
m.Fu.12a7 <- update(m.Fu.12a6, 
                    control = list(adapt_delta=0.99))

pp_check(m.Fu.12a7, ndraws=100)

