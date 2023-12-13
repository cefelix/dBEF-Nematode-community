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

#convert characters to factors:
  dBEF_nem$treatment <- as.factor(dBEF_nem$treatment)
  dBEF_nem$plot <- as.factor(dBEF_nem$plot)
  dBEF_nem <- dBEF_nem %>% # a colorcode of sowndiv for plotting
    mutate(col.sowndiv = as.factor(sowndiv), .after=sowndiv)

#adding "yearblock" which could be a random factor:
dBEF_nem <- dBEF_nem %>%
  mutate(yearblock = as.factor(paste(dBEF_nem$year, dBEF_nem$block, sep="")),
         .after = block)

#rename func.group to funcdiv for naming consistency:
dBEF_nem <- dBEF_nem %>%
  rename(funcdiv = func.group)

#### log transformation ####
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
         .after=Om_per100g) %>%
  mutate(Om_per100gLog = log(Om_per100gLog)) %>%
  #log ind_per100g
  mutate(ind_per100gLog = log(ind_per100g),
         .after = ind_per100g)

#log transformation for each cp-group:
dBEF_nem <- dBEF_nem %>% 
  #cp-1:
  mutate(cp1_per100gLog.hurdle = ifelse(cp1_per100g == 0, 
                                        0, log(cp1_per100g)),
         .after = cp1_per100g) %>%
  #cp-2:
  mutate(cp2_per100gLog.hurdle = ifelse(cp2_per100g == 0, 
                                        0, log(cp2_per100g)),
         .after = cp2_per100g) %>%
  #cp-3:
  mutate(cp3_per100gLog.hurdle = ifelse(cp3_per100g == 0, 
                                        0, log(cp3_per100g)),
         .after = cp3_per100g) %>%
  #cp-4:
  mutate(cp4_per100gLog.hurdle = ifelse(cp4_per100g == 0, 
                                        0, log(cp4_per100g)),
         .after = cp4_per100g) %>%
  #cp-5:
  mutate(cp5_per100gLog.hurdle = ifelse(cp5_per100g == 0, 
                                        0, log(cp5_per100g)),
         .after = cp5_per100g) 

  

#densities
dBEF_nem$Fu_per100gLog %>% density() %>% plot() #normal
dBEF_nem$Ba_per100gLog %>% density() %>% plot() #much probability on the right
dBEF_nem$Pl_per100gLog %>% density() %>% plot() #normal
dBEF_nem$Pr_per100gLog %>% density() %>% plot() #many zeros
dBEF_nem$Om_per100gLog %>% density() %>% plot() #way too many zeros, bimodal?
dBEF_nem$ind_per100gLog %>% density() %>% plot()

#densities non transformed:
dBEF_nem$Fu_per100g %>% density() %>% plot() #normal
dBEF_nem$Ba_per100g %>% density() %>% plot() #much probability on the right
dBEF_nem$Pl_per100g %>% density() %>% plot() #normal
dBEF_nem$Pr_per100g %>% density() %>% plot() #many zeros
dBEF_nem$Om_per100g %>% density() %>% plot() #way too many zeros, bimodal?
dBEF_nem$ind_per100g %>% density() %>% plot()

#### response standardization ####
#standardize each response density:
  scale_Ba <- (dBEF_nem$Ba_per100g %>% scale())[,1] 
  scale_Fu <- (dBEF_nem$Fu_per100g %>% scale())[,1] 
  scale_Pr <- (dBEF_nem$Pr_per100g %>% scale())[,1] 
  scale_Pl <- (dBEF_nem$Pl_per100g %>% scale())[,1] 
  scale_Om <- (dBEF_nem$Om_per100g %>% scale())[,1] 

#append to dBEF_nem:
dBEF_nem <- dBEF_nem %>%
  mutate(Ba_per100gStd = scale_Ba, 
         .after = Ba_per100g) %>%
  mutate(Fu_per100gStd = scale_Fu,
         .after = Fu_per100g) %>%
  mutate(Pl_per100gStd = scale_Pl, 
         .after = Pl_per100g) %>%
  mutate(Pr_per100gStd = scale_Pr, 
         .after = Pr_per100g) %>%
  mutate(Om_per100gStd = scale_Om, 
         .after = Om_per100g) 
rm(scale_Ba, scale_Fu,  scale_Pl, scale_Pr, scale_Om)

#density plots:
dBEF_nem$Fu_per100gStd %>% density() %>% plot()
dBEF_nem$Ba_per100gStd %>% density() %>% plot()
dBEF_nem$Pl_per100gStd %>% density() %>% plot()
dBEF_nem$Pr_per100gStd %>% density() %>% plot()
dBEF_nem$Om_per100gStd %>% density() %>% plot()

#histograms:
dBEF_nem$Ba_per100gStd %>% hist(breaks=30)

#pairs:
pairs(~Ba_per100g +  Fu_per100g + Pl_per100g + Om_per100g + Pr_per100g, data=dBEF_nem)
pairs(~Ba_per100g + log(sowndiv) +  SWC_gravimetric, data=dBEF_nem)

#### predictor standardization and log-transformation #####
dBEF_nem <- dBEF_nem %>%
  mutate(sowndivLog = log(sowndiv, base = 2),
         .after = sowndiv) %>%
  mutate(realdivLog = log(realdiv, base = 2),
         .after = realdiv)

dBEF_nem <- dBEF_nem %>%
  mutate(soilDW.Log = log(soilDW),
         .after = soilDW)


#### re-level treatment so that 1 = +SH+PH, 2 = +SH-PH, 3 = -SH-PH ####
for (i in 1:nrow(dBEF_nem)) {
  if(dBEF_nem$treatment[i] == 1) {
    dBEF_nem$treatment[i] = 3
  }
  else if(dBEF_nem$treatment[i] == 2) {
    dBEF_nem$treatment[i] = 2
  }
  else if(dBEF_nem$treatment[i] == 3) {
    dBEF_nem$treatment[i] = 1
  }
}
dBEF_nem$treatment <- factor(dBEF_nem$treatment, levels = c("1", "2", "3"))
dBEF_nem$treatment %>% str()

#### 1e-3 constant addition to zeros: ####

dBEF_nem <- dBEF_nem %>%
  mutate(Ba_per100gZeroC = ifelse(Ba_per100g == 0, 0.001, Ba_per100g), 
         .after=Ba_per100g) %>% 
  mutate(Fu_per100gZeroC = ifelse(Fu_per100g == 0, 0.001, Fu_per100g), 
         .after=Fu_per100g) %>%
  mutate(Pr_per100gZeroC = ifelse(Pr_per100g == 0, 0.001, Pr_per100g), 
         .after=Pr_per100g) %>%
  mutate(Pl_per100gZeroC = ifelse(Pl_per100g == 0, 0.001, Pl_per100g), 
         .after=Pl_per100g) %>%
  mutate(Om_per100gZeroC = ifelse(Om_per100g == 0, 0.001, Om_per100g), 
         .after=Om_per100g) 

####adding a week variable in '21 data####
dBEF_nem <- dBEF_nem %>%
  mutate(week = ifelse(year == 2021,
                       ifelse(block == "B1"| block == "B2", "W1", "W2"),
                       NA), .after=block)

#### creating a custom brms hurdle_gaussian family ####
  #this is the family: hurdle_gaussian as coded by Andrew Heiss:
  #https://www.andrewheiss.com/blog/2022/05/09/hurdle-lognormal-gaussian-brms/#exponentially-distributed-outcomes-with-zeros

#1 - create a custom brms family:
hurdle_gaussian <- 
  # Create a custom family that is logit if y = 0, normal/gaussian if not
  custom_family("hurdle_gaussian", 
                dpars = c("mu", "sigma", "hu"),
                links = c("identity", "log", "logit"),
                lb = c(NA, 0, NA),
                type = "real")

#2 - stan code to handle the sampling (idk what exactly it does):

stan_funs <- "
    real hurdle_gaussian_lpdf(real y, real mu, real sigma, real hu) { 
      if (y == 0) { 
        return bernoulli_lpmf(1 | hu); 
      } else { 
        return bernoulli_lpmf(0 | hu) +  
               normal_lpdf(y | mu, sigma); 
      } 
    }
  "
# Prepare Stan code for use in brm()
stanvars <- stanvar(scode = stan_funs, block = "functions")

#3 - post processing functions:  
posterior_predict_hurdle_gaussian <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  theta <- brms::get_dpar(prep, "hu", i = i)
  
  hu <- runif(prep$ndraws, 0, 1)
  ifelse(hu < theta, 0, rnorm(prep$ndraws, mu,sigma))
}

posterior_epred_hurdle_gaussian <- function(prep) {
  with(prep$dpars, mu * (1 - hu))
}

#4 - to use code, pass it to brm() using brm(... , stanvars = stanvars)  
ggplot(dBEF_nem, aes(x=block, y=SWC_gravimetric, col = col.sowndiv))+
  geom_jitter()

####creating subset df's for '21 and '17 and no 60 species ####
dBEF_nem21 <- subset(dBEF_nem, year==2021)
dBEF_nem17 <- subset(dBEF_nem, year==2017)



