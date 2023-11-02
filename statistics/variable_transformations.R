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
dBEF_nem <- dBEF_nem %>% # a colorcode of sowndiv for plotting
  mutate(col.sowndiv = as.factor(sowndiv), .after=sowndiv)

#adding "yearblock" which could be a random factor:
dBEF_nem <- dBEF_nem %>%
  mutate(yearblock = as.factor(paste(dBEF_nem$year, dBEF_nem$block, sep="")),
         .after = block)

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
  

#densities
dBEF_nem$Fu_per100gLog %>% density() %>% plot() #normal
dBEF_nem$Ba_per100gLog %>% density() %>% plot() #much probability on the right
dBEF_nem$Pl_per100gLog %>% density() %>% plot() #normal
dBEF_nem$Pr_per100gLog %>% density() %>% plot() #many zeros
dBEF_nem$Om_per100gLog %>% density() %>% plot() #way too many zeros, bimodal?
dBEF_nem$ind_per100gLog %>% density() %>% plot()

#### standardization ####
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
pairs(~Ba_per100g +  Fu_per100g + Pl_per100g + Om_per100g + Pr_per100g, data=dBEF_nem)


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



