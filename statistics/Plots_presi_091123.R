library(brms)
library(tidyr)
library(ggplot2)
library(tidybayes)
library(modelr)

#run variable transformations.R before this script!

#how to extract data from conditonal_effects():
  #https://discourse.mc-stan.org/t/change-aesthetics-conditional-effects/13500/5

####Predators 2021####
load("./statistics/brms/231108_Pr_hurdle.RData")
#best: m.Pr.hurdle51, log(pr)~log(sowndiv)

#get model predictions for regression plot   
  predictions <- conditional_effects(m.Pr.hurdle51)[[3]]
  predictions$estimate__

  ggplot(dBEF_nem21, aes(x=sowndivLog, y=Pr_per100g, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=exp(estimate__), col=treatment))
  
  rm(m.Pr.hurdle11, m.Pr.hurdle21, m.Pr.hurdle22, m.Pr.hurdle31, 
     m.Pr.hurdle42, m.Pr.hurdle51, m.Pr.hurdle61)

  
####Predators 2017####
load("./statistics/17_remodelling/brms/231108_Pr17_hurdle.RData")
#the good one is 51, log(Pr~log(sowndiv))
  predictions <- conditional_effects(m17.Pr.hurdle51)[[3]]

  ggplot(dBEF_nem17, aes(x=sowndivLog, y=Pr_per100g, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=exp(estimate__), col=treatment))
    #thats a poor fit as it seems!
    (dBEF_nem17$Pr_per100g ==0) %>% sum() # 25 zeros could be responsible...
  
  #lets plot the abundance on the log scale:
  ggplot(dBEF_nem17, aes(x=sowndivLog, y=Pr_per100gLog.hurdle, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=estimate__, col=treatment))

rm(m17.Pr.hurdle21, m17.Pr.hurdle51)


#### abundance 2021 ####



