####loading data and packages ####

library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
dBEF_nem %>%
  str()

####examining distributions####

dBEF_nem #%>%#
hist(dBEF_nem$ind_per100g, breaks = seq(min(dBEF_nem$ind_per100g), max(dBEF_nem$ind_per100g), length.out=30))

a <- dBEF_nem$ind_per100g %>%
  log()
hist(a)

dBEF_nem$ind_per100g %>%
  hist()

hist(dBEF_nem$ind_per100g, breaks = seq(min(dBEF_nem$ind_per100g), max(dBEF_nem$ind_per100g), length.out=30))



####simple model on overall density####


m1.1 <- brm(ind_per100g ~ sowndiv,
          data = dBEF_nem, family = "lognormal",
          chains = 3,
          cores = 3,
          iter = 2000, warmup = 1000)
summary(m1.1)


m1.2 <- brm(ind_per100g ~ sowndiv + (1|block),
         data = dBEF_nem, family = "lognormal",
         chains = 3,
         cores = 3,
         iter = 2000, warmup = 1000)
summary(m1.2)
#

m2.1 <- brm(ind_per100g ~ sowndiv,
          data = dBEF_nem, family = "gaussian",
          chains = 3,
          cores = 3,
          iter = 2000, warmup = 1000)
summary(m2.1)

m2.2 <- brm(ind_per100g ~ sowndiv + (1|block),
          data = dBEF_nem, family = "gaussian",
          chains = 3,
          cores = 3,
          iter = 2000, warmup = 1000)
summary(m2.2)


####saving brm outputs####
save(m1.1, m1.2, m2.1, m2.2,
     file = "./statistics/brms/230928_Lm.RData")





####old####


parallel::detectCores()



m <- brm(bf(CR ~ 1+ treatment + sowndiv +treatment:sowndiv + (1 + treatment|plot)),
         family = zero_one_inflated_beta(),
         chains = 3,
         cores = 3,
         iter = 2000,
         #backend = 'cmdstanr',
         data = data.analysis[is.na(data.analysis$CR)==F,])
summary(m)

pp_check(m)
plot(m)


conditional_effects(m)

plot(conditional_effects(m), ask = F) #plots predictions

#exponent of vegan::diversity
#paper/blog by jost

#abundances of whole commnity
#diversity of complete nematodes
#abundances of feeding groups
#mi 


#package : tidybayes

#https://ben18785.shinyapps.io/distribution-zoo/