####loading data and packages ####

library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
dBEF_nem %>%
  str()
load(file = "./statistics/brms/230929_Lm.RData")

dBEF_nem$treatment <- dBEF_nem$treatment %>%
  as.factor()

####examining distributions####

dBEF_nem #%>%#
hist(dBEF_nem$ind_per100g, breaks = seq(min(dBEF_nem$ind_per100g), max(dBEF_nem$ind_per100g), length.out=30))

#lets log transform this
dBEF_nem <- dBEF_nem %>%
  mutate(Log_ind_per100g = log(ind_per100g), .after = ind_per100g)

hist(dBEF_nem$Log_ind_per100g)




####simple model on overall density####

#m1: abun ~ sowndiv+(1|block) fam=lognormal
#m2: log(abun) ~ sowndiv+(1|block), fam=gaussian
#m3: log(abun) ~ sowndiv+treatment+(1|block), fam=gaussian
#m4: log(abun) ~ sowndiv+SH+PH+(1|block), fam= gaussian, SH+PH numeric
#m5: log(abun) ~ sowndiv+SH+(1|block), fam= gaussian, SH factor
#m6: log(abun) ~ sowndiv*SH+(1|block), fam= gaussian, SH numeric

m1.1 <- brm(ind_per100g ~ sowndiv + (1|block),
         data = dBEF_nem, family = "lognormal",
         chains = 3,
         cores = 3,
         iter = 2000, warmup = 1000)
summary(m1.1) #5 divergent transitions

m1.2 <- brm(ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "lognormal",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta=0.99))
summary(m1.2) #still 1 divergent transition

m1.3 <- brm(ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "lognormal",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta=0.9999))
summary(m1.3) # no divergent transitions
#
pp_check(m1.3)


m2.1 <- brm(Log_ind_per100g ~ sowndiv + (1|block),
          data = dBEF_nem, family = "gaussian",
          chains = 3,
          cores = 3,
          iter = 2000, warmup = 1000)
summary(m2.1) # 3 divergent transitions

m2.2 <- brm(Log_ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.99))
summary(m2.2) # 1 divergent transition

m2.3 <- brm(Log_ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.9999))

summary(m2.3) #184 transitions after warmup exceeding max. treedepth h -> increase h above 10

m2.4 <- brm(Log_ind_per100g ~ sowndiv + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.9999,
                           max_treedepth = 11))
summary(m2.4)

pp_check(m2.4)

m3.1 <- brm(Log_ind_per100g ~ sowndiv + treatment + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000)
summary(m3.1) #6 divergent transitions

m3.2 <- brm(Log_ind_per100g ~ sowndiv + treatment + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.99))
summary(m3.2)

pp_check(m3.2)


#include SH and PH as numeric variables
m4.1 <- brm(Log_ind_per100g ~ sowndiv + SH + PH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000)
summary(m4.1)
pp_check(m4.1)

m4.2 <- brm(Log_ind_per100g ~ sowndiv + SH + PH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.99))
summary(m4.2)

m4.3 <- brm(Log_ind_per100g ~ sowndiv + SH + PH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.9999))
summary(m4.3)

pp_check(m4.2)



#include SH and PH as a factors
dBEF_nem$SH <- dBEF_nem$SH %>% as.factor()
dBEF_nem$PH <- dBEF_nem$PH %>% as.factor()

m5.1 <- brm(Log_ind_per100g ~ sowndiv + SH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 10))
summary(m5.1)

#but actually, we want to see whether the effect of sowndiv depends on SH
#thus: interaction
dBEF_nem$SH <- as.numeric(dBEF_nem$SH)

m6.1 <- brm(Log_ind_per100g ~ sowndiv * SH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000)
summary(m6.1) # 4 divergent transitions

m6.2 <- brm(Log_ind_per100g ~ sowndiv * SH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.99))
summary(m6.2) #1 divergent transition

m6.3 <- brm(Log_ind_per100g ~ sowndiv * SH + (1|block),
            data = dBEF_nem, family = "gaussian",
            chains = 3,
            cores = 3,
            iter = 2000, warmup = 1000,
            control = list(adapt_delta = 0.9999))
summary(m6.3)

pp_check(m6.3)

           



####compare models####

#gaussian (logtransformed data)
summary(m2.4)

#lognormal(untransformed data)
summary(m1.3)


####saving brm outputs####
save(m1.1, m1.2,m1.3, 
     m2.1, m2.2, m2.3, m2.4,
     m3.1, m3.2,
     m4.1, m4.2, m4.3,
     m5.1,
     m6.1, m6.2, m6.3,
     file = "./statistics/brms/231004_Lm.RData")





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