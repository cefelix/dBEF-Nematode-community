#packages:
library(brms)
library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)

SEED = 19111996



#check scetchy values like 0 or NA:
  #(this means that numerator or denominator, respectively, are zero):
summary((dBEF_nem21$EI == 0)) #7 zeros, 1 NAs
  summary((dBEF_nem21$CI == 0)) #5 zeros, 5 NAs
  summary((dBEF_nem21$Fu..Fu.Ba. == 0)) #4 zeros, 0 NAs
  summary((dBEF_nem21$MI == 0)) #0 zeros, 1 NAs
  summary((dBEF_nem21$SI == 0)) #10 zeros, 1 NAs


#complicated way:
dBEF_nem21_EI <- dBEF_nem21 %>% 
  filter(is.na(EI) == FALSE)
dBEF_nem21_EI$EI %>% density() %>% plot() #could use a hurdle model here as well!

#easy way:
#### exploration: Enrichment index #### 
dBEF_nem21 %>% filter(is.na(EI)==FALSE) %>% pull(EI) %>% density %>% plot() #hurdle_gaussian
ggplot(data = filter(dBEF_nem21, is.na(EI) == FALSE), 
       aes(x= sowndivLog, y= EI , color = treatment))+
  geom_jitter(width=0.2)

#### exploration: Channel index (0 to 100): 100 * (0.8*Fu2 / (3.2*Ba1 + 0.8*Fu2)) ####
dBEF_nem21 %>% filter(is.na(CI)==FALSE) %>% pull(CI) %>% density %>% plot() #weird peak at 100
  #a CI of 100 means no Ba-1 in the sample!
  #Thus: plot a density distribution where zero Ba-1 and zero Fu-2 samples are excluded:
  dBEF_nem21 %>% filter(is.na(CI)==FALSE) %>%
    filter(CI!=0 & CI !=100) %>%
    pull(CI) %>% density %>% plot() #only 119 samples left after this restriction
ggplot(data = filter(dBEF_nem21, is.na(CI) == FALSE), 
       aes(x= sowndivLog, y= CI , color = treatment))+
  geom_jitter(width=0.2)  
    #follow up: model "chance of of having no Ba-1" and 
    #"chance of of having no Fu-2" as a Bernoulli process?
  
#### exploration: Channel ratio (0 to 1): Fu / (Fu+Ba) ####
dBEF_nem21$Fu..Fu.Ba. %>% density() %>% plot()  #hurdle_gaussian
  dBEF_nem21 %>% filter(Fu..Fu.Ba. == 1) %>% pull(Fu..Fu.Ba.) %>% str() 
    #11 samples with zero Fungivores
  dBEF_nem21 %>% filter(Fu..Fu.Ba. == 0) %>% pull(Fu..Fu.Ba.) %>% str()
    #4 samples with zero Bacterivores
ggplot(data = filter(dBEF_nem21, is.na(Fu..Fu.Ba) == FALSE), 
       aes(x= sowndivLog, y=Fu..Fu.Ba , color = treatment))+
  geom_jitter(width=0.2)  
  
  
#### exploration: Structure index #### 
dBEF_nem21%>% filter(is.na(SI)==FALSE) %>% pull(SI) %>% density %>% plot() #hurdle_gaussian
ggplot(data = filter(dBEF_nem21, is.na(SI) == FALSE), 
       aes(x = sowndivLog, y = SI, color = treatment))+
  geom_jitter(width=0.2)  

#### exploration: Maturity index ####
dBEF_nem21%>% filter(is.na(MI)==FALSE) %>% pull(MI) %>% density %>% plot() #hurdle_gaussian
ggplot(data = filter(dBEF_nem21, is.na(MI) == FALSE), 
       aes(x = sowndivLog, y = MI , color = treatment))+
  geom_jitter(width=0.2)  

#### EI ~ sowndivLogStd*treatment + (1|block/plot), fam=hurdle_gamma ####
m.EI.21 <- brm(
  bf(EI ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_gamma,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)
) # 1 divergent transition
  #tail ESS too low
pp_check(m.EI.21, ndraws=100)

m.EI.22 <- update(m.EI.21, 
                  iter=3000, warmup=1500,
                  control = list(adapt_delta=0.999))
pp_check(m.EI.22, ndraws=100)

#hurdle against treatment:
m.EI.31 <- brm(
  bf(EI ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dBEF_nem21, 
  family = hurdle_gamma,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)
)
pp_check(m.EI.31, ndraws=100)
#still having probability mass at EI>100 -.-

#### EI ~ sowndivLogStd*treatment + (1|block/plot), fam=zero_one_inflated_beta ####
dBEF_nem21 <- dBEF_nem21 %>% 
  mutate(EI_ZeroOne = EI/100,
         .after = EI)

m.EI.41 <- brm(EI_ZeroOne ~ sowndivLogStd*treatment + (1|block/plot),
  data = dBEF_nem21, 
  family = zero_one_inflated_beta,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)
) #2 divergent transitions

m.EI.42 <- update(m.EI.41, 
                  control=list(adapt_delta=0.999))

pp_check(m.EI.42, ndraws=100)

blob <- dBEF_nem21_EI %>% 
  filter(EI_ZeroOne > 0.49 & EI_ZeroOne < 0.51) #%>%
  pull(EI_ZeroOne) %>%
  density() %>%
  plot()

nrow(dBEF_nem21)  
(dBEF_nem21$Ba1 == 0) %>% sum() #116 out of 240
(dBEF_nem21$Ba2 == 0) #22 out of 240   

#only the zero Ba1 subset:
dBEF_nem21_EI %>% 
  filter(Ba1 == 0) %>%
pull(EI_ZeroOne) %>%
  density() %>%
  plot()

#only the non-zero Ba1 subset:
dBEF_nem21_EI %>% 
  filter(Ba1 != 0) %>%
  pull(EI_ZeroOne) %>%
  density() %>%
  plot()


#### EI ~ sowndivLogStd*treatment + (1|block/plot), fam=Beta ####
dBEF_nem21 <- dBEF_nem21 %>% 
  mutate(EI_ZeroOne = EI/100,
         .after = EI)

dBEF_nem21_EI <- subset(dBEF_nem21, 
                        EI_ZeroOne != 0 & EI_ZeroOne != 1)
dBEF_nem21_EI$EI_ZeroOne %>% max()

m.EI.51 <- brm(EI_ZeroOne ~ sowndivLogStd*treatment + (1|block/plot),
               data = dBEF_nem21_EI, 
               family = Beta,
               stanvars = stanvars, #necessary to use custom brms families!
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               seed = SEED,
               control = list(adapt_delta=0.99)
) 
pp_check(m.EI.51, ndraws=100)





#### WRONG: EI ~ sowndivLogStd*treatment + (1|block/plot), fam=hurdle_gaussian ####
#problem: EI is bounded between zero and 1


m.EI.11 <-  brm(
  bf(EI ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.EI.11, ndraws=100) # a little to peaky

#another try with log-transformed EI
dBEF_nem21 <- dBEF_nem21 %>% 
  mutate(EI.Log = ifelse(EI == 0, 0, log(EI)),
         .after = EI)
  dBEF_nem21 %>% filter(is.na(EI.Log)==FALSE) %>% pull(EI.Log) %>% density %>% plot()
#the same model with log-transformed response:
  m.EI.21 <-  brm(
    bf(EI.Log ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ 1),
    data = dBEF_nem21, 
    family = hurdle_gaussian,
    stanvars = stanvars, #necessary to use custom brms families!
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) # 1 divergent transition
  
  m.EI.22 <- update(m.EI.21, 
                    control = list(adapt_delta=0.999))
  

  pp_check(m.EI.22, ndraws=100)  #still fit of same quality

#the same model with exp-transformed response    
  dBEF_nem21 <- dBEF_nem21 %>% 
    mutate(EI.Exp = ifelse(EI == 0, 0, exp(EI)),
           .after = EI)
  dBEF_nem21 %>% filter(is.na(EI.Exp)==FALSE) %>% 
    pull(EI.Log) %>% 
    density %>% plot() 
  
  m.EI.31 <-  brm(
    bf(EI.Exp ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ 1),
    data = dBEF_nem21, 
    family = hurdle_gaussian,
    stanvars = stanvars, #necessary to use custom brms families!
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) 
      #3 transitions exceeded max_treedepth
      #R-hat is Inf -> more iterations
      #Bulk ESS too low, tail ESS -> more iter
  
  m.EI.32 <- update(m.EI.31,
                    iter=4000, warmup=2000,
                    control = list(max_treedepth=15))#this does not work!
  
#the same model with square root transformed response:
  dBEF_nem21 <- dBEF_nem21 %>% 
    mutate(EI.Sqrt = ifelse(EI == 0, 0, sqrt(EI)),
           .after = EI)
  m.EI.41 <-  brm(
    bf(EI.Sqrt ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ 1),
    data = dBEF_nem21, 
    family = hurdle_gaussian,
    stanvars = stanvars, #necessary to use custom brms families!
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) 
  pp_check(m.EI.41, ndraws = 100) #same fit as non sqrt-transformed response
  
  ####zero one inflated bet for fu ba proportions####
  SEED = 123
  m.CR_zeroone <- brm(
    bf(Fu..Fu.Ba. ~ sowndivLogStd*treatment + (1|block/plot),
       zoi ~ 1,
       coi ~ 1),
    data = dBEF_nem21, 
    family = zero_one_inflated_beta,
    stanvars = stanvars, #necessary to use custom brms families!
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) 
    
  
  
#### hurdle channel ratio:  CR ~ sowndivLogStd*treatment + (1|block/plot), fam=hurdle_gaussian ####
  #brm can't read Fu..Fu.Ba, so:
  (dBEF_nem21$Fu..Fu.Ba. ==0) %>% sum()
  #try: zero-1-inflated beta distribution
  dBEF_nem21 <- dBEF_nem21 %>% 
    mutate(CR = Fu..Fu.Ba.,
           .before = Fu..Fu.Ba.)
  
  
  
  m.CR.11 <-  brm(
    bf(CR ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ 1),
    data = dBEF_nem21, 
    family = hurdle_gaussian,
    stanvars = stanvars, #necessary to use custom brms families!
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) 
      #3 divergent transitions
  
  m.CR.12 <- update(m.CR.11, 
                    control = list(adapt_delta=0.999))
      #3 transitions exceeded max_treedepth
  
  m.CR.13 <- update(m.CR.12, 
                    control = list(max_treedepth=12))
  
  
  pp_check(m.CR.13, ndraws=100) 
    #overestimating the amount of values > 0 
    #these are implausible values, CR is bounded between [0,1]!
  
  #use a prior that restricts CR to be between [0,1]:
  
####SI 11 hurdle: SI ~ sowndivLogStd*treatment + (1|block/plot), fam=hurdle_gaussian ####
  m.SI.11 <-  brm(
    bf(SI ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ 1),
    data = dBEF_nem21, 
    family = hurdle_gaussian,
    stanvars = stanvars, #necessary to use custom brms families!
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99))
  
  pp_check(m.SI.11, ndraws = 100)
  summary(m.SI.11)
  
    
#### hurdle MI: MI ~ sowndivLogStd*treatment + (1|block/plot), fam=hurdle_gaussian ####
  m.MI.11 <-  brm(
    bf(MI ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ 1),
    data = dBEF_nem21, 
    family = hurdle_gaussian,
    stanvars = stanvars, #necessary to use custom brms families!
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #1 divergent transition
  
  m.MI.12 <- update(m.MI.11, 
                    control= list(adapt_delta = 0.9999))
  
  pp_check(m.MI.13, ndraws = 100) #dont need a hurdle here, as no zeros in the data
  
#### MI ~ sowndivLogStd*treatment + (1|block/plot),   
  m.MI.21 <- brm(MI ~ sowndivLogStd*treatment + (1|block/plot),
    data = dBEF_nem21, 
    family = gaussian,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) 
  
  pp_check(m.MI.21, ndraws = 100)
  summary(m.MI.21)  
  
#### save models ####
  save(m.EI.11, #EI non-transformed
       m.EI.21, m.EI.22, #EI log-transformed
       m.EI.31, m.EI.32, #EI exp-transformed, DOES NOT WORK!
       m.EI.41, #EI sqrt-transformed, same fit as untransformed
       m.CR.11, m.CR.12, m.CR.13,
       m.SI.11, #SI, hurdled
       m.MI.11, m.MI.12, #hurdled MI models, but unnecessary as no zeros in MI
       m.MI.21, #MI without hurdle term
       file = "./statistics/brms/231107_NematodeIndices.RData")
load(file = "./statistics/brms/231107_NematodeIndices.RData")

conditional_effects(m.SI.11)


conditional_effects(m.MI.21)
conditional_effects(m.MI.12)



