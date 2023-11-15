
#### 11 hurdle: Pr_per100gLog ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Pr.hurdle11 <- brm(
  bf(Pr_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle11, ndraws=100)


#### 21 hurdle: Pr_per100gLog ~ sowndiv*treatment + (1|block/plot), fam=hurdle_gaussian ####
SEED = 19111996
m.Pr.hurdle21 <- brm(
  bf(Pr_per100gLog.hurdle ~ sowndiv*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

m.Pr.hurdle22 <- update(m.Pr.hurdle21,
                        control = list(adapt_delta = 0.999))

pp_check(m.Pr.hurdle22, ndraws=100)

#### 31 hurdle: Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot), fam=hurdle_lognormal ####
SEED = 19111996
m.Pr.hurdle31 <- brm(
  bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle31, ndraws=100)+
  xlim(0,300)

#### 41 hurdle: Pr_per100g ~ sowndiv*treatment + (1|block/plot), fam=hurdle_lognormal ####
SEED = 19111996
m.Pr.hurdle41 <- brm(
  bf(Pr_per100g ~ sowndiv*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

m.Pr.hurdle42 <- update(m.Pr.hurdle41,
                        control=list(adapt_delta=0.999)
)

pp_check(m.Pr.hurdle42, ndraws=100)+
  xlim(0,300)

#### 51 hurdle: Pr_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot), fam=hurdle_gaussian ####
SEED = 19111996
m17.Pr.hurdle51 <- brm(
  bf(Pr_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m17.Pr.hurdle51, ndraws=100)

conditional_effects(Pr_per100gLog.hurdle)

#### 61 hurdle: Pr_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot), fam=hurdle_gaussian ####
SEED = 19111996
m17.Pr.hurdle61 <- brm(
  bf(Pr_per100g ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Pr.hurdle61, ndraws=100)+
  xlim(0,300)



####comparison - waic ####

waic(m.Pr.hurdle31, m.Pr.hurdle42)
loo(m.Pr.hurdle31, m.Pr.hurdle42)

waic(m.Pr.hurdle31)
loo(m.Pr.hurdle31)
#### save hurdle models ####
save(#m.Pr.hurdle11, #Pr.Log ~ sowndivLogStd, fam=hurdle_gaussian
     m17.Pr.hurdle21, #Pr.Log ~ sowndiv, fam=hurdle_gaussian
     #m.Pr.hurdle31, #Pr ~ sowndivLogStd, fam=hurdle_lognormal
     #m.Pr.hurdle42, #Pr ~ sowndiv, fam=hurdle_lognormal
     m17.Pr.hurdle51, #Pr.Log ~ sowndivLog, fam=hurdle_gaussian
     #.Pr.hurdle61, #Pr ~ sowndivLog, fam=hurdle_lognormal
     file = "./statistics/17_remodelling/brms/231108_Pr17_hurdle.RData")


load(file = "./statistics/brms/231108_Pr17_hurdle.RData")
pp_check(m.Pr.hurdle21, ndraws=100)
conditional_effects(m17.Pr.hurdle51)

