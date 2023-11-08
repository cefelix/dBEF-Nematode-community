#### hurdle11: Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot), family = hurdle_lognormal ####
m.Fu.hurdle11 <- brm(
  bf(Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
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


dBEF_nem17$Fu_per100gLog.hurdle %>% density() %>% plot()


m.Fu.hurdle21 <- brm(
  bf(Fu_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
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
  data = dBEF_nem17, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Fu.hurdle31, ndraws=100)


#### hurdle41: Fu_per100gLog.hurdle ~ sowndiv*treatment + (1|block/plot), family = hurdle_gaussian ####
m.Fu.hurdle41 <- brm(
  bf(Fu_per100gLog.hurdle ~ sowndiv*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m.Fu.hurdle41, ndraws=100)

#### hurdle51: Fu_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot), family = hurdle_gaussian ####
SEED = 22061996
m17.Fu.hurdle51 <- brm(
  bf(Fu_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem17, 
  family = hurdle_gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m17.Fu.hurdle51, ndraws=100)

#### 61: Fu_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot), family = gaussian ####
SEED = 22061996
m17.Fu61 <- brm(
  Fu_per100gLog.hurdle ~ sowndivLog*treatment + (1|block/plot),
  data = dBEF_nem17, 
  family = gaussian,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))

pp_check(m17.Fu61, ndraws=100)

####saving Fu hurdle models####
save(#m.Fu.hurdle11,
     #m.Fu.hurdle21,
     #m.Fu.hurdle31,
     #m.Fu.hurdle41,
     m17.Fu.hurdle51,
     m17.Fu61,
     file="./statistics/17_remodelling/brms/231108_Fu17_hurdle.RData")

#load(file="./statistics/brms/231107_Fu_hurdle.RData")
conditional_effects(m17.Fu.hurdle51)
conditional_effects(m17.Fu61)
summary(m17.Fu61)

m17.Fu61 <- m17.Fu.hurdle61



