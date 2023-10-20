#load data:
library(rethinking)
data(milk)

d <- 
  milk %>%
  drop_na(ends_with("_s"))
rm(milk)

d <-
  d %>%
  mutate(neocortex = neocortex.perc / 100)
dim(d)

#load brms:
detach(package:rethinking, unload = T)
library(brms)


#model fitting
inits <- list(Intercept = mean(d$kcal.per.g),
              sigma     = sd(d$kcal.per.g))

inits_list <-list(inits, inits, inits, inits)

b6.11 <- 
  brm(data = d, family = gaussian,
      kcal.per.g ~ 1,
      prior = c(prior(uniform(-1000, 1000), class = Intercept),
                prior(uniform(0, 100), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      inits = inits_list,
      seed = 6)

inits <- list(Intercept = mean(d$kcal.per.g),
              neocortex = 0,
              sigma     = sd(d$kcal.per.g))
inits_list <-list(inits, inits, inits, inits)

b6.12 <- 
  brm(data = d, family = gaussian,
      kcal.per.g ~ 1 + neocortex,
      prior = c(prior(uniform(-1000, 1000), class = Intercept),
                prior(uniform(-1000, 1000), class = b),
                prior(uniform(0, 100), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      inits = inits_list,
      seed = 6)

inits <- list(Intercept   = mean(d$kcal.per.g),
              `log(mass)` = 0,
              sigma       = sd(d$kcal.per.g))
inits_list <-list(inits, inits, inits, inits)

b6.13 <-
  update(b6.12, 
         newdata = d,
         formula = kcal.per.g ~ 1 + log(mass),
         inits   = inits_list)

inits <- list(Intercept   = mean(d$kcal.per.g),
              neocortex   = 0,
              `log(mass)` = 0,
              sigma       = sd(d$kcal.per.g))
inits_list <-list(inits, inits, inits, inits)

b6.14 <- 
  update(b6.13, 
         newdata = d,
         formula = kcal.per.g ~ 1 + neocortex + log(mass),
         inits   = inits_list)


####comparing WAIC values####
#model terms:
#11: kcal ~ 1
#12: kcal ~ 1+neocortex
#13: kcal ~ 1+log(mass)
#14: kcal ~ 1+neocortex+log(mass)


waic(b6.14)
waic(b6.11)

b6.11 <- add_criterion(b6.11, "waic")
b6.11$waic
