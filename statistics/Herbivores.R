
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
  dBEF_nemSH1 <- subset(dBEF_nem, SH == 1)
  dBEF_nemSH5 <- subset(dBEF_nem, SH == 5)
  dBEF_nemSH15 <- subset(dBEF_nem, SH == 15)
  dBEF_nemSH19 <- subset(dBEF_nem, SH == 19)
  
  dBEF_nem21 <- subset(dBEF_nem, year == 2021)


#### exploration ####
p.1 <- ggplot(dBEF_nemSH1, aes(y = Pl_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  labs(title = "SH 1")+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")

p.5 <- ggplot(dBEF_nemSH5, aes(y = Pl_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")+
  labs(title = "SH 5")

p.15 <- ggplot(dBEF_nemSH15, aes(y = Pl_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 15")+
  geom_smooth(method="lm")

p.19 <- ggplot(dBEF_nemSH19, aes(y = Pl_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 19")+
  geom_smooth(method="lm")

grid.arrange(p.1, p.5, p.15, p.19)
rm(p.1, p.5, p.15, p.19)

dBEF_nem$





#### Pl_per100gLog ~ sowndiv*treatment + (1|block) ####
m.Pl.31 <- brm(Pl_per100gLog ~ sowndiv*treatment + (1|block),
               data = dBEF_nem21, family = "gaussian",
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               control = list(adapt_delta=0.9))
pp_check(m.Pl.31, ndraws = 100)


#### Pl_per100gLog ~ sowndiv*treatment + (1|block/plot) ####
m.Pl.31b <- brm(Pl_per100gLog ~ sowndiv*treatment + (1|block/plot),
               data = dBEF_nem21, family = "gaussian",
               chains = 3,
               cores = 3,
               iter = 2000, warmup = 1000,
               control = list(adapt_delta=0.9))
update(m.Pl.31b,
       control=)
pp_check(m.Pl.31b, ndraws = 100)


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


#### Pl_per100g ~ sowndiv * treatment + (1|block), fam=lognormal ####
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

####Pl_per100g ~ sowndivLogStd * treatment + (1|block/plot) ####
m.Pl.nest11 <- brm(Pl_per100gLog ~ sowndivLogStd*treatment + (1|block/plot),
                   data = dBEF_nem21, family = "gaussian",
                   chains = 3,
                   cores = 3,
                   iter = 2000, warmup = 1000,
                   control = list(adapt_delta=0.9)) #7 divergent transitions

m.Pl.nest12 <- update(m.Pl.nest11,
                      control=list(adapt_delta=0.99))

pp_check(m.Pl.nest12, ndraws = 100)



#### hurdle: Pl_per100g ~ sowndivLogStd * treatment + (1|block/plot) ####
  #lets have a quick look at the data:
  dBEF_nem21$Pl_per100gLog %>% density() %>% plot()
  sum(dBEF_nem21$Pl_per100g == 0) #2 samples have zero herbivores

SEED <- 22061996

dBEF_nem21$Pl_per100g %>% density %>% plot()
dBEF_nem21$Pl_per100gLog %>% density %>% plot()

dBEF_nem21 <- dBEF_nem21 %>%
  mutate(Pl_per100gLog.hurdle = ifelse(Pl_per100g == 0, 
                                       0, log(Pl_per100g)),
         .after = Pl_per100gLog)

m.Pl.hurdle11 <- brm(
  bf(Pl_per100gLog.hurdle ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ 1),
  data = dBEF_nem21, 
  family = "hurdle_lognormal",
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99))


m.Pl.hurdle12 <- update(m.Pl.hurdle11, 
                             seed = SEED,
                             iter = 3000, warmup=1500,
                             control = list(adapt_delta=0.99))

m.Pl.hurdle13 <- update(m.Pl.hurdle12, 
                             control = list(adapt_delta = 0.999))

m.Pl.hurdle14 <- update(m.Pl.hurdle13,
                             iter = 4000, warmup=2000,
                             control = list(adapt_delta = 0.9999))
                             # still 4 divergent transitions -.-
                             # maybe specify more restictive priors


pp_check(m.Pl.hurdle14, ndraws = 100)

#### saving Pl hurdle models ####
save(m.Pl.hurdle11,
     m.Pl.hurdle12,
     m.Pl.hurdle13,
     m.Pl.hurdle14,
     file="./statistics/brms/231103_Pl_hurdle.RData")



