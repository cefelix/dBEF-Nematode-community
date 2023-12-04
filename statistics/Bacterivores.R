####data and packages####
  library(brms)
  library(dplyr)
  library(tidyr)
  library(bayesplot)
  library(ggplot2)
  library(gridExtra)
  library(hexbin)
  library(GGally)
  library(emmeans)

# a seed:
SEED = 22061996


#### exploration ####
  dBEF_nemSH1 <- subset(dBEF_nem, SH == 1)
  dBEF_nemSH5 <- subset(dBEF_nem, SH == 5)
  dBEF_nemSH15 <- subset(dBEF_nem, SH == 15)
  dBEF_nemSH19 <- subset(dBEF_nem, SH == 19)
  
p.1 <- ggplot(dBEF_nemSH1, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  labs(title = "SH 1")+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")

p.5 <- ggplot(dBEF_nemSH5, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  geom_smooth(method="lm")+
  labs(title = "SH 5")

p.15 <- ggplot(dBEF_nemSH15, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 15")+
  geom_smooth(method="lm")

p.19 <- ggplot(dBEF_nemSH19, aes(y = Ba_per100g, x=log(sowndiv)) )+
  geom_jitter(aes(col=col.sowndiv))+
  scale_y_continuous(limits= c(0, 2600))+
  labs(title= "SH 19")+
  geom_smooth(method="lm")

grid.arrange(p.1, p.5, p.15, p.19)
rm(p.1, p.5, p.15, p.19,
   dBEF_nemSH19, dBEF_nemSH15, dBEF_nemSH5, dBEF_nemSH1)

#### 31b hurdle: Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot), hu~term, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 22061996
dat <- subset(dBEF_nem21, sowndiv != 60)
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog)

m.Ba.hurdle31b <- brm(
  bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~ sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  stanvars = stanvars, #necessary to use custom brms families!
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all fine

pp_check(m.Ba.hurdle31b, ndraws=100)+
  xlim(0,1000)

#### 31b_Ba1-4 hurdle: Ba1-4 ~ sowndivLogStd*treatment + (1|block/plot), hu~term, family = hurdle_lognormal, no 60 sp. plots ####
SEED = 22061996
dat <- subset(dBEF_nem21, sowndiv != 60)
#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog)

sum(dat$Ba1 == 0) #109
sum(dat$Ba2 == 0) #22
sum(dat$Ba3 == 0) #218 of 228
sum(dat$Ba4 == 0) #139

m.Ba1.hurdle31b <- brm(
  bf(Ba1 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~  sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all fine

pp_check(m.Ba1.hurdle31b, ndraws=100)+
  xlim(0,200) 

m.Ba2.hurdle31b <- brm(
  bf(Ba2 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~  sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good

pp_check(m.Ba2.hurdle31b, ndraws=100)+
  xlim(0,700) 

#skipping Ba-cp3 as there are only 10 samples with a non-zero density

m.Ba4.hurdle31b <- brm(
  bf(Ba4 ~ sowndivLogStd*treatment + (1|block/plot),
     hu ~  sowndivLogStd*treatment + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) 
  #all good

pp_check(m.Ba4.hurdle31b, ndraws=100)+
  xlim(0,100) 


###31ab Ba model comparisons####
m.Ba.hurdle31a <- add_criterion(m.Ba.hurdle31a, "loo")
m.Ba.hurdle31b <- add_criterion(m.Ba.hurdle31b, "loo")

m.Ba.hurdle31a$criteria  #looic: 2502.1 (for looic: smaller = better)
m.Ba.hurdle31b$criteria  #looic: 2500

#high pareto-k values: look whether p_loo is bigger/smaller than true number of parameters
#as described here: https://mc-stan.org/loo/reference/loo-glossary.html
#to get the true number of paremeters, follow this:
#https://discourse.mc-stan.org/t/determine-number-of-parameters-in-brms-gamm-to-compare-to-p-loo-value/23804/3

#for this you have to re-compile the model: https://github.com/quentingronau/bridgesampling/issues/7
rstan::get_num_upars(m.Ba.hurdle31b$fit) #177, thus p > N/5
m.Ba.hurdle31b$criteria  #p_loo = 38.1
  
pp_check(m.Ba.hurdle31b, ndraws=100)+
  xlim(0,1200)

#for elpd: higher = better  
loo_compare(m.Ba.hurdle31a, m.Ba.hurdle31b) 




####51b hurdle: Ba ~ sowndivLogStd*treatment*week + (1|block/plot), hu~term, family=hurdle_lognormal, no 60 sp. plots####
m.Ba.hurdle51b <- brm(
  bf(Ba_per100g ~ sowndivLogStd + treatment + week + 
       sowndivLogStd:treatment + sowndivLogStd:week + treatment:week +
       sowndivLogStd:treatment:week + (1|block/plot),
     hu ~ sowndivLogStd + treatment + week + 
       sowndivLogStd:treatment + sowndivLogStd:week + treatment:week +
       sowndivLogStd:treatment:week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #all good

pp_check(m.Ba.hurdle51b, ndraws=100)+
  xlim(0,1000)

summary(m.Ba.hurdle51b)

####all 2way interactions####
m.Ba.hurdle51c <- brm(
  bf(Ba_per100g ~ sowndivLogStd + treatment + week + 
       sowndivLogStd:treatment + sowndivLogStd:week + treatment:week + (1|block/plot),
     hu ~ week + (1|block/plot)),
  data = dat, 
  family = hurdle_lognormal,
  chains = 3,
  cores = 3,
  iter = 2000, warmup = 1000,
  seed = SEED,
  control = list(adapt_delta=0.99)) #1 diverg

m.Ba.hurdle52c <- update(m.Ba.hurdle51c ,
                         control=list(adapt_delta=0.999)) #986 exceeded max_treedepth 

m.Ba.hurdle52c <- update(m.Ba.hurdle51c ,
                         control=list(adapt_delta=0.999,
                                      max_treedepth=12)) #4divergent transitions
pp_check(m.Ba.hurdle52c , ndraws=100)+ 
  xlim(0,1000)


#look at predictions
  #make a 3 way interaction
  conditions <- make_conditions(m.Ba.hurdle51d, "week")
  predictions.Ba<- conditional_effects(m.Ba.hurdle51d, "sowndivLogStd:treatment", conditions = conditions)[[1]]
  
  predictions = predictions.Ba
  treatments = c("+SH+PH", "+SH-PH", "-SH-PH")
  cols=c("brown2", "darkolivegreen", "dodgerblue3")
  BREAKS <- unique(dat$sowndivLogStd) %>% 
    sort()
  
  p = ggplot(dat, aes(x=sowndivLogStd, y=Ba_per100g) )+
    #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, http://127.0.0.1:35231/graphics/plot_zoom_png?width=1200&height=900
    #                                  fill=treatment), 
    #            alpha=0.2, show.legend=FALSE)+
    geom_jitter(width=0.2,  alpha=0.4, shape=19,
                aes(col=treatment))+
    geom_line(data=predictions, aes(x= sowndivLogStd, y=estimate__, 
                                    linetype=week, col=treatment),
              linewidth= 1, show.legend = FALSE)+
    #scale_color_manual(labels=treatments, values = cols)+
    scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                       labels = c("1", "2", "4", "8", "16"))+ 
    scale_y_continuous(name = "bacterivores per 100g DW")+
    theme_bw()+
    theme(legend.position ="bottom")
  
  p
  
  #compare the weeks
  dat2 <- dat %>% 
    mutate(sowndivLogStd = as.factor(sowndivLogStd))
  ggplot(dat2, aes(x=sowndivLogStd, y=Ba_per100g, col=treatment, shape=week) )+
    #geom_jitter(width = 0.2)+
    geom_point(position = position_dodge(width=0.7), alpha=0.9)+
    geom_boxplot(alpha=0.4, outlier.shape = NA)
    
  
  rm(dat2)

sum(dat$Ba_per100g ==0)
#### savemodels ####

save(  m.Ba.hurdle31b,
       m.Ba1.hurdle31b, 
       m.Ba2.hurdle31b,
       #Ba-3 skipped, as 10 of 228 samples have a non-zero density
       m.Ba4.hurdle31b, 
     
     #m.Ba.hurdle41a,
     #m.Ba.hurdle41b,
     
     file="./statistics/brms/231204_Ba.RData")

load(file="./statistics/brms/231127_Ba.RData")
conditional_effects(m.Ba.hurdle11)
