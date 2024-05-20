


#data:
#exclude 60 sp.:
dat <- subset(dBEF_nem21, func.group != 60) 

datW1 <- subset(dat, week=="W1")
datW2 <- subset(dat, week=="W2")



treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")


#### Ba ~ funcdiv, W1-p, W2-, both- ####
  #W1-p: 0 div, W1-d: 2 div
  #W1-p: 0 div, W1-d: 3 div
  #both-p: 0 div, both-d: 0 div

  #different slopes between week 1 and 2 for treatment 1 and 2

    SEED = 22061996
    beta_coeff_priors <- prior(normal(0,20), class = "b")  
    sum(subset(dat, week=="W1")$Ba_per100g == 0) #9
    sum(subset(dat, week=="W2")$Ba_per100g == 0) #2
    
    
    #for week 1:
    m.Ba_func.groupW1_p <- brm(
      bf(Ba_per100g ~ func.group*treatment + (1|block/plot),
         hu ~ func.group*treatment + (1|block/plot)),
      data = subset(dat, week=="W1"), 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #all good
    
    summary(m.Ba_func.groupW1_p) #remove hurdle interaction
    pp_check(m.Ba_func.groupW1_p, ndraws=100)+
      xlim(0,2000)
    
        m.Ba_func.groupW1_p2 <- brm(
          bf(Ba_per100g ~ func.group*treatment + (1|block/plot),
             hu ~ func.group + treatment + (1|block/plot)),
          data = subset(dat, week=="W1"), 
          prior = beta_coeff_priors,
          family = hurdle_lognormal,
          stanvars = stanvars, #necessary to use custom brms families!
          chains = 3,
          cores = 3,
          iter = 2000, warmup = 1000,
          seed = SEED,
          control = list(adapt_delta=0.99)) #1 div
        
        summary(m.Ba_func.groupW1_p2) #remove hu~treatment
        
        m.Ba_func.groupW1_p3 <- brm(
          bf(Ba_per100g ~ func.group*treatment + (1|block/plot),
             hu ~ func.group + (1|block/plot)),
          data = subset(dat, week=="W1"), 
          prior = beta_coeff_priors,
          family = hurdle_lognormal,
          stanvars = stanvars, #necessary to use custom brms families!
          chains = 3,
          cores = 3,
          iter = 2000, warmup = 1000,
          seed = SEED,
          control = list(adapt_delta=0.99)) #7 div transitions
    
        summary(m.Ba_func.groupW1_p3) #remove hu~funcdiv
        
        m.Ba_func.groupW1_p4 <- brm(
          bf(Ba_per100g ~ func.group*treatment + (1|block/plot),
             hu ~ 1),
          data = subset(dat, week=="W1"), 
          prior = beta_coeff_priors,
          family = hurdle_lognormal,
          stanvars = stanvars, #necessary to use custom brms families!
          chains = 3,
          cores = 3,
          iter = 2000, warmup = 1000,
          seed = SEED,
          control = list(adapt_delta=0.99)) #3 div
        pp_check(m.Ba_func.groupW1_p4, ndraws=100)+
          xlim(0,1800)
        summary(m.Ba_func.groupW1_p4, prob=0.9) #narrower prior
        
        
        beta_coeff_priors2 <- prior(normal(0,10), class = "b")  
        m.Ba_func.groupW1_p5 <- update(m.Ba_func.groupW1_p4,
                                       prior = beta_coeff_priors2,
                                       seed = SEED) #7 div 
        pp_check(m.Ba_func.groupW1_p5, ndraws=100)+
          xlim(0,1800)
        summary(m.Ba_func.groupW1_p5, prob=0.9)
        
        beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
        m.Ba_func.groupW1_p6 <- update(m.Ba_func.groupW1_p4,
                                       prior = beta_coeff_priors3,
                                       seed = SEED) #7 div 
        pp_check(m.Ba_func.groupW1_p6, ndraws=100)+
          xlim(0,1800)
        summary(m.Ba_func.groupW1_p6, prob=0.9)
        
        
        
        loo(m.Ba_func.groupW1_p4, m.Ba_func.groupW1_p5, m.Ba_func.groupW1_p6)
    
    #with default prior:
    m.Ba_func.groupW1_d <- brm(
      bf(Ba_per100g ~ func.group*treatment + (1|block/plot),
         hu ~ func.group*treatment + (1|block/plot)),
      data = subset(dat, week=="W1"), 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99))  #2 div
    
    summary(m.Ba_func.group_d)
    pp_check(m.Ba_func.groupW1_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Ba_func.groupW1_p, m.Ba_func.groupW1_d)
    print(rstan::get_elapsed_time(m.Ba_func.groupW1_p$fit))
    print(rstan::get_elapsed_time(m.Ba_func.groupW1_d$fit))
    
    
    
#for week 2:
    m.Ba_func.groupW2_p <- brm(
      bf(Ba_per100g ~ func.group*treatment + (1|block/plot),
         hu ~ 1),
      data = subset(dat, week=="W2"), 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #all good
    
    summary(m.Ba_func.groupW2_p) 
    pp_check(m.Ba_func.groupW2_p, ndraws=100)
    
    #with default priors
    m.Ba_func.groupW2_d <- brm(
      bf(Ba_per100g ~ func.group*treatment + (1|block/plot),
         hu ~ 1),
      data = subset(dat, week=="W2"), 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99))  #3 div
    
    summary(m.Ba_func.groupW2_d)
    pp_check(m.Ba_func.groupW2_d, ndraws=100)
    
    #compare them  
    loo(m.Ba_func.groupW2_p, m.Ba_func.groupW2_d)
    print(rstan::get_elapsed_time(m.Ba_func.groupW2_p$fit))
    print(rstan::get_elapsed_time(m.Ba_func.groupW2_d$fit))
    loo(m.Ba_func.groupW2_d, m.Ba_func.groupW2_p)
    
    #for both weeks  
    m.Ba_func.group_p <- brm(
      bf(Ba_per100g ~ func.group*treatment + (1|block/plot),
         hu ~ func.group*treatment + (1|block/plot)),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #all good
    
    summary(m.Ba_func.group_p, prob=0.9)
    pp_check(m.Ba_func.group_p, ndraws=100)+
      xlim(0,2000)
    
    #with default priors
    m.Ba_func.group_d <- brm(
      bf(Ba_per100g ~ func.group*treatment + (1|block/plot),
         hu ~ func.group*treatment + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #all good
    
    summary(m.Ba_func.group_d, prob=0.9)
    pp_check(m.Ba_func.group_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Ba_func.group_p, m.Ba_func.group_d)
    print(rstan::get_elapsed_time(m.Ba_func.group_p$fit))
    print(rstan::get_elapsed_time(m.Ba_func.group_d$fit))  
    
  #save them:
    save(m.Ba_func.group_d, m.Ba_func.group_p,
         m.Ba_func.groupW1_d, m.Ba_func.groupW1_p,
          m.Ba_func.groupW1_p4, m.Ba_func.groupW1_p5, m.Ba_func.groupW1_p6,
         m.Ba_func.groupW2_d, m.Ba_func.groupW2_p,
         file = "./statistics/brms/231208_Ba_funcdiv_priors.RData")

#### Ba plot ####
    pred.Ba_prior1 <- conditional_effects(m.Ba_func.groupW1_p)[[3]]
    pred.Ba_def1 <- conditional_effects(m.Ba_func.groupW1_d)[[3]]
    summary(m.Ba_func.groupW1_p)
    
    pred.Ba_prior2 <- conditional_effects(m.Ba_func.groupW2_p)[[3]]
    pred.Ba_def2 <- conditional_effects(m.Ba_func.groupW2_d)[[3]]
    summary(m.Ba_func.groupW2_p)
    
    pred.Ba_prior  <- conditional_effects(m.Ba_func.group_p)[[3]]
    pred.Ba_def  <- conditional_effects(m.Ba_func.group_d)[[3]]
    summary(m.Ba_func.group_d)
    
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    BREAKS = unique(dat$func.group)    
    
    
    
    
ggplot(data = dat, aes(x= func.group, y=Ba_per100g))+
  #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
  #                                 fill=treatment), 
  #           alpha=0.2, show.legend=FALSE)+
  geom_jitter(data =datW1,
              width=0.2, shape=19, alpha=0.3, 
              aes(col=treatment))+
  geom_jitter(data=datW2,
              width=0.2, shape=18, alpha=0.9, 
              aes(col=treatment))+
  #predictions week1
  geom_line(data=pred.Ba_prior1, aes(x= func.group, y=estimate__, 
                                     linetype="dashed", col=treatment),
            linewidth= 0.5, show.legend = FALSE)+
  geom_line(data=pred.Ba_def1, aes(x= func.group, y=estimate__, 
                                   linetype="solid", col=treatment),
            linewidth= 0.5, show.legend = FALSE)+
  #predictions week 2:
  geom_line(data=pred.Ba_prior2, aes(x= func.group, y=estimate__, 
                                     linetype="dashed", col=treatment),
            linewidth= 0.5, show.legend = FALSE)+
  geom_line(data=pred.Ba_def2, aes(x= func.group, y=estimate__, 
                                   linetype="solid", col=treatment),
            linewidth=0.5, show.legend = FALSE)+
  #models for both weeks:
  geom_line(data=pred.Ba_prior, aes(x= func.group, y=estimate__, 
                                    linetype="dashed", col=treatment),
            linewidth= 0.7, show.legend = FALSE)+
  geom_line(data=pred.Ba_def, aes(x= func.group, y=estimate__, 
                                  linetype="solid", col=treatment),
            linewidth=0.7, show.legend = FALSE)+
  scale_color_manual(labels=treatments2, values = cols)+
  scale_x_continuous(name = "functional group richness")+
  scale_y_continuous(name = "Ba per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")  

#### Fu ~ funcdiv plot ####
ggplot(data = dat, aes(x= func.group, y=Fu_per100g))+
  #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
  #                                 fill=treatment), 
  #           alpha=0.2, show.legend=FALSE)+
  geom_jitter(data =datW1,
              width=0.2, shape=19, alpha=0.3, 
              aes(col=treatment))+
  geom_jitter(data=datW2,
              width=0.2, shape=18, alpha=0.9, 
              aes(col=treatment))+
  scale_color_manual(labels=treatments2, values = cols)+
  scale_x_continuous(name = "functional group richness")+
  scale_y_continuous(name = "Fu per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")  

#### Pl ~ funcdiv plot ####
ggplot(data = dat, aes(x= func.group, y=Pl_per100g))+
  #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
  #                                 fill=treatment), 
  #           alpha=0.2, show.legend=FALSE)+
  geom_jitter(data =datW1,
              width=0.2, shape=19, alpha=0.3, 
              aes(col=treatment))+
  geom_jitter(data=datW2,
              width=0.2, shape=18, alpha=0.9, 
              aes(col=treatment))+
  scale_color_manual(labels=treatments2, values = cols)+
  scale_x_continuous(name = "functional group richness")+
  scale_y_continuous(name = "Pl per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")  


#### Pr ~ funcdiv plot ####
ggplot(data = dat, aes(x= func.group, y=Pr_per100g))+
  #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
  #                                 fill=treatment), 
  #           alpha=0.2, show.legend=FALSE)+
  geom_jitter(data =datW1,
              width=0.2, shape=19, alpha=0.3, 
              aes(col=treatment))+
  geom_jitter(data=datW2,
              width=0.2, shape=18, alpha=0.9, 
              aes(col=treatment))+
  scale_color_manual(labels=treatments2, values = cols)+
  scale_x_continuous(name = "functional group richness")+
  scale_y_continuous(name = "Fu per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")  


#### Om ~ funcdiv plot ####
ggplot(data = dat, aes(x= func.group, y=Om_per100g))+
  #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
  #                                 fill=treatment), 
  #           alpha=0.2, show.legend=FALSE)+
  geom_jitter(data =datW1,
              width=0.2, shape=19, alpha=0.3, 
              aes(col=treatment))+
  geom_jitter(data=datW2,
              width=0.2, shape=18, alpha=0.9, 
              aes(col=treatment))+
  scale_color_manual(labels=treatments2, values = cols)+
  scale_x_continuous(name = "functional group richness")+
  scale_y_continuous(name = "Fu per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")  
