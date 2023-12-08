library(brms)
library(rstan)

# fitting trophic group densities ~ sowndiv

#data:
    #exclude 60 sp.:
    dat <- subset(dBEF_nem21, sowndiv != 60) 
    #standardize:  
    dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                          .after = sowndivLog) %>%
      mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
             .after = realdivLog)
    
    datW1 <- subset(dat, week=="W1")
    datW2 <- subset(dat, week=="W2")
    


####Ba~sowndiv: W1-d, W2-d , both- ####
  SEED = 22061996
  beta_coeff_priors <- prior(normal(0,20), class = "b")  
  sum(subset(dat, week=="W1")$Ba_per100g == 0) #9
  sum(subset(dat, week=="W2")$Ba_per100g == 0) #2
  

#for week 1:
  m.Ba_sowndivW1_p <- brm(
    bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ sowndivLogStd*treatment + (1|block/plot)),
    data = subset(dat, week=="W1"), 
    prior = beta_coeff_priors,
    family = hurdle_lognormal,
    stanvars = stanvars, #necessary to use custom brms families!
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #9div 
  
  summary(m.Ba_sowndivW1_p) 
  pp_check(m.Ba_sowndivW1_p, ndraws=100)+
    xlim(0,2000)

  
  #with default prior:
    m.Ba_sowndivW1_d <- brm(
      bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = subset(dat, week=="W1"), 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99))  #0 div
    
    summary(m.Ba_sowndiv_d)
      pp_check(m.Ba_sowndivW1_d, ndraws=100)+
        xlim(0,2000)
      
    #compare them  
      loo(m.Ba_sowndivW1_p, m.Ba_sowndivW1_d)
      print(rstan::get_elapsed_time(m.Ba_sowndivW1_p$fit))
      print(rstan::get_elapsed_time(m.Ba_sowndivW1_d$fit))
      

      
#for week 2:
m.Ba_sowndivW2_p <- brm(
  bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
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

summary(m.Ba_sowndivW2_p) 
pp_check(m.Ba_sowndivW2_p, ndraws=100)

  #with default priors
  m.Ba_sowndivW2_d <- brm(
    bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ 1),
    data = subset(dat, week=="W2"), 
    family = hurdle_lognormal,
    stanvars = stanvars, #necessary to use custom brms families!
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99))  #all good
  
  summary(m.Ba_sowndivW2_d)
  pp_check(m.Ba_sowndivW2_d, ndraws=100)
  
  #compare them  
  loo(m.Ba_sowndivW2_p, m.Ba_sowndivW2_d)
  print(rstan::get_elapsed_time(m.Ba_sowndivW2_p$fit))
  print(rstan::get_elapsed_time(m.Ba_sowndivW2_d$fit))
  loo(m.Ba_sowndivW2_d, m.Ba_sowndivW2_p)
  
#for both weeks  
  m.Ba_sowndiv_p <- brm(
    bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ sowndivLogStd*treatment + (1|block/plot)),
    data = dat, 
    prior = beta_coeff_priors,
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) #all good
  
  summary(m.Ba_sowndiv_p)
  pp_check(m.Ba_sowndiv_p, ndraws=100)+
    xlim(0,2000)
  
  #with default priors
  m.Ba_sowndiv_d <- brm(
    bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ sowndivLogStd*treatment + (1|block/plot)),
    data = dat, 
    family = hurdle_lognormal,
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) 
  
  summary(m.Ba_sowndiv_d)
  pp_check(m.Ba_sowndiv_d, ndraws=100)+
    xlim(0,2000)
  
  #compare them  
  loo(m.Ba_sowndiv_p, m.Ba_sowndiv_d)
  print(rstan::get_elapsed_time(m.Ba_sowndiv_p$fit))
  print(rstan::get_elapsed_time(m.Ba_sowndiv_d$fit))  



####Ba plot predictions for week1####
  pred.Ba_prior <- conditional_effects(m.Ba_sowndivW1_p)[[3]]
  pred.Ba_def <- conditional_effects(m.Ba_sowndivW1_d)[[3]]
  
  treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
  cols=c("brown2", "darkolivegreen", "dodgerblue3")
  BREAKS = unique(dat$sowndivLogStd)
  
  
  ggplot(subset(dat, week == "W1"), aes(x= sowndivLogStd, y=Ba_per100g))+
    #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
    #                                 fill=treatment), 
    #           alpha=0.2, show.legend=FALSE)+
    geom_jitter(width=0.2, shape=19, alpha=0.4, 
                aes(col=treatment))+
    geom_line(data=pred.Ba_prior, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype="dashed", col=treatment),
              linewidth= 1, show.legend = FALSE)+
    geom_line(data=pred.Ba_def, aes(x= sowndivLogStd, y=estimate__, 
                                    linetype="solid", col=treatment),
              linewidth= 1, show.legend = FALSE)+
    scale_color_manual(labels=treatments2, values = cols)+
    scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                       labels = c("1", "2", "4", "8", "16"))+
    scale_y_continuous(name = "Ba per 100g DW")+
    theme_bw()+
    theme(legend.position ="bottom")  
  
####Ba plot predictions for week2####
pred.Ba_prior <- conditional_effects(m.Ba_sowndivW2_p)[[3]]
pred.Ba_def <- conditional_effects(m.Ba_sowndivW2_d)[[3]]

treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")
BREAKS = unique(dat$sowndivLogStd)

  
    ggplot(subset(dat, week == "W2"), aes(x= sowndivLogStd, y=Ba_per100g))+
      #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
      #                                 fill=treatment), 
      #           alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.2, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=pred.Ba_prior, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype="dashed", col=treatment),
                linewidth= 1, show.legend = FALSE)+
      geom_line(data=pred.Ba_def, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype="solid", col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "Ba per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")
    
####Ba plot predictions for both weeks ####
    pred.Ba_prior <- conditional_effects(m.Ba_sowndiv_p)[[3]]
    pred.Ba_def <- conditional_effects(m.Ba_sowndiv_d)[[3]]
    
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    BREAKS = unique(dat$sowndivLogStd)
    
    
    ggplot(dat, aes(x= sowndivLogStd, y=Ba_per100g))+
      #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
      #                                 fill=treatment), 
      #           alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.2, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=pred.Ba_prior, aes(x= sowndivLogStd, y=estimate__, 
                                        linetype="dashed", col=treatment),
                linewidth= 1, show.legend = FALSE)+
      geom_line(data=pred.Ba_def, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype="solid", col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "Ba per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")    

#### Fu ~ sowndiv:W1-d, W2-d, both-d####
 
    SEED = 22061996
    beta_coeff_priors <- prior(normal(0,20), class = "b")  
    sum(subset(dat, week=="W1")$Fu_per100g == 0) #3
    sum(subset(dat, week=="W2")$Fu_per100g == 0) #1
    
    
    #for week 1:
    m.Fu_sowndivW1_p <- brm(
      bf(Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = subset(dat, week=="W1"), 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #11 div
    
    summary(m.Fu_sowndivW1_p)
    pp_check(m.Fu_sowndivW1_p, ndraws=100)+
      xlim(0,2000)
    
    #with default prior:
    m.Fu_sowndivW1_d <- brm(
      bf(Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = subset(dat, week=="W1"), 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED, 
      control = list(adapt_delta=0.99)) #1 div
    
    summary(m.Fu_sowndivW1_d)
    pp_check(m.Fu_sowndivW1_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Fu_sowndivW1_p, m.Fu_sowndivW1_d)
    print(rstan::get_elapsed_time(m.Fu_sowndivW1_p$fit))
    print(rstan::get_elapsed_time(m.Fu_sowndivW1_d$fit))
    
    
  #for week 2:
    m.Fu_sowndivW2_p <- brm(
      bf(Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = subset(dat, week=="W2"), 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #2 div
    
    summary(m.Fu_sowndivW2_p)
    pp_check(m.Fu_sowndivW2_p, ndraws=100)
    
    #with default priors
    m.Fu_sowndivW2_d <- brm(
      bf(Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = subset(dat, week=="W2"), 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #all good
    
    summary(m.Fu_sowndivW2_d)
    pp_check(m.Fu_sowndivW2_d, ndraws=100)
    
    #compare them  
    loo(m.Fu_sowndivW2_p, m.Fu_sowndivW2_d)
    print(rstan::get_elapsed_time(m.Fu_sowndivW2_p$fit))
    print(rstan::get_elapsed_time(m.Fu_sowndivW2_d$fit))
    
    #for both weeks  
    m.Fu_sowndiv_p <- brm(
      bf(Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99))
    
    summary(m.Fu_sowndiv_p, prob=0.9)
    pp_check(m.Fu_sowndiv_p, ndraws=100)+
      xlim(0,2000)
    
    #with default priors
    m.Fu_sowndiv_d <- brm(
      bf(Fu_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) 
    
    summary(m.Fu_sowndiv_d, prob=0.9)
    pp_check(m.Fu_sowndiv_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Fu_sowndiv_p, m.Fu_sowndiv_d)
    print(rstan::get_elapsed_time(m.Fu_sowndiv_p$fit))
    print(rstan::get_elapsed_time(m.Fu_sowndiv_d$fit))  
    
####Fu plot predictions ####
    pred.Fu_prior1 <- conditional_effects(m.Fu_sowndivW1_p)[[3]]
    pred.Fu_def1 <- conditional_effects(m.Fu_sowndivW1_d)[[3]]
    summary(m.Fu_sowndivW1_p)
    
    pred.Fu_prior2 <- conditional_effects(m.Fu_sowndivW2_p)[[3]]
    pred.Fu_def2 <- conditional_effects(m.Fu_sowndivW2_d)[[3]]
    summary(m.Fu_sowndivW2_p)
    
    pred.Fu_prior  <- conditional_effects(m.Fu_sowndiv_p)[[3]]
    pred.Fu_def  <- conditional_effects(m.Fu_sowndiv_d)[[3]]
    summary(m.Fu_sowndiv_d)
    
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    BREAKS = unique(dat$sowndivLogStd)
    
    
    ggplot(data = dat, aes(x= sowndivLogStd, y=Fu_per100g))+
      #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
      #                                 fill=treatment), 
      #           alpha=0.2, show.legend=FALSE)+
      geom_jitter(data =datW1,
                  width=0.2, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_jitter(data=datW2,
                  width=0.2, shape=18, alpha=0.4, 
                  aes(col=treatment))+
      #predictions week1
      geom_line(data=pred.Fu_prior1, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype="dashed", col=treatment),
              linewidth= 0.5, show.legend = FALSE)+
      geom_line(data=pred.Fu_def1, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype="solid", col=treatment),
              linewidth= 0.5, show.legend = FALSE)+
      #predictions week 2:
      geom_line(data=pred.Fu_prior2, aes(x= sowndivLogStd, y=estimate__, 
                                         linetype="dashed", col=treatment),
                linewidth= 0.5, show.legend = FALSE)+
      geom_line(data=pred.Fu_def2, aes(x= sowndivLogStd, y=estimate__, 
                                       linetype="solid", col=treatment),
                linewidth=0.5, show.legend = FALSE)+
      #models for both weeks:
      geom_line(data=pred.Fu_prior, aes(x= sowndivLogStd, y=estimate__, 
                                         linetype="dashed", col=treatment),
                linewidth= 0.7, show.legend = FALSE)+
      geom_line(data=pred.Fu_def, aes(x= sowndivLogStd, y=estimate__, 
                                       linetype="solid", col=treatment),
                linewidth=0.7, show.legend = FALSE)+
      
    scale_color_manual(labels=treatments2, values = cols)+
    scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                       labels = c("16", "8", "4", "2", "1"))+
    scale_y_continuous(name = "Fu per 100g DW")+
    theme_bw()+
    theme(legend.position ="bottom")  
    
    plot(x = subset(dat, week=="W1")$sowndivLogStd, y = subset(dat, week=="W1")$Fu_per100g )

    
    
#### Pl ~ sowndiv: W1-d, W2-d, both-p ####
    
    SEED = 22061996
    beta_coeff_priors <- prior(normal(0,20), class = "b")  
    sum(subset(dat, week=="W1")$Pl_per100g == 0) #2
    sum(subset(dat, week=="W2")$Pl_per100g == 0) #0
    
#for week 1:
    m.Pl_sowndivW1_p <- brm(
      bf(Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = subset(dat, week=="W1"), 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #30 div
    
    summary(m.Pl_sowndivW1_p, prob=0.9)
    pp_check(m.Pl_sowndivW1_p, ndraws=100)+
      xlim(0,2000)
    
    #with default prior:
    m.Pl_sowndivW1_d <- brm(
      bf(Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = subset(dat, week=="W1"), 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #8 div

    summary(m.Pl_sowndivW1_d, prob=0.9)
    pp_check(m.Pl_sowndivW1_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Pl_sowndivW1_p, m.Pl_sowndivW1_d)
    print(rstan::get_elapsed_time(m.Pl_sowndivW1_p$fit))
    print(rstan::get_elapsed_time(m.Pl_sowndivW1_d$fit))
    
    
  #for week 2:
    m.Pl_sowndivW2_p <- brm(
      bf(Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = subset(dat, week=="W2"), 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #8 div
    
    summary(m.Pl_sowndivW2_p, prob=0.9)
    pp_check(m.Pl_sowndivW2_p, ndraws=100)
    
    #with default priors
    m.Pl_sowndivW2_d <- brm(
      bf(Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = subset(dat, week=="W2"), 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED, 
      control = list(adapt_delta=0.99))  #5 div
    
    summary(m.Pl_sowndivW2_d, prob=0.9)
    pp_check(m.Pl_sowndivW2_d, ndraws=100)
    
    #compare them  
    loo(m.Pl_sowndivW2_p, m.Pl_sowndivW2_d)
    print(rstan::get_elapsed_time(m.Pl_sowndivW2_p$fit))
    print(rstan::get_elapsed_time(m.Pl_sowndivW2_d$fit))
    loo(m.Pl_sowndivW2_d, m.Pl_sowndivW2_p)
    
    #for both weeks  
    m.Pl_sowndiv_p <- brm(
      bf(Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED, 
      control = list(adapt_delta=0.99)) #all good
    
    summary(m.Pl_sowndiv_p, prob=0.9)
    pp_check(m.Pl_sowndiv_p, ndraws=100)+
      xlim(0,2000)
    
    #with default priors
    m.Pl_sowndiv_d <- brm(
      bf(Pl_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ 1),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) 
    
    
    summary(m.Pl_sowndiv_d, prob=0.9)
    pp_check(m.Pl_sowndiv_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Pl_sowndiv_p, m.Pl_sowndiv_d)
    print(rstan::get_elapsed_time(m.Pl_sowndiv_p$fit))
    print(rstan::get_elapsed_time(m.Pl_sowndiv_d$fit))  
    

#### Pr ~ sowndiv: W1-p, W2-p, both-p  ####
    #in W1: p has less divergent transitions, but slightly wors ELPD
    
    SEED = 22061996
    beta_coeff_priors <- prior(normal(0,20), class = "b")  
    sum(subset(dat, week=="W1")$Pr_per100g == 0) #44
    sum(subset(dat, week=="W2")$Pr_per100g == 0) #13
    
    #for week 1:
    m.Pr_sowndivW1_p <- brm(
      bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = subset(dat, week=="W1"), 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #3 div
    
    summary(m.Pr_sowndivW1_p, prob=0.9)
    pp_check(m.Pr_sowndivW1_p, ndraws=100)+
      xlim(0,2000)
    
    #with default prior:
    m.Pr_sowndivW1_d <- brm(
      bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = subset(dat, week=="W1"), 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #42 div
    
    summary(m.Pr_sowndivW1_d, prob=0.9)
    pp_check(m.Pr_sowndivW1_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Pr_sowndivW1_p, m.Pr_sowndivW1_d)
    print(rstan::get_elapsed_time(m.Pr_sowndivW1_p$fit))
    print(rstan::get_elapsed_time(m.Pr_sowndivW1_d$fit))
    
    
  #for week 2:
    m.Pr_sowndivW2_p <- brm(
      bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = subset(dat, week=="W2"), 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #6 div trans
    
    summary(m.Pr_sowndivW2_p, prob=0.9)
    pp_check(m.Pr_sowndivW2_p, ndraws=100)
    
    #with default priors
    m.Pr_sowndivW2_d <- brm(
      bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = subset(dat, week=="W2"), 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99))  #9 div trans
    
    summary(m.Pr_sowndivW2_d, prob=0.9)
    pp_check(m.Pr_sowndivW2_d, ndraws=100)
    
    #compare them  
    loo(m.Pr_sowndivW2_p, m.Pr_sowndivW2_d)
    print(rstan::get_elapsed_time(m.Pr_sowndivW2_p$fit))
    print(rstan::get_elapsed_time(m.Pr_sowndivW2_d$fit))
    loo(m.Pr_sowndivW2_d, m.Pr_sowndivW2_p)
    
  #for both weeks  
    m.Pr_sowndiv_p <- brm(
      bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #all good
    
    summary(m.Pr_sowndiv_p, prob=0.9)
    pp_check(m.Pr_sowndiv_p, ndraws=100)+
      xlim(0,2000)
    
    #with default priors
    m.Pr_sowndiv_d <- brm(
      bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #1 div trans
    
    summary(m.Pr_sowndiv_d, prob=0.9)
    pp_check(m.Pr_sowndiv_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Pr_sowndiv_p, m.Pr_sowndiv_d)
    print(rstan::get_elapsed_time(m.Pr_sowndiv_p$fit))
    print(rstan::get_elapsed_time(m.Pr_sowndiv_d$fit))  
    
#### Om ~ sowndiv: W1-p, W2-p, both-d ####
    #W1: d has better elpd, but p less divergent transitions
    #W2: 11 vs 12 divergent transitions, elpd basically same
    #both: 10 div in p, zero in d, elpd slightly better in p
    
    SEED = 22061996
    beta_coeff_priors <- prior(normal(0,20), class = "b")  
    sum(subset(dat, week=="W1")$Om_per100g == 0) #89
    sum(subset(dat, week=="W2")$Om_per100g == 0) #42
    
    #for week 1:
    m.Om_sowndivW1_p <- brm(
      bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = subset(dat, week=="W1"), 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99))  #15 div
    
    summary(m.Om_sowndivW1_p, prob=0.9)
    pp_check(m.Om_sowndivW1_p, ndraws=100)+
      xlim(0,2000)
    
    #with default prior:
    m.Om_sowndivW1_d <- brm(
      bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = subset(dat, week=="W1"), 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #25 div
    
    summary(m.Om_sowndivW1_d, prob=0.9)
    pp_check(m.Om_sowndivW1_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Om_sowndivW1_p, m.Om_sowndivW1_d)
    print(rstan::get_elapsed_time(m.Om_sowndivW1_p$fit))
    print(rstan::get_elapsed_time(m.Om_sowndivW1_d$fit))
    
    
    #for week 2:
    m.Om_sowndivW2_p <- brm(
      bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = subset(dat, week=="W2"), 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #11 div
    
    summary(m.Om_sowndivW2_p, prob=0.9)
    pp_check(m.Om_sowndivW2_p, ndraws=100)
    
    #with default priors
    m.Om_sowndivW2_d <- brm(
      bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = subset(dat, week=="W2"), 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #12 div 
    
    summary(m.Om_sowndivW2_d, prob=0.9)
    pp_check(m.Om_sowndivW2_d, ndraws=100)
    
    #compare them  
    loo(m.Om_sowndivW2_p, m.Om_sowndivW2_d)
    print(rstan::get_elapsed_time(m.Om_sowndivW2_p$fit))
    print(rstan::get_elapsed_time(m.Om_sowndivW2_d$fit))
    loo(m.Om_sowndivW2_d, m.Om_sowndivW2_p)
    
    #for both weeks  
    m.Om_sowndiv_p <- brm(
      bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #10 div
    
    summary(m.Om_sowndiv_p, prob=0.9)
    pp_check(m.Om_sowndiv_p, ndraws=100)+
      xlim(0,220)
    
    #with default priors
    m.Om_sowndiv_d <- brm(
      bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #all good
    
    summary(m.Om_sowndiv_d, prob=0.9)
    pp_check(m.Om_sowndiv_d, ndraws=100)+
      xlim(0,220)
    
    #compare them  
    loo(m.Om_sowndiv_p, m.Om_sowndiv_d)
    print(rstan::get_elapsed_time(m.Om_sowndiv_p$fit))
    print(rstan::get_elapsed_time(m.Om_sowndiv_d$fit))  
    
        
    
#### save models ####
    
    save(m.Om_sowndivW1_d, m.Om_sowndivW1_p,
         m.Om_sowndivW2_d, m.Om_sowndivW2_p,
         m.Om_sowndiv_d, m.Om_sowndiv_p,
         file = "./statistics/brms/231207_Om_sowndiv_priors.RData")
    
    
    save(m.Pr_sowndivW1_d, m.Pr_sowndivW1_p,
         m.Pr_sowndivW2_d, m.Pr_sowndivW2_p,
         m.Pr_sowndiv_d, m.Pr_sowndiv_p,
         file = "./statistics/brms/231207_Pr_sowndiv_priors.RData")
    
    
     save(m.Pl_sowndivW1_d, m.Pl_sowndivW1_p,
         m.Pl_sowndivW2_d, m.Pl_sowndivW2_p,
         m.Pl_sowndiv_d, m.Pl_sowndiv_p,
         file = "./statistics/brms/231207_Pl_sowndiv_priors.RData")
    
    
    save(m.Fu_sowndivW1_d, m.Fu_sowndivW1_p,
         m.Fu_sowndivW2_d, m.Fu_sowndivW2_p,
         m.Fu_sowndiv_d, m.Fu_sowndiv_p,
         file = "./statistics/brms/231207_Fu_sowndiv_priors.RData")
    
    save(m.Ba_sowndivW1_d, m.Ba_sowndivW1_p,
         m.Ba_sowndivW2_d, m.Ba_sowndivW2_p,
         m.Ba_sowndiv_d, m.Ba_sowndiv_p,
         file = "./statistics/brms/231207_Ba_sowndiv_priors.RData")

load(file = "./statistics/brms/231206_troph_sowndiv_priors.RData")
