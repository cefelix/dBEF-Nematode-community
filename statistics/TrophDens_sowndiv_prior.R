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



####Ba~sowndiv####
SEED = 22061996
beta_coeff_priors <- prior(normal(0,20), class = "b")  
      
#for week 1:
  m.Ba_sowndivW1_p <- brm(
    bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
       hu ~ sowndivLogStd*treatment + (1|block/plot)),
    data = dat[dat$week=="W1"], 
    prior = beta_coeff_priors,
    family = hurdle_lognormal,
    stanvars = stanvars, #necessary to use custom brms families!
    chains = 3,
    cores = 3,
    iter = 2000, warmup = 1000,
    seed = SEED,
    control = list(adapt_delta=0.99)) 
  
  pp_check(m.Ba_sowndivW1_p, ndraws=100)+
    xlim(0,2000)

  #with default prior:
    m.Ba_sowndivW1_d <- brm(
      bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
         hu ~ sowndivLogStd*treatment + (1|block/plot)),
      data = dat[dat$week=="W1"], 
      family = hurdle_lognormal,
      stanvars = stanvars, #necessary to use custom brms families!
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) 
    
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
  control = list(adapt_delta=0.99))
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
    control = list(adapt_delta=0.99)) 
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
    control = list(adapt_delta=0.99))
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

    #### save models ####

    save(m.Ba_sowndivW1_d, m.Ba_sowndivW1_p,
         m.Ba_sowndivW2_d, m.Ba_sowndivW2_p,
         m.Ba_sowndiv_d, m.Ba_sowndiv_p,
         file = "./statistics/brms/231206_troph_sowndiv_priors.RData")

load(file = "./statistics/brms/231206_troph_sowndiv_priors.RData")
