library(brms)
library(rstan)
library(ggplot2)
library(emmeans)
library(bayestestR)

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
    
    #a seed:
    SEED = 22061996
    
#priors    
    beta_coeff_priors <- prior(normal(0,20), class = "b")  
    beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
    beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
    
    #the narrowest prior is still basically flat (at realistic values): 
    ggplot(data.frame(density(rnorm(1e5, 0, 2)) ), 
           aes(x=x, y=y))+
      geom_point()+
      xlim(-.5,.5)
    
    exp(0.5) 
    #this would mean that increasing plant diversity by 1 SD would lead to 64% more individuals per dry weight of soil
    
####Ba~sowndiv: W1-d4, W2-p , both- ####
    
    
    #W2: _p has 1 div, _d has 4 div --> despite slightly better elpd (less than 2 SE difference), choose _p 
    load(file = "./statistics/brms/231213_Ba_sowndiv_priors.RData")
    
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
    control = list(adapt_delta=0.99)) #2 div 
  summary(m.Ba_sowndivW1_p, prob=0.9) 
  
    #removing hu-interaction:
    m.Ba_sowndivW1_p2 <- update(m.Ba_sowndivW1_p,
                                bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                   hu ~ sowndivLogStd + treatment + (1|block/plot)),
                                seed = SEED) #2 div
    summary(m.Ba_sowndivW1_p2, prob=0.9) 
    
    #removing hu-sowndivLogStd:
    m.Ba_sowndivW1_p3 <- update(m.Ba_sowndivW1_p2,
                                bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                   hu ~  treatment + (1|block/plot)),
                                seed = SEED) # 0 div
    summary(m.Ba_sowndivW1_p3, prob=0.9) 
    
    #removing hu-treatment:
    m.Ba_sowndivW1_p4 <- update(m.Ba_sowndivW1_p3,
                                bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                   hu ~  1),
                                seed = SEED)
    summary(m.Ba_sowndivW1_p4, prob=0.9) #12 div
  
      #using a narrower prior:
      beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
      m.Ba_sowndivW1_p5 <- update(m.Ba_sowndivW1_p4,
                                  prior = beta_coeff_priors2,
                                  seed = SEED) #8 div
      summary(m.Ba_sowndivW1_p5, prob=0.9) 
      
      #using a even narrower prior:
      beta_coeff_priors <- prior(normal(0,2), class = "b")  
      m.Ba_sowndivW1_p6 <- update(m.Ba_sowndivW1_p5,
                                  prior = beta_coeff_priors3,
                                  seed = SEED) #11 div
      summary(m.Ba_sowndivW1_p6, prob=0.9) 
  
      #  _p5 best
  summary(m.Ba_sowndivW1_p5, prob=0.9) 
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
      control = list(adapt_delta=0.99))  #1 div
    summary(m.Ba_sowndivW1_d, prob=0.9) 
    
    #removing hu-interaction
    m.Ba_sowndivW1_d2 <- update(m.Ba_sowndivW1_d,
                                bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                   hu ~ sowndivLogStd + treatment + (1|block/plot)),
                                seed = SEED) #2 div
    summary(m.Ba_sowndivW1_d2, prob=0.9) 
    
    #removing hu-sowndivLogStd
    m.Ba_sowndivW1_d3 <- update(m.Ba_sowndivW1_d2,
                                bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                   hu ~ treatment + (1|block/plot)),
                                seed = SEED) # 4div
    summary(m.Ba_sowndivW1_d3, prob=0.9) 
    
    
    #removing hu-sowndivLogStd
    m.Ba_sowndivW1_d4 <- update(m.Ba_sowndivW1_d3,
                                bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                   hu ~ 1),
                                seed = SEED) # 7div
    summary(m.Ba_sowndivW1_d4, prob=0.9) 
    
    #_d4 best
    
    summary(m.Ba_sowndivW1_d4, prob=0.9)
      pp_check(m.Ba_sowndivW1_d4, ndraws=100)+
        xlim(0,2000)
      
    #compare them  
      loo(m.Ba_sowndivW1_p5, m.Ba_sowndivW1_d4)
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
  control = list(adapt_delta=0.99)) #1 div
  summary(m.Ba_sowndivW2_p, prob=0.9) 


    #using a narrower prior:
    beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
    m.Ba_sowndivW2_p2 <- update(m.Ba_sowndivW1_p,
                                prior = beta_coeff_priors2,
                                seed = SEED) # 3 div
    summary(m.Ba_sowndivW2_p2, prob=0.9)
    
    #using a even narrower prior:
    beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
    m.Ba_sowndivW2_p3 <- update(m.Ba_sowndivW1_p2,
                                prior = beta_coeff_priors3,
                                seed = SEED) # 6div
    summary(m.Ba_sowndivW2_p3, prob=0.9)
    
    #best _p

summary(m.Ba_sowndivW2_p, prob=0.9) 
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
    control = list(adapt_delta=0.99))  #4 div
  summary(m.Ba_sowndivW2_d, prob=0.9)
  
  summary(m.Ba_sowndivW2_d)
  pp_check(m.Ba_sowndivW2_d, ndraws=100)
  
  #compare them  
  loo(m.Ba_sowndivW2_p, m.Ba_sowndivW2_d)
  print(rstan::get_elapsed_time(m.Ba_sowndivW2_p$fit))
  print(rstan::get_elapsed_time(m.Ba_sowndivW2_d$fit))
  loo(m.Ba_sowndivW2_d, m.Ba_sowndivW2_p)

####Ba plot predictions ####
  # no substantial differences in slope and intercept --> don't include week into model
  
  #pred.Ba_prior1 <- conditional_effects(m.Ba_sowndivW1_p)[[3]]
  pred.Ba_def1 <- conditional_effects(m.Ba_sowndivW1_d4)[[3]]
  summary(m.Ba_sowndivW1_d4)
  
  pred.Ba_prior2 <- conditional_effects(m.Ba_sowndivW2_p)[[3]]
  #pred.Ba_def2 <- conditional_effects(m.Ba_sowndivW2_d)[[3]]
  summary(m.Ba_sowndivW2_p)
  
  pred.Ba_prior  <- conditional_effects(m.Ba_sowndiv_p)[[3]]
  pred.Ba_def  <- conditional_effects(m.Ba_sowndiv_d)[[3]]
  summary(m.Ba_sowndiv_d, prob=0.9)
  
  treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
  cols=c("brown2", "darkolivegreen", "dodgerblue3")
  BREAKS = unique(dat$sowndivLogStd)
  
  
  ggplot(data = dat, aes(x= sowndivLogStd, y=Ba_per100g))+
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
    geom_line(data=pred.Ba_def1, aes(x= sowndivLogStd, y=estimate__, 
                                     linetype="solid", col=treatment),
              linewidth= 0.5, linetype=3,
              show.legend = FALSE)+
    #predictions week 2:
    geom_line(data=pred.Ba_prior2, aes(x= sowndivLogStd, y=estimate__, 
                                       col=treatment),
              linewidth= 0.5, linetype=6,
              show.legend = FALSE)+
    #models for both weeks:
    #geom_line(data=pred.Ba_prior, aes(x= sowndivLogStd, y=estimate__, 
    #                                   col=treatment),
    #          linewidth= 0.7, linetype="solid",
    #          show.legend = FALSE)+
    scale_color_manual(labels=treatments2, values = cols)+
    scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                       labels = c("16", "8", "4", "2", "1"))+
    scale_y_continuous(name = "Ba per 100g DW", limits = c(-1, 500))+
    theme_bw()+
    theme(legend.position ="bottom")  
  
  #no substantial difference in slope or intercept --> dont include week into model!

####Ba ~ sowndiv, both weeks: _p5  ####
  #best is _p5, as it is the only one with 0 divergent transitions and elpd doesnt differ by more than 2 SE
  
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
    control = list(adapt_delta=0.99))  #0 div
  summary(m.Ba_sowndiv_d)
  
    #remove hurdle interaction
    m.Ba_sowndiv_d2 <- update(m.Ba_sowndiv_d,
                              bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                 hu ~ sowndivLogStd + treatment + (1|block/plot))) #0 div
    summary(m.Ba_sowndiv_d2, prob=0.9)
    
    #remove hurdle ~ treatment
    m.Ba_sowndiv_d3 <- update(m.Ba_sowndiv_d2,
                              bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                 hu ~ sowndivLogStd + (1|block/plot))) #0 div
    summary(m.Ba_sowndiv_d3, prob=0.9)
    
    #remove hurdle ~ sowndivLogStd
    m.Ba_sowndiv_d4 <- update(m.Ba_sowndiv_d3,
                              bf(Ba_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                 hu ~ 1)) #2 div
    summary(m.Ba_sowndiv_d4, prob=0.9)
    
    #add a normal(0,20) prior:
    beta_coeff_priors <- prior(normal(0,20), class = "b")  
    m.Ba_sowndiv_p5 <- update(m.Ba_sowndiv_d4,
                                prior = beta_coeff_priors,
                                seed = SEED)  #0 div
    summary(m.Ba_sowndiv_p5, prob=0.9)
    
    #add a normal(0,5) prior:
    beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
    m.Ba_sowndiv_p6 <- update(m.Ba_sowndiv_p5,
                              prior = beta_coeff_priors2,
                              seed = SEED) 
    summary(m.Ba_sowndiv_p6, prob=0.9) #1 div
    
    #add a normal(0,2) prior:
    beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
    m.Ba_sowndiv_p7 <- update(m.Ba_sowndiv_p6,
                              prior = beta_coeff_priors2,
                              seed = SEED)
    summary(m.Ba_sowndiv_p7, prob=0.9) #1 div
    pp_check(m.Ba_sowndiv_p7, ndraws=100)+
      xlim(0,2000)
    
    loo(m.Ba_sowndiv_p5,m.Ba_sowndiv_p7)
    #best _p7
  
  summary(m.Ba_sowndiv_p5, prob=0.9)
  pp_check(m.Ba_sowndiv_p5, ndraws=100)+
    xlim(0,1000)
  
  #compare them  
  loo(m.Ba_sowndiv_p5, m.Ba_sowndiv_p6, m.Ba_sowndiv_p7, m.Ba_sowndiv_d4) 
  #best is _p5, as it has the highest elpd at 0 divergent transitions 
  print(rstan::get_elapsed_time(m.Ba_sowndiv_p$fit))
  print(rstan::get_elapsed_time(m.Ba_sowndiv_d$fit))  
  
  emt = emtrends(m.Ba_sowndiv_p5, specs = c("treatment"), var="sowndivLogStd")
  summary(emt, point.est=mean, level = .9) 
  emt.pairs <- pairs(emt)
  summary(emt.pairs, point.est=mean, level = .9)
  bayestestR::p_direction(emt.pairs)
  p_direction.brmsfit(m.Ba_sowndiv_d4)
  
#save them:
  save(m.Ba_sowndivW1_d, m.Ba_sowndivW1_p,
        m.Ba_sowndivW1_p2, m.Ba_sowndivW1_p3, m.Ba_sowndivW1_p4, m.Ba_sowndivW1_p5, m.Ba_sowndivW1_p6,
        m.Ba_sowndivW1_d2, m.Ba_sowndivW1_d3, m.Ba_sowndivW1_d4,
       m.Ba_sowndivW2_d, m.Ba_sowndivW2_p,
        m.Ba_sowndivW2_p2, m.Ba_sowndivW2_p3,
       m.Ba_sowndiv_d, m.Ba_sowndiv_p,
        m.Ba_sowndiv_d2, m.Ba_sowndiv_d3, m.Ba_sowndiv_d4,
        m.Ba_sowndiv_p5, m.Ba_sowndiv_p6, m.Ba_sowndiv_p7,
       file = "./statistics/brms/231213_Ba_sowndiv_priors.RData")
  
#remove from workspace to prevent crashing
  rm(m.Ba_sowndivW1_d, m.Ba_sowndivW1_p,
    m.Ba_sowndivW1_p2, m.Ba_sowndivW1_p3, m.Ba_sowndivW1_p4, m.Ba_sowndivW1_p5, m.Ba_sowndivW1_p6,
    m.Ba_sowndivW1_d2, m.Ba_sowndivW1_d3, m.Ba_sowndivW1_d4,
    m.Ba_sowndivW2_d, m.Ba_sowndivW2_p,
    m.Ba_sowndivW2_p2, m.Ba_sowndivW2_p3,
    m.Ba_sowndiv_d, m.Ba_sowndiv_p,
    m.Ba_sowndiv_d2, m.Ba_sowndiv_d3, m.Ba_sowndiv_d4, #keeping the final model loaded
    #m.Ba_sowndiv_p5, #keeping the final model loaded
    m.Ba_sowndiv_p6, m.Ba_sowndiv_p7)
  
  

#### Fu ~ sowndiv:W1-p2, W2-d, both- ####
    #W1: _p2: 1 div, _d 7 div, elpd for p2 lower --> choose p2 
    #W2: _p: 1 div, _d 0 div, elpd for p slightly lower (less than 2 se) --> choose d
    load(file = "./statistics/brms/231213_Fu_sowndiv_priors.RData")
  
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
      control = list(adapt_delta=0.99)) #8 div
    summary(m.Fu_sowndivW1_p, prob=0.9)
    
        #using a narrower prior:
        beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
        m.Fu_sowndivW1_p2 <- update(m.Fu_sowndivW1_p,
                                    prior = beta_coeff_priors2,
                                    seed = SEED) #1 div
        summary(m.Fu_sowndivW1_p2, prob=0.9)
        
        #using a even narrower prior:
        beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
        m.Fu_sowndivW1_p3 <- update(m.Fu_sowndivW1_p2,
                                    prior = beta_coeff_priors3,
                                    seed = SEED) #4 div
        summary(m.Fu_sowndivW1_p3, prob=0.9)
    
    #best _p2    
    summary(m.Fu_sowndivW1_p2, prob=0.9)
    pp_check(m.Fu_sowndivW1_p2, ndraws=100)+
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
      control = list(adapt_delta=0.99)) #7 div
    summary(m.Fu_sowndivW1_d, prob=0.9)
    
    summary(m.Fu_sowndivW1_d, prob=0.9)
    pp_check(m.Fu_sowndivW1_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Fu_sowndivW1_p2, m.Fu_sowndivW1_d) #choose _p2
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
      control = list(adapt_delta=0.99)) #1 div
    summary(m.Fu_sowndivW2_p, prob=0.9)
    
        #using a narrower prior:
        beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
        m.Fu_sowndivW2_p2 <- update(m.Fu_sowndivW2_p,
                                    prior = beta_coeff_priors2,
                                    seed = SEED) # 5div
        summary(m.Fu_sowndivW2_p2, prob=0.9)
        
        #using a even narrower prior:
        beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
        m.Fu_sowndivW2_p3 <- update(m.Fu_sowndivW2_p2,
                                    prior = beta_coeff_priors3,
                                    seed = SEED) # 2 div
        summary(m.Fu_sowndivW2_p3, prob=0.9)
        
    #best _p
    summary(m.Fu_sowndivW2_p)
    pp_check(m.Fu_sowndivW2_p3, ndraws=100)
    
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
    summary(m.Fu_sowndivW2_d, prob=0.9)
    
    summary(m.Fu_sowndivW2_d)
    pp_check(m.Fu_sowndivW2_d, ndraws=100)
    
    #compare them  
    loo(m.Fu_sowndivW2_p, m.Fu_sowndivW2_d)
    print(rstan::get_elapsed_time(m.Fu_sowndivW2_p$fit))
    print(rstan::get_elapsed_time(m.Fu_sowndivW2_d$fit))
    
    
####Fu plot predictions ####
    #strong difference in intercept/proportion of zeros between week 1 and week 2 --> 
    #include week as additive term in hurdle and non-hurdle part of the model
    
    pred.Fu_prior1 <- conditional_effects(m.Fu_sowndivW1_p2)[[3]]
    #pred.Fu_def1 <- conditional_effects(m.Fu_sowndivW1_d)[[3]]
    summary(m.Fu_sowndivW1_p)
    
    #pred.Fu_prior2 <- conditional_effects(m.Fu_sowndivW2_p)[[3]]
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
                                       col=treatment),
              linewidth= 0.5, linetype=1,
              show.legend = FALSE)+
      #predictions week 2:
      geom_line(data=pred.Fu_def2, aes(x= sowndivLogStd, y=estimate__, 
                                       col=treatment),
                linewidth=0.5, linetype=2, 
                show.legend = FALSE)+
    scale_color_manual(labels=treatments2, values = cols)+
    scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                       labels = c("16", "8", "4", "2", "1"))+
    scale_y_continuous(name = "Fu per 100g DW")+
    theme_bw()+
    theme(legend.position ="bottom")  
    
    #strong difference in intercept --> 
    #include week as additive term in lognormal part of the model
    plot(x = subset(dat, week=="W1")$sowndivLogStd, y = subset(dat, week=="W1")$Fu_per100g )

    
    
  #### Fu ~ sowndivLogStd, both weeks: _p3 ####
    #_p3 2 div, _d 2 div, elpd in _p3 slightly better but less than 2 SE --> choose _p3
    
    #for both weeks  
    m.Fu_sowndiv_p <- brm(
      bf(Fu_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
         hu ~ 1),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99))  #2 div
    summary(m.Fu_sowndiv_p, prob=0.9)
    
    #add a normal(0,5) prior:
    beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
    m.Fu_sowndiv_p2 <- update(m.Fu_sowndiv_p,
                              prior = beta_coeff_priors2, #all good
                              seed = SEED) #0 div
    summary(m.Fu_sowndiv_p2, prob=0.9)
    
    #add a normal(0,2) prior:
    beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
    m.Fu_sowndiv_p3 <- update(m.Fu_sowndiv_p2,
                              prior = beta_coeff_priors3,
                              seed = SEED)
    summary(m.Fu_sowndiv_p3, prob=0.9) #2 div trans
    
    loo(m.Fu_sowndiv_p, m.Fu_sowndiv_p2) #same, use _p as it has a wider prior
    
    
    summary(m.Fu_sowndiv_p, prob=0.9)
    pp_check(m.Fu_sowndiv_p, ndraws=100)+
      xlim(0,2000)
    
    #with default priors
    m.Fu_sowndiv_d <- brm(
      bf(Fu_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
         hu ~ 1),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #2 div
    summary(m.Fu_sowndiv_d, prob=0.9)
    
    summary(m.Fu_sowndiv_d, prob=0.9)
    pp_check(m.Fu_sowndiv_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Fu_sowndiv_p, m.Fu_sowndiv_p2, m.Fu_sowndiv_p3, m.Fu_sowndiv_d)  
      #--> best is _p3, as _p2 without any divergent transitions does have significantly worse elpd (>2 SE)
    print(rstan::get_elapsed_time(m.Fu_sowndiv_p$fit))
    print(rstan::get_elapsed_time(m.Fu_sowndiv_d$fit))  
    
#save models:
    save(m.Fu_sowndivW1_d, m.Fu_sowndivW1_p,
          m.Fu_sowndivW1_p2, m.Fu_sowndivW1_p3,
         m.Fu_sowndivW2_d, m.Fu_sowndivW2_p,
          m.Fu_sowndivW2_p2, m.Fu_sowndivW2_p3,
         m.Fu_sowndiv_d, m.Fu_sowndiv_p,
          m.Fu_sowndiv_p2, m.Fu_sowndiv_p3,
         file = "./statistics/brms/231213_Fu_sowndiv_priors.RData")
    
#remove models to prevent crashes:
    rm(m.Fu_sowndivW1_d, m.Fu_sowndivW1_p,
      m.Fu_sowndivW1_p2, m.Fu_sowndivW1_p3,
      m.Fu_sowndivW2_d, m.Fu_sowndivW2_p,
      m.Fu_sowndivW2_p2, m.Fu_sowndivW2_p3,
      m.Fu_sowndiv_d, m.Fu_sowndiv_p,
      m.Fu_sowndiv_p2#, m.Fu_sowndiv_p3
      )
    
    
#### Pl ~ sowndiv: W1-p2, W2-p, both-p ####
    #W1: _p2: 1 div, _d 19 div, elpd for _p2 equal --> choose _p2
    #W2: _p: 2 div, _d 11 div, elpd for _p slightly worse (but less than 2 SE) --> choose _p
    load(file = "./statistics/brms/231213_Pl_sowndiv_priors.RData")
    
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
      control = list(adapt_delta=0.99)) #6 div
      
      summary(m.Pl_sowndivW1_p, prob=0.9)
      
      #using a narrower prior:
      beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
      m.Pl_sowndivW1_p2 <- update(m.Pl_sowndivW1_p,
                                 prior = beta_coeff_priors2,
                                 seed = SEED) #1 div
      summary(m.Pl_sowndivW1_p2, prob=0.9)
      
      #using a even narrower prior:
      beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
      m.Pl_sowndivW1_p3 <- update(m.Pl_sowndivW1_p2,
                                  prior = beta_coeff_priors3,
                                  seed = SEED)
      summary(m.Pl_sowndivW1_p3, prob=0.9) #2 div
        
    #best _p2
    summary(m.Pl_sowndivW1_p2, prob=0.9)
    pp_check(m.Pl_sowndivW1_p, ndraws=100)+
      xlim(0,2500)
    
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
      control = list(adapt_delta=0.99)) #19 div

    summary(m.Pl_sowndivW1_d, prob=0.9)
    pp_check(m.Pl_sowndivW1_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Pl_sowndivW1_p, m.Pl_sowndivW1_p2, m.Pl_sowndivW1_p3, m.Pl_sowndivW1_d)
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
      control = list(adapt_delta=0.99)) #2 div
      summary(m.Pl_sowndivW2_p, prob=0.9)
    
        #using a narrower prior:
        beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
        m.Pl_sowndivW2_p2 <- update(m.Pl_sowndivW2_p,
                                    prior = beta_coeff_priors2,
                                    seed = SEED) #11 div
        summary(m.Pl_sowndivW2_p2, prob=0.9)
        
        #using a even narrower prior:
        beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
        m.Pl_sowndivW2_p3 <- update(m.Pl_sowndivW2_p2,
                                    prior = beta_coeff_priors3,
                                    seed = SEED) #14 div
        summary(m.Pl_sowndivW2_p3, prob=0.9)
    
        #best: _p
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
      control = list(adapt_delta=0.99))  #11 div
    
    summary(m.Pl_sowndivW2_d, prob=0.9)
    pp_check(m.Pl_sowndivW2_d, ndraws=100)
    
    #compare them  
    loo(m.Pl_sowndivW2_p, m.Pl_sowndivW2_p2, m.Pl_sowndivW2_p3, m.Pl_sowndivW2_d) 
      #elpd of _d slightly worse, but less than 2 SE
    print(rstan::get_elapsed_time(m.Pl_sowndivW2_p$fit))
    print(rstan::get_elapsed_time(m.Pl_sowndivW2_d$fit))
    loo(m.Pl_sowndivW2_d, m.Pl_sowndivW2_p)
    
    
    
####Pl plot predictions ####
    #different intercept  between weeks 
    #--> include week as additive term in lognormal part of the model
    
    pred.Pl_prior1 <- conditional_effects(m.Pl_sowndivW1_p2)[[3]]
    #pred.Pl_def1 <- conditional_effects(m.Pl_sowndivW1_d)[[3]]
    summary(m.Pl_sowndivW1_p2, prob =0.9)
    
    pred.Pl_prior2 <- conditional_effects(m.Pl_sowndivW2_p)[[3]]
    #pred.Pl_def2 <- conditional_effects(m.Pl_sowndivW2_d)[[3]]
    summary(m.Pl_sowndivW2_p, prob=0.9)
    
    pred.Pl_prior  <- conditional_effects(m.Pl_sowndiv_p)[[3]]
    pred.Pl_def  <- conditional_effects(m.Pl_sowndiv_d)[[3]]
    summary(m.Pl_sowndiv_d)
    
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    BREAKS = unique(dat$sowndivLogStd)
    
    
    ggplot(data = dat, aes(x= sowndivLogStd, y=Pl_per100g))+
      #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
      #                                 fill=treatment), 
      #           alpha=0.2, show.legend=FALSE)+
      geom_jitter(data =datW1,
                  width=0.2, shape=15, alpha=0.4, 
                  aes(col=treatment))+
      geom_jitter(data=datW2,
                  width=0.2, shape=17, alpha=0.8, 
                  aes(col=treatment))+
      #predictions week1
      geom_line(data=pred.Pl_prior1, aes(x= sowndivLogStd, y=estimate__, 
                                          col=treatment),
                linewidth= 0.5, linetype=1,
                show.legend = FALSE)+
      #predictions week 2:
      geom_line(data=pred.Pl_prior2, aes(x= sowndivLogStd, y=estimate__, 
                                         col=treatment),
                linewidth= 0.5, linetype=2,
                show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("16", "8", "4", "2", "1"))+
      scale_y_continuous(name = "Pl per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")  

#### Pl ~ sowndiv, both weeks: _p3 ####
    #_p3 has 1 div transition (least) and best elpd
    
    #for both weeks  
    m.Pl_sowndiv_p <- brm(
      bf(Pl_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
         hu ~ 1),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED, 
      control = list(adapt_delta=0.99))  #10 div
    summary(m.Pl_sowndiv_p, prob=0.9) #week is significant --> keep it
    
    #using a narrower prior:
    beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
    m.Pl_sowndiv_p2 <- update(m.Pl_sowndiv_p,
                                prior = beta_coeff_priors2,
                                seed = SEED) #0 div
    summary(m.Pl_sowndiv_p2, prob=0.9) #56 div 
    
    #using an even narrower prior:
    beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
    m.Pl_sowndiv_p3 <- update(m.Pl_sowndiv_p,
                              prior = beta_coeff_priors3,
                              seed = SEED) #0 div
    summary(m.Pl_sowndiv_p3, prob=0.9) #1 div
    
    loo(m.Pl_sowndiv_p, m.Pl_sowndiv_p2, m.Pl_sowndiv_p3) #no substantial difference (less than 2 SE) -->
      # best: _p2
    
    
    pp_check(m.Pl_sowndiv_p2, ndraws=100)+
      xlim(0,2000)
    
    #with default priors
    m.Pl_sowndiv_d <- brm(
      bf(Pl_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
         hu ~ 1),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #19 div
    summary(m.Pl_sowndiv_d, prob=0.9)
    
    #best: _d2
    summary(m.Pl_sowndiv_d2, prob=0.9)
    pp_check(m.Pl_sowndiv_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Pl_sowndiv_p, m.Pl_sowndiv_p2, m.Pl_sowndiv_p3, m.Pl_sowndiv_d)
      #choose _p3, as it has least divergent transitions
    
    print(rstan::get_elapsed_time(m.Pl_sowndiv_p$fit))
    print(rstan::get_elapsed_time(m.Pl_sowndiv_d$fit))  
    
    
#save the models
    save(m.Pl_sowndivW1_d, m.Pl_sowndivW1_p,
          m.Pl_sowndivW1_p2, m.Pl_sowndivW1_p3,
         m.Pl_sowndivW2_d, m.Pl_sowndivW2_p,
          m.Pl_sowndivW2_p2, m.Pl_sowndivW2_p3,
         m.Pl_sowndiv_d, m.Pl_sowndiv_p,
          m.Pl_sowndiv_p2, m.Pl_sowndiv_p3,
         file = "./statistics/brms/231213_Pl_sowndiv_priors.RData")
    
#remove to prevent crashes:
    rm(m.Pl_sowndivW1_d, m.Pl_sowndivW1_p,
      m.Pl_sowndivW1_p2, m.Pl_sowndivW1_p3,
      m.Pl_sowndivW2_d, m.Pl_sowndivW2_p,
      m.Pl_sowndivW2_p2, m.Pl_sowndivW2_p3,
      m.Pl_sowndiv_d, m.Pl_sowndiv_p,
      m.Pl_sowndiv_p2#, m.Pl_sowndiv_p3
      )
    
   
    
    
    
#### Pr ~ sowndiv: W1-p5, W2-p5, both-  ####
    #W1: _p5 6 div, _d3 9 div, loo is equal --> choose _p5
    #W2: _p5 5 div, _d4 13 div, loo is slightly better in _p5 (less than 2 SE) --> choose _p5
    load(file = "./statistics/brms/231213_Pr_sowndiv_priors.RData")
    
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
      control = list(adapt_delta=0.99)) #4 div
    summary(m.Pr_sowndivW1_p, prob=0.9)
    
    
        #remove hu~sowndivLogStd:treatment
        m.Pr_sowndivW1_p2 <-update(m.Pr_sowndivW1_p,
                                   bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                   hu ~ sowndivLogStd + treatment + (1|block/plot)),
                                   seed = SEED) # 18 div, exceeded max_treedepth, tail ESS too low(?)
        summary(m.Pr_sowndivW1_p2, prob=0.9)
        
        #remove hu~sowndivLogStd
        m.Pr_sowndivW1_p3 <-update(m.Pr_sowndivW1_p2,
                                   bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                      hu ~ treatment + (1|block/plot)),
                                   seed = SEED) #10 div, exceeded max_treedepth
        summary(m.Pr_sowndivW1_p3, prob=0.9)
        
        #increase max_treedepth
        m.Pr_sowndivW1_p4 <-update(m.Pr_sowndivW1_p3,
                                   control=list(max_treedepth=12),
                                   seed = SEED) #7 div
        summary(m.Pr_sowndivW1_p4, prob=0.9)
        
        
        #narrower priors: 
        beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
        m.Pr_sowndivW1_p5 <-update(m.Pr_sowndivW1_p4,
                                   prior=beta_coeff_priors2,
                                   seed=SEED) #6 div
        summary(m.Pr_sowndivW1_p5, prob=0.9)
        
        #even narrower priors: 
        beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
        m.Pr_sowndivW1_p6 <-update(m.Pr_sowndivW1_p5,
                                   prior=beta_coeff_priors3,
                                   seed=SEED) #31 div
        summary(m.Pr_sowndivW1_p6, prob=0.9)
        
        loo(m.Pr_sowndivW1_p5, m.Pr_sowndivW1_p4)
        #no elpd difference bigger than 2 SE --> 
        #best fit: _p5
    
    summary(m.Pr_sowndivW1_p5, prob=0.9)
    pp_check(m.Pr_sowndivW1_p5, ndraws=100)+
      xlim(0,200)
    
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
      control = list(adapt_delta=0.99)) #8 div
    summary(m.Pr_sowndivW1_d, prob=0.9)
    
    
        #remove hu~sowndivLogStd:treatment
        m.Pr_sowndivW1_d2 <-update(m.Pr_sowndivW1_d,
                                   bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                      hu ~ sowndivLogStd + treatment + (1|block/plot)),
                                   seed = SEED) #0 div
        summary(m.Pr_sowndivW1_d2, prob=0.9)
        
        #remove hu~sowndivLogStd
        m.Pr_sowndivW1_d3 <-update(m.Pr_sowndivW1_d2,
                                   bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                      hu ~ treatment + (1|block/plot)),
                                   seed = SEED) #9 div
        summary(m.Pr_sowndivW1_d3, prob=0.9)
        
        #best: _d3
    
    summary(m.Pr_sowndivW1_d3, prob=0.9)
    pp_check(m.Pr_sowndivW1_d, ndraws=100)+
      xlim(0,2000)
    
    #compare them  
    loo(m.Pr_sowndivW1_p4, m.Pr_sowndivW1_p5, m.Pr_sowndivW1_p6, m.Pr_sowndivW1_d3) #d_3 is best
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
      
      #remove hu ~ sowndivLogStd:treatment
      m.Pr_sowndivW2_p2 <- update(m.Pr_sowndivW2_p,
                                  bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                     hu ~ sowndivLogStd + treatment + (1|block/plot)),
                                  seed = SEED) #2 div
      summary(m.Pr_sowndivW2_p2, prob=0.9)
      
      #remove hu ~ sowndivLogStd
      m.Pr_sowndivW2_p3 <- update(m.Pr_sowndivW2_p2,
                                  bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                     hu ~ treatment + (1|block/plot)),
                                  seed = SEED)  #8 div
      summary(m.Pr_sowndivW2_p3, prob=0.9)
      
      #remove hu ~ treatment
      m.Pr_sowndivW2_p4 <- update(m.Pr_sowndivW2_p3,
                                  bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                     hu ~ 1),
                                  seed = SEED) #15 div
      summary(m.Pr_sowndivW2_p4, prob=0.9)
      
      #narrower prior:
      beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
      m.Pr_sowndivW2_p5 <-update(m.Pr_sowndivW2_p4,
                                 prior=beta_coeff_priors2,
                                 seed=SEED) #5 div
      summary(m.Pr_sowndivW2_p5, prob=0.9)
      
      #even narrower priors: 
      beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
      m.Pr_sowndivW2_p6 <-update(m.Pr_sowndivW2_p5,
                                 prior=beta_coeff_priors3,
                                 seed=SEED) #17 div
      summary(m.Pr_sowndivW2_p6, prob=0.9)
      
      loo(m.Pr_sowndivW2_p4, m.Pr_sowndivW2_p5 )
      #best: _p5
      
      summary(m.Pr_sowndivW2_p4, prob=0.9)
      pp_check(m.Pr_sowndivW2_p4, ndraws=100)
    
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
        control = list(adapt_delta=0.99))  #7 div
      summary(m.Pr_sowndivW2_d, prob=0.9)
      
      #remove hu ~ sowndivLogStd:treatment
      m.Pr_sowndivW2_d2 <- update(m.Pr_sowndivW2_d,
                                  bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                     hu ~ sowndivLogStd + treatment + (1|block/plot)),
                                  seed = SEED)
      summary(m.Pr_sowndivW2_d2, prob=0.9) #11 div
      
      #remove hu ~ sowndivLogStd
      m.Pr_sowndivW2_d3 <- update(m.Pr_sowndivW2_d2,
                                  bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                     hu ~ treatment + (1|block/plot)),
                                  seed = SEED) # 9 div
      summary(m.Pr_sowndivW2_d3, prob=0.9)
      
      #remove hu ~ treatment
      m.Pr_sowndivW2_d4 <- update(m.Pr_sowndivW2_d3,
                                  bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                     hu ~ 1),
                                  seed = SEED) #13 div
      summary(m.Pr_sowndivW2_d4, prob=0.9)
      
      
      #best _d4
      summary(m.Pr_sowndivW2_d2, prob=0.9) 
      pp_check(m.Pr_sowndivW2_d2, ndraws=100)
      
      #compare them  
      loo(m.Pr_sowndivW2_p5, m.Pr_sowndivW2_d4)
      print(rstan::get_elapsed_time(m.Pr_sowndivW2_p$fit))
      print(rstan::get_elapsed_time(m.Pr_sowndivW2_d$fit))
      loo(m.Pr_sowndivW2_d, m.Pr_sowndivW2_p)
    
    
####Pr plot predictions ####
    #difference in proportion of zeros/intercept and slopes 
    # --> include week as an 3- fold interaction term for the lognormal part of the model, 
    #and as an additive term for the hurdle part
    
    pred.Pr_prior1 <- conditional_effects(m.Pr_sowndivW1_p5)[[3]]
    #pred.Pr_def1 <- conditional_effects(m.Pr_sowndivW1_d)[[3]]
    summary(m.Pr_sowndivW1_p)
    
    pred.Pr_prior2 <- conditional_effects(m.Pr_sowndivW2_p5)[[3]]
    #pred.Pr_def2 <- conditional_effects(m.Pr_sowndivW2_d)[[3]]
    summary(m.Pr_sowndivW2_p)
    
    pred.Pr_prior  <- conditional_effects(m.Pr_sowndiv_p)[[3]]
    pred.Pr_def  <- conditional_effects(m.Pr_sowndiv_d)[[3]]
    summary(m.Pr_sowndiv_d)
    
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    BREAKS = unique(dat$sowndivLogStd)
    
    
    ggplot(data = dat, aes(x= sowndivLogStd, y=Pr_per100g))+
      #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
      #                                 fill=treatment), 
      #           alpha=0.2, show.legend=FALSE)+
      geom_jitter(data =datW1,
                  width=0.2, shape=15, alpha=0.4, 
                  aes(col=treatment))+
      geom_jitter(data=datW2,
                  width=0.2, shape=17, alpha=0.8, 
                  aes(col=treatment))+
      #predictions week1
      geom_line(data=pred.Pr_prior1, aes(x= sowndivLogStd, y=estimate__, 
                                        col=treatment),
                linewidth= 0.5, linetype=6,
                show.legend = FALSE)+
      #predictions week 2:
      geom_line(data=pred.Pr_prior2, aes(x= sowndivLogStd, y=estimate__, 
                                       col=treatment),
                linewidth=0.5, linetype=1,
                show.legend = FALSE)+
      
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("16", "8", "4", "2", "1"))+
      scale_y_continuous(name = "Pr per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")  
    
#### Pr ~ sowndiv, both weeks: _p7 ####
    #for both weeks  
    m.Pr_sowndiv_p <- brm(
      bf(Pr_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
         hu ~ sowndivLogStd + treatment + week + (1|block/plot)),
      data = dat, 
      prior = beta_coeff_priors,
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #0 div
    summary(m.Pr_sowndiv_p, prob=0.9)
    
    #with default priors
    m.Pr_sowndiv_d <- brm(
      bf(Pr_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
         hu ~ sowndivLogStd + treatment + week + (1|block/plot)),
      data = dat, 
      family = hurdle_lognormal,
      chains = 3,
      cores = 3,
      iter = 2000, warmup = 1000,
      seed = SEED,
      control = list(adapt_delta=0.99)) #1 div trans
    summary(m.Pr_sowndiv_d, prob=0.9)
      #simplifying with the default priors from now on
    
    #remove 3 fold interaction:
    m.Pr_sowndiv_d2 <- update(m.Pr_sowndiv_d,
                             bf(Pr_per100g ~ sowndivLogStd*treatment + week*treatment + week*sowndivLogStd + (1|block/plot),
                                hu ~ sowndivLogStd + treatment + week + (1|block/plot)),
                             seed = SEED) #13 div
    summary(m.Pr_sowndiv_d2, prob=0.9)
    
    #remove week*treatment:
    m.Pr_sowndiv_d3 <- update(m.Pr_sowndiv_d2,
                              bf(Pr_per100g ~ sowndivLogStd*treatment + week*sowndivLogStd + (1|block/plot),
                                 hu ~ sowndivLogStd + treatment + week + (1|block/plot)),
                              seed = SEED) #6 div
    summary(m.Pr_sowndiv_d3, prob=0.9)
    
    #remove week*sowndivLogStd:
    m.Pr_sowndiv_d4 <- update(m.Pr_sowndiv_d3,
                              bf(Pr_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                                 hu ~ sowndivLogStd + treatment + week + (1|block/plot)),
                              seed = SEED) #3 div, tail ESS too low
    summary(m.Pr_sowndiv_d4, prob=0.9) 
    
    #remove hu~sowndivLogStd:
    m.Pr_sowndiv_d5 <- update(m.Pr_sowndiv_d4,
                              bf(Pr_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                                 hu ~ treatment + week + (1|block/plot)),
                              seed = SEED) #19 div
    summary(m.Pr_sowndiv_d5, prob=0.9) 
    
    #remove week:
    m.Pr_sowndiv_d6 <- update(m.Pr_sowndiv_d5,
                              bf(Pr_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                 hu ~ treatment + week + (1|block/plot)),
                              seed = SEED) #1 div
    summary(m.Pr_sowndiv_d6, prob=0.9) 
        #NOTE: removing the very marginally sigifnificant hu~treatment 
        #does not substantially change anything  in the other terms
    
    #add a wide normal(0,20) prior
    beta_coeff_priors <- prior(normal(0,20), class = "b")  
    m.Pr_sowndiv_p7 <- update(m.Pr_sowndiv_d6,
                              prior = beta_coeff_priors,
                              seed = SEED) #0 div
    summary(m.Pr_sowndiv_p7, prob=0.9) 
    
    #add a normal(0,5) prior
    beta_coeff_priors2 <- prior(normal(0,5), class = "b")  
    m.Pr_sowndiv_p8 <- update(m.Pr_sowndiv_d6,
                              prior = beta_coeff_priors2,
                              seed = SEED) #1 div
    summary(m.Pr_sowndiv_p8, prob=0.9) 
    
    #add a normal(0,2) prior
    beta_coeff_priors3 <- prior(normal(0,2), class = "b")  
    m.Pr_sowndiv_p9 <- update(m.Pr_sowndiv_d6,
                              prior = beta_coeff_priors3,
                              seed = SEED)
    summary(m.Pr_sowndiv_p9, prob=0.9) 
    
    summary(m.Pr_sowndiv_p7, prob=0.90) 
    pp_check(m.Pr_sowndiv_p7, ndraws=100)+
      xlim(0,200)
    
    #compare them  
    loo(m.Pr_sowndiv_d6, m.Pr_sowndiv_p7, m.Pr_sowndiv_p8, m.Pr_sowndiv_p9)
    #choosing _p7, as it is the model with highest elpd and 0 divergent transitions
    print(rstan::get_elapsed_time(m.Pr_sowndiv_p$fit))
    print(rstan::get_elapsed_time(m.Pr_sowndiv_d$fit))  
    
#save the models: 
    save(m.Pr_sowndivW1_d, m.Pr_sowndivW1_p,
          m.Pr_sowndivW1_d2, m.Pr_sowndivW1_d3,
          m.Pr_sowndivW1_p2, m.Pr_sowndivW1_p3, m.Pr_sowndivW1_p4, m.Pr_sowndivW1_p5, m.Pr_sowndivW1_p6,
         m.Pr_sowndivW2_d, m.Pr_sowndivW2_p,
          m.Pr_sowndivW2_d2, m.Pr_sowndivW2_d3, m.Pr_sowndivW2_d4,
          m.Pr_sowndivW2_p2, m.Pr_sowndivW2_p3, m.Pr_sowndivW2_p4, m.Pr_sowndivW2_p5, m.Pr_sowndivW2_p6,
         m.Pr_sowndiv_d, m.Pr_sowndiv_p,
          m.Pr_sowndiv_d2,m.Pr_sowndiv_d3, m.Pr_sowndiv_d4,m.Pr_sowndiv_d5, m.Pr_sowndiv_d6,
          m.Pr_sowndiv_p7, m.Pr_sowndiv_p8, m.Pr_sowndiv_p9, 
         file = "./statistics/brms/231213_Pr_sowndiv_priors.RData")
    
#unload the models to prevent crashes:
    rm(m.Pr_sowndivW1_d, m.Pr_sowndivW1_p,
       m.Pr_sowndivW1_d2, m.Pr_sowndivW1_d3,
       m.Pr_sowndivW1_p2, m.Pr_sowndivW1_p3, m.Pr_sowndivW1_p4, m.Pr_sowndivW1_p5, m.Pr_sowndivW1_p6,
       m.Pr_sowndivW2_d, m.Pr_sowndivW2_p,
       m.Pr_sowndivW2_d2, m.Pr_sowndivW2_d3, m.Pr_sowndivW2_d4,
       m.Pr_sowndivW2_p2, m.Pr_sowndivW2_p3, m.Pr_sowndivW2_p4, m.Pr_sowndivW2_p5, m.Pr_sowndivW2_p6,
       m.Pr_sowndiv_d, m.Pr_sowndiv_p,
       m.Pr_sowndiv_d2,m.Pr_sowndiv_d3, m.Pr_sowndiv_d4,m.Pr_sowndiv_d5, m.Pr_sowndiv_d6,
       #m.Pr_sowndiv_p7, #keep the best
       m.Pr_sowndiv_p8, m.Pr_sowndiv_p9)
    
    
    
    
    
#### Om ~ sowndiv: W1-d3, W2-p, both-d4 ####
    #W1: _p3 13 div, _d3 10 div, p has slightly better elpd (less than 2 SE) -> d3 is best
    #W2: _p 1 div, _d 6 div, _p has slightly better elpd (less than 2 SE) --> _p is best
    #both: 10 div in p, zero in d, elpd slightly better in p
    load(file = "./statistics/brms/231213_Om_sowndiv_priors.RData")
    
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
      control = list(adapt_delta=0.99))  #0 div
    summary(m.Om_sowndivW1_p, prob=0.9)
    
        #remove hu~sowndivLogStd:treatment
        m.Om_sowndivW1_p2 <- update(m.Om_sowndivW1_p,
                                    bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                       hu ~ sowndivLogStd + treatment + (1|block/plot)),
                                    newdata = subset(dat, week=="W1"),
                                    seed = SEED) #28 div
        summary(m.Om_sowndivW1_p2, prob=0.9)
        
        #remove hu~sowndivLogStd
        m.Om_sowndivW1_p3 <- update(m.Om_sowndivW1_p,
                                    bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                       hu ~ treatment + (1|block/plot)),
                                    newdata = subset(dat, week=="W1"),
                                    seed = SEED) #13 div, bulk too low
        summary(m.Om_sowndivW1_p3, prob=0.9)
        
        
        #try a narrower prior:
        beta_coeff_priors2<- prior(normal(0,5), class = "b")
        m.Om_sowndivW1_p4 <- update(m.Om_sowndivW1_p3,
                                    prior=beta_coeff_priors2,
                                    newdata = subset(dat, week=="W1"),
                                    seed = SEED) #30 div, bulk ESS too low
        summary(m.Om_sowndivW1_p4, prob=0.9)
        
        #try an even narrower prior:
        beta_coeff_priors3<- prior(normal(0,2), class = "b")
        m.Om_sowndivW1_p5 <- update(m.Om_sowndivW1_p3,
                                    prior=beta_coeff_priors3,
                                    newdata = subset(dat, week=="W1"),
                                    seed = SEED) #probably worse, re-run to check as fit is lost
        summary(m.Om_sowndivW1_p5, prob=0.9)
        
        loo(m.Om_sowndivW1_p3,m.Om_sowndivW1_p4) #less than 2 SE of elpd difference -> choose model with les div transitions
        #decision: _p3 is best
      
      summary(m.Om_sowndivW1_p3, prob=0.9)
      pp_check(m.Om_sowndivW1_p3, ndraws=100)+
        xlim(0,50)
    
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
      control = list(adapt_delta=0.99)) #107 div
    summary(m.Om_sowndivW1_d, prob=0.9)
    
      #priciple of parsimony: remove hurdle interaction
        m.Om_sowndivW1_d2 <- update(m.Om_sowndivW1_d,
                                    bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                        hu ~ sowndivLogStd + treatment + (1|block/plot)),
                                    newdata = subset(dat, week=="W1"),
                                    seed = SEED)#65 div
        summary(m.Om_sowndivW1_d2, prob=0.9)
        
      #parsimony: remove hu~sowndivLogStd
        m.Om_sowndivW1_d3 <- update(m.Om_sowndivW1_d,
                                    bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                       hu ~ treatment + (1|block/plot)),
                                    newdata = subset(dat, week=="W1"),
                                    seed = SEED) #10 div
        summary(m.Om_sowndivW1_d3, prob=0.9) 
        
        #_d3 is best
      summary(m.Om_sowndivW1_d3, prob=0.9)
      pp_check(m.Om_sowndivW1_d3, ndraws=100)+
        xlim(0,100)
    
    #compare them  
    loo(m.Om_sowndivW1_p3, m.Om_sowndivW1_d3) #p3 has slightly better elpd (less than 2 SE)
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
      control = list(adapt_delta=0.99)) #1 div
    summary(m.Om_sowndivW2_p, prob=0.9)
    
     #with a narrower prior:
      beta_coeff_priors2<- prior(normal(0,5), class = "b")
      m.Om_sowndivW2_p2 <- update(m.Om_sowndivW2_p, 
                                  prior = beta_coeff_priors2,
                                  seed = SEED) #2 div
      summary(m.Om_sowndivW2_p2, prob=0.9)
      
      #with an even narrower prior:
      beta_coeff_priors3<- prior(normal(0,2), class = "b")
      m.Om_sowndivW2_p3 <- update(m.Om_sowndivW2_p, 
                                  prior = beta_coeff_priors3,
                                  seed = SEED) #9 div
      summary(m.Om_sowndivW2_p3, prob=0.9)
      
      loo(m.Om_sowndivW2_p, m.Om_sowndivW2_p2) # slightly better (less than 2 SE)
      #choose: _p
    
    summary(m.Om_sowndivW2_p2, prob=0.9)
    pp_check(m.Om_sowndivW2_p2, ndraws=100)
    
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
      control = list(adapt_delta=0.99)) #6 div 
    
    summary(m.Om_sowndivW2_d, prob=0.9)
    pp_check(m.Om_sowndivW2_d, ndraws=100)
    
    #compare them  
    loo(m.Om_sowndivW2_p, m.Om_sowndivW2_d) #_p has slightly better elpd (less than 2 SE)
    print(rstan::get_elapsed_time(m.Om_sowndivW2_p$fit))
    print(rstan::get_elapsed_time(m.Om_sowndivW2_d$fit))
    loo(m.Om_sowndivW2_d, m.Om_sowndivW2_p)
    

    
####Om plot predictions ####
    #different intercepts by week --> include week as additive term in hurdle and lognormal part of the model
    #different slopes by week: fit 3way interaction week:treatment:sowndivLogStd
    
    
    #pred.Om_prior1 <- conditional_effects(m.Om_sowndivW1_p4)[[3]]
    pred.Om_def1 <- conditional_effects(m.Om_sowndivW1_d3)[[3]]
    summary(m.Om_sowndivW1_d3)
    
    pred.Om_prior2 <- conditional_effects(m.Om_sowndivW2_p)[[3]]
    #pred.Om_def2 <- conditional_effects(m.Om_sowndivW2_d)[[3]]
    summary(m.Om_sowndivW2_p)
    
    pred.Om_prior  <- conditional_effects(m.Om_sowndiv_p4)[[3]]
    pred.Om_def  <- conditional_effects(m.Om_sowndiv_d4)[[3]]
    summary(m.Om_sowndiv_d4)
    
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    BREAKS = unique(dat$sowndivLogStd)
    
    
   p.Om <- ggplot(data = dat, aes(x= sowndivLogStd, y=Om_per100g))+
      #geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
      #                                 fill=treatment), 
      #           alpha=0.2, show.legend=FALSE)+
      geom_jitter(data =datW1,
                  width=0.2, shape=15, alpha=0.4, 
                  aes(col=treatment))+
      geom_jitter(data=datW2,
                  width=0.2, shape=17, alpha=0.8, 
                  aes(col=treatment))+
      #predictions week1
      geom_line(data=pred.Om_def1, aes(x= sowndivLogStd, y=estimate__, 
                                         col=treatment),
                linetype=6, 
                linewidth= 0.5, show.legend = FALSE)+
      #predictions week 2:
      geom_line(data=pred.Om_prior2, aes(x= sowndivLogStd, y=estimate__, 
                                          col=treatment),
                linetype=1,
                linewidth= 0.5, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("16", "8", "4", "2", "1"))+
      scale_y_continuous(name = "Om per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")  
   
   p.Om
   
   ggsave(filename = "./plots model comparison/Om_21.png" ,plot = p.Om,
          height = 4,
          width = 4)
   
   
#### Om ~ sowndiv, both weeks: _d7 ####
   #_d7: 1 div transition, _p8 and _p9 1 div, elpd in _d7 marginally better (less than 2SE)
   # --> choose d7
   
   #for both weeks 
   #with prior
   m.Om_sowndiv_p <- brm(
     bf(Om_per100g ~ sowndivLogStd*treatment*week +(1|block/plot),
        hu ~ sowndivLogStd + treatment + week + (1|block/plot)),
     data = dat, 
     prior = beta_coeff_priors,
     family = hurdle_lognormal,
     chains = 3,
     cores = 3,
     iter = 2000, warmup = 1000,
     seed = SEED,
     control = list(adapt_delta=0.99)) #10 div
   summary(m.Om_sowndiv_p, prob =0.9)
   
   #with default priors
   m.Om_sowndiv_d <- brm(
     bf(Om_per100g ~ sowndivLogStd*treatment*week + (1|block/plot),
        hu ~ sowndivLogStd + treatment + week +  (1|block/plot)),
     data = dat, 
     family = hurdle_lognormal,
     chains = 3,
     cores = 3,
     iter = 2000, warmup = 1000,
     seed = SEED,
     control = list(adapt_delta=0.99)) #all good
   summary(m.Om_sowndiv_d, prob =0.9)
   
   #parsimony: remove week:sowndivLogStd:treatment:
   m.Om_sowndiv_d2 <- update(m.Om_sowndiv_d,
                             bf(Om_per100g ~ sowndivLogStd*treatment+ treatment*week + sowndivLogStd*week + (1|block/plot),
                                hu ~ sowndivLogStd + treatment + (1|block/plot)),
                             seed=SEED)
   summary(m.Om_sowndiv_d2, prob = 0.9)
   emt.pairs <- emtrends(m.Om_sowndiv_d2, specs = c("treatment", "week"), var="sowndivLogStd") %>%
     pairs()
   bayestestR::p_direction(emt.pairs) #t1W1 - t1W2 : 71%; t2W1 -t2W2
   
   #parsimony: remove treatment*week
   m.Om_sowndiv_d3 <- update(m.Om_sowndiv_d,
                             bf(Om_per100g ~ sowndivLogStd*treatment + sowndivLogStd*week + (1|block/plot),
                                hu ~ sowndivLogStd + treatment + (1|block/plot)),
                             seed=SEED) #1 div
   summary(m.Om_sowndiv_d3, prob = 0.9)
   
   #parsimony: remove sowndivLogStd*week
   m.Om_sowndiv_d4 <- update(m.Om_sowndiv_d,
                             bf(Om_per100g ~ sowndivLogStd*treatment + week + (1|block/plot),
                                hu ~ sowndivLogStd + treatment + (1|block/plot)),
                             seed=SEED) #39 div
   summary(m.Om_sowndiv_d4, prob = 0.9)
   
   #parsimony: remove week
   m.Om_sowndiv_d5 <- update(m.Om_sowndiv_d,
                             bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                hu ~ sowndivLogStd + treatment + (1|block/plot)),
                             seed=SEED)
   summary(m.Om_sowndiv_d5, prob = 0.9)
   
   #parsimony: remove hurdle~sowndivLogStd
   m.Om_sowndiv_d6 <- update(m.Om_sowndiv_d,
                             bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                hu ~ treatment + (1|block/plot)),
                             seed=SEED) #0 div
   summary(m.Om_sowndiv_d6, prob = 0.9)
   
   #parsimony: remove hurdle~treatment
   m.Om_sowndiv_d7 <- update(m.Om_sowndiv_d,
                             bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                hu ~ 1),
                             seed=SEED) #1 div
   summary(m.Om_sowndiv_d7, prob = 0.9) 
   
   #add a normal(0,20) prior for beta coefficients
   m.Om_sowndiv_p7 <- update(m.Om_sowndiv_d,
                             bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                hu ~ 1),
                             prior = beta_coeff_priors,
                             seed=SEED) #7 div
   summary(m.Om_sowndiv_p7, prob = 0.9)
   
   #add a normal(0,5) prior for beta coefficients
   m.Om_sowndiv_p8 <- update(m.Om_sowndiv_d,
                             bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                hu ~ 1),
                             prior = beta_coeff_priors2,
                             seed=SEED) #1 div
   summary(m.Om_sowndiv_p8, prob = 0.9)
   
   #add a normal(0,2) prior for beta coefficients
   m.Om_sowndiv_p9 <- update(m.Om_sowndiv_d,
                             bf(Om_per100g ~ sowndivLogStd*treatment + (1|block/plot),
                                hu ~ 1),
                             prior = beta_coeff_priors3,
                             seed=SEED) #1div
   summary(m.Om_sowndiv_p9, prob = 0.9)
   
   pp_check(m.Om_sowndiv_p9, ndraws=100)+
     xlim(0,200)
   loo(m.Om_sowndiv_p9, m.Om_sowndiv_p8, m.Om_sowndiv_p7, m.Om_sowndiv_d7) 
   #best elpd _d7, 1 divergent transition like _p8 and _p9 --> choose d7
   
   
   

#save them:
   save(m.Om_sowndivW1_d, m.Om_sowndivW1_p, #basic
        m.Om_sowndivW1_d2 , m.Om_sowndivW1_d3, 
        m.Om_sowndivW1_p2, m.Om_sowndivW1_p3, m.Om_sowndivW1_p4, #simplified according to Ockham's razor
        m.Om_sowndivW2_d, m.Om_sowndivW2_p,
        m.Om_sowndivW2_p2 ,m.Om_sowndivW2_p3, # _d is most simplified already
        m.Om_sowndiv_d, m.Om_sowndiv_p,
        m.Om_sowndiv_p7, m.Om_sowndiv_p8, m.Om_sowndiv_p9, 
        m.Om_sowndiv_d2, m.Om_sowndiv_d3, m.Om_sowndiv_d4,m.Om_sowndiv_d5, m.Om_sowndiv_d6, m.Om_sowndiv_d7,
        file = "./statistics/brms/231213_Om_sowndiv_priors.RData")
   
#unload to prevent crashes:    
   rm(m.Om_sowndivW1_d, m.Om_sowndivW1_p, #basic
      m.Om_sowndivW1_d2 , m.Om_sowndivW1_d3, 
      m.Om_sowndivW1_p2, m.Om_sowndivW1_p3, m.Om_sowndivW1_p4, 
      m.Om_sowndivW2_d, m.Om_sowndivW2_p,
      m.Om_sowndivW2_p2 ,m.Om_sowndivW2_p3, 
      m.Om_sowndiv_d, m.Om_sowndiv_p,
      m.Om_sowndiv_p7, m.Om_sowndiv_p8, m.Om_sowndiv_p9, 
      m.Om_sowndiv_d2, m.Om_sowndiv_d3, m.Om_sowndiv_d4,m.Om_sowndiv_d5, m.Om_sowndiv_d6#, m.Om_sowndiv_d7
      )

####save best fit models####
   save(m.Ba_sowndiv_p5, m.Fu_sowndiv_p3, m.Om_sowndiv_d7, m.Pl_sowndiv_p3, m.Pr_sowndiv_p7,
        file = "./statistics/brms/231215_Trophic_sowndiv_BestFits.RData")
   
####emmeans####
   emt = emtrends(m.Om_sowndiv_d7, specs = c("treatment"), var="sowndivLogStd")
   summary(emt, point.est=mean, level = .9) 
   emt.pairs <- pairs(emt)
   summary(emt.pairs, point.est=mean, level = .9)
   bayestestR::p_direction(emt.pairs)
   