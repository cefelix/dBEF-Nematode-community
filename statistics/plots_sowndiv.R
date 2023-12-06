library(brms)
library(bayesplot)
library(tidybayes)
library(dplyr)
library(ggplot2)


#our data as in the models
dat <- subset(dBEF_nem21, sowndiv != 60) #60 sp. plots
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog) %>%
  mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
         .after = realdivLog)


####Load models: ####
  #best fit models: trophic ~ sowndiv
  load(file="./statistics/brms/231205_densBEST_sowndiv.RData")
  #easier model names:
  m.Ba.sowndiv <- m.Ba.2way_b2
    m.Fu.sowndiv <- m.Fu.2way_b
    m.Pl.sowndiv <- m.Pl.2way_d
    m.Pr.sowndiv <- m.Pr.2way_d2
    m.Om.sowndiv <- m.Om.2way_d
  rm(m.Ba.2way_b2, m.Fu.2way_b, m.Pl.2way_d, m.Pr.2way_d2, m.Om.2way_d)
  
  
#### Ba ~ sowndiv####

#using conditional effects:
    predictions <- conditional_effects(m.Ba.sowndiv, re_formula = NULL,
                                       prob = 0.9)[[4]]
    #when using conditional effects, to account for random effects set re_formula=NUll]
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    BREAKS = unique(dat$sowndivLogStd)
    ggplot(dat, aes(x=sowndivLogStd, y=Ba_per100g) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.2, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=predictions, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "Ba per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")



#taking draws from posterior distribution
  #add predictions as column to df:
  dat2 <- dat %>% add_predicted_draws(m.Ba.sowndiv, ndraws = 20)
  
  ggplot(data = dat2, aes(y=Ba_per100g, x=sowndivLogStd))+
    stat_dist_lineribbon(aes(y = .prediction, col = treatment), .width = c(.90),
                         alpha = 0.2)+
    geom_jitter(data = dat, aes(y=Ba_per100g, x=sowndivLogStd, col=treatment),
                width=0.2)

#### Fu ~ sowndiv####

  #using conditional effects:
    predictions <- conditional_effects(m.Fu.sowndiv, re_formula = NULL,
                                       prob = 0.9)[[4]]
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    BREAKS = unique(dat$sowndivLogStd)
    ggplot(dat, aes(x=sowndivLogStd, y=Fu_per100g) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.2, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=predictions, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("16", "8", "4", "2", "1"))+
      scale_y_continuous(name = "Fu per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")
  

#### Pl ~ sowndiv ####

#using conditional effects:
    predictions <- conditional_effects(m.Pl.sowndiv, re_formula = NULL,
                                       prob = 0.9)[[3]]
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    BREAKS = unique(dat$sowndivLogStd)
    ggplot(dat, aes(x=sowndivLogStd, y=Pl_per100g) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.2, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=predictions, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("16", "8", "4", "2", "1"))+
      scale_y_continuous(name = "Pl per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")
    
#### Pr ~ sowndiv ####

#using conditional effects:
    predictions <- conditional_effects(m.Pr.sowndiv, re_formula = NULL, 
                                       prob = 0.9)[[3]]
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    BREAKS = unique(dat$sowndivLogStd)
    ggplot(dat, aes(x=sowndivLogStd, y=Pr_per100g) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.2, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=predictions, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("16", "8", "4", "2", "1"))+
      scale_y_continuous(name = "Pr per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")   
    
#### Om ~ sowndiv ####
    
    #using conditional effects:
    predictions <- conditional_effects(m.Om.sowndiv, re_formula = NULL, 
                                       prob =0.9)[[3]]
    treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    BREAKS = unique(dat$sowndivLogStd)
    ggplot(dat, aes(x=sowndivLogStd, y=Om_per100g) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_jitter(width=0.2, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      geom_line(data=predictions, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      scale_color_manual(labels=treatments2, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("16", "8", "4", "2", "1"))+
      scale_y_continuous(name = "Om per 100g DW")+
      theme_bw()+
      theme(legend.position ="bottom")   
    
    
####summarise model results in a table####
    summary(m.Ba.sowndiv, prob=0.9) #check again model selection, non-significant terms!
    summary(m.Fu.sowndiv, prob=0.9) #check again, non-signif.
    summary(m.Pl.sowndiv, prob=0.9)
    summary(m.Pr.sowndiv, prob=0.9) #moist week marginally decreases prob(zero)
    summary(m.Om.sowndiv, prob=0.95) #moist week significantly decreases prob(zero) 
    
    
    
    

