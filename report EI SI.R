#load the selected models:
load("./statistics/brms/240215_EI_mselect.RData")
load("./statistics/brms/240215_SI_mselect.RData")



####  reporting SI  ####
min(dat$sowndivLogStd)

#these three predictions result in the same estimate and CI:
    conditional_effects(m.SI.sowndiv_gaus_p5, robust = TRUE,
                        method = "posterior_epred") 
    emmeans(m.SI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd",
            at = list(sowndivLogStd = -1))
    emtrends(m.SI.sowndiv_gaus_p5, specs = c("treatment", "sowndivLogStd"), var = "sowndivLogStd",
            at = list(sowndivLogStd = unique(dat$sowndivLogStd))) 
            #so the slope is the same everywhere, because its a gaussian family
    
#lets check pd
    emtrends(m.SI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      bayestestR::p_direction() #effectof sowndiv at treatment 3 is marginally significant at pd=96.97%
    
    emtrends(m.SI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>%
      pairs() %>% bayestestR::p_direction() #no different slopes between treatments
    
    emtrends(m.SI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      bayestestR::p_direction() #effectof sowndiv at treatment 3 is marginally significant at pd=97.40%
    
    emtrends(m.SI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd") %>%
      pairs() %>% bayestestR::p_direction() #no different slopes between treatments
    
# and the slopes:
    #increase in SI by increasing sowndiv by 1
    emtrends(m.SI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") 
    
    #icrease in SI by doubling SOWNdiv
    emtrends(m.SI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>%
      summary(., point.est = mean) %>% #reporting mean, not median
      as.data.frame() %>%
      mutate_at((vars(sowndivLogStd.trend, lower.HPD, upper.HPD) ), ~ .*0.7249)  %>% #multplies every column in vars with 0.7249
      mutate_at((vars(sowndivLogStd.trend, lower.HPD, upper.HPD) ), round, digits=2) 

#SI estimates at min/max sown diversity    
  emmeans(m.SI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd",
          at = list(sowndivLogStd = min(dat$sowndivLogStd)))
  emmeans(m.SI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd",
          at = list(sowndivLogStd = max(dat$sowndivLogStd)))
  
#SI estimates at min max realdiv:
  emmeans(m.SI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd",
          at = list(realivLogStd = min(dat$realdivLogStd)))
  emmeans(m.SI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd",
          at = list(realdivLogStd = max(dat$realdivLogStd)))

####reporting EI ####
conditional_effects(m.EI.sowndiv_gaus_p5, robust = FALSE, spaghetti = TRUE, ndraws = 50) #robust=FALSE: shows mean as point estimate
conditional_effects(m.EI.realdiv_gaus_p5)


?marginal_effects()
marginaleffects::avg_predictions(m.EI.realdiv_gaus_p5, newdata = dat)


emmeans(m.EI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd",
        at = list(sowndivLogStd = min(dat$sowndivLogStd)))
emmeans(m.EI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd",
        at = list(sowndivLogStd = max(dat$sowndivLogStd)))

#pd:
emtrends(m.EI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
  bayestestR::p_direction()
emtrends(m.EI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
  bayestestR::p_direction()

#slopes:
emtrends(m.EI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>%
  summary(., point.est = mean) %>% #reporting mean, not median
  as.data.frame() %>%
  mutate_at((vars(sowndivLogStd.trend, lower.HPD, upper.HPD) ), ~ .*0.7249)  %>% #multplies every column in vars with 0.7249 (1/sd(sowndivLog)
  mutate_at((vars(sowndivLogStd.trend, lower.HPD, upper.HPD) ), round, digits=2) 

emtrends(m.EI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd") %>%
  summary(., point.est = mean) %>% #reporting mean, not median
  as.data.frame() %>%
  mutate_at((vars(realdivLogStd.trend, lower.HPD, upper.HPD) ), ~ .*0.8283)  %>% #multplies every column in vars with 1/sd(realdivLog)
  mutate_at((vars(realdivLogStd.trend, lower.HPD, upper.HPD) ), round, digits=2) 

#### plotting ####

#EI vs SI:
ggplot(dat, aes(x=SI, y=EI, col=treatment))+
  geom_point(alpha = 0.5, aes(size = sowndivLog))+
  scale_size(range = c(1,2) )+
  facet_wrap(~sowndivLog)



#traditional scatter regression SI:
    predictions <- conditional_effects(m.SI.realdiv_gaus_p5, method = "posterior_epred")[[3]]
    predictions.emm <- emmeans(m.SI.realdiv_gaus_p5, specs = c("treatment", "realdivLogStd"), var = "realdivLogStd",
                               #re_formula = NULL,
                               at = list(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 100))) %>% 
      as.data.frame() %>%
      mutate(SI = mean(dat$SI, na.rm = TRUE)) #need this to plot with geom_ribbon()
    
    treatments = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    BREAKS <- c
    
  #estimated marginal means:  
    ggplot(dat, aes(x=realdivLogStd, y=SI) )+
     geom_ribbon(data=predictions.emm, aes(ymin= lower.HPD, ymax=upper.HPD, 
                                            fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_line(data=predictions.emm, aes(x= realdivLogStd, y=emmean, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      geom_jitter(width=0.02, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      scale_color_manual(labels=treatments, values = cols)+
      scale_x_continuous(name = "sown plant diversity", #breaks = BREAKS,
                         labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "SI")+
      theme_classic()+
      theme(legend.position ="bottom")
    
#conditional effects:  
    ggplot(dat, aes(x=realdivLogStd, y=SI) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                 alpha=0.2, show.legend=FALSE)+
      geom_line(data=predictions, aes(x= realdivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      
      geom_jitter(width=0.02, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      scale_color_manual(labels=treatments, values = cols)+
      scale_x_continuous(name = "sown plant diversity", #breaks = BREAKS,
                         labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "SI")+
      theme_classic()+
      theme(legend.position ="bottom")
    
####WORK HERE####
   grand_mean.EI <- m.EI.realdiv_gaus_p5 %>% 
      epred_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                        treatment = c("1", "2", "3")),
                  re_formula = NA) #expectation of the predictive posterior distribution
    
    post_predict.EI <- m.EI.realdiv_gaus_p5 %>% 
      predicted_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                            treatment = c("1", "2", "3")),
                      re_formula = NA) #predictive posterior distribution
    #aesthetics:
    treat_labs = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    
    ggplot(grand_mean.EI, aes(x=realdivLogStd, y=.epred, col=treatment))+
      stat_lineribbon(geom = "ribbon", data = post_predict.EI, 
                      aes(x=realdivLogStd, y=.prediction, fill=treatment), 
                      .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment), 
                      .width = .95, alpha=0.3, colour=NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", 
                      aes( fill=NA, col=treatment#, linetype =c(rep("dashed", 6e4), rep("solid", 6e4), rep("dotted", 6e4) )
                           ), show.legend = FALSE)+
      geom_point(data = dat, aes(x=realdivLogStd, y=EI))+
      scale_y_continuous(name = "EI")+
      scale_x_continuous(name = "realized plant richness",
                         breaks = seq(min(dat$realdivLogStd), 1.916, length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$realdivLogStd), 1.92) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      facet_wrap(~treatment, labeller = labeller(treatment = c("1" = "+SH+PH", "2" = "+SH-PH", "3" = "-SH-SH")))+
      theme_bw()+
      ggtitle("EI ~ realdiv*treatment + (1|block/plot)")
      theme(legend.position = "none")
    
    emmeans(m.EI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd",
            at = list(realdivLogStd = min(dat$realdivLogStd)))
    
#### end work here ####    
#scatter EI:    
    predictions <- conditional_effects(m.EI.realdiv_gaus_p5, method = "posterior_epred")[[3]]
    treatments = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    BREAKS <- c
    
    ggplot(dat, aes(x=realdivLogStd, y=EI) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_line(data=predictions, aes(x= realdivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      
      geom_jitter(width=0.02, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      scale_color_manual(labels=treatments, values = cols)+
      scale_x_continuous(name = "sown plant diversity", #breaks = BREAKS,
                         labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "EI")+
      #theme_classic()+
      theme(legend.position ="bottom")+
      facet_wrap(~treatment)


ggplot(dat, aes(x=sowndivLog, y=EI, col = treatment))+
  geom_point()
ggplot(dat, aes(x=sowndivLog, y=SI, col = treatment))+
  geom_point()
