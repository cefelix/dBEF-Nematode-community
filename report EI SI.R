#load the selected models:
load("./statistics/brms/240215_EI_mselect.RData")
load("./statistics/brms/240215_SI_mselect.RData")

library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(brms)
library(tidybayes)
library(bayestestR)
library(emmeans)


####  reporting SI  ####
min(dat$sowndivLogStd)
emmeans(m.SI.sowndiv_gaus_p5, specs = c("treatment", "sowndivLogStd")) %>% bayestestR::describe_posterior()

#these three predictions result in the same estimate and CI:
    conditional_effects(m.SI.sowndiv_gaus_p5, robust = TRUE,
                        method = "posterior_epred") 
    emmeans(m.SI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd",
            at = list(sowndivLogStd = -1)) 
    emtrends(m.SI.sowndiv_gaus_p5, specs = c("treatment", "sowndivLogStd"), var = "sowndivLogStd",
            at = list(sowndivLogStd = unique(dat$sowndivLogStd))) 
            #so the slope is the same everywhere, because its a gaussian family

# report the posterior:
    emtrends(m.SI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>% bayestestR::describe_posterior() 
    bayestestR::describe_posterior(m.SI.sowndiv_gaus_p5) %>% plot()
    bayestestR::describe_posterior(m.SI.sowndiv_gaus_p5)
    
#lets check pd
    emtrends(m.SI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      bayestestR::p_direction() #effect of sowndiv at treatment 3 is marginally significant at pd=96.97%
    
    emtrends(m.SI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>%
      pairs() %>% bayestestR::p_direction() #no different slopes between treatments
    
    emtrends(m.SI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      bayestestR::p_direction() #effect of sowndiv at treatment 3 is marginally significant at pd=97.40%
     
    emtrends(m.SI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd") %>%
      pairs() %>% bayestestR::p_direction() #no different slopes between treatments
    
# and the slopes:
    #increase in SI by increasing sowndiv by 1 SD
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

emmeans(m.EI.sowndiv_gaus_p5, specs = c("treatment"))


?marginal_effects()
marginaleffects::avg_predictions(m.EI.realdiv_gaus_p5, newdata = dat)

# emmeans sown
emmeans(m.EI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd",
        at = list(sowndivLogStd = min(dat$sowndivLogStd)),
        epred=TRUE)
emmeans(m.EI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd",
        at = list(sowndivLogStd = max(dat$sowndivLogStd)),
        epred=TRUE)
emmeans(m.EI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd",
        at = list(sowndivLogStd = 0),
        epred=TRUE)
emmeans(m.EI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>% pairs() %>% p_direction()

#emmeans realized 
emmeans(m.EI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd",
        at = list(sowndivLogStd = min(dat$realdivLogStd)))
emmeans(m.EI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd",
        at = list(sowndivLogStd = max(dat$realdivLogStd)))
emmeans(m.EI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd") %>% pairs() %>% p_direction()

#pd:
emtrends(m.EI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
  bayestestR::p_direction()
emtrends(m.EI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
  bayestestR::p_direction()

  #differences between slopes
  emtrends(m.EI.sowndiv_gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>%
    pairs() %>% bayestestR::p_direction() #no different slopes between treatments

  emtrends(m.EI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd") %>%
    pairs() %>% bayestestR::p_direction() #no different slopes between treatments
  
  

emmeans(m.EI.realdiv_gaus_p5, pairwise ~ treatment|realdivLogStd, cov.reduce = range) 
emmeans(m.EI.realdiv_gaus_p5, pairwise ~ treatment|realdivLogStd, cov.reduce = range) %>% 
  bayestestR::p_direction()



emmeans(m.EI.sowndiv_gaus_p5, specs = c("treatment")) 
emmeans(m.EI.sowndiv_gaus_p5, specs = c("treatment")) %>% 
  pairs() %>% p_direction()
  m.EI.sowndiv_gaus_p5 %>% p_direction() #same as running it on the emmeans output

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



#### plotting EI vs SI (4 quadrants plot) ####
#all in one: 
cols=c("brown2", "darkolivegreen", "dodgerblue3")

#entire dataset in one plot:
p.EI_Si <- ggplot(dat, aes(x=SI, y=EI, col=treatment))+
  geom_hline(yintercept = 50 )+
  geom_vline(xintercept = 50 )+
  geom_point(alpha = 0.6 )+
  scale_size(1 )+
  theme_bw()+
  ggtitle("")+
  scale_color_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
  theme(legend.position = "bottom")
p.EI_Si

#wrapped against sown plant richness
p.EI_SI.wrap <- ggplot(dat, aes(x=SI, y=EI, col=treatment))+
  geom_hline(yintercept = 50 )+
  geom_vline(xintercept = 50 )+
  geom_point(alpha = 0.6 )+
  scale_size(1 )+
  theme_bw()+
  ggtitle("a",
          subtitle = "Faunal profile at each sown plant richness" )+
  scale_color_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
  facet_wrap(~sowndiv, ncol=5, nrow =1)+
  #theme(aspect.ratio = 1)
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        title = element_text(size = 9.5))

  theme(legend.position = c(0.85, 0.25),
        legend.key.size = unit(1.2, "cm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size= 12),
        aspect.ratio = 1,
        title = element_text(size = 10.5))
p.EI_SI.wrap

    
####plotting EI ~ realdiv ####

emmip(m.EI.realdiv_gaus_p5, treatment ~ realdivLogStd, CIs = TRUE,
      cov.reduce = FALSE, ylab = "predicted EI")

   grand_mean.EI <- m.EI.realdiv_gaus_p5 %>% 
      epred_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                        treatment = c("1", "2", "3")),
                  re_formula = NA)  #expectation of the predictive posterior distribution // 
                                    #same as: estimated marginal means and 95% HPDI

    post_predict.EI <- m.EI.realdiv_gaus_p5 %>% 
      predicted_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                            treatment = c("1", "2", "3")),
                      re_formula = NA) #predictive posterior distribution
    #aesthetics:
    treat_labs = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    
  p.EI.realdiv <-  ggplot(grand_mean.EI, aes(x=realdivLogStd, y=.epred, col=treatment))+
      stat_lineribbon(geom = "ribbon", data = post_predict.EI, 
                      aes(x=realdivLogStd, y=.prediction, fill=treatment), 
                      .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment), 
                      .width = .95, alpha=0.3, colour=NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", 
                      aes( fill=NA, col=treatment#, linetype =c(rep("dashed", 6e4), rep("solid", 6e4), rep("dotted", 6e4) )
                           ), show.legend = FALSE)+
      geom_point(data = dat, aes(x=realdivLogStd, y=EI), alpha = 0.6 )+
      scale_y_continuous(name = "EI")+
      scale_x_continuous(name = "realized plant richness",
                         breaks = seq(min(dat$realdivLogStd), 1.916, length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$realdivLogStd), 1.92) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      facet_wrap(~treatment, labeller = labeller(treatment = c("1" = "+SH+PH", "2" = "+SH-PH", "3" = "-SH-SH")))+
      theme_bw()+
      ggtitle("b",
              "EI ~ realdiv*treatment + (1|block/plot)")+
      theme(legend.position = "none",
            title = element_text(size = 9.5))
    
    p.EI.realdiv  
      
    emmeans(m.EI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd",
            at = list(realdivLogStd = min(dat$realdivLogStd)))
    
####plotting EI ~ sowndiv ####
    
    emmip(m.EI.realdiv_gaus_p5, treatment ~ realdivLogStd, CIs = TRUE,
          cov.reduce = FALSE, ylab = "predicted EI")
    
    grand_mean.EI <- m.EI.sowndiv_gaus_p5 %>% 
      epred_draws(newdata = expand.grid(sowndivLogStd = seq(min(dat$sowndivLogStd), max(dat$sowndivLogStd), length.out = 20),
                                        treatment = c("1", "2", "3")),
                  re_formula = NA)  #expectation of the predictive posterior distribution // 
    #same as: estimated marginal means and 95% HPDI
    
    post_predict.EI <- m.EI.sowndiv_gaus_p5 %>% 
      predicted_draws(newdata = expand.grid(sowndivLogStd = seq(min(dat$sowndivLogStd), max(dat$sowndivLogStd), length.out = 20),
                                            treatment = c("1", "2", "3")),
                      re_formula = NA) #predictive posterior distribution
    #aesthetics:
    treat_labs = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    
    p.EI.sowndiv <-  ggplot(grand_mean.EI, aes(x=sowndivLogStd, y=.epred, col=treatment))+
      stat_lineribbon(geom = "ribbon", data = post_predict.EI, 
                      aes(x=sowndivLogStd, y=.prediction, fill=treatment), 
                      .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment), 
                      .width = .95, alpha=0.3, colour=NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", 
                      aes( fill=NA, col=treatment#, linetype =c(rep("dashed", 6e4), rep("solid", 6e4), rep("dotted", 6e4) )
                      ), show.legend = FALSE)+
      geom_point(data = dat, aes(x=sowndivLogStd, y=EI), alpha = 0.6 )+
      scale_y_continuous(name = "EI")+
      scale_x_continuous(name = "sown plant richness",
                         breaks = seq(min(dat$sowndivLogStd), 1.916, length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$sowndivLogStd), 1.92) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      facet_wrap(~treatment, labeller = labeller(treatment = c("1" = "+SH+PH", "2" = "+SH-PH", "3" = "-SH-SH")))+
      theme_bw()+
      ggtitle("b",
              "EI ~ sowndiv*treatment + (1|block/plot)")+
      theme(legend.position = "none",
            title = element_text(size = 9.5))
    
    p.EI.sowndiv  
    
    emmeans(m.EI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd",
            at = list(realdivLogStd = min(dat$realdivLogStd)))
    

#### plotting SI ~ realdiv ####    
    grand_mean.SI <- m.SI.realdiv_gaus_p5 %>% 
      epred_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                        treatment = c("1", "2", "3")),
                  re_formula = NA) #expectation of the predictive posterior distribution
    grand_mean.SI %>% filter()
    
    post_predict.SI <- m.SI.realdiv_gaus_p5 %>% 
      predicted_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                            treatment = c("1", "2", "3")),
                      re_formula = NA) #predictive posterior distribution
    
    #aesthetics:#aesthetics:mean()
    treat_labs = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    
    p.SI.realdiv <-  ggplot(grand_mean.SI, aes(x=realdivLogStd, y=.epred, col=treatment))+
      stat_lineribbon(geom = "ribbon", data = post_predict.SI, 
                      aes(x=realdivLogStd, y=.prediction, fill=treatment), 
                      .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment), 
                      .width = .95, alpha=0.3, colour=NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", 
                      aes( fill=NA, col=treatment#, linetype =c(rep("dashed", 6e4), rep("solid", 6e4), rep("dotted", 6e4) )
                      ), show.legend = FALSE)+
      geom_point(data = dat, aes(x=realdivLogStd, y=SI), alpha = 0.6)+
      scale_y_continuous(name = "SI", 
                         breaks = c(0, 25, 50, 75, 100))+
      scale_x_continuous(name = "realized plant richness",
                         breaks = seq(min(dat$realdivLogStd), 1.916, length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$realdivLogStd), 1.92) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      facet_wrap(~treatment, labeller = labeller(treatment = c("1" = "+SH+PH", "2" = "+SH-PH", "3" = "-SH-SH")))+
      theme_bw()+
      ggtitle("c",
              subtitle = "SI ~ realdiv*treatment + (1|block/plot)")+
    theme(legend.position = "none",
          title = element_text(size = 9.5))
    
    p.SI.realdiv  
    
    emmeans(m.SI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd",
            at = list(realdivLogStd = min(dat$realdivLogStd)))
    
    
    
    
#### plotting SI ~ sowniv ####    
    grand_mean.SI <- m.SI.sowndiv_gaus_p5 %>% 
      epred_draws(newdata = expand.grid(sowndivLogStd = seq(min(dat$sowndivLogStd), max(dat$sowndivLogStd), length.out = 20),
                                        treatment = c("1", "2", "3")),
                  re_formula = NA) # draws from the expectation of the predictive posterior distribution
                                   # corresponds to brms::posterior_epred()
    
    post_predict.SI <- m.SI.sowndiv_gaus_p5 %>% 
      predicted_draws(newdata = expand.grid(sowndivLogStd = seq(min(dat$sowndivLogStd), max(dat$sowndivLogStd), length.out = 20),
                                            treatment = c("1", "2", "3")),
                      re_formula = NA) # draws from the predictive posterior distribution
                                       # corresponds to brms::posterior_predict
    
    #aesthetics:#aesthetics:mean()
    treat_labs = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    
    p.SI.sowndiv <-  ggplot(grand_mean.SI, aes(x=sowndivLogStd, y=.epred, col=treatment))+
      stat_lineribbon(geom = "ribbon", data = post_predict.SI, 
                      aes(x=sowndivLogStd, y=.prediction, fill=treatment), 
                      .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment), 
                      .width = .95, alpha=0.3, colour=NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", 
                      aes( fill=NA, col=treatment#, linetype =c(rep("dashed", 6e4), rep("solid", 6e4), rep("dotted", 6e4) )
                      ), show.legend = FALSE)+
      geom_point(data = dat, aes(x=sowndivLogStd, y=SI), alpha = 0.6)+
      scale_y_continuous(name = "SI", 
                         breaks = c(0, 25, 50, 75, 100))+
      scale_x_continuous(name = "sown plant richness",
                         breaks = seq(min(dat$sowndivLogStd), max(dat$sowndivLogStd), length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$sowndivLogStd), max(dat$sowndivLogStd)) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      facet_wrap(~treatment, labeller = labeller(treatment = c("1" = "+SH+PH", "2" = "+SH-PH", "3" = "-SH-SH")))+
      theme_bw()+
      ggtitle("c",
              subtitle = "SI ~ sowndiv*treatment + (1|block/plot)")+
      theme(legend.position = "none",
            title = element_text(size = 9.5))
    
    p.SI.sowndiv  
    
    emmeans(m.SI.realdiv_gaus_p5, specs = c("treatment"), var = "realdivLogStd",
            at = list(realdivLogStd = min(dat$realdivLogStd)))
    
    
    
#### combine plots, save plots ####
  #regression against realized SR
    p.EI_SI.combined  <- ggarrange(p.EI_SI.wrap, ggarrange(p.EI.realdiv, p.SI.realdiv, nrow = 1, ncol = 2), 
              ncol = 1, nrow = 2)
    p.EI_SI.combined
    
    ggsave(p.EI_SI.combined, filename = "./plots/240224_EI_SI_real_combined.png",
           dpi=300, width = 15, height = 15, units = "cm")
    
    
  #regression against sown SR
    p.EI_SI.combined  <- ggarrange(p.EI_SI.wrap, ggarrange(p.EI.sowndiv, p.SI.sowndiv, nrow = 1, ncol = 2), 
                                   ncol = 1, nrow = 2)
    p.EI_SI.combined
    
    ggsave(p.EI_SI.combined, filename = "./plots/240224_EI_SI_sown_combined.png",
           dpi=300, width = 15, height = 15, units = "cm")
    
### pp_check plots EI SI ####
    load("./statistics/brms/240215_EI_mselect.RData")
    load("./statistics/brms/240215_SI_mselect.RData")
    
    pp.EI.r <- m.EI.realdiv_gaus_p5 %>%
      pp_check(ndraws = 100)+
      ylim(0,.045)+
      ggtitle("a: EI ~ realdivLogStd * treatment + (1|block/plot)")+
      #labs(caption = "SI ~ realdivLogStd * treatment + (1|block/plot)" )+
      theme(text = element_text(size=8))

    pp.SI.r <- m.SI.realdiv_gaus_p5 %>%
      pp_check(ndraws = 100)+
      ylim(0,.03)+
      ggtitle("b: SI ~ realdivLogStd * treatment + (1|block/plot)")+
      labs(caption = "SI ~ realdivLogStd * treatment + (1|block/plot)" )+
      theme(text = element_text(size=8))
    
    pp.EI.s <- m.EI.sowndiv_gaus_p5 %>%
      pp_check(ndraws = 100)+
      ylim(0,.045)+
      ggtitle("c: EI ~ sowndivLogStd * treatment + (1|block/plot)")+
      #labs(caption = "SI ~ realdivLogStd * treatment + (1|block/plot)" )+
      theme(text = element_text(size=8))

    pp.SI.s <- m.SI.sowndiv_gaus_p5 %>%
      pp_check(ndraws = 100)+
      ylim(0,.03)+
      ggtitle("d: SI ~ sowndivLogStd * treatment + (1|block/plot)")+
      labs(caption = "SI ~ sowndivLogStd * treatment + (1|block/plot)" )+
      theme(text = element_text(size=8))
    
    p.pp_check_EI.SI <- ggarrange(plotlist= list(pp.EI.r, pp.EI.s, pp.SI.r, pp.SI.s),
              ncol=2, nrow = 2)
    p.pp_check_EI.SI
    
    ggsave("./plots/240216_EI_SI_pp_checks.png", plot= p.pp_check_EI.SI,
           bg = "white",
           dpi=300, width = 16, height = 8, units = "cm") 
    