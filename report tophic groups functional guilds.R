#calculations to report the results on the densities of trophic groups and functional guilds:
library(tidybayes)

#### emmeans dens_all ####
  #abun_all
  emmeans(m.abun_all.sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=1)
          ) 
  
  #dens all, EMM: 
  emmeans(m.dens_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=min(dat$sowndivLogStd)),
          #epred = TRUE,
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp) 
  
  emmeans(m.dens_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=max(dat$sowndivLogStd)),
          #epred = TRUE,
          re_formula = NA) %>% as.data.frame() %>%170-
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  emtrends(m.dens_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>%
    pairs() %>% p_direction
  
  #epred = TRUE
  emmeans(m.dens_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=min(dat$sowndivLogStd)),
          epred = TRUE,
          re_formula = NA)
 
#### emmeans dens Om, Pr ####   
  #dens Om, EMM: 
  emmeans(m.Om_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=min(dat$sowndivLogStd)),
          #epred = TRUE,
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  emmeans(m.Om_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=max(dat$sowndivLogStd)),
          #epred = TRUE,
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
        
        #descriptive means and sd Om       
        mean(subset(dat, treatment == 1 & sowndiv == 1)$Om_per100g)
        sd(subset(dat, treatment == 3 & sowndiv == 1)$Om_per100g)
        #without zeros:
        mean(subset(dat, treatment == 1 & sowndiv == 1 & Om_per100g!=0)$Om_per100g)
        sd(subset(dat, treatment == 1 & sowndiv == 1 & Om_per100g!=0)$Om_per100g)
        
        mean(subset(dat, treatment == 1 & sowndiv == 16)$Om_per100g)
        sd(subset(dat, treatment == 3 & sowndiv == 16)$Om_per100g)
        #without zeros:
        mean(subset(dat, treatment == 1 & sowndiv == 16 & Om_per100g!=0)$Om_per100g)
        sd(subset(dat, treatment == 1 & sowndiv == 16 & Om_per100g!=0)$Om_per100g)
  
  #Pr
  predictions <- conditional_effects(m.Pr_sowndiv_p7, method = "posterior_epred")[[3]]
  predictions %>% 
    filter(effect1__ == min(dat$sowndivLogStd) | 
             effect1__ == max(dat$sowndivLogStd))
  
  predictions <- conditional_effects(m.Pr_realdiv_p7, method = "posterior_epred")[[3]]
  predictions %>% 
    filter(effect1__ == min(dat$realdivLogStd) | 
             effect1__ == max(dat$realdivLogStd))
  
####emmeans Ba ####
  
  emmeans(m.Ba_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=min(dat$sowndivLogStd)),
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  emmeans(m.Ba_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=max(dat$sowndivLogStd)),
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  #Ba1
  emmeans(m.Ba1_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=min(dat$sowndivLogStd)),
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  emmeans(m.Ba1_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=max(dat$sowndivLogStd)),
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  #Ba2
  emmeans(m.Ba2_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=min(dat$sowndivLogStd)),
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  emmeans(m.Ba2_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=max(dat$sowndivLogStd)),
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  #Ba4
  emmeans(m.Ba4_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=mean(dat$sowndivLogStd)),
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  emmeans(m.Ba4_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=max(dat$sowndivLogStd)),
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  #so a problem is that sum(Ba1, Ba2, Ba4) is bigger than overall Ba density --> 
  #models predict too high values for each functional guild, or too low values for the density of all bacterivores
  
  pp_check(m.Ba_sowndiv_p5, ndraws = 100) #predicted density of bacterivores is too low at high predicted Ba densities (at realistic values) 
  
  
  
####emmeans Fu ####  
  
  emmeans(m.Fu_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=mean(dat$sowndivLogStd)),
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  emmeans(m.Fu_sowndiv_p5, specs = c("treatment"),
          at = list(sowndivLogStd=max(dat$sowndivLogStd)),
          re_formula = NA) %>% as.data.frame() %>%
    mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
  
  #Fu2
    emmeans(m.Fu2_sowndiv_p5, specs = c("treatment"),
            at = list(sowndivLogStd=mean(dat$sowndivLogStd)),
            re_formula = NA) %>% as.data.frame() %>%
      mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
    
    emmeans(m.Fu2_sowndiv_p5, specs = c("treatment"),
            at = list(sowndivLogStd=max(dat$sowndivLogStd)),
            re_formula = NA) %>% as.data.frame() %>%
      mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
    
    #Fu3
    emmeans(m.Fu3_sowndiv_p5, specs = c("treatment"),
            at = list(sowndivLogStd=mean(dat$sowndivLogStd)),
            re_formula = NA) %>% as.data.frame() %>%
      mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
    
    emmeans(m.Fu3_sowndiv_p5, specs = c("treatment"),
            at = list(sowndivLogStd=max(dat$sowndivLogStd)),
            re_formula = NA) %>% as.data.frame() %>%
      mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
    
    #Fu4
    emmeans(m.Fu4_sowndiv_p5, specs = c("treatment"),
            at = list(sowndivLogStd=mean(dat$sowndivLogStd)),
            re_formula = NA) %>% as.data.frame() %>%
      mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
    
    emmeans(m.Fu4_sowndiv_p5, specs = c("treatment"),
            at = list(sowndivLogStd=max(dat$sowndivLogStd)),
            re_formula = NA) %>% as.data.frame() %>%
      mutate_at(vars(emmean, lower.HPD, upper.HPD), exp)
    

#### regression plot density all #####
    load("./statistics/brms/240225_abun_dens_all_mselect.RData")
    
    emm_draws.dens <- emmeans(m.dens_sowndiv_p5, specs = c("treatment", "sowndivLogStd"), var = "sowndivLogStd",
             at = list(sowndivLogStd = unique(dat$sowndivLogStd)), 
             #at = list(realdivLogStd = 0),
             re_formula = NA ) %>% 
      gather_emmeans_draws(level =0.95) 
      
    #aesthetics:
    treat_labs = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    #wrapped against treatment
    #p.dens.sown <-  ggplot(epred_draws, aes(x=sowndivLogStd, y=.epred, col=treatment))+
    p.dens.sown <-  ggplot(emm_draws.dens, aes(x=sowndivLogStd, y=exp(.value), col=treatment))+
      #stat_lineribbon(geom = "ribbon", data = post_predict_dens.sown, 
      #                aes(x=sowndivLogStd, y=.prediction, fill=treatment), 
      #                .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment), 
                      .width = .95, alpha=0.3, colour=NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", 
                      aes( fill=NA, col=treatment#, linetype =c(rep("dashed", 6e4), rep("solid", 6e4), rep("dotted", 6e4) )
                      ), show.legend = FALSE)+
      geom_point(data = dat, shape=19, alpha=0.6, aes(x=sowndivLogStd, y=ind_per100g))+
      scale_y_continuous(name = "nematodes per 100g ")+
      scale_x_continuous(name = "sown plant richness",
                         breaks = seq(min(dat$sowndivLogStd), max(dat$sowndivLogStd), length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$sowndivLogStd), max(dat$sowndivLogStd)) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      facet_wrap(~treatment, labeller = labeller(treatment = c("1" = "+SH+PH", "2" = "+SH-PH", "3" = "-SH-SH")))+
      theme_bw()+
      ggtitle("", subtitle = "density ~ sowndivLogStd * treatment + (1 | block/plot)")+
      theme(legend.position = "none")
    
    p.dens.sown
    
    ggsave(p.dens.sown, file = "./plots/240228_regression_all.png",
         dpi=300, width = 15, height = 8, units = "cm")
    
    
####average marginal effect plots for abundance of all nematodes as on andrew heiss' blogpost ####
  #https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/#continuous-effect
  
# emtrends()-based AMEs
  load("./statistics/brms/240225_abun_dens_all_mselect.RData")
  
  #this is not used:
  all_sown_ame.1unit <- emtrends(m.dens_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                                at = list(sowndivLogStd = min(dat$sowndivLogStd)),
                                re_formula = NA ) %>% 
    gather_emmeans_draws() %>%
    mutate(.value = 1/sd(dat$sowndivLog)*.value)
  
  cols=c("brown2", "darkolivegreen", "dodgerblue3")
  
  p.all_sown.ame <- ggplot(all_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
    ggtitle("", subtitle = "dens ~ sowndiv*treatment + (1|block/plot)")+
    stat_halfeye(slab_alpha = 0.3) +
    labs(x = NULL,
         y = "Density", fill = "treatment") +
    #xlim(-0.6, 0.6)+
    theme_bw() + 
    scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
    theme(legend.position = "none")
  p.all_sown.ame
  

  
  
  
####average marginal effect plots for densities of Pr and Om ~ realdiv ####
  load("./statistics/brms/240221_TrophDens_realdiv_mselect.RData")
  
  Pr_real_ame.1unit <- emtrends(m.Pr_realdiv_p5, specs = c("treatment"), var = "realdivLogStd",
                                  #at = list(sowndivLogStd = min(dat$realdivLogStd)),
                                  at = list(realdivLogStd = 0),
                                  re_formula = NA ) %>% 
    gather_emmeans_draws() %>%
    mutate(.value = 1/sd(dat$realdivLog)*.value) #this gives us the marginal effect of an increase by 1 unit of sowndiv 
 
      Om_real_ame.1unit <- emtrends(m.Om_realdiv_p5, specs = c("treatment"), var = "realdivLogStd",
                                    at = list(realdivLogStd = 0),
                                    re_formula = NA ) %>% 
        gather_emmeans_draws() %>%
        mutate(.value = 1/sd(dat$realdivLog)*.value)
          
          # what does epred=TRUE do in emtrends()?:
          m.Om_realdiv_p5 %>% emtrends(~ realdivLogStd*treatment ,
                                       var = "realdivLogStd", 
                                       at = list(realdivLogStd=c(-1.5,-0.5, 0.5,1.5),
                                                 treatment = c("1", "2", "3")),
                                        epred=TRUE, re_formula = NA) #this gives the absolute increase (on the response scale)
                                          #meaning: at realdivLogStd = -1.5 under treatment 1 the density of omnivores would increase by 
                                          #2.45 for a 1 SD increase in realdivLogStd
          
          m.Om_realdiv_p5 %>% emtrends(~ realdivLogStd*treatment ,
                                       var = "realdivLogStd", 
                                       at = list(realdivLogStd=c(-1.5,-0.5, 0.5,1.5),
                                                 treatment = c("1", "2", "3")),
                                       epred=FALSE, re_formula = NA) #this gives the percentual increase:
                                        #meaning: a 1 SD increase in realdivLogStd leads to a 16.2% increase 
                                        #in omnivore density (under treatment 1)  
          
  #plots:    
      cols=c("brown2", "darkolivegreen", "dodgerblue3")
      
      p.Pr_real.ame <- ggplot(Pr_real_ame.1unit, aes(x = .value, fill = factor(treatment))) +
        ggtitle("a", subtitle = "Pr ~ realdiv*treatment + (1|block/plot)")+
        stat_halfeye(slab_alpha = 0.3) +
        labs(x = NULL,
             y = "Density", fill = "treatment") +
        xlim(-0.6, 0.6)+
        theme_bw() + 
        scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
        theme(legend.position = "none")
      p.Pr_real.ame  
      
      p.Om_real.ame <- ggplot(Om_real_ame.1unit, aes(x = .value, fill = factor(treatment))) +
        ggtitle("b", subtitle = "Om ~ realdiv*treatment + (1|block/plot)")+
        stat_halfeye(slab_alpha = 0.3) +
        labs(x = NULL,
             y = "Density", fill = "treatment") +
        xlim(-0.6, 0.6)+
        theme_bw() + 
        scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
        labs(x = "Marginal effect of a duplication of sown plant richness \n on predatory (a) / omnivorous (b) nematode density",
             y = "Density", fill = "treatment") +
        theme(legend.position = "bottom")
      p.Om_real.ame
      
      
      
    ggarrange(p.Pr_real.ame, p.Om_real.ame, nrow = 2, ncol = 1, heights = c(1, 1.4))
      
#### average marginal effect plots for densities of Ba ~ sowndiv ####
load("./statistics/brms/240221_TrophDens_sowndiv_mselect.RData")

Ba_sown_ame.1unit <- emtrends(m.Ba_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                              at = list(sowndivLogStd = min(dat$sowndivLogStd)),
                              re_formula = NA ) %>% 
  gather_emmeans_draws() %>%
  mutate(.value = 1/sd(dat$sowndivLog)*.value) #this gives us the marginal effect of an increase by 1 unit of sowndiv 

    Ba1_sown_ame.1unit <- emtrends(m.Ba1_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                                  at = list(sowndivLogStd = min(dat$sowndivLogStd)),
                                  re_formula = NA ) %>% 
      gather_emmeans_draws() %>%
      mutate(.value = 1/sd(dat$sowndivLog)*.value)
    
    Ba2_sown_ame.1unit <- emtrends(m.Ba2_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                                  at = list(sowndivLogStd = min(dat$sowndivLogStd)),
                                  re_formula = NA ) %>% 
      gather_emmeans_draws() %>%
      mutate(.value = 1/sd(dat$sowndivLog)*.value)
    
    Ba4_sown_ame.1unit <- emtrends(m.Ba4_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                                  at = list(sowndivLogStd = min(dat$sowndivLogStd)),
                                  re_formula = NA ) %>% 
      gather_emmeans_draws() %>%
      mutate(.value = 1/sd(dat$sowndivLog)*.value)

#plots: 
  cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
  p.Ba_sown.ame <- ggplot(Ba_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
    ggtitle("a", subtitle = "Ba ~ sowndiv*treatment + (1|block/plot)")+
    stat_halfeye(slab_alpha = 0.3) +
    labs(x = NULL,
         y = "Density", fill = "treatment") +
    xlim(-0.6, 0.6)+
    theme_bw() + 
    scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
    theme(legend.position = "none")
  p.Ba_sown.ame  
    
    p.Ba1_sown.ame <- ggplot(Ba1_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
      ggtitle("b", subtitle = "Ba1 ~ sowndiv*treatment + (1|block/plot)")+
      stat_halfeye(slab_alpha = 0.3) +
      labs(x = NULL,
           y = "Density", fill = "treatment") +
      xlim(-0.6, 0.6)+
      theme_bw() + 
      scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
      theme(legend.position = "none")
    
    p.Ba1_sown.ame  
    
    p.Ba2_sown.ame <- ggplot(Ba2_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
      stat_halfeye(slab_alpha = 0.3) +
      ggtitle("c", subtitle = "Ba2 ~ sowndiv*treatment + (1|block/plot)")+
      labs(x = NULL,
           y = "Density", fill = "treatment") +
      xlim(-0.6, 0.6)+
      theme_bw() + 
      scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
      theme(legend.position = "none")
    p.Ba2_sown.ame  
    
    p.Ba4_sown.ame <- ggplot(Ba4_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
      stat_halfeye(slab_alpha = 0.3) +
      ggtitle("d", subtitle = "Ba4 ~ sowndiv*treatment + (1|block/plot)")+
      labs(x = "Marginal effect of a duplication of sown plant richness \n on bacterivorous nematode density",
           y = "Density", fill = "treatment") +
      xlim(-0.6, 0.6)+
      theme_bw() + 
      scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
      scale_color_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
      theme(legend.position = "bottom")
    p.Ba4_sown.ame  

    ggarrange(p.Ba_sown.ame, p.Ba1_sown.ame, p.Ba2_sown.ame, p.Ba4_sown.ame, ncol = 1, nrow = 4, heights = c(1,1,1,1.5))    

    
#### average marginal effect plots for densities of Fu ~ sowndiv####
  #gather draws from emtrends
    Fu_sown_ame.1unit <- emtrends(m.Fu_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                                  #at = list(sowndivLogStd = min(dat$sowndivLogStd)),
                                  at = list(sowndivLogStd = 0),
                                  re_formula = NA ) %>% 
      gather_emmeans_draws() %>%
      mutate(.value = 1/sd(dat$sowndivLog)*.value) #this gives us the marginal effect of an increase by 1 unit of sowndiv 
    
        Fu2_sown_ame.1unit <- emtrends(m.Fu2_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                                      #at = list(sowndivLogStd = min(dat$sowndivLogStd)),
                                      at = list(sowndivLogStd = 0),
                                      re_formula = NA ) %>% 
          gather_emmeans_draws() %>%
          mutate(.value = 1/sd(dat$sowndivLog)*.value) 
        
        Fu3_sown_ame.1unit <- emtrends(m.Fu3_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                                       #at = list(sowndivLogStd = min(dat$sowndivLogStd)),
                                       at = list(sowndivLogStd = 0),
                                       re_formula = NA ) %>% 
          gather_emmeans_draws() %>%
          mutate(.value = 1/sd(dat$sowndivLog)*.value) 
        
        Fu4_sown_ame.1unit <- emtrends(m.Fu4_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                                       #at = list(sowndivLogStd = min(dat$sowndivLogStd)),
                                       at = list(sowndivLogStd = 0),
                                       re_formula = NA ) %>% 
          gather_emmeans_draws() %>%
          mutate(.value = 1/sd(dat$sowndivLog)*.value) 
    
  #single density plots:  
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
        
    p.Fu_sown.ame <- ggplot(Fu_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
      stat_halfeye(slab_alpha = 0.3) +
      ggtitle("a", subtitle = "Fu ~ sowndiv*treatment + (1|block/plot)")+
      labs(x = NULL,
           y = "Density", fill = "treatment") +
      xlim(-0.6, 0.6)+
      theme_bw() + 
      scale_size(1)+
      scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
      theme(legend.position = "none")
    p.Fu_sown.ame  
    
        p.Fu2_sown.ame <- ggplot(Fu2_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
          stat_halfeye(slab_alpha = 0.3) +
          ggtitle("b", subtitle = "Fu2 ~ sowndiv*treatment + (1|block/plot)")+
          labs(x = NULL,
               y = "Density", fill = "treatment") +
          xlim(-0.6, 0.6)+
          theme_bw() + 
          scale_size(1)+ 
          scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
          theme(legend.position = "none")
        p.Fu2_sown.ame  
        
        p.Fu3_sown.ame <- ggplot(Fu3_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
          stat_halfeye(slab_alpha = 0.3) +
          ggtitle("c", subtitle = "Fu3 ~ sowndiv*treatment + (1|block/plot)")+
          labs(x = NULL,
               y = "Density", fill = "treatment") +
          xlim(-0.6, 0.6)+
          theme_bw() + 
          scale_size(1)+
          scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
          theme(legend.position = "none")
        p.Fu3_sown.ame  
        
        p.Fu4_sown.ame <- ggplot(Fu4_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
          stat_halfeye(slab_alpha = 0.3) +
          ggtitle("d", "Fu4 ~ sowndiv*treatment + (1|block/plot)")+
          labs(x = "Marginal effect of a duplication of sown plant richness \n on fungivorous nematode density",
               y = "Density", fill = "treatment") +
          xlim(-0.6, 0.6)+
          theme_bw() + 
          scale_size(1)+
          scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
          scale_color_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
          theme(legend.position = "bottom")
        p.Fu4_sown.ame  

        ggarrange(p.Fu_sown.ame, p.Fu2_sown.ame, p.Fu3_sown.ame, p.Fu4_sown.ame, 
                  ncol=1, nrow=4, heights = c(1,1,1,1.45))    
        

#### average marginal effects of all trophic groups ~ sowndiv ####
  #use ame of Ba and Fu from before
  #calculate ame of Pl:
    Pl_sown_ame.1unit <- emtrends(m.Pl_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                                  #at = list(sowndivLogStd = min(dat$sowndivLogStd)),
                                  at = list(sowndivLogStd = 0),
                                  re_formula = NA ) %>% 
      gather_emmeans_draws() %>%
      mutate(.value = 1/sd(dat$sowndivLog)*.value) #this gives us the marginal effect of an increase by 1 unit of sowndiv 
    
    Pr_sown_ame.1unit <- emtrends(m.Pr_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                                  #at = list(sowndivLogStd = min(dat$realdivLogStd)),
                                  at = list(sowndivLogStd = 0),
                                  re_formula = NA ) %>% 
      gather_emmeans_draws() %>%
      mutate(.value = 1/sd(dat$realdivLog)*.value) #this gives us the marginal effect of an increase by 1 unit of sowndiv 
    
    Om_sown_ame.1unit <- emtrends(m.Om_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd",
                                  at = list(sowndivLogStd = 0),
                                  re_formula = NA ) %>% 
      gather_emmeans_draws() %>%
      mutate(.value = 1/sd(dat$sowndivLog)*.value)    
    
        
    
  #the plots:
    p.BaALL_sown.ame <- ggplot(Ba_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
      ggtitle("a", subtitle = "Ba ~ sowndiv*treatment + (1|block/plot)")+
      stat_halfeye(slab_alpha = 0.3) +
      labs(x = NULL,
           y = "Density", fill = "treatment") +
      xlim(-0.6, 0.6)+
      theme_bw() + 
      scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
      theme(legend.position = "none")
    p.BaALL_sown.ame      
        
        
    p.FuALL_sown.ame <- ggplot(Fu_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
      stat_halfeye(slab_alpha = 0.3) +
      ggtitle("b", subtitle = "Fu ~ sowndiv*treatment + (1|block/plot)")+
      labs(x = NULL,
           y = "Density", fill = "treatment") +
      xlim(-0.6, 0.6)+
      theme_bw() + 
      scale_size(1)+
      scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
      theme(legend.position = "none")
    p.FuALL_sown.ame  
    
    
    p.PlALL_sown.ame <- ggplot(Pl_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
      ggtitle("c", subtitle = "Pl ~ sowndiv*treatment + (1|block/plot)")+
      stat_halfeye(slab_alpha = 0.3) +
      labs(x = NULL,
           y = "Density", fill = "treatment") +
      xlim(-0.6, 0.6)+
      theme_bw() + 
      scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
      theme(legend.position = "none")
    p.PlALL_sown.ame  
    
    
    p.PrALL_sown.ame <- ggplot(Pr_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
      ggtitle("d", subtitle = "Pr ~ sowndiv*treatment + (1|block/plot)")+
      stat_halfeye(slab_alpha = 0.3) +
      labs(x = NULL,
           y = "Density", fill = "treatment") +
      xlim(-0.6, 0.6)+
      theme_bw() + 
      scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
      theme(legend.position = "none")
    p.PrALL_sown.ame  
        
    p.OmALL_sown.ame <- ggplot(Om_sown_ame.1unit, aes(x = .value, fill = factor(treatment))) +
      ggtitle("e", subtitle = "Om ~ sowndiv*treatment + (1|block/plot)")+
      stat_halfeye(slab_alpha = 0.3) +
      labs(x = NULL,
           y = "Density", fill = "treatment") +
      xlim(-0.6, 0.6)+
      theme_bw() + 
      scale_fill_manual(labels = c("+SH+PH", "+SH-PH", "-SH-PH" ), values = cols)+
      labs(x = "Marginal effect of a doubling of sown plant richness on nematode density",
           y = "Density", fill = "treatment") +
      theme(legend.position = "bottom")
    p.OmALL_sown.ame
    
    p.legend.ameALL <- ggplot()+
      labs(x = "Marginal effect of a duplication of sown plant richness on the \ndensity of bacterivorous (a), fungivorous (b), herbivorous (c), \npredatory (d), and omnivorous (e) nematodes")
      
        
    
    
          
        
#### arrange density plots ####
     p.ame.PrOm_real <- ggarrange(p.Pr_real.ame, p.Om_real.ame, 
                                 nrow = 2, ncol = 1, heights = c(1, 1.4))  
        p.ame.PrOm_real
        
     p.ame.Fu_sown <-   ggarrange(p.Fu_sown.ame, p.Fu2_sown.ame, p.Fu3_sown.ame, p.Fu4_sown.ame, 
                                  ncol=1, nrow=4, heights = c(1,1,1,1.5))        
        p.ame.Fu_sown 
        
     p.ame.Ba_sown <-   ggarrange(p.Ba_sown.ame, p.Ba1_sown.ame, p.Ba2_sown.ame, p.Ba4_sown.ame, 
                                  ncol=1, nrow=4, heights = c(1,1,1,1.45))        
        p.ame.Ba_sown 
        
      p.ame.troph_sown <- ggarrange(p.BaALL_sown.ame, p.FuALL_sown.ame, p.PlALL_sown.ame, p.PrALL_sown.ame, p.OmALL_sown.ame,
                                    ncol=1, nrow=5, heights = c(1,1,1,1,1.6))  
      p.ame.troph_sown 
        
    ggsave(p.ame.Fu_sown, filename = "./plots/240226_ame_Fu_sown.png",
           dpi=300, width = 15, height = 25, units = "cm")    
    
    ggsave(p.ame.Ba_sown, filename = "./plots/240226_ame_Ba_sown.png",
           dpi=300, width = 15, height = 25, units = "cm") 
    
    ggsave(p.ame.PrOm_real, filename = "./plots/240226_ame_PrOm_real.png",
           dpi=300, width = 15, height = 12.5, units = "cm")   
    
    ggsave(p.ame.troph_sown, filename = "./plots/240327_ame_trophdens_sown.png",
           dpi=300, width = 15, height = 22, units = "cm") 
    
    
####OLD: regression plots ~sowndiv ####
    #Pr
    predictions <- conditional_effects(m.Pr_realdiv_p7)[[3]]
    treatments = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    ggplot(dat, aes(x=realdivLogStd, y=Pr_per100g) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_line(data=predictions, aes(x= realdivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      
      geom_jitter(width=0.02, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      scale_color_manual(labels=treatments, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "Pr per 100g")+
      theme_classic()+
      theme(legend.position ="bottom")
    
    #Ba
    predictions <- conditional_effects(m.Ba4_sowndiv_p5)[[3]]
    treatments = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    ggplot(dat, aes(x=sowndivLogStd, y=Ba4) )+
      geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                        fill=treatment), 
                  alpha=0.2, show.legend=FALSE)+
      geom_line(data=predictions, aes(x= sowndivLogStd, y=estimate__, 
                                      linetype=treatment, col=treatment),
                linewidth= 1, show.legend = FALSE)+
      
      geom_jitter(width=0.2, shape=19, alpha=0.4, 
                  aes(col=treatment))+
      scale_color_manual(labels=treatments, values = cols)+
      scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                         labels = c("1", "2", "4", "8", "16"))+
      scale_y_continuous(name = "Ba4 per 100g")+
      theme_classic()+
      theme(legend.position ="bottom")
    
    conditional_effects(m.Ba4_sowndiv_p5)
    
    
    
    
    
    
#### OLD: pd between treatments, ~sowndiv ####
    #trophic groups:
    emtrends(m.Ba_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none
    
    emtrends(m.Fu_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none (t1-t3: 94.07%)
    
    emtrends(m.Pl_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none
    
    emtrends(m.Pr_sowndiv_p7, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #t1-t2: 97.57%, t1-t3: 99.43%
    conditional_effects(m.Pr_sowndiv_p7)
    
    emtrends(m.Om_sowndiv_p7, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #t2-t3: 97.37%,  (t1-t2: 94.87%)
    
    #functional guilds:
    emtrends(m.Ba1_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none
    emtrends(m.Ba2_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none
    emtrends(m.Ba4_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #t1:t3 99.57%, t2-t3: 99.43%
    
    emtrends(m.Fu2_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none (t1-t3 92.60%)
    emtrends(m.Fu3_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none 
    emtrends(m.Fu4_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #t1-t2: 97.03% (t1-t3 92.30%)
    
    emtrends(m.Pl2_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none (t1-t3 95.00%)
    emtrends(m.Pl3_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none
    emtrends(m.Pl4_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none 
    
    emtrends(m.Pr4_Om4_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #t1-t2: 99.60%, t1-t3: 96.20%
    
#### OLD: pd between treatments, ~realdiv ####
    #trophic groups:
    emtrends(m.Ba_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none
    
    emtrends(m.Fu_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none (t1-t3: 92.87%)
    
    emtrends(m.Pl_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none
    
    emtrends(m.Pr_realdiv_p7, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #t1-t2: 97.77%, t1-t3: 99.43%
    conditional_effects(m.Pr_realdiv_p7)
    
    emtrends(m.Om_realdiv_p7, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #t2-t3: 99.13%,  t1-t2: 98.23%
    conditional_effects(m.Om_realdiv_p7)
    
    #functional guilds:
    emtrends(m.Ba1_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none
    emtrends(m.Ba2_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none
    emtrends(m.Ba4_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #t1:t3 99.73%, t2-t3: 99.37%
    
    emtrends(m.Fu2_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none (t1-t3 91.07%)
    emtrends(m.Fu3_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none 
    emtrends(m.Fu4_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #t1-t2: 96.13% 
    
    emtrends(m.Pl2_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none (t1-t3: 77.63%, compared to 95% with sowndiV!)
    emtrends(m.Pl3_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none
    emtrends(m.Pl4_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #none 
    
    emtrends(m.Pr4_Om4_realdiv_p5, specs = c("treatment"), var = "realdivLogStd") %>% 
      pairs() %>%
      bayestestR::p_direction() #t1-t2: 99.77%, t1-t3: 97.20%  
    