#report Hill numbers
#load("./statistics/brms/240214_Hill_sowndiv_mselect.RData")
library(ggplot2)
library(viridis)
library(ggpubr)
  
  Hill.sown.summary <- read.xlsx("./statistics/240214_Model_summaries.xlsx", sheet = 'Hill numbers ~ sowndiv')
  Hill.real.summary <- read.xlsx("./statistics/240214_Model_summaries.xlsx", sheet = 'Hill numbers ~ realdiv')
  
  Hill.sown.summary %>% select(response, treatment) %>% mutate(diff.real_sown = Hill.real.summary$mean.trend - Hill.sown.summary$mean.trend,
                                                               pd.sown = Hill.sown.summary$pd,
                                                               pd.real = Hill.real.summary$pd,
                                                               diff.pd = Hill.real.summary$pd - Hill.sown.summary$pd)
  Hill.real.summary$mean.trend - Hill.sown.summary$mean.trend

#data preparation:
  dat <- subset(dBEF_nem21, sowndiv != 60) 
  dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                        .after = sowndivLog) %>%
    mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
           .after = realdivLog) %>%
    mutate(funcdivStd = ( (funcdiv - mean(funcdiv)) / sd(funcdiv) ),
           .after = funcdiv)
#### emmmeans for reporting, sowndiv####
  load("./statistics/brms/240214_Hill_sowndiv_mselect.RData")
  emtrends(m.all.Shannon.gaus_p5, specs = c("treatment"), var = "sowndivLogStd") %>% pairs() %>% p_direction()
  
  emmeans(m.all.Shannon.gaus_p5, specs = c("treatment"),
          at = list(sowndivLogStd=min(dat$sowndivLogStd)),
          #epred = TRUE,
          re_formula = NA) 
  
  emmeans(m.all.Shannon.gaus_p5, specs = c("treatment"),
          at = list(sowndivLogStd=max(dat$sowndivLogStd)),
          #epred = TRUE,
          re_formula = NA) 
  
#### emmmeans for reporting, realdiv####
  load("./statistics/brms/240214_Hill_realdiv_mselect.RData")
  emtrends(m.all.Shannon.gaus_p5, specs = c("treatment"), var = "realdivLogStd") %>% pairs() %>% p_direction()

  
#### plotting exp(H') ~ realdiv ####
load("./statistics/brms/240214_Hill_realdiv_mselect.RData")

grand_mean.Hill <- m.all.Shannon.gaus_p5 %>% 
  epred_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                    treatment = c("1", "2", "3")),
              re_formula = NA) #expectation of the predictive posterior distribution

post_predict.Hill <- m.all.Shannon.gaus_p5 %>% 
  predicted_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                        treatment = c("1", "2", "3")),
                  re_formula = NA) #predictive posterior distribution
#aesthetics:
treat_labs = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")

#wrapped against treatment
    p.Hill.realdiv.wrap <-  ggplot(grand_mean.Hill, aes(x=realdivLogStd, y=.epred, col=treatment))+
      stat_lineribbon(geom = "ribbon", data = post_predict.Hill, 
                      aes(x=realdivLogStd, y=.prediction, fill=treatment), 
                      .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment), 
                      .width = .95, alpha=0.3, colour=NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", 
                      aes( fill=NA, col=treatment#, linetype =c(rep("dashed", 6e4), rep("solid", 6e4), rep("dotted", 6e4) )
                      ), show.legend = FALSE)+
      geom_point(data = dat, shape=19, alpha=0.6, aes(x=realdivLogStd, y=Hill_q1))+
      scale_y_continuous(name = "exp(H')")+
      scale_x_continuous(name = "realized plant richness",
                         breaks = seq(min(dat$realdivLogStd), 1.916, length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$realdivLogStd), 1.92) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      facet_wrap(~treatment, labeller = labeller(treatment = c("1" = "+SH+PH", "2" = "+SH-PH", "3" = "-SH-SH")))+
      theme_bw()+
      ggtitle("Effective number of families")+
      theme(legend.position = "none")
    
    p.Hill.realdiv.wrap
    
    ggsave("./plots/240216_Hill_realdiv_wrapped.png", plot=p.Hill.realdiv.wrap,
           dpi=300, width = 15, height = 8, units = "cm")

#overlapped:   
    
    p.Hill.realdiv.overl <-  ggplot(grand_mean.Hill, aes(x=realdivLogStd, y=.epred, col=treatment))+
      #stat_lineribbon(geom = "ribbon", data = post_predict.Hill, 
       #               aes(x=realdivLogStd, y=.prediction, fill=treatment), 
       #               .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment, colour=NA), 
                      .width = .95, alpha=0.2, colour = NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", data = subset(grand_mean.Hill, treatment == 1),
                      aes(fill=NA, colour=treatment), show.legend = FALSE)+
      stat_lineribbon(geom = "line", data = subset(grand_mean.Hill, treatment == 2),
                      aes(fill=NA, colour=treatment), show.legend = FALSE)+
      stat_lineribbon(geom = "line", data = subset(grand_mean.Hill, treatment == 3),
                      aes(fill=NA, colour=treatment), show.legend = FALSE)+
      geom_point(data = dat, shape=19, alpha=0.6, aes(x=realdivLogStd, y=Hill_q1))+
      scale_y_continuous(name = "exp(H')")+
      scale_x_continuous(name = "realized plant richness",
                         breaks = seq(min(dat$realdivLogStd), 1.916, length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$realdivLogStd), 1.92) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      theme(legend.position ="bottom")+
      theme_bw()+
      ggtitle("Effective number of families")+
      theme(legend.position = "none")
    
    p.Hill.realdiv.overl 
    ggsave("./plots/240216_Hill_realdiv_overlapped.png", plot= p.Hill.realdiv.overl,
           dpi=300, width = 15, height = 8, units = "cm")
    

#### plotting exp(H') ~ sowndiv ####
    load("./statistics/brms/240214_Hill_sowndiv_mselect.RData")
    
    grand_mean.Hill <- m.all.Shannon.gaus_p5 %>% 
      epred_draws(newdata = expand.grid(sowndivLogStd = seq(min(dat$sowndivLogStd), max(dat$sowndivLogStd), length.out = 20),
                                        treatment = c("1", "2", "3")),
                  re_formula = NA) #expectation of the predictive posterior distribution
    
    post_predict.Hill <- m.all.Shannon.gaus_p5 %>% 
      predicted_draws(newdata = expand.grid(sowndivLogStd = seq(min(dat$sowndivLogStd), max(dat$sowndivLogStd), length.out = 20),
                                            treatment = c("1", "2", "3")),
                      re_formula = NA) #predictive posterior distribution
    
    emm_draws.Hill.sown <- emmeans(m.all.Shannon.gaus_p5, specs = c("treatment", "sowndivLogStd"), var = "sowndivLogStd",
                              at = list(sowndivLogStd = unique(dat$sowndivLogStd)),
                              re_formula = NA ) %>% 
      gather_emmeans_draws() #draws from the marginal posterior distribution (see documentation of gather_emmeans_draws)
    
    #aesthetics:
    treat_labs = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    #wrapped against treatment
    p.Hill.sowndiv.wrap <-  ggplot(emm_draws.Hill.sown , aes(x=sowndivLogStd, y=.value, col=treatment))+
      #stat_lineribbon(geom = "ribbon", data =  grand_mean.Hill, 
      #                aes(x=sowndivLogStd, y=.epred, fill=treatment), 
      #                .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment), 
                      .width = .95, alpha=0.3, colour=NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", 
                      aes( fill=NA, col=treatment#, linetype =c(rep("dashed", 6e4), rep("solid", 6e4), rep("dotted", 6e4) )
                      ), show.legend = FALSE)+
      geom_point(data = dat, shape=19, alpha=0.6, aes(x=sowndivLogStd, y=Hill_q1))+
      scale_y_continuous(name = "exp(H')")+
      scale_x_continuous(name = "sown plant richness",
                         breaks = seq(min(dat$sowndivLogStd), max(dat$sowndivLogStd), length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$sowndivLogStd), max(dat$sowndivLogStd)) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      facet_wrap(~treatment, labeller = labeller(treatment = c("1" = "+SH+PH", "2" = "+SH-PH", "3" = "-SH-SH")))+
      theme_bw()+
      ggtitle("Effective number of families", subtitle = "exp(H') ~ sowndivLogStd * treatment + (1 | block/plot)")+
      theme(legend.position = "none")
    
    p.Hill.sowndiv.wrap
    
    ggsave("./plots/240229_Hill_sowndiv_wrapped.png", plot=p.Hill.sowndiv.wrap,
           dpi=300, width = 15, height = 8, units = "cm")
    
    #overlapped:   
    
    p.Hill.sowndiv.overl <-  ggplot(grand_mean.Hill, aes(x=sowndivLogStd, y=.epred, col=treatment))+
      #stat_lineribbon(geom = "ribbon", data = post_predict.Hill, 
      #               aes(x=sowndivLogStd, y=.prediction, fill=treatment), 
      #               .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment, colour=NA), 
                      .width = .95, alpha=0.2, colour = NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", data = subset(grand_mean.Hill, treatment == 1),
                      aes(fill=NA, colour=treatment), show.legend = FALSE)+
      stat_lineribbon(geom = "line", data = subset(grand_mean.Hill, treatment == 2),
                      aes(fill=NA, colour=treatment), show.legend = FALSE)+
      stat_lineribbon(geom = "line", data = subset(grand_mean.Hill, treatment == 3),
                      aes(fill=NA, colour=treatment), show.legend = FALSE)+
      geom_point(data = dat, shape=19, alpha=0.6, aes(x=sowndivLogStd, y=Hill_q1))+
      scale_y_continuous(name = "exp(H')")+
      scale_x_continuous(name = "sownized plant richness",
                         breaks = seq(min(dat$sowndivLogStd), max(dat$sowndivLogStd), length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$sowndivLogStd), max(dat$sowndivLogStd)) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      theme(legend.position ="bottom")+
      theme_bw()+
      ggtitle("Effective number of families")+
      theme(legend.position = "none")
    
    p.Hill.sowndiv.overl 
    ggsave("./plots/240216_Hill_sowndiv_overlapped.png", plot= p.Hill.sowndiv.overl,
           dpi=300, width = 15, height = 8, units = "cm")
    
#### number of effective families in Ba ####
    load("./statistics/brms/240214_Hill_realdiv_mselect.RData")
    
    grand_mean.Hill.Ba <- m.Ba.Shannon.gamma_p5 %>% 
      epred_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                        treatment = c("1", "2", "3")),
                  re_formula = NA) #expectation of the predictive posterior distribution
    
    post_predict.Hill.Ba <- m.Ba.Shannon.gamma_p5 %>% 
      predicted_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                            treatment = c("1", "2", "3")),
                      re_formula = NA) #predictive posterior distribution
    #aesthetics:
    treat_labs = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    #wrapped against treatment
    p.Hill.Ba.realdiv.wrap <-  ggplot(grand_mean.Hill.Ba, aes(x=realdivLogStd, y=.epred, col=treatment))+
      stat_lineribbon(geom = "ribbon", data = post_predict.Hill.Ba, 
                      aes(x=realdivLogStd, y=.prediction, fill=treatment), 
                      .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment), 
                      .width = .95, alpha=0.3, colour=NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", 
                      aes( fill=NA, col=treatment#, linetype =c(rep("dashed", 6e4), rep("solid", 6e4), rep("dotted", 6e4) )
                      ), show.legend = FALSE)+
      geom_point(data = dat, shape=19, alpha=0.6, aes(x=realdivLogStd, y=Hill_q1.Ba))+
      scale_y_continuous(name = "exp(H')")+
      scale_x_continuous(name = "realized plant richness",
                         breaks = seq(min(dat$realdivLogStd), 1.916, length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$realdivLogStd), 1.92) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      facet_wrap(~treatment, labeller = labeller(treatment = c("1" = "+SH+PH", "2" = "+SH-PH", "3" = "-SH-SH")))+
      theme_bw()+
      ggtitle("Effective number of bacterivorous families")+
      theme(legend.position = "none")
    
    p.Hill.Ba.realdiv.wrap
    
    ggsave("./plots/240216_Hill_Ba_realdiv_wrapped.png", plot=p.Hill.Ba.realdiv.wrap,
           dpi=300, width = 15, height = 8, units = "cm")
    
#### number of effective families in Pr and Om ####    
    load("./statistics/brms/240214_Hill_realdiv_mselect.RData")
    
    grand_mean.Hill.PrOm <- m.PrOm.Shannon.gamma_p5 %>% 
      epred_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                        treatment = c("1", "2", "3")),
                  re_formula = NA) #expectation of the predictive posterior distribution
    
    post_predict.Hill.PrOm <- m.PrOm.Shannon.gamma_p5 %>% 
      predicted_draws(newdata = expand.grid(realdivLogStd = seq(min(dat$realdivLogStd), max(dat$realdivLogStd), length.out = 20),
                                            treatment = c("1", "2", "3")),
                      re_formula = NA) #predictive posterior distribution
    #aesthetics:
    treat_labs = c("+SH+PH", "+SH-PH", "-SH-PH")
    cols=c("brown2", "darkolivegreen", "dodgerblue3")
    
    #wrapped against treatment
    p.Hill.PrOm.realdiv.wrap <-  ggplot(grand_mean.Hill.PrOm, aes(x=realdivLogStd, y=.epred, col=treatment))+
      stat_lineribbon(geom = "ribbon", data = post_predict.Hill.PrOm, 
                      aes(x=realdivLogStd, y=.prediction, fill=treatment), 
                      .width=.95, alpha =0.1, colour = NA, show.legend = FALSE)+ 
      stat_lineribbon(geom = "lineribbon", 
                      aes( fill= treatment), 
                      .width = .95, alpha=0.3, colour=NA, show.legend = FALSE)+
      stat_lineribbon(geom = "line", 
                      aes( fill=NA, col=treatment#, linetype =c(rep("dashed", 6e4), rep("solid", 6e4), rep("dotted", 6e4) )
                      ), show.legend = FALSE)+
      geom_point(data = dat, shape=19, alpha=0.6, aes(x=realdivLogStd, y=Hill_q1.PrOm))+
      scale_y_continuous(name = "exp(H')")+
      scale_x_continuous(name = "realized plant richness",
                         breaks = seq(min(dat$realdivLogStd), 1.916, length.out=5) %>% round(3),
                         labels = c("1", "2", "4", "8", "16"),
                         limits = c(min(dat$realdivLogStd), 1.92) )+
      scale_color_manual(labels=treat_labs, values = cols)+
      facet_wrap(~treatment, labeller = labeller(treatment = c("1" = "+SH+PH", "2" = "+SH-PH", "3" = "-SH-SH")))+
      theme_bw()+
      ggtitle("Effective number of families in the predatory and omnivorous trophic groups")+
      theme(legend.position = "none")
    
    p.Hill.PrOm.realdiv.wrap
    
    ggsave("./plots/240216_Hill_PrOm_realdiv_wrapped.png", plot=p.Hill.PrOm.realdiv.wrap,
           dpi=300, width = 15, height = 8, units = "cm")

####pp_check plots####
    load("./statistics/brms/240214_Hill_realdiv_mselect.RData")
    
    pp.all.r <- m.all.Shannon.gaus_p5 %>%
      pp_check(ndraws = 100)+
      ggtitle("a: exp(H') ~ realdivLogStd * treatment + (1|block/plot)")+
      theme(text = element_text(size=6))
    
    pp.Ba.r <- m.Ba.Shannon.gamma_p5 %>%
      pp_check(ndraws = 100)+
      ggtitle(expression("b: exp(H')"[Ba]*" ~ realdivLogStd * treatment + (1|block/plot)"))+
      theme(text = element_text(size=6))
    
    pp.Fu.r <- m.Fu.Shannon.gaus_p5 %>%
      pp_check(ndraws = 100)+
      ggtitle(expression("c: exp(H')"[Fu]*" ~ realdivLogStd * treatment + (1|block/plot)"))+
      theme(text = element_text(size=6))
    
    pp.Pl.r <- m.Pl.Shannon.gamma_p5 %>%
      pp_check(ndraws = 100)+
      ggtitle(expression("d: exp(H')"[Pl]*" ~ realdivLogStd * treatment + (1|block/plot)"))+
      theme(text = element_text(size=6))
    
    pp.PrOm.r <- m.PrOm.Shannon.gamma_p5 %>%
      pp_check(ndraws = 100)+
      ggtitle(expression("e: exp(H')"[Pr+Om]*" ~ realdivLogStd * treatment + (1|block/plot)"))+
      theme(text = element_text(size=6))
    
    load("./statistics/brms/240214_Hill_sowndiv_mselect.RData")
    #this overwrites every loaded realdiv model with its sowndiv counterpart, so be careful when re-running code from above!
    
    pp.all.s <- m.all.Shannon.gaus_p5 %>%
      pp_check(ndraws = 100)+
      ggtitle("f: exp(H') ~ sowndivLogStd * treatment + (1|block/plot)")+
      theme(text = element_text(size=6))
    
    pp.Ba.s <- m.Ba.Shannon.gamma_p5 %>%
      pp_check(ndraws = 100)+
      ggtitle(expression("g: exp(H')"[Ba]*" ~ sowndivLogStd * treatment + (1|block/plot)"))+
      theme(text = element_text(size=6))
    
    pp.Fu.s <- m.Fu.Shannon.gaus_p5 %>%
      pp_check(ndraws = 100)+
      ggtitle(expression("h: exp(H')"[Fu]*" ~ sowndivLogStd * treatment + (1|block/plot)"))+
      theme(text = element_text(size=6))
    
    pp.Pl.s <- m.Pl.Shannon.gamma_p5 %>%
      pp_check(ndraws = 100)+
      ggtitle(expression("i: exp(H')"[Pl]*" ~ sowndivLogStd * treatment + (1|block/plot)"))+
      theme(text = element_text(size=6))
    
    pp.PrOm.s <- m.PrOm.Shannon.gamma_p5 %>%
      pp_check(ndraws = 100)+
      ggtitle(expression("j: exp(H')"[Pr+Om]*" ~ sowndivLogStd * treatment + (1|block/plot)"))+
      theme(text = element_text(size=6))

p.pp_check_Shannon <- ggarrange(plotlist= list(pp.all.r, pp.all.s, pp.Ba.r, pp.Ba.s, pp.Fu.r, 
                                               pp.Fu.s, pp.Pl.r, pp.Pl.s, pp.PrOm.r, pp.PrOm.s),
                                ncol=2, nrow = 5)

p.pp_check_Shannon    

ggsave("./plots/240216_Hill_pp_checks.png", plot= p.pp_check_Shannon ,
       bg = "white",
       dpi=300, width = 16, height = 16, units = "cm") 
        
    
    