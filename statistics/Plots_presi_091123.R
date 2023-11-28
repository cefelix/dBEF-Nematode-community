library(brms)
library(tidyr)
library(ggplot2)
library(tidybayes)
library(modelr)
library(emmeans)

#run variable transformations.R before this script!

#how to extract data from conditonal_effects():
  #https://discourse.mc-stan.org/t/change-aesthetics-conditional-effects/13500/5

#apparently, conditional_effects "computes the mean of posterior predictive distribution", 
#thus it "includes both hurdle and pure log normal part" (Paul Buerkner):
#https://discourse.mc-stan.org/t/brms-does-the-lognormal-part-of-the-hurdle-lognormal-regression-include-zeros-into-analysis/18136/3


####Predators 2021####
load("./statistics/brms/231117_Pr.RData")

#linear regression from 1-60 species:  
  predictions <- conditional_effects(m.Pr.hurdle21b)[[3]]
  predictions$estimate__

  ggplot(dBEF_nem21, aes(x=sowndivLog, y=Pr_per100g, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=estimate__, col=treatment))
  
  
#from 1-16 species  
  predictions <- conditional_effects(m.Fu.hurdle21a)[[3]]
  predictions$estimate__
  
  ggplot(dBEF_nem21, aes(x=sowndivLog, y=Fu_per100g, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=estimate__, col=treatment))

  
####Predators 2017####
load("./statistics/17_remodelling/brms/231108_Pr17_hurdle.RData")
#the good one is 51, log(Pr~log(sowndiv))
  predictions <- conditional_effects(m17.Pr.hurdle51)[[3]]
  pp_check(m17.Pr.hurdle51, ndraws=100)
  
  dBEF_nem17 %>% filter(treatment==3)%>%
    ggplot(aes(x=sowndivLog, y=Pr_per100g))+
    geom_jitter()+
    geom_smooth(method="lm")
  #This might mean that i should use treatment~hurdle instead of 1~hurdle

  ggplot(dBEF_nem17, aes(x=sowndivLog, y=Pr_per100gLog, col=treatment) )+
    geom_jitter(position = position_dodge(width=0.5))+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=estimate__, col=treatment))
    #thats a poor fit as it seems!
    (dBEF_nem17$Pr_per100g ==0) %>% sum() # 25 zeros could be responsible...
  
  #lets plot the abundance on the log scale:
  ggplot(dBEF_nem17, aes(x=sowndivLog, y=Pr_per100gLog.hurdle, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=estimate__, col=treatment))

rm(m17.Pr.hurdle21, m17.Pr.hurdle51)


#### abundances (21 and 17), no 60sp. ####
load("./statistics/brms/231127_abundance_OffSet.RData")

predictions.abun21 <- conditional_effects(m.abun.nBinomOffS.31)[[3]]
predictions.abun17 <- conditional_effects(m.abun.nBinomOffS.31vog)[[3]]

mcmc21 <- mcmc_plot(m.abun.nBinomOffS.31, prob=0.85, prob_outer=0.95, 
          variable = c("b_sowndivLogStd:treatment2", "b_sowndivLogStd:treatment3"))+
  scale_y_discrete(labels=c("b_sowndivLog:treatment2" = "SH:19y | PH:5y", "b_sowndivLog:treatment3"= "SH:19y | PH:19y"))+
  xlim(-0.4,0.4)

mcmc21
ggsave(filename = "./plots presi/MCMC_abun21.png",
       height = 4,
       width = 4)

BREAKS <- unique(dat$sowndivLogStd) %>% 
  sort()

#for 21 data
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60) %>% 
  mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
         .after = sowndivLog)

predictions = predictions.abun21
treatments = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")

  
p.abun21 <- ggplot(dat, aes(x=sowndivLogStd, y=total_nematodes) )+
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
  scale_y_continuous(name = "No. of nematodes")+
  theme_classic()+
  theme(legend.position ="bottom")
  

p.abun21
summary(m.abun.nBinomOffS.31)

ggsave(filename = "./plots presi/abun21_withCI02.png" ,plot = p.abun21,
       height = 4,
       width = 4)

#for 17 data
dat = dBEF_nem17 %>% 
  filter(sowndiv != 60) %>% 
  mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
         .after = sowndivLog)

predictions = predictions.abun17
treatments = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")


p.abun17 <- ggplot(dat, aes(x=sowndivLogStd, y=total_nematodes) )+
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
  scale_y_continuous(name = "No. of nematodes")+
  theme_classic()+
  theme(legend.position ="bottom")


p.abun17
summary(m.abun.nBinomOffS.31vog)

ggsave(filename = "./plots presi/abun17.png" ,plot = p.abun17,
       height = 4,
       width = 4)



####Hill numbers 2021, no 60 sp plots####
load("./statistics/brms/231127_Hill.RData")


#all: m.41.Hill_q1 
#Ba: m.41b.Hill_q1.Ba
#Fu: m.41b.Hill_q1.Fu

predictions.Hill_q1.21 <- conditional_effects(m.41.Hill_q1)[[3]]
predictions.Hill_q1.21.Ba <- conditional_effects(m.41.Hill_q1.Ba)[[3]]
predictions.Hill_q1.21.Fu <- conditional_effects(m.41.Hill_q1.Fu)[[3]]

#for all (21 data)
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60) %>% 
  mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
         .after = sowndivLog)

predictions = predictions.Hill_q1.21
treatments = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")
BREAKS = c(0,1,2,3,4)


p.Hill1.21 <- ggplot(dat, aes(x=sowndivLog, y=Hill_q1) )+
  geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                    fill=treatment), 
              alpha=0.2, show.legend=FALSE)+
  geom_line(data=predictions, aes(x= sowndivLog, y=estimate__, 
                                  linetype=treatment, col=treatment),
            linewidth= 1, show.legend = FALSE)+
  
  geom_jitter(width=0.2, shape=19, alpha=0.4, 
              aes(col=treatment))+
  scale_color_manual(labels=treatments, values = cols)+
  scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                     labels = c("1", "2", "4", "8", "16"))+
  scale_y_continuous(name = "exp(H')  for all Nematodes")+
  theme_classic()+
  theme(legend.position ="bottom")


p.Hill1.21
summary(m.41.Hill_q1)

ggsave(filename = "./plots presi/Hill1_21.png" ,plot = p.Hill1.21,
       height = 4,
       width = 4)


#for Bacterivores (21 data)
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60) %>% 
  mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
         .after = sowndivLog)

predictions = predictions.Hill_q1.21.Ba
treatments = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")
BREAKS = c(0,1,2,3,4)


p.Hill1.21.Ba <- ggplot(dat, aes(x=sowndivLog, y=Hill_q1.Ba) )+
  geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                    fill=treatment), 
              alpha=0.2, show.legend=FALSE)+
  geom_line(data=predictions, aes(x= sowndivLog, y=estimate__, 
                                  linetype=treatment, col=treatment),
            linewidth= 1, show.legend = FALSE)+
  
  geom_jitter(width=0.2, shape=19, alpha=0.4, 
              aes(col=treatment))+
  scale_color_manual(labels=treatments, values = cols)+
  scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                     labels = c("1", "2", "4", "8", "16"))+
  scale_y_continuous(name = "exp(H')  for Bacterivores")+
  theme_classic()+
  theme(legend.position ="bottom")


p.Hill1.21.Ba
summary(m.41.Hill_q1.Ba)

ggsave(filename = "./plots presi/Hill1Ba_21.png" ,plot = p.Hill1.21.Ba,
       height = 4,
       width = 4)


#for Fungivores (21 data)
dat = dBEF_nem21 %>% 
  filter(sowndiv != 60) %>% 
  mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
         .after = sowndivLog)

predictions = predictions.Hill_q1.21.Fu
treatments = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")
BREAKS = c(0,1,2,3,4)


p.Hill1.21.Fu <- ggplot(dat, aes(x=sowndivLog, y=Hill_q1.Fu) )+
  geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                    fill=treatment), 
              alpha=0.2, show.legend=FALSE)+
  geom_line(data=predictions, aes(x= sowndivLog, y=estimate__, 
                                  linetype=treatment, col=treatment),
            linewidth= 1, show.legend = FALSE)+
  
  geom_jitter(width=0.2, shape=19, alpha=0.4, 
              aes(col=treatment))+
  scale_color_manual(labels=treatments, values = cols)+
  scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                     labels = c("1", "2", "4", "8", "16"))+
  scale_y_continuous(name = "exp(H')  for Fungivores")+
  theme_classic()+
  theme(legend.position ="bottom")


p.Hill1.21.Fu
summary(m.41.Hill_q1.Fu)

ggsave(filename = "./plots presi/Hill1Fu_21.png" ,plot = p.Hill1.21.Fu,
       height = 4,
       width = 4)

####bacterivores 2021####
load(file="./statistics/brms/231127_Ba.RData")
loo(m.Ba.hurdle31a, m.Ba.hurdle31b)

predictions.Ba <- conditional_effects(m.Ba.hurdle31b)[[3]]
predictions.Ba1 <- conditional_effects(m.Ba1.hurdle31b)[[3]]
predictions.Ba2 <- conditional_effects(m.Ba2.hurdle31b)[[3]]
predictions.Ba4 <- conditional_effects(m.Ba4.hurdle31b)[[3]]

#our data as in the models
dat <- subset(dBEF_nem21, sowndiv!=60)
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog)

#
BREAKS <- unique(dat$sowndivLogStd) %>% 
  sort()

#for all bacterivores:
predictions = predictions.Ba
treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")

p.Ba <- ggplot(dat, aes(x=sowndivLogStd, y=Ba_per100g) )+
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
  scale_y_continuous(name = "Bacterivores per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")

p.Ba
ggsave(filename = "./plots presi/Ba_21.png" ,plot = p.Ba,
       height = 4,
       width = 4)


#for Ba-cp1:
predictions = predictions.Ba1
treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")

p.Ba1 <- ggplot(dat, aes(x=sowndivLogStd, y=Ba1) )+
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
  scale_y_continuous(name = "Ba (cp1) per 100g DW", limits = c(-1, 200))+
  theme_bw()+
  theme(legend.position ="bottom")

p.Ba1

ggsave(filename = "./plots presi/Ba1_21.png" ,plot = p.Ba1,
       height = 4,
       width = 4)


#for Ba-cp2:
predictions = predictions.Ba2
treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")

p.Ba2 <- ggplot(dat, aes(x=sowndivLogStd, y=Ba2) )+
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
  scale_y_continuous(name = "Ba (cp2) per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")

p.Ba2

ggsave(filename = "./plots presi/Ba2_21.png" ,plot = p.Ba2,
       height = 4,
       width = 4)


#for Ba-cp4:
predictions = predictions.Ba4
treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")

p.Ba4 <- ggplot(dat, aes(x=sowndivLogStd, y=Ba4) )+
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
  scale_y_continuous(name = "Ba (cp4) per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")

p.Ba4

ggsave(filename = "./plots presi/Ba4_21.png" ,plot = p.Ba4,
       height = 4,
       width = 4)


summary(m.Ba.hurdle31b)




#### Fungivores 2021 ####
load(file="./statistics/brms/231127_Fu.RData")
loo(m.Fu.hurdle31a, m.Fu.hurdle34b)

predictions.Fua  <- conditional_effects(m.Fu.hurdle31a)[[3]]
predictions.Fu2a <- conditional_effects(m.Fu2.hurdle31a)[[3]]
predictions.Fu3a <- conditional_effects(m.Fu3.hurdle32a)[[3]]
predictions.Fu4a <- conditional_effects(m.Fu4.hurdle32a)[[3]]

predictions.Fub  <- conditional_effects(m.Fu.hurdle34b)[[3]]
predictions.Fu2b <- conditional_effects(m.Fu2.hurdle31b)[[3]]
predictions.Fu3b <- conditional_effects(m.Fu3.hurdle31b)[[3]]
predictions.Fu4b <- conditional_effects(m.Fu4.hurdle33b)[[3]]

#our data as in the models
dat <- subset(dBEF_nem21, sowndiv!=60)
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog)
BREAKS <- unique(dat$sowndivLogStd) %>% 
  sort()

#for all Fungivores:
predictions = predictions.Fua
treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")

p.Fu <- ggplot(dat, aes(x=sowndivLogStd, y=Fu_per100g) )+
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
  scale_y_continuous(name = "Fungivores per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")

p.Fu
ggsave(filename = "./plots presi/Fu_21.png" ,plot = p.Fu,
       height = 4,
       width = 4)

#for Fu-2:
predictions = predictions.Fu2b
treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")

p.Fu2 <- ggplot(dat, aes(x=sowndivLogStd, y=Fu2) )+
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
  scale_y_continuous(name = "Fu (cp2) per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")

p.Fu2
ggsave(filename = "./plots presi/Fu2_21.png" ,plot = p.Fu2,
       height = 4,
       width = 4)


#for Fu-3:
predictions = predictions.Fu3b
treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")

p.Fu3 <- ggplot(dat, aes(x=sowndivLogStd, y=Fu3) )+
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
  scale_y_continuous(name = "Fu (cp3) per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")

p.Fu3

ggsave(filename = "./plots presi/Fu3_21.png" ,plot = p.Fu3,
       height = 4,
       width = 4)

#for Fu-4:
predictions = predictions.Fu4b
treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")

p.Fu4 <- ggplot(dat, aes(x=sowndivLogStd, y=Fu4) )+
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
  scale_y_continuous(name = "Fu (cp4) per 100g DW")+
  theme_bw()+
  theme(legend.position ="bottom")

p.Fu4

ggsave(filename = "./plots presi/Fu4_21.png" ,plot = p.Fu4,
       height = 4,
       width = 4)


####probability of direction, slopes ####
#abun 21
m.abun.nBinomOffS.31 %>% summary()

emt = emtrends(m.abun.nBinomOffS.31, "treatment", var="sowndivLogStd")
#get slopes:
  summary(emt, point.est=mean, level = .95) 
emt.pairs <- pairs(emt)
  summary(emt.pairs, point.est=mean, level = .95)
  bayestestR::p_direction(emt.pairs) #probability of difference between the treatments in the output to be negative
  
#abun 17
m.abun.nBinomOffS.31vog

emt = emtrends(m.abun.nBinomOffS.31vog, "treatment", var="sowndivLogStd")
summary(emt, point.est=mean, level = .95) 
emt.pairs <- pairs(emt)
summary(emt.pairs, point.est=mean, level = .95)
bayestestR::p_direction(emt.pairs) #probability of difference between the treatments in the output to be negative

#H' all
m.41.Hill_q1

emt = emtrends(m.41.Hill_q1, "treatment", var="sowndivLog")
summary(emt, point.est=mean, level = .95) 
emt.pairs <- pairs(emt)
summary(emt.pairs, point.est=mean, level = .95)
bayestestR::p_direction(emt.pairs) #probability of difference between the treatments in the output to be negative

#H' Fu
m.41.Hill_q1.Fu
emt = emtrends(m.41.Hill_q1.Fu, "treatment", var="sowndivLog")
summary(emt, point.est=mean, level = .95) 
emt.pairs <- pairs(emt)
summary(emt.pairs, point.est=mean, level = .95)
bayestestR::p_direction(emt.pairs) #probability of difference between the treatments in the output to be negative

#H' Ba
m.41.Hill_q1.Ba

emt = emtrends(m.41.Hill_q1.Ba, "treatment", var="sowndivLog")
summary(emt, point.est=mean, level = .95) 
emt.pairs <- pairs(emt)
summary(emt.pairs, point.est=mean, level = .95)
bayestestR::p_direction(emt.pairs) #probability of difference between the treatments in the output to be negative

#dens Fu
m.Fu.hurdle31a

emt = emtrends(m.Fu.hurdle31a, "treatment", var="sowndivLogStd")
summary(emt, point.est=mean, level = .95)
emt.pairs <- pairs(emt)
summary(emt.pairs, point.est=mean, level = .95)
bayestestR::p_direction(emt.pairs) #probability of difference between the treatments in the output to be negative

#dens Ba
m.Ba.hurdle31b

emt = emtrends(m.Ba.hurdle31b, "treatment", var="sowndivLogStd")
summary(emt, point.est=mean, level = .95)
emt.pairs <- pairs(emt)
summary(emt.pairs, point.est=mean, level = .95)
bayestestR::p_direction(emt.pairs) #probability of difference between the treatments in the output to be negative



#### to do####
pp_check(m.Pr.hurdle51, ndraws=100)


emt = emtrends(m.Pr.hurdle21b, "treatment", var="sowndivLog")
summary(emt, point.est=mean)
summary(emt, point.est=mean, level = .9) #get slopes

emt.pairs <- pairs(emt)
summary(emt.pairs, point.est=mean, level = .95) #get differences in slopes betweem treatments
bayestestR::p_direction(emt.pairs) #probability of difference between the treatments in the output to be negative
bayestestR::p_direction(m.Pr.hurdle51)

#try beta distribution to fit indices

#dont log transform only part of data
# -> use either hurdle_lognormal or
#hurdle_gamma

ggplot(data = dBEF_nem21, aes(x=sowndivLog, y=Hill_q0, col=treatment))+
  geom_jitter(width=0.3)

dat=subset(dBEF_nem21, sowndiv!=60)
dat$sowndivLogStd %>% density() %>% plot()

