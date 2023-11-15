library(brms)
library(tidyr)
library(ggplot2)
library(tidybayes)
library(modelr)

#run variable transformations.R before this script!

#how to extract data from conditonal_effects():
  #https://discourse.mc-stan.org/t/change-aesthetics-conditional-effects/13500/5

####Predators 2021####
load("./statistics/brms/231108_Pr_hurdle.RData")
#best: m.Pr.hurdle51, log(pr)~log(sowndiv)

#get model predictions for regression plot   
  predictions <- conditional_effects(m.Pr.hurdle51)[[3]]
  predictions$estimate__

  ggplot(dBEF_nem21, aes(x=sowndivLog, y=Pr_per100g, col=treatment) )+
    geom_jitter(width=0.3)+
    geom_smooth(data=predictions, aes(x= sowndivLog, y=exp(estimate__), col=treatment))
  
  rm(m.Pr.hurdle11, m.Pr.hurdle21, m.Pr.hurdle22, m.Pr.hurdle31, 
     m.Pr.hurdle42, m.Pr.hurdle51, m.Pr.hurdle61)

  
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


#### abundance 2021 ####
load("./statistics/brms/231108_abundance_OffSet.RData")
# good one: m.abun.nBinomOffS.22
summary(m.abun.nBinomOffS.22)
mcmc21 <- mcmc_plot(m.abun.nBinomOffS.22, prob=0.85, prob_outer=0.95, 
          variable = c("b_sowndivLog:treatment2", "b_sowndivLog:treatment3"))+
  scale_y_discrete(labels=c("b_sowndivLog:treatment2" = "SH:19y | PH:5y", "b_sowndivLog:treatment3"= "SH:19y | PH:19y"))+
  xlim(-0.4,0.4)

mcmc21
ggsave(filename = "./plots presi/MCMC_abun21.png",
       height = 4,
       width = 4)

  predictions <- conditional_effects(m.abun.nBinomOffS.22)[[3]]
  cols = c("darkolivegreen4", "coral3", "cornflowerblue")
  #colsCI = c("chartreuse3", "brown1", "cyan3")
  treatments=c("-SH -PH", "+SH -PH", "+SH +PH")
  treatments2=c(" SH:5y | PH:5y", "SH:19y | PH:5y", "SH:19y | PH:19y")
  
p <- ggplot(dBEF_nem21, aes(x=sowndivLog, y=total_nematodes) )+
  geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                    fill=treatment), 
              alpha=0.2, show.legend=FALSE)+
  geom_line(data=predictions, aes(x= sowndivLog, y=estimate__, 
                                  linetype=treatment, col=treatment),
            linewidth= 1, show.legend = FALSE)+
  
  geom_jitter(width=0.3, shape=19, alpha=0.4, 
              aes(col=treatment))+
  #ggtitle("brm(total_nematodes ~ log(sowndiv) * treatment\n+ offset(log(soilDW)) + (1 | block/plot), \nfam=negbinomial)")+
  scale_color_manual(labels=treatments, values = cols)+
  scale_x_continuous(name = "sown plant diversity", breaks = c(0, 1, 2, 3, 4, 6),
                     labels = c("1", "2", "4", "8", "16", "60"))+
  scale_y_continuous(name = "No. of nematodes")+
  theme_classic()+
  theme(legend.position ="bottom")
  

p

ggsave(filename = "./plots presi/abun21_withCI02.png" ,plot = p,
       height = 4,
       width = 4)


#lets look at the 3 outliers on the left:
outl.abun21 <- subset(dBEF_nem21, treatment == 1 & sowndiv==1 & total_nematodes > 400)
outl.abun21$Pl_per100g
outl.abun21$Pr_per100g
outl.abun21$Ba_per100g

####abundance 2017####
load("./statistics/17_remodelling/brms/231108_abundance17_OffSet.RData")
# good one: m.abun.nBinomOffS.22
summary(m.abun17.nBinomOffS.22)
mcmc_plot(m.abun17.nBinomOffS.22)

predictions <- conditional_effects(m.abun17.nBinomOffS.22)[[3]]
cols = c("darkolivegreen4", "coral3", "cornflowerblue")
#colsCI = c("chartreuse3", "brown1", "cyan3")
treatments=c("-SH -PH", "+SH -PH", "+SH +PH")
treatments2=c(" SH:1y | PH:1y", "SH:15y | PH:1y", "SH:15y | PH:15y")

p <- ggplot(dBEF_nem17, aes(x=sowndivLog, y=total_nematodes) )+
  geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                    fill=treatment), 
              alpha=0.2, show.legend=FALSE)+
  geom_jitter(width=0.3, shape=19, alpha=0.4, 
              aes(col=treatment))+
  geom_line(data=predictions, aes(x= sowndivLog, y=estimate__, 
                                  linetype=treatment, col=treatment),
            linewidth= 1, show.legend = FALSE)+
  #ggtitle("brm(total_nematodes ~ log(sowndiv) * treatment\n+ offset(log(soilDW)) + (1 | block/plot), \nfam=negbinomial)")+
  
  scale_color_manual(labels=treatments, values = cols)+
  scale_x_continuous(name = "sown plant diversity", breaks = c(0, 1, 2, 3, 4, 6),
                     labels = c("1", "2", "4", "8", "16", "60"))+
  scale_y_continuous(name = "No. of nematodes")+
  theme_classic()+
  theme(legend.position ="bottom")


p

ggsave(filename = "./plots presi/abun17_withCI.png" ,plot = p,
       width=6)

####hill numbers 2021####
load(file = "./statistics/brms/231108_hillNumbers.RData" )
#m.Hill_q1.12

summary(m.Hill_q1.12)

predictions <- conditional_effects(m.Hill_q1.12)[[3]]
cols = c("darkolivegreen4", "coral3", "cornflowerblue")
#colsCI = c("chartreuse3", "brown1", "cyan3")
treatments=c("-SH -PH", "+SH -PH", "+SH +PH")
treatments2=c(" SH:5y | PH:5y", "SH:19y | PH:5y", "SH:19y | PH:19y")

p <- ggplot(dBEF_nem21, aes(x=sowndivLog, y=Hill_q1) )+
  geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                    fill=treatment), 
              alpha=0.2, show.legend=FALSE)+
  geom_jitter(width=0.3, shape=19, alpha=0.4, 
              aes(col=treatment))+
  geom_line(data=predictions, aes(x= sowndivLog, y=estimate__, 
                                    linetype=treatment, col=treatment),
            linewidth= 1, show.legend = FALSE)+
  
  
  #ggtitle("brm(Hill_q1 ~ log(sowndiv)*treatment + (1|block/plot), fam=gaussian)")+
  
  scale_color_manual(labels=treatments, values = cols)+
  scale_x_continuous(name = "sown plant diversity", breaks = c(0, 1, 2, 3, 4, 6),
                     labels = c("1", "2", "4", "8", "16", "60"))+
  scale_y_continuous(name = "exp(Shannon's H')")+
  theme_classic()+
  theme(legend.position ="bottom")

p

ggsave(filename = "./plots presi/Hillq1_21_withCI.png" ,plot = p,
       width = 6)

mcmc_plot(m.Hill_q1.12, prob=0.85, prob_outer=0.95) 

MCMC_hillq1_21 <- mcmc_plot(m.Hill_q1.12, prob=0.85, prob_outer=0.95,
                            variable = c("b_sowndivLog:treatment2", "b_sowndivLog:treatment3"))+
  scale_y_discrete(labels=c("b_sowndivLog:treatment2" = "SH:19y | PH:5y", "b_sowndivLog:treatment3"= "SH:19y | PH:19y"))+
  xlim(-1.2, 1.2)

MCMC_hillq1_21

ggsave(filename = "./plots presi/Hillq1_21_MCMC.png" , plot = MCMC_hillq1_21,
       width = 6)


####Hill numbers 2021 without 60 sp plots####
dBEF_nem21.no60 <- subset(dBEF_nem21, sowndiv != 60)
load(file = "./statistics/brms/231108_hillNumbers.RData" )
#m.no60.Hill_q1.12

summary(m.no60.Hill_q1.12)

predictions <- conditional_effects(m.no60.Hill_q1.12)[[3]]
cols = c("darkolivegreen4", "coral3", "cornflowerblue")
#colsCI = c("chartreuse3", "brown1", "cyan3")
treatments=c("-SH -PH", "+SH -PH", "+SH +PH")
treatments2=c(" SH:5y | PH:5y", "SH:19y | PH:5y", "SH:19y | PH:19y")

p <- ggplot(dBEF_nem21.no60, aes(x=sowndivLog, y=Hill_q1) )+
  geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                    fill=treatment), 
              alpha=0.2, show.legend=FALSE)+
  geom_jitter(width=0.3, shape=19, alpha=0.4, 
              aes(col=treatment))+
  geom_line(data=predictions, aes(x= sowndivLog, y=estimate__, 
                                  linetype=treatment, col=treatment),
            linewidth= 1, show.legend = FALSE)+
  
  
  ggtitle("brm(Hill_q1 ~ log(sowndiv)*treatment + (1|block/plot),\nfam=gaussian)")+
  
  scale_color_manual(labels=treatments2, values = cols)+
  scale_x_continuous(name = "sown plant diversity", breaks = c(0, 1, 2, 3, 4, 6),
                     labels = c("1", "2", "4", "8", "16", "60"))+
  scale_y_continuous(name = "Hill q1")+
  theme_bw()

p
ggsave(filename = "./plots presi/Hillq1_withou60sp1_21_withCI.png" ,plot = p)

#### to do####
pp_check(m.Pr.hurdle51, ndraws=100)


emt = emtrends(m.Pr.hurdle51, "treatment", var="sowndivLog")
summary(emt, point.est=mean)
summary(emt, point.est=mean, level = .9) #get slopes

emt.pairs <- pairs(emt)
summary(emt.pairs, point.est=mean, level = .95) #get differences in slopes betweem treatments
bayestestR::p_direction(emt.pairs) #probability of difference between the treatments in the output to be negative
bayestestR::p_direction(m.Pr.hurdle51)

#dont log transform only part of data
# -> use either hurdle_lognormal or
#hurdle_gamma



