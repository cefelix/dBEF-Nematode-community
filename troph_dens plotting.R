library(brms)
library(ggplot2)
library(tidybayes)
library(emmeans)
library(tidyr)

#loading the selected models:
  load("./statistics/brms/240109_TrophDens_sowndiv_mselect.RData")
  load("./statistics/brms/240109_TrophDens_realdiv_mselect.RData")
  load("./statistics/brms/240110_TrophDens_funcdiv_mselect.RData")

#data:
  #exclude 60 sp.:
  dat <- subset(dBEF_nem21, sowndiv != 60) 
  #standardize:  
  dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                        .after = sowndivLog) %>%
    mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
           .after = realdivLog)
  
#a guide to "correctly calculating posterior predictions and average marginal effects with multilievel Bayesian models":
  #https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/

#### troph dens ~ sowndiv: regression scatter plot ####
#using conditional effects:
predictions <- conditional_effects(m.Ba_sowndiv_p5, 
                                   re_formula = NULL, # this includes random effects when predicting
                                   prob = 0.9)[[3]] #90 % CI
  
str(predictions)
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
  dat2 <- dat %>% add_predicted_draws(m.Ba_sowndiv_p5, ndraws = 20)

    ggplot(data = dat2, aes(y=Ba_per100g, x=sowndivLogStd))+
      stat_dist_lineribbon(aes(y = .prediction, col = treatment), .width = c(.90),
                           alpha = 0.2)+
      geom_jitter(data = dat, aes(y=Ba_per100g, x=sowndivLogStd, col=treatment))
      

#conditional effects for each treatment                                    
    posterior_epred(model,
                    newdata = ,
                    re_formula = NULL)    
    
        
#using emmeans
    summary(m.Ba_sowndiv_p5)
    Ba.sowndiv_95CI <- emmeans(m.Ba_sowndiv_p5,
            ~sowndivLogStd + treatment,
            at = list(sowndivLogStd = seq(from = min(dat$sowndivLogStd), to = max(dat$sowndivLogStd),
                                          length.out = 100)),
            epred = TRUE) #gives us the expectation of predicted draws
    
    
#angelos' code for probability of direction:    
    emt = emtrends( m.Ba_sowndiv_p5, specs = c("treatment"), var="sowndivLogStd")
    summary(emt, point.est=mean, level = .9) 
    emt.pairs <- pairs(emt)
    summary(emt.pairs, point.est=mean, level = .9)
    bayestestR::p_direction(emt.pairs)
   
    
    
    
