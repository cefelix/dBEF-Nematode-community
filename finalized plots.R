#the final plots:
library(ggplot2)
library(bayesplot)
library(dplyr)

####regression plots: Trophic group densities####



#### posterior forrest plots with ggplot like in heiss' blog ####
#mcmc_areas(m.Ba_sowndiv_p5, prob=0.9 , prob_outer = 0.95, pars=c("b_sowndivLogStd:treatment"))
at_div <- 


p.Ba.ame <- m.Ba_sowndiv_p5 %>%
  emtrends(~ sowndivLogStd+treatment,
           var = "sowndivLogStd",
           at = list(#sowndivLogStd = c(0),
                     treatment = levels(dat$treatment)),
           epred = TRUE, 
           re_formula = NA) %>% 
  gather_emmeans_draws() %>%
  ggplot(., aes(x = .value, fill = factor(treatment))) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect of a\n 1 SD increase in sown plant diversity at Sown diversity = 4 species",
       y = "Density", fill = "History treatment") +
  xlim(-250, 250)+
  #theme(legend.position = "bottom")
  theme(legend.position = "NONE")



p.Fu.ame <- m.Fu_sowndiv_p5 %>%
  emtrends(~ sowndivLogStd+treatment,
           var = "sowndivLogStd",
           at = list(sowndivLogStd = c(-1.5, 1.5),
                     treatment = levels(dat$treatment)),
           epred = TRUE, 
           re_formula = NA) %>% 
  gather_emmeans_draws() %>%
  ggplot(., aes(x = .value, fill = factor(treatment))) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect of a\n 1 SD increase in sown plant diversity at Sown diversity = 4 species",
       y = "Density", fill = "History treatment") +
  xlim(-250, 250)+
  #theme(legend.position = "bottom")
  theme(legend.position = "NONE")

p.Pl.ame <- m.Pl_sowndiv_p5 %>%
  emtrends(~ sowndivLogStd+treatment,
           var = "sowndivLogStd",
           at = list(sowndivLogStd = c(-1.5, 1.5),
                     treatment = levels(dat$treatment)),
           epred = TRUE, 
           re_formula = NA) %>% 
  gather_emmeans_draws() %>%
  ggplot(., aes(x = .value, fill = factor(treatment))) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect of a\n 1 SD increase in sown plant diversity at Sown diversity = 4 species",
       y = "Density", fill = "History treatment") +
  xlim(-250, 250)+
  theme(legend.position = "bottom")

p.Pr.ame <- m.Pr_sowndiv_p7 %>%
  emtrends(~ sowndivLogStd+treatment,
           var = "sowndivLogStd",
           at = list(sowndivLogStd = c(-1.5, 1.5),
                     treatment = levels(dat$treatment)),
           #epred = TRUE, 
           re_formula = NA) %>% 
  gather_emmeans_draws() %>%
  ggplot(., aes(x = .value, fill = factor(treatment))) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect of a\n 1 SD increase in sown plant diversity at Sown diversity = 4 species",
       y = "Density", fill = "History treatment") +
  #xlim(-50, 50)+
  theme(legend.position = "bottom")

p.Om.ame <- m.Om_sowndiv_p7 %>%
  emtrends(~ sowndivLogStd+treatment,
           var = "sowndivLogStd",
           at = list(sowndivLogStd = c(0),
                     treatment = levels(dat$treatment)),
           #epred = TRUE, 
           re_formula = NA) %>% 
  gather_emmeans_draws() %>%
  ggplot(., aes(x = .value, fill = factor(treatment))) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect of a\n 1 SD increase in sown plant diversity at Sown diversity = 4 species",
       y = "Density", fill = "History treatment") +
  #xlim(-50, 50)+
  theme(legend.position = "bottom")



p.Ba.ame 
p.Fu.ame
p.Pl.ame 
p.Pr.ame 
p.Om.ame

library(gridExtra)
grid.arrange(p.Ba.ame, p.Fu.ame, p.Pl.ame) #epred=TRUE
grid.arrange(p.Pr.ame, p.Om.ame) #this is with epred=FALSE

#### the effect of epred = TRUE on the same scale as original: #####

log((lower.HPD+100)/100)

m.Ba_sowndiv_p5 %>%
  emtrends(~ sowndivLogStd+treatment,
           var = "sowndivLogStd",
           at = list(#sowndivLogStd = c(0),
             treatment = levels(dat$treatment)),
           epred = TRUE, 
           re_formula = NA) %>% 
  gather_emmeans_draws() %>%
  mutate(value = log((.value+100)/100) ) %>%
  ggplot(., aes(x = .value, fill = factor(treatment))) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect of a\n 1 SD increase in sown plant diversity at Sown diversity = 4 species",
       y = "Density", fill = "History treatment") +
  #xlim(-250, 250)+
  theme(legend.position = "bottom")


#### a scatter plot showing the abundance of nematodes for each week ####
ggplot(dat, aes(x=week, y = total_nematodes))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.3, alpha=0.7, aes(col=treatment))+
  theme_bw()




#### a regression of the density of Predators ####
predictions <- conditional_effects(m.Pr_sowndiv_p7, re_formula = NA,
                                   prob = 0.9)[[3]]
#when using conditional effects, to account for random effects set re_formula=NUll]
treatments2 = c("+SH+PH", "+SH-PH", "-SH-PH")
cols=c("brown2", "darkolivegreen", "dodgerblue3")

BREAKS = unique(dat$sowndivLogStd)
plot.regression.Pr <- ggplot(dat, aes(x=sowndivLogStd, y=Pr_per100g) )+
  geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                    fill=treatment), 
              alpha=0.2, show.legend=FALSE)+
  geom_jitter(width=0.05, shape=19, alpha=0.4, 
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

plot.regression.Pr  
conditional_effects(m.Pr_sowndiv_p7, ask=FALSE)


