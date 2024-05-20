#

#posterior predictions across treatments 
load("./statistics/brms/240109_TrophDens_sowndiv_mselect.RData")

grand_mean_treat <- m.Ba_sowndiv_p5 %>%
  epred_draws(newdata = expand.grid(treatment = c("1", "2", "3"),
                                    sowndivLogStd = unique(dat$sowndivLogStd)),
              re_formula = NA)

plot_grand_mean_treat <- ggplot(grand_mean_treat, 
                                  aes(x = .epred, y = "Grand mean", 
                                      fill = treatment))+
  stat_halfeye()

plot_grand_mean_treat


# emtrends based average marginal effects
  #avoiding hardcoding:
  div <- "sowndivLogStd"

  generate_ame <- function(model, response, predictor, epred = FALSE ){
    ame <- model %>%
      emtrends( var = predictor, specs = c("treatment"),
                at = list(predictor = 0),
                epred = epred, 
                re_formula = NA ) %>%
      gather_emmeans_draws() %>% #here the original code ends
      mutate(response = response)
    return(ame)
  }
  
  ame.Ba <- m.Ba_sowndiv_p5 %>% 
    generate_ame(response = "Ba", predictor = "sowndivLogStd")
    ame.Ba1 <- m.Ba1_sowndiv_p5 %>% 
      generate_ame(response = "Ba1", predictor = "sowndivLogStd")
    ame.Ba2 <- m.Ba2_sowndiv_p5 %>% 
      generate_ame(response = "Ba2", predictor = "sowndivLogStd")
    ame.Ba4 <- m.Ba4_sowndiv_p5 %>% 
      generate_ame(response = "Ba4", predictor = "sowndivLogStd")
    
  ame.Fu <- m.Fu_sowndiv_p5 %>% 
    generate_ame(response = "Fu", predictor = "sowndivLogStd")
  ame.Pl <- m.Pl_sowndiv_p5 %>% 
    generate_ame(response = "Pl", predictor = "sowndivLogStd")
  ame.Pr <- m.Pr_sowndiv_p7 %>% 
    generate_ame(response = "Pr", predictor = "sowndivLogStd")
  ame.Om <- m.Om_sowndiv_p7 %>% 
    generate_ame(response = "Om", predictor = "sowndivLogStd")

  ame.Ba.cp <- rbind(ame.Ba, ame.Ba1,ame.Ba2, ame.Ba4)
  
  
plot_grand_mean_treat_AME <- ggplot(ame.Pr,
                                    aes(x = .value, color= treatment)
                                    )+
  #xlim(-100, 100)+
  geom_density(alpha=0.5)#+
  #stat_halfeye(slab_alpha=0.55)+
  #facet_wrap(~treatment)
plot_grand_mean_treat_AME


