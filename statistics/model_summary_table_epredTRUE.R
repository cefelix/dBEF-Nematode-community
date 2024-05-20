#creating a neat summary table for all models:
#packages:
library(brms)
library(tidyr)
library(emmeans)
library(dplyr)
library(tidybayes)
library(ggplot2)
library(modelbased)
library(bayestestR)

#exclude 60 sp.:
dat <- subset(dBEF_nem21, sowndiv != 60) 

#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog) %>%
  mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
         .after = realdivLog)

#Explanation how to un-standardize coefficients here:
#https://discourse.mc-stan.org/t/how-to-rescale-coefficients-of-lognormal-regression-to-be-change-in-y-for-unit-increase-of-x/22472

# a function which creates a table with point estimates (mean of posterior) and HPDI (95% CI):
summarise_models <- function(brmsfit, predictor, unstandardized_predictor, level=.95) {
  response <- deparse(brmsfit$formula)[1] %>% #this row extracts the model term and converts it to a string
    gsub(".*= (.+) ~.*", "\\1", .) #and this excludes everything which is not the response
  
  
  #extracting Bulk/Tail ESS:
  ms <- summary(brmsfit) #extract summary from model
  fixed <- ms$fixed
  fixed$predictor <- rownames(fixed)
  ESS <- cbind( #a 2x3 matrix with bulk/tail ESS for sowndiv:treatment 
    rbind(fixed$Bulk_ESS[fixed$predictor == predictor] %>% round(2),                  
          fixed$Bulk_ESS[fixed$predictor == paste(predictor, "treatment2", sep = ":")] %>% round(2),
          fixed$Bulk_ESS[fixed$predictor == paste(predictor, "treatment3", sep = ":")] %>% round(2)
    ), 
    rbind(fixed$Tail_ESS[fixed$predictor == predictor] %>% round(2),
          fixed$Tail_ESS[fixed$predictor == paste(predictor, "treatment2", sep = ":")] %>% round(2),
          fixed$Tail_ESS[fixed$predictor == paste(predictor, "treatment3", sep = ":")] %>% round(2)
    )
  )
  colnames(ESS) <- c("bulk ESS", "tail ESS")
  
  family <- family(brmsfit)$family
  
  #the mean of the posterior predictive distribution, HPDIs, and PDs
  emt.s <- emtrends( brmsfit, specs = c("treatment"), var = predictor,
                     epred = TRUE, 
                     re_formula = NA) %>% 
    summary(., point.est = mean, level = level) %>% #get slope estimates mean and HPDI
    mutate(pd = estimate_slopes(brmsfit, trend = predictor, at = "treatment", ci= level)$pd, .after = "upper.HPD") %>%
    mutate( mean.trend = eval(as.symbol(paste(predictor, "trend", sep="."))), .before = "lower.HPD" ) %>%
    select( -as.symbol(paste(predictor, "trend", sep=".")) )
  
  #adding significance, based on probability of direction
  support <- rep(NA, nrow(emt.s)) 
  for (i in 1:nrow(emt.s) ) {
    if (emt.s$pd[i] < 0.95 ) {
      support[i] <- " "
    }
    else if (emt.s$pd[i] >= 0.95 & emt.s$pd[i] < 0.975)  {
      support[i] <- "."
    }
    else if (emt.s$pd[i] >= 0.975 & emt.s$pd[i] < 0.995)  {
      support[i] <- "*"
    }
    else if (emt.s$pd[i] >= 0.995 & emt.s$pd[i] < 0.9995)  {
      support[i] <- "**"
    }
    else if (emt.s$pd[i] >= 0.9995 )  {
      support[i] <- "***"
    }
  }
  
  divs <- 1 #Placeholder: add divergent transitions
  row <- cbind(response, predictor, family, emt.s[], support, ESS#, divs, pd_t1, pd_t2, pd_t3
  )
  return(row)
}

#### summarise_models() for abundance/density ~ sowndiv  #####
#load models
load("./statistics/brms/240221_TrophDens_sowndiv_mselect.RData")
load("./statistics/brms/240205_Ba_cp_sowndiv_mselect.RData")
load("./statistics/brms/240205_Fu_cp_sowndiv_mselect.RData") 
load("./statistics/brms/240205_Pl_cp_sowndiv_mselect.RData")
load("./statistics/brms/240221_Pr_Om_cp_sowndiv_mselect.RData")
load("./statistics/brms/240225_abun_dens_all_mselect.RData")

#~sowndivLogStd
#a list with all final models:
sown.list <- list(m.abun_all.sowndiv_p5, m.dens_sowndiv_p5,
                  m.Ba_sowndiv_p5, m.Ba1_sowndiv_p5, m.Ba2_sowndiv_p5, m.Ba4_sowndiv_p5,
                  m.Fu_sowndiv_p5, m.Fu2_sowndiv_p5, m.Fu3_sowndiv_p5, m.Fu4_sowndiv_p5,
                  m.Pl_sowndiv_p5, m.Pl2_sowndiv_p5, m.Pl3_sowndiv_p5, m.Pl4_sowndiv_p5,
                  m.Pr_sowndiv_p5, m.Om_sowndiv_p5, m.Pr_Om_sowndiv_p5)

# a df of the model summaries:
sowndiv.summary <- lapply(sown.list, summarise_models, predictor = "sowndivLogStd", level =.95) %>% 
  bind_rows()
sowndiv.summary




#un-standardize beta coefficients:
#explanation of the un-standardizing is here:
#https://discourse.mc-stan.org/t/how-to-rescale-coefficients-of-lognormal-regression-to-be-change-in-y-for-unit-increase-of-x/22472
sowndiv.summary <- sowndiv.summary %>%
  mutate(., mean.trend_destand = (sowndiv.summary$mean.trend / sd(dat$sowndivLog)) %>% round(., 3), .before = "pd") %>%
  mutate(lower.HPD.2 = ((sowndiv.summary$lower.HPD) / sd(dat$sowndivLog) ) %>% round(., 3), .before = "pd") %>%
  mutate(upper.HPD.2 = ((sowndiv.summary$upper.HPD) / sd(dat$sowndivLog) ) %>% round(., 3), .before = "pd") %>%
  mutate_at((vars(mean.trend, lower.HPD, upper.HPD, pd) ), round, digits=3)  

sowndiv.summary <- sowndiv.summary %>% 
  mutate(mean.trend = round(mean.trend, 3),
         lower.HPD = round(lower.HPD, 3),
         upper.HPD = round(upper.HPD, 3),
         pd = round(pd, 4))
sowndiv.summary 


#### summarise_models() for abundance/density ~ realdiv  #####
load("./statistics/brms/240221_TrophDens_realdiv_mselect.RData")
load("./statistics/brms/240205_Ba_cp_realdiv_mselect.RData")
load("./statistics/brms/240205_Fu_cp_realdiv_mselect.RData") 
load("./statistics/brms/240205_Pl_cp_realdiv_mselect.RData")
load("./statistics/brms/240221_Pr_Om_cp_realdiv_mselect.RData")
load("./statistics/brms/240225_abun_dens_all_mselect.RData")

#list of models:
real.list <- list(m.abun_all.realdiv_p5, m.dens_realdiv_p5,
                  m.Ba_realdiv_p5, m.Ba1_realdiv_p5, m.Ba2_realdiv_p5, m.Ba4_realdiv_p5,
                  m.Fu_realdiv_p5, m.Fu2_realdiv_p5, m.Fu3_realdiv_p5, m.Fu4_realdiv_p5,
                  m.Pl_realdiv_p5, m.Pl2_realdiv_p5, m.Pl3_realdiv_p5, m.Pl4_realdiv_p5,
                  m.Pr_realdiv_p5, m.Om_realdiv_p5, m.Pr_Om_realdiv_p5)

# a df of the model summaries:
realdiv.summary <- lapply(real.list, summarise_models, predictor = "realdivLogStd", level =.95) %>% 
  bind_rows()
realdiv.summary

#un-standardize beta coefficients:
realdiv.summary <- realdiv.summary %>%
  mutate(., realdivLog.trend = (realdiv.summary$mean.trend / sd(dat$realdivLog)) %>% round(., 3), .before = "bulk ESS") %>%
  mutate(lower.HPD.2 = ((realdiv.summary$lower.HPD) / sd(dat$realdivLog)) %>% round(., 3), .before = "bulk ESS") %>% 
  mutate(upper.HPD.2 = ((realdiv.summary$upper.HPD) / sd(dat$realdivLog)) %>% round(., 3), .before = "bulk ESS") %>%
  mutate_at((vars(mean.trend, lower.HPD, upper.HPD, pd) ), round, digits=3)  

realdiv.summary <- realdiv.summary %>%
  mutate(mean.trend = round(mean.trend, 3),
         lower.HPD = round(lower.HPD, 3),
         upper.HPD = round(upper.HPD, 3),
         pd = round(pd, 4))
realdiv.summary

#### summarise_models() for Effective number of families ####  
#for ~ sowndiv:
load("./statistics/brms/240214_Hill_sowndiv_mselect.RData")
Hill.sown.list <- list(m.all.Shannon.gaus_p5, m.Ba.Shannon.gamma_p5, m.Fu.Shannon.gaus_p5, m.Pl.Shannon.gamma_p5, m.PrOm.Shannon.gamma_p5)

Hill.sown.summary <- lapply(Hill.sown.list, summarise_models, predictor = "sowndivLogStd", level =.95) %>%
  bind_rows()

#un-standardize beta coefficients:
Hill.sown.summary <- Hill.sown.summary %>%
  mutate(., sowndivLog.trend = Hill.sown.summary$mean.trend / sd(dat$sowndivLog) , .before = "bulk ESS") %>%
  mutate(lower.HPD.2 = (Hill.sown.summary$lower.HPD) / sd(dat$sowndivLog) , .before = "bulk ESS") %>% 
  mutate(upper.HPD.2 = (Hill.sown.summary$upper.HPD) / sd(dat$sowndivLog) %>% round(., 3), .before = "bulk ESS") %>%
  #and a transform for families with a log-link:
  mutate(sowndivLog.trend = ifelse(family == "gamma", exp(sowndivLog.trend), sowndivLog.trend) %>% round(., 3)) %>%
  mutate(lower.HPD.2 = ifelse(family == "gamma", exp(lower.HPD.2), lower.HPD.2) %>% round(., 3)) %>%
  mutate(upper.HPD.2 = ifelse(family == "gamma", exp(upper.HPD.2), upper.HPD.2) %>% round(., 3)) %>%
  mutate_at((vars(mean.trend, lower.HPD, upper.HPD, pd) ), round, digits=3)  

Hill.sown.summary <- Hill.sown.summary %>%
  mutate(mean.trend = round(mean.trend, 3),
         lower.HPD = round(lower.HPD, 3),
         upper.HPD = round(upper.HPD, 3),
         pd = round(pd, 4))

Hill.sown.summary 

#for ~ realdiv
load("./statistics/brms/240214_Hill_realdiv_mselect.RData")
Hill.real.list <- list(m.all.Shannon.gaus_p5, m.Ba.Shannon.gamma_p5, m.Fu.Shannon.gaus_p5, m.Pl.Shannon.gamma_p5, m.PrOm.Shannon.gamma_p5)

Hill.real.summary <- lapply(Hill.real.list, summarise_models, predictor = "realdivLogStd", level =.95) %>%
  bind_rows()

#un-standardize beta coefficients:
Hill.real.summary <- Hill.real.summary %>%
  mutate(., realdivLog.trend = Hill.real.summary$mean.trend / sd(dat$realdivLog) , .before = "bulk ESS") %>%
  mutate(lower.HPD.2 = (Hill.real.summary$lower.HPD) / sd(dat$realdivLog) , .before = "bulk ESS") %>% 
  mutate(upper.HPD.2 = (Hill.real.summary$upper.HPD) / sd(dat$realdivLog) , .before = "bulk ESS") %>%
  #and a transform for families with a log-link:
  mutate(realdivLog.trend = ifelse(family == "gamma", exp(realdivLog.trend), realdivLog.trend) %>% round(., 3)) %>%
  mutate(lower.HPD.2 = ifelse(family == "gamma", exp(lower.HPD.2), lower.HPD.2) %>% round(., 3)) %>%
  mutate(upper.HPD.2 = ifelse(family == "gamma", exp(upper.HPD.2), upper.HPD.2) %>% round(., 3)) %>%
  mutate_at((vars(mean.trend, lower.HPD, upper.HPD, pd) ), round, digits=3)  

Hill.real.summary <- Hill.real.summary %>%
  mutate(mean.trend = round(mean.trend, 3),
         lower.HPD = round(lower.HPD, 3),
         upper.HPD = round(upper.HPD, 3),
         pd = round(pd, 4))

Hill.real.summary 


#### summarise_models() for indices ####
load("./statistics/brms/240215_SI_mselect.RData")
load("./statistics/brms/240215_EI_mselect.RData")

indices.sown.list <- list(m.EI.sowndiv_gaus_p5, m.SI.sowndiv_gaus_p5)
indices.sown.summary <- lapply(indices.sown.list, summarise_models, predictor = "sowndivLogStd", level =.95) %>%
  bind_rows()

indices.sown.summary <- indices.sown.summary %>%
  mutate(., sowndivLog.trend = indices.sown.summary$mean.trend / sd(dat$sowndivLog) , .before = "bulk ESS") %>%
  mutate(lower.HPD.2 = (indices.sown.summary$lower.HPD) / sd(dat$sowndivLog) , .before = "bulk ESS") %>% 
  mutate(upper.HPD.2 = (indices.sown.summary$upper.HPD) / sd(dat$sowndivLog) %>% round(., 3), .before = "bulk ESS")

#for real:
indices.real.list <- list(m.EI.realdiv_gaus_p5, m.SI.realdiv_gaus_p5)
indices.real.summary <- lapply(indices.real.list, summarise_models, predictor = "realdivLogStd", level =.95) %>%
  bind_rows()

indices.real.summary <- indices.real.summary %>%
  mutate(., sowndivLog.trend = indices.real.summary$mean.trend / sd(dat$realdivLog) , .before = "bulk ESS") %>%
  mutate(lower.HPD.2 = (indices.real.summary$lower.HPD) / sd(dat$realdivLog) , .before = "bulk ESS") %>% 
  mutate(upper.HPD.2 = (indices.real.summary$upper.HPD) / sd(dat$realdivLog) %>% round(., 3), .before = "bulk ESS")

indices.summary <- rbind(indices.real.summary, indices.sown.summary)

#### write each model summary into a .xlsx sheet  ####
library(openxlsx)
list.msummaries <- list('Densities ~ sowndiv' = sowndiv.summary, 
                        'Densities ~ realdiv' = realdiv.summary ,
                        'Hill numbers ~ sowndiv' = Hill.sown.summary,
                        'Hill numbers ~ realdiv' = Hill.real.summary,
                        'EI SI ~ sowndiv realdiv' = indices.summary)

write.xlsx(list.msummaries, 
           file = "./statistics/240226_Model_summaries_epredT.xlsx")

####OLD: the difference between emtrends(), emtrends(epred=TRUE) and posterior_epred() ####
emtrends(m.Ba4_sowndiv_p5, var="sowndivLogStd") %>%
  summary(., point.est = mean, level = .95) %>%
  mutate(sowndivLogStd.trend = exp(sowndivLogStd.trend / sd(dat$sowndivLog))-1,
         lower.HPD = exp(lower.HPD / sd(dat$sowndivLog))-1, 
         upper.HPD = exp(upper.HPD / sd(dat$sowndivLog))-1)

ce <- conditional_effects(m.Pr4_Om4_sowndiv_p5)$`sowndivLogStd:treatment` 
ce
cet1 <- subset(ce, treatment == 1) 
cet1$estimate__ %>% mean()
cet1$lower__ %>% mean()
cet1$upper__ %>% mean()

posterior_summary(m.Ba4_sowndiv_p5, variable = "sowndivLogStd")
library(posterior)
draws <- as_draws_array(m.Ba4_sowndiv_p5)
summarise_draws(draws)


emtrends(m.Ba4_sowndiv_p5, var = "sowndivLogStd", specs = c("treatment"),
         epred = TRUE) %>%
  summary(., point.est = mean, level = .95)

posterior_epred(m.Ba4_sowndiv_p5) #%>% summary()

spread_draws(m.Ba4_sowndiv_p5)



a <- pd_table(m.Ba_sowndiv_p5, predictor = "sowndivLogStd")

####NOT USED: custom forest ####
custom_forest_plot <- function(mod, variable) {
  post <- as_draws_array(mod) 
  ci.95 <- ci(post, method = 'ETI') %>% filter(Parameter %in% variable) %>% rowwise() %>% mutate(mu = (CI_low+CI_high)/2)
  ci.85 <- ci(post, method = 'ETI',ci = 0.85) %>% filter(Parameter %in%  variable) %>% rowwise() %>% mutate(mu = (CI_low+CI_high)/2)
  ci.50 <- ci(post, method = 'ETI', ci = 0.5) %>% filter(Parameter %in%  variable) %>% rowwise() %>% mutate(mu = (CI_low+CI_high)/2)  
  plot <- ggplot(ci.95, aes(x = mu, y = Parameter))+
    geom_vline(xintercept = 0, color = 'grey', linewidth = 0.8)+
    geom_linerange(aes(xmin = CI_low, xmax = CI_high), linewidth = 1, color = '#464645')+
    geom_linerange(data = ci.85, aes(xmin = CI_low, xmax = CI_high), linewidth = 2.6, color = '#a40b0b')+
    geom_point(shape = 21, size = 5, color = 'black', fill = '#dc3a3a')+
    theme_bw(base_size = 20, base_line_size = 20/44)+
    theme(axis.text = element_text(color = 'black'))+
    labs(x = NULL, y = NULL)
  return(plot)
}

custom_forest_plot(m.Ba4_sowndiv_p5, variable = "treatment")
custom_forest_plot(m.Pr4_Om4_sowndiv_p5, variable = "treatment")

post <- as_draws_array(m.Ba4_sowndiv_p5) 
ci.95 <- ci(post, method = 'ETI')

summary(m.Pr4_Om4_sowndiv_p5)
