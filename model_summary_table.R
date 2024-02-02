#creating a neat summary table for all models:
#packages:
library(brms)
library(tidyr)
library(emmeans)
library(dplyr)
library(tidybayes)

#exclude 60 sp.:
dat <- subset(dBEF_nem21, sowndiv != 60) 

#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog) %>%
  mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
         .after = realdivLog)


####just trash to try out code ####
ggplot(dat, aes(x= sowndivLogStd, y=Coverage, col= treatment))+
  geom_boxplot()

ggplot(dat, aes(x= sowndivLogStd, y=Coverage, col= week))+
  geom_boxplot()


#probability of direction:
Ba.prob_dir <- emtrends(m.Ba_sowndiv_p5, specs = c("treatment"), var="sowndivLogStd") %>%
  pairs() %>% bayestestR::p_direction()
Ba.prob_dir %>% str()
Ba.prob_dir


#this at least works:
m.Pr_sowndiv_p7 %>%
  emtrends(~ sowndivLogStd + treatment, var = "sowndivLogStd", epred = TRUE,
           at = list(sowndivLogStd = seq(-1.5, 1.5, 0.3),
                     treatment = c("1", "2", "3")),
           level = 0.9
           ) #%>% summary(point.est=mean, level=0.9)

m.Pr_sowndiv_p7 %>%
  emtrends(~ sowndivLogStd + treatment, var = "sowndivLogStd", epred = TRUE,
           at = list(treatment = c("1", "2", "3")),
           level = 0.9
  ) #%>% summary(point.est=mean, level=0.9)


emtrends(m.Pr_sowndiv_p7, specs = c("treatment"), var="sowndivLogStd") %>%
  summary(point.est = mean, level =.9)

conditional_effects(m.Om_sowndiv_p7, method = "posterior_epred") #%>% summary()
posterior_epred(m.Pr_sowndiv_p7) #the brms version
epred_draws(m.Pr_sowndiv_p7, newdata = dat) #the tidybayes version


summary(m.Pr_sowndiv_p7)
posterior_epred(m.Pr_sowndiv_p7)


#### summarise_models() for TrophDens ~sowndiv, ~realdiv, ~funcdiv #####
  #load models
  load("./statistics/brms/240109_TrophDens_sowndiv_mselect.RData")
  load("./statistics/brms/240109_TrophDens_realdiv_mselect.RData")
  load("./statistics/brms/240131_TrophDens_funcdiv_mselect.RData")


#a function which creates a table with point estimates and HPDI for regression coefficients for each treatment:
  summarise_models <- function(brmsfit, predictor, unstandardized_predictor, level=.95) {
  response <- deparse(brmsfit$formula)[1] %>% #this row extracts the model term and converts it to a string
    gsub(".*= (.+) ~.*", "\\1", .) #and this excludes everything which is not the response

  
  #extracting Bulk/Tail ESS:
  ms <- summary(brmsfit) #extract summary from model
  fixed <- ms$fixed
  fixed$predictor <- rownames(fixed)
  ESS <- cbind( #a 2x3 matrix with bulk/tail ESS for sowndiv:treatment 
      rbind(fixed$Bulk_ESS[fixed$predictor == predictor],                  
            fixed$Bulk_ESS[fixed$predictor == paste(predictor, "treatment2", sep = ":")],
            fixed$Bulk_ESS[fixed$predictor == paste(predictor, "treatment3", sep = ":")]
      ), 
      rbind(fixed$Tail_ESS[fixed$predictor == predictor],
            fixed$Tail_ESS[fixed$predictor == paste(predictor, "treatment2", sep = ":")],
            fixed$Tail_ESS[fixed$predictor == paste(predictor, "treatment3", sep = ":")]
      )
    )
  colnames(ESS) <- c("bulk ESS", "tail ESS")
  
  family <- family(brmsfit)$family
  
  #the mean of the expected value of the posterior predictive distribution, HPDIs, and PDs
  emt.s <- emtrends( brmsfit, specs = c("treatment"), var = predictor,
                     epred = TRUE) %>% 
    summary(., point.est = mean, level = level) %>% #get slope estimates mean and HPDI
    mutate(mean.trend = log((eval(as.symbol(paste(predictor, "trend", sep=".")))+100)/100), .before = "lower.HPD",
              #https://stackoverflow.com/questions/9057006/getting-strings-recognized-as-variable-names-in-r (2nd answer)
              #eval(as.symbol(paste(predictor, "trend", sep="."))) necessary to extract any predictor 
              #back-transforms to the same scale as emtrends(epred=FALS) yields
           lower.HPD = log((lower.HPD+100)/100),
           upper.HPD = log((upper.HPD+100)/100)) %>%
    mutate(pd = modelbased::estimate_slopes(brmsfit, trend = predictor, at = "treatment", ci= level)$pd, .after = "upper.HPD") %>%
    select( -as.symbol(paste(predictor, "trend", sep=".")) )
  
  #### start exclude: significance based on HPDIs ####
  emt.s90 <- emtrends( brmsfit, specs = c("treatment"), var = predictor,
                       epred = TRUE) %>% 
    summary(., point.est = mean, level = .9) %>%
  mutate(mean.trend = log((eval(as.symbol(paste(predictor, "trend", sep=".")))+100)/100), .before = "lower.HPD",
         lower.HPD = log((lower.HPD+100)/100),
         upper.HPD = log((upper.HPD+100)/100)) %>%
  select( -as.symbol(paste(predictor, "trend", sep=".")) )
  
  emt.s95 <- emtrends( brmsfit, specs = c("treatment"), var = predictor,
                       epred = TRUE) %>% 
    summary(., point.est = mean, level = .9) %>%
    mutate(mean.trend = log((eval(as.symbol(paste(predictor, "trend", sep=".")))+100)/100), .before = "lower.HPD",
           lower.HPD = log((lower.HPD+100)/100),
           upper.HPD = log((upper.HPD+100)/100)) %>%
    select( -as.symbol(paste(predictor, "trend", sep=".")) )
  

  #adding letters, indicating different significance levels based on HPDs (in relation having no effect) : 
  significance.HPD <- rep(NA, nrow(emt.s))
  for ( i in 1 : nrow(emt.s95)){
    if( (emt.s95$upper.HPD[i] > 0 & emt.s95$lower.HPD[i] > 0 ) |  #positive relationship (95% CI)
        (emt.s95$upper.HPD[i] < 0 & emt.s95$lower.HPD[i] < 0 )) { #negative (95% CI)
      significance.HPD[i] = "a"
    }
    else if( (emt.s90$upper.HPD[i] > 0 & emt.s90$lower.HPD[i] > 0 ) |  #positive (90% CI)
             (emt.s90$upper.HPD[i] < 0 & emt.s90$lower.HPD[i] < 0 )) { #negative (90% CI)
      significance.HPD[i] = "b"
    }
    
    
    else {
      significance.HPD[i] = ""
    }
  }
  #### end exclude: significance based on HPDIs ####
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
  
  ####work here####
  modelbased <- modelbased::estimate_slopes(brmsfit, trend = predictor, at = "treatment", ci= level)
  modelbased::estimate_contrasts(brmsfit, trend = predictor, contrast = "treatment", ci= level)
  #here it is old
  
  
#### exclude: p_direction between treatments ####   
  #adding probability of direction to each treatment:
    pd <- emtrends(brmsfit, specs = c("treatment"), var = predictor) %>%
      pairs() %>%
      bayestestR::p_direction() #a 2x3 pd matrix
  #re-arranging the 2x3 matrix into three columns:
    pd_t1 <- c("", round(pd$pd[pd$Parameter == "treatment1 - treatment2"], 3),  round(pd$pd[pd$Parameter == "treatment1 - treatment3"], 3)) # a column showing pd of each treatment against treatment 1
    pd_t2 <- c(round(pd$pd[pd$Parameter == "treatment1 - treatment2"], 3), "",  round(pd$pd[pd$Parameter == "treatment2 - treatment3"], 3)) # a column showing pd of each treatment against treatment 2
    pd_t3 <- c(round(pd$pd[pd$Parameter == "treatment1 - treatment3"], 3),  round(pd$pd[pd$Parameter == "treatment2 - treatment3"], 3), "") # a column showing pd of each treatment against treatment 3
  
#adding significance codes to pd: . (>.9), * (>.95), ** (>.99)
    for (i in 1:length(pd_t1)) {
      if (pd_t1[i] > 0.995 ) {
        pd_t1[i] <- paste(pd_t1[i], "**")
      }
      else if (pd_t1[i] > 0.975 ) {
        pd_t1[i] <- paste(pd_t1[i], "*")
      }
      else if (pd_t1[i] > 0.95 ) {
        pd_t1[i] <- paste(pd_t1[i], ".")
      }
      else {
        pd_t1[i] <- paste(pd_t1[i], " ")
      }
    }
    
    for (i in 1:length(pd_t2)) {
      if (pd_t2[i] > 0.995 ) {
        pd_t2[i] <- paste(pd_t2[i], "**")
      }
      else if (pd_t2[i] > 0.975 ) {
        pd_t2[i] <- paste(pd_t2[i], "*")
      }
      else if (pd_t2[i] > 0.95 ) {
        pd_t2[i] <- paste(pd_t2[i], ".")
      }
      else {
        pd_t2[i] <- paste(pd_t2[i], " ")
      }
    }
    
    for (i in 1:length(pd_t3)) {
      if (pd_t3[i] > 0.995 ) {
        pd_t3[i] <- paste(pd_t3[i], "**")
      }
      else if (pd_t3[i] > 0.975 ) {
        pd_t3[i] <- paste(pd_t3[i], "*")
      }
      else if (pd_t3[i] > 0.95 ) {
        pd_t3[i] <- paste(pd_t3[i], ".")
      }
      else {
        pd_t3[i] <- paste(pd_t3[i], " ")
      }
    }
    
       
  divs <- 1 #Placeholder: add divergent transitions
  row <- cbind(response, predictor, family, emt.s[], support, ESS#, divs, pd_t1, pd_t2, pd_t3
               )
  return(row)
}

 
#~sowndivLogStd
  #a list with all final models:
  sown.list <- list(m.Ba_sowndiv_p5, m.Fu_sowndiv_p5, m.Pl_sowndiv_p5, m.Pr_sowndiv_p7, m.Om_sowndiv_p7)
  
  # a df of the model summaries:
  sowndiv.summary <- lapply(sown.list, summarise_models, predictor = "sowndivLogStd", level =.9) %>% 
    bind_rows()
  sowndiv.summary
  
  
#un-standardize beta coefficients:
  sowndiv.summary <- sowndiv.summary %>%
    mutate(., mean.trend_destand = exp(sowndiv.summary$mean.trend / sd(dat$sowndivLog)), .before = "pd") %>%
    mutate(lower.HPD.2 = exp((sowndiv.summary$lower.HPD) / sd(dat$sowndivLog)), .before = "pd") %>%
    mutate(upper.HPD.2 = exp((sowndiv.summary$upper.HPD) / sd(dat$sowndivLog)), .before = "pd") 
  
  sowndiv.summary
 
  
#~realdivLogStd    
  #list of models:
  real.list <- list(m.Ba_realdiv_p5, m.Fu_realdiv_p5, m.Pl_realdiv_p5, m.Pr_realdiv_p7, m.Om_realdiv_p7)

  # a df of the model summaries:
  realdiv.summary <- lapply(real.list, summarise_models, predictor = "realdivLogStd", level =.9) %>% 
    bind_rows()
  realdiv.summary
  
  #un-standardize beta coefficients:
  realdiv.summary <- realdiv.summary %>%
    mutate(., realdivLog.trend = exp(realdiv.summary$LogStd.trend / sd(dat$realdivLog)), .before = "bulk ESS") %>%
    mutate(lower.HPD.2 = exp((realdiv.summary$lower.HPD) / sd(dat$realdivLog)), .before = "bulk ESS") %>%
    mutate(upper.HPD.2 = exp((realdiv.summary$upper.HPD) / sd(dat$realdivLog)), .before = "bulk ESS") 
    
  realdiv.summary

#~funcdivLogStd    
  #list of models:
  func.list <- list(m.Ba_funcdiv_p5, m.Fu_funcdiv_p5, m.Pl_funcdiv_p5, m.Pr_funcdiv_p7, m.Om_funcdiv_p7)
  
  # a df of the model summaries:
  funcdiv.summary <- lapply(func.list, summarise_models, predictor = "funcdivStd", level =.9) %>% 
    bind_rows()
  funcdiv.summary
  
  #un-standardize beta coefficients:
  funcdiv.summary <- funcdiv.summary %>%
    mutate(., funcdiv.trend = funcdiv.summary$LogStd.trend / sd(dat$funcdiv), .before = "bulk ESS") %>%
    mutate(lower.HPD.2 = (funcdiv.summary$lower.HPD) / sd(dat$funcdiv), .before = "bulk ESS") %>%
    mutate(upper.HPD.2 = (funcdiv.summary$upper.HPD) / sd(dat$funcdiv), .before = "bulk ESS") 
  
  funcdiv.summary
  
  
  
#### summarise_models() for Shannon H' ~sowndiv, ~realdiv, ~funcdiv #####
  
#sowndiv:
  load("./statistics/brms/240116_Hill_sowndiv_mselect.RData")

  Shannon.sown.list <-  list(m.all.Shannon.gaus_p5, m.Ba.Shannon.gaus_p5, m.Fu.Shannon.gamma_p5,
                       m.Pl.Shannon.gaus_p5, m.Pr.Shannon.gamma_p5)  
  
  Shannon.sowndiv.summary <- lapply(Shannon.sown.list, summarise_models, predictor = "sowndivLogStd", level =.9) %>% 
    bind_rows()
  Shannon.sowndiv.summary
  
  ####work here
  bayestestR::rope(Shannon.sowndiv.summary$lower.HPD)
  
  #un-standardize beta coefficients:
  Shannon.sowndiv.summary <- Shannon.sowndiv.summary %>%
    mutate(., sowndivLog.trend = Shannon.sowndiv.summary$sowndivLogStd.trend / sd(dat$sowndivLog), .before = "bulk ESS") %>%
    mutate(lower.HPD.2 = (Shannon.sowndiv.summary$lower.HPD) / sd(dat$sowndivLog), .before = "bulk ESS") %>%
    mutate(upper.HPD.2 = (Shannon.sowndiv.summary$upper.HPD) / sd(dat$sowndivLog), .before = "bulk ESS")
    
  Shannon.sowndiv.summary

#realdiv:
  load("./statistics/brms/240116_Hill_realdiv_mselect.RData")
  
  Shannon.real.list <-  list(m.all.Shannon.gaus_p5, m.Ba.Shannon.gaus_p5, m.Fu.Shannon.gamma_p5,
                             m.Pl.Shannon.gaus_p5, m.Pr.Shannon.gamma_p5)  
  
  Shannon.realdiv.summary <- lapply(Shannon.real.list, summarise_models, predictor = "realdivLogStd", level =.9) %>% 
    bind_rows()
  Shannon.realdiv.summary

  #un-standardize beta coefficients:
  Shannon.realdiv.summary <- Shannon.realdiv.summary %>%
    mutate(., realdivLog.trend = Shannon.realdiv.summary$realdivLogStd.trend / sd(dat$realdivLog), .before = "bulk ESS") %>%
    mutate(lower.HPD.2 = (Shannon.realdiv.summary$lower.HPD) / sd(dat$realdivLog), .before = "bulk ESS") %>%
    mutate(upper.HPD.2 = (Shannon.realdiv.summary$upper.HPD) / sd(dat$realdivLog), .before = "bulk ESS")
  
  Shannon.realdiv.summary  

  
#funcdiv: YET TO COME
  
  
#### write each model summary into a .xlsx sheet  ####
  library(openxlsx)
  list.msummaries <- list('TrophDens ~ sowndiv' = sowndiv.summary, 
                          'TrophDens ~ realdiv' = realdiv.summary, 
                          'TrophDens ~ funcdiv' = funcdiv.summary,
                          'Shannon H ~ sowndiv' = Shannon.sowndiv.summary, 
                          'Shannon H ~ realdiv' = Shannon.realdiv.summary)
  
  write.xlsx(list.msummaries, 
        file = "./statistics/240116_Model_summaries.xlsx")

    
####plot the most significant relationships: Trophic group densities####
  library(ggplot2)
  
  
  ggplot(data = dat, aes(x=sowndivLog, y=Pr_per100g))+
    geom_jitter(width = 0.1, shape =19, alpha = 0.4,
                aes(col=treatment))+
    facet_wrap(~treatment)
  
  predictions <- conditional_effects(m.Pr_sowndiv_p7, 
                                    #re_formula = NULL, # this includes random effects when predicting
                                     prob = 0.9,
                                     method = "posterior_epred")[[3]] 
  predictions <- conditional_effects(m.Pr_sowndiv_p7, 
                                     re_formula = NULL, # this includes random effects when predicting
                                     prob = 0.9,
                                     method = "posterior_predict")[[3]] 
  predictions
  
  emt.s <- emtrends(m.Pr_sowndiv_p7, specs = c("treatment"), var = predictor) %>% 
    summary(., point.est = mean, level = 0.9) #get slope estimates mean and HPDI
  emt.s
  
  
  
  
  emtrends(m.Pr_sowndiv_p7, specs = c("treatment"), var = "sowndivLogStd")
  
  ggplot(dat, aes(x=sowndivLogStd, y=Pr_per100g) )+
    geom_ribbon(data=predictions, aes(ymin= lower__, ymax=upper__, 
                                      fill=treatment), 
                alpha=0.2, show.legend=FALSE)+
    geom_jitter(width=0.1, shape=19, alpha=0.4, 
                aes(col=treatment))+
    geom_line(data=predictions, aes(x= sowndivLogStd, y=estimate__, 
                                    linetype=treatment, col=treatment),
              linewidth= 1, show.legend = FALSE)+
    facet_wrap(~treatment)+
    scale_color_manual(labels=treatments2, values = cols)+
    scale_x_continuous(name = "sown plant diversity", breaks = BREAKS,
                       labels = c("1", "2", "4", "8", "16"))+
    scale_y_continuous(name = "Ba per 100g DW")+
    theme_bw()+
    theme(legend.position ="bottom")
  
  
####the probability of zero predators: ####
  plot_emmeans1 <- m.Pr_sowndiv_p7 |> 
    emmeans(~ sowndivLogStd*treatment, var = "sowndivLogStd", dpar = "hu",
            #regrid = "response", 
            #tran = "log", type = "response",
            at = list(sowndivLogStd = seq(min(dat$sowndivLogStd), max(dat$sowndivLogStd), 0.01),
                      treatment = unique(dat$treatment))) |> 
    gather_emmeans_draws() |> 
    mutate(.value = exp(.value)) |> 
    
    
    ggplot(aes(x = sowndivLogStd, y = .value)) +
    stat_lineribbon(size = 1) +
    labs(x = "sowndiv (standardized)", y = "Predicted density of Pr per 100g SDW",
         subtitle = "Hurdle part of the model (dpar = \"hu\")",
         fill = "Credible interval") +
    theme(legend.position = "bottom")
  
  plot_emmeans1
  
  ref_grid(m.Pr_sowndiv_p7)
  
  



#### old ####
  
  emt = emtrends( m.Pr_sowndiv_p, specs = c("treatment", "week"), var="sowndivLogStd")
  summary(emt, point.est=mean, level = .9) 
  emt.pairs <- pairs(emt)
  summary(emt.pairs, point.est=mean, level = .9)
  bayestestR::p_direction(emt.pairs)  
  
  


predictor = "sowndivLogStd"
brmsfit = m.Ba_sowndiv_p5

summary(m.Pr_sowndiv_p7) %>% str()  

library(rstan)
foo <- get_sampler_params(m.Pr_sowndiv_p7$fit)
foo[[1]]

rstan::get_num_divergent(m.Ba_sowndiv_p5)

emt.s <- emtrends( brmsfit, specs = c("treatment"), var = predictor) %>% 
  summary(., point.est = mean, level = 0.95)
emt.s[1,]


m.list <- list(m.Ba_sowndiv_p5, m.Fu_sowndiv_p5)

#test:
    m.list[[2]] %>% summary() #works
    m.Ba_sowndiv_p5 %>% summary() #works
     
    i = 1
    model = m.list[[i]]
    model #works
    
    
for (i in m.list) {
  model = m.list[[i]]
  model %>% summary %>%
    print()
}


for (i in m.list) {
  i[[spec_pars]]
}

Ba.ms$spec_pars$Bulk_ESS


ms <- m.Ba_sowndiv_p5 %>%
  summary()



