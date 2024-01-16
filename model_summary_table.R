#creating a neat summary table for all models:
#packages:
library(brms)
library(tidyr)
library(emmeans)
library(dplyr)




#probability of direction:
Ba.prob_dir <- emtrends(m.Ba_sowndiv_p5, specs = c("treatment"), var="sowndivLogStd") %>%
  pairs() %>% bayestestR::p_direction()
Ba.prob_dir


#### summarise_models() for TrophDens ~sowndiv, ~realdiv, ~funcdiv #####
  #load models
  load("./statistics/brms/240109_TrophDens_sowndiv_mselect.RData")
  load("./statistics/brms/240109_TrophDens_realdiv_mselect.RData")
  load("./statistics/brms/240110_TrophDens_funcdiv_mselect.RData")




#a function which creates a table with point estimates and HPDI for regression coefficients for each treatment:
summarise_models <- function(brmsfit, predictor, level=.95) {
  response <- deparse(brmsfit$formula)[1] %>% #this row extracts the model term and converts it to a string
    gsub(".*= (.+) ~.*", "\\1", .) #and this excludes everything which is not the response

  
  ms <- summary(brmsfit) #extract summary from model
  ESS <- cbind(ms$spec_pars$Bulk_ESS, #bulk ESS
               ms$spec_pars$Tail_ESS) #tail ESS
  colnames(ESS) <- c("bulk ESS", "tail ESS")
  
  emt.s <- emtrends( brmsfit, specs = c("treatment"), var = predictor) %>% 
    summary(., point.est = mean, level = level) #get slope estimates mean and HPDI
  emt.s90 <- emtrends( brmsfit, specs = c("treatment"), var = predictor) %>% 
    summary(., point.est = mean, level = .9)
  emt.s95 <- emtrends( brmsfit, specs = c("treatment"), var = predictor) %>% 
    summary(., point.est = mean, level = .9)

  #adding letters, indicating different significance levels: 
  significance <- rep(NA, nrow(emt.s))
  for ( i in 1 : nrow(emt.s95)){
    if( (emt.s95$upper.HPD[i] > 0 & emt.s95$lower.HPD[i] > 0 ) |  #positive relationship (95% CI)
        (emt.s95$upper.HPD[i] < 0 & emt.s95$lower.HPD[i] < 0 )) { #negative (95% CI)
      significance[i] = "a"
    }
    else if( (emt.s90$upper.HPD[i] > 0 & emt.s90$lower.HPD[i] > 0 ) |  #positive (90% CI)
             (emt.s90$upper.HPD[i] < 0 & emt.s90$lower.HPD[i] < 0 )) { #negative (90% CI)
      significance[i] = "b"
    }
    
    
    else {
      significance[i] = ""
    }
  }

  
  divs <- 1 #Placeholder: add divergent transitions
  row <- cbind(response, predictor, emt.s[], ESS, divs, significance)
  return(row)
}

 
#~sowndivLogStd
  #a list with all final models:
  sown.list <- list(m.Ba_sowndiv_p5, m.Fu_sowndiv_p5, m.Pl_sowndiv_p5, m.Pr_sowndiv_p7, m.Om_sowndiv_p7)
  
  # a df of the model summaries:
  sowndiv.summary <- lapply(sown.list, summarise_models, predictor = "sowndivLogStd", level =.9) %>% 
    bind_rows()
  sowndiv.summary
  
  
#~realdivLogStd    
  #list of models:
  real.list <- list(m.Ba_realdiv_p5, m.Fu_realdiv_p5, m.Pl_realdiv_p5, m.Pr_realdiv_p7, m.Om_realdiv_p7)

  # a df of the model summaries:
  realdiv.summary <- lapply(real.list, summarise_models, predictor = "realdivLogStd", level =.9) %>% 
    bind_rows()
  realdiv.summary

#~funcdivLogStd    
  #list of models:
  func.list <- list(m.Ba_funcdiv_p5, m.Fu_funcdiv_p5, m.Pl_funcdiv_p5, m.Pr_funcdiv_p7, m.Om_funcdiv_p7)
  
  # a df of the model summaries:
  funcdiv.summary <- lapply(func.list, summarise_models, predictor = "funcdivStd", level =.9) %>% 
    bind_rows()
  funcdiv.summary
  
  
#### summarise_models() for Shannon H' ~sowndiv, ~realdiv, ~funcdiv #####
  #load models
  load("./statistics/brms/240109_TrophDens_sowndiv_mselect.RData")
  
  




#### old ####


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



