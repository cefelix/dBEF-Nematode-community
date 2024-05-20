#exclude 60 sp.:
dat <- subset(dBEF_nem21, sowndiv != 60) 

#standardize:  
dat <- dat %>% mutate(sowndivLogStd = ( (sowndivLog - mean(sowndivLog)) / sd(sowndivLog) ),
                      .after = sowndivLog) %>%
  mutate(realdivLogStd = ( (realdivLog - mean(realdivLog)) / sd(realdivLog) ),
         .after = realdivLog)


#### a function that creates a table with probabilities of direction between treatments ####
pd_table <- function(model, predictor, lvl, epred = FALSE) {
  
  #get the response
  response <- deparse(model$formula)[1] %>% #this row extracts the model term and converts it to a string
    gsub(".*= (.+) ~.*", "\\1", .) #and this excludes everything which is not the response
  
  #get probabilities of direction between treatments
  emt.p <- emtrends(model, specs = c("treatment"), var = predictor, epred = epred) %>%
    pairs() %>% p_direction
  
    pd_t1_t2 <- emt.p$pd[1] %>% round(.,3) 
    pd_t1_t3 <- emt.p$pd[2] %>% round(.,3)
    pd_t2_t3 <- emt.p$pd[3] %>% round(.,3)
  
  #get the trend which is bigger     
  emt <- emtrends(model, specs = c("treatment"), var = predictor, epred = epred) %>% 
    as.data.frame() %>%
    mutate(mean.trend = eval(as.symbol(paste(predictor, "trend", sep="."))), .before = "lower.HPD") 

  #only show the relation between trends if pd > 0.95
  trend_12 <- ifelse(pd_t1_t2 > 0.95, 
                    yes = ifelse(emt$mean.trend[1] > emt$mean.trend[2], yes = "t1 > t2", no = "t2 > t1" ), #if pd > 0.95, show direction
                    no ="") #if pd < 0.95, dont show it
  trend_13 <- ifelse(pd_t1_t3 > 0.95, 
                    yes = ifelse(emt$mean.trend[1] > emt$mean.trend[3], yes = "t1 > t3", no = "t3 > t1" ), #if pd > 0.95, show direction
                    no ="") #if pd < 0.95, dont show it
  trend_23 <- ifelse(pd_t2_t3 > 0.95, 
                    yes = ifelse(emt$mean.trend[2] > emt$mean.trend[3], yes = "t2 > t3", no = "t3 > t2" ), #if pd > 0.95, show direction
                    no ="") #if pd < 0.95, dont show it
  
  ####ugly code to add significance asterisk####
  pd_t1_t2 <- ifelse(pd_t1_t2 > 0.995, 
         yes = paste(pd_t1_t2, "**"),
         no  = ifelse(pd_t1_t2 > .975,
                     yes = paste(pd_t1_t2, " *"),
                     no  = ifelse(pd_t1_t2 > .95,
                                  yes = paste(pd_t1_t2, " ."),
                                  no = paste(pd_t1_t2, "  ")))
  )
  
  pd_t1_t3 <- ifelse(pd_t1_t3 > 0.995, 
                    yes = paste(pd_t1_t3, "**"),
                    no  = ifelse(pd_t1_t3 > .975,
                                 yes = paste(pd_t1_t3, " *"),
                                 no  = ifelse(trend_13 > .95,
                                              yes = paste(pd_t1_t3, " ."),
                                              no = paste(pd_t1_t3, "  ")))
  )
  
  pd_t2_t3 <- ifelse(pd_t2_t3 > 0.995, 
                    yes = paste(pd_t2_t3, "**"),
                    no  = ifelse(pd_t2_t3 > .975,
                                 yes = paste(pd_t2_t3, " *"),
                                 no  = ifelse(pd_t2_t3 > .95,
                                              yes = paste(pd_t2_t3, " ."),
                                              no = paste(pd_t2_t3, "  ")))
  )
  
  #percent of the posterior laying in range of practical equivalence (ROPE)
  #https://easystats.github.io/bayestestR/articles/region_of_practical_equivalence.html
  #https://easystats.github.io/bayestestR/articles/guidelines.html
  
  emt.rope <- emtrends(model, specs = c("treatment"), var = predictor, epred = epred) %>%
    pairs() %>% rope()
  
  ROPE_t12 <- (emt.rope$ROPE_Percentage[1] * 100 ) %>% round(3)
  ROPE_t13 <- (emt.rope$ROPE_Percentage[2] * 100 ) %>% round(3)
  ROPE_t23 <- (emt.rope$ROPE_Percentage[3] * 100 ) %>% round(3)
  
  
  #return relevant output
  cbind(response, predictor, pd_t1_t2, pd_t1_t3, pd_t2_t3, trend_12, trend_13, trend_23, ROPE_t12, ROPE_t13, ROPE_t23) %>%
    as.data.frame() %>%
    return()
  
}

  

#test the function:
pd_table(m.Ba_sowndiv_p5, predictor = "sowndivLogStd")

emtrends(m.Pr4_Om4_sowndiv_p5, specs = c("treatment"), var = "sowndivLogStd") %>%
  pairs() %>% p_direction

m.Pr4_Om4_sowndiv_p5 %>% summary()

####apply it on all density ~ sowndiv models####
load("./statistics/brms/240109_TrophDens_sowndiv_mselect.RData")
load("./statistics/brms/240205_Ba_cp_sowndiv_mselect.RData")
load("./statistics/brms/240205_Fu_cp_sowndiv_mselect.RData") 
load("./statistics/brms/240205_Pl_cp_sowndiv_mselect.RData")
load("./statistics/brms/240205_Pr_Om_cp_sowndiv_mselect.RData")

m.sowndiv.list <- list(m.Ba_sowndiv_p5, m.Ba1_sowndiv_p5, m.Ba2_sowndiv_p5, m.Ba4_sowndiv_p5,
                       m.Fu_sowndiv_p5, m.Fu2_sowndiv_p5, m.Fu3_sowndiv_p5, m.Fu4_sowndiv_p5,
                       m.Pl_sowndiv_p5, m.Pl2_sowndiv_p5, m.Pl3_sowndiv_p5, m.Pl4_sowndiv_p5,
                       m.Om_sowndiv_p7,
                       m.Pr_sowndiv_p7,
                       m.Pr4_Om4_sowndiv_p5)

pd_summary_sowndiv <- lapply(m.sowndiv.list, pd_table, predictor = "sowndivLogStd") %>% 
  bind_rows()
pd_summary_sowndiv

####apply it on all density ~ realdiv models####
load("./statistics/brms/240109_TrophDens_realdiv_mselect.RData")
load("./statistics/brms/240205_Ba_cp_realdiv_mselect.RData")
load("./statistics/brms/240205_Fu_cp_realdiv_mselect.RData") 
load("./statistics/brms/240205_Pl_cp_realdiv_mselect.RData")
load("./statistics/brms/240205_Pr_Om_cp_realdiv_mselect.RData")

m.realdiv.list <- list(m.Ba_realdiv_p5, m.Ba1_realdiv_p5, m.Ba2_realdiv_p5, m.Ba4_realdiv_p5,
                       m.Fu_realdiv_p5, m.Fu2_realdiv_p5, m.Fu3_realdiv_p5, m.Fu4_realdiv_p5,
                       m.Pl_realdiv_p5, m.Pl2_realdiv_p5, m.Pl3_realdiv_p5, m.Pl4_realdiv_p5,
                       m.Om_realdiv_p7,
                       m.Pr_realdiv_p7,
                       m.Pr4_Om4_realdiv_p5)

pd_summary_realdiv <- lapply(m.realdiv.list, pd_table, predictor = "realdivLogStd") %>% 
  bind_rows()
pd_summary_realdiv

rope(m.Pr_realdiv_p7, effects = "fixed")
m.Pr_realdiv_p7 %>% emtrends(specs = c("treatment"), var = "realdivLogStd", epred = epred) %>%
  pairs() %>%
  rope()

#### apply it on all indices ####
load("./statistics/brms/240130_MI_sowndiv.RData")
load("./statistics/brms/240130_MI_sowndiv.RData")

#### save in a xlsx file ####
library(openxlsx)
list.pd_summaries <- list('pd densities ~ sowndiv' = pd_summary_sowndiv,
                        'pd densities ~ realdiv' = pd_summary_realdiv)

write.xlsx(list.pd_summaries, 
      file = "./statistics/240207_pd_summaries.xlsx")

#save some RAM
    rm(m.Ba_realdiv_p5, m.Ba1_realdiv_p5, m.Ba2_realdiv_p5, m.Ba4_realdiv_p5,
       m.Fu_realdiv_p5, m.Fu2_realdiv_p5, m.Fu3_realdiv_p5, m.Fu4_realdiv_p5,
       m.Pl_realdiv_p5, m.Pl2_realdiv_p5, m.Pl3_realdiv_p5, m.Pl4_realdiv_p5,
       m.Om_realdiv_p7,
       m.Pr_realdiv_p7,
       m.Pr4_Om4_realdiv_p5)
    
    rm(m.Ba_sowndiv_p5, m.Ba1_sowndiv_p5, m.Ba2_sowndiv_p5, m.Ba4_sowndiv_p5,
       m.Fu_sowndiv_p5, m.Fu2_sowndiv_p5, m.Fu3_sowndiv_p5, m.Fu4_sowndiv_p5,
       m.Pl_sowndiv_p5, m.Pl2_sowndiv_p5, m.Pl3_sowndiv_p5, m.Pl4_sowndiv_p5,
       m.Om_sowndiv_p7,
       m.Pr_sowndiv_p7,
       m.Pr4_Om4_sowndiv_p5)


    
    
    
    