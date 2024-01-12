#creating a neat summary table for all models:
#packages:
library(brms)


#load models
load("./statistics/brms/240109_TrophDens_sowndiv_mselect.RData")
load("./statistics/brms/240109_TrophDens_realdiv_mselect.RData")
load("./statistics/brms/240110_TrophDens_funcdiv_mselect.RData")


#Ba ~ sowndiv
  Ba.ms <- summary(m.Ba_sowndiv_p5)
  #extract bulk/tail ESS
    ESS <- cbind(Ba.ms$spec_pars$Bulk_ESS, 
               Ba.ms$spec_pars$Tail_ESS)
    ESS
    
  #estimated marginal means for each treatment:
    Ba.emt.s <- emtrends( m.Ba_sowndiv_p5, specs = c("treatment"), var="sowndivLogStd") %>% 
      summary(., point.est = mean, level = .95)
    Ba.emt.s
    cbind(Ba.emt.s, ESS)

#Ba ~ realdiv
    Ba.ms <- summary(m.Ba_realdiv_p5)
    #extract bulk/tail ESS
    ESS <- cbind(Ba.ms$spec_pars$Bulk_ESS, 
                 Ba.ms$spec_pars$Tail_ESS)
    ESS
    
    #estimated marginal means for each treatment:
    Ba.emt.s <- emtrends( m.Ba_realdiv_p5, specs = c("treatment"), var="realdivLogStd") %>% 
      summary(., point.est = mean, level = .95)
    Ba.emt.s
    cbind(Ba.emt.s, ESS)
    
#Ba ~ funcdiv
    Ba.ms <- summary(m.Ba_funcdiv_p5)
    #extract bulk/tail ESS
    ESS <- cbind(Ba.ms$spec_pars$Bulk_ESS, 
                 Ba.ms$spec_pars$Tail_ESS)
    ESS
    
    #estimated marginal means for each treatment:
    Ba.emt.s <- emtrends( m.Ba_funcdiv_p5, specs = c("treatment"), var="funcdivStd") %>% 
      summary(., point.est = mean, level = .95)
    Ba.emt.s
    cbind(Ba.emt.s, ESS)
    
        
    

#probability of direction:
Ba.prob_dir <- emtrends(m.Ba_sowndiv_p5, specs = c("treatment"), var="sowndivLogStd") %>%
  pairs() %>% bayestestR::p_direction()
Ba.prob_dir








#### under construction: piping through it #####
#check out: https://www.youtube.com/watch?v=DZN5nXtEHcg

give_me_a_row <- function(brmsfit, predictor, level=.95) {
  ms <- summary(brmsfit)
  ESS <- cbind(ms$spec_pars$Bulk_ESS, 
               ms$spec_pars$Tail_ESS) 
  emt.s <- emtrends( brmsfit, specs = c("treatment"), var = predictor) %>% 
    summary(., point.est = mean, level = level)
  row <- cbind(emt.s[], ESS, divs)
  return(row)
}
  
give_me_a_row(m.Ba_sowndiv_p5, 
              predictor = "sowndivLogStd", 
              level = 0.9)

sown.list <- list(m.Ba_sowndiv_p5, m.Fu_sowndiv_p5, m.Pl_sowndiv_p5, m.Pr_sowndiv_p7, m.Om_sowndiv_p7)
lapply(sown.list, give_me_a_row, predictor = "sowndivLogStd")


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



