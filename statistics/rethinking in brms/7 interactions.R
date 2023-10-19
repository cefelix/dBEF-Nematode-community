#get data
library(rethinking)
data(rugged)
  d <- rugged
  detach(package:rethinking, unload=T)
library(brms)
  rm(rugged)

  
#### 7.1 - building an interaction #####  
#wrangle data
  library(tidyverse)
  
  # make the log version of criterion
  d <- 
    d %>%
    mutate(log_gdp = log(rgdppc_2000))
  
  # extract countries with GDP data
  dd <-
    d %>%
    filter(complete.cases(rgdppc_2000))
  
  # split the data into countries in Africa and not in Africa
  d.A1 <-
    dd %>%
    filter(cont_africa == 1) #africa
  
  d.A0 <-
    dd %>%
    filter(cont_africa == 0) #not africa
  
#predict log_gdp:
  # a univariable model for africa
  b7.1 <-
    brm(data = d.A1, family = gaussian,
        log_gdp ~ 1 + rugged,
        prior = c(prior(normal(8, 100), class = Intercept),
                  prior(normal(0, 1), class = b),
                  prior(uniform(0, 10), class = sigma)),
        iter = 2000, warmup = 1000, chains = 4, cores = 4,
        seed = 7)
  # a univariable model for non-africa
  b7.2 <-
    update(b7.1, 
           newdata = d.A0)
  # wrangle data for figure 7.2:
  nd <- 
    tibble(rugged = seq(from = 0, to = 6.3, length.out = 30))
  
  f_b7.1 <-
    fitted(b7.1, newdata = nd) %>%
    as_tibble() %>%
    bind_cols(nd)
  
  f_b7.2 <-
    fitted(b7.2, newdata = nd) %>%
    as_tibble() %>%
    bind_cols(nd)
  
  # here we'll put both in a single data object, with `f_b7.1` stacked atop `f_b7.2`
  f <-
    full_join(f_b7.1, f_b7.2) %>%
    mutate(cont_africa = rep(c("Africa", "not Africa"), each = 30))
  
  # plot figure 7.2
  library(ggthemes)
  dd %>%
    mutate(cont_africa = ifelse(cont_africa == 1, "Africa", "not Africa")) %>%
    
    ggplot(aes(x = rugged)) +
    geom_smooth(data = f,
                aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                    fill = cont_africa, color = cont_africa),
                stat = "identity", 
                alpha = 1/4, size = 1/2) +
    geom_point(aes(y = log_gdp, color = cont_africa),
               size = 2/3) +
    scale_x_continuous("Terrain Ruggedness Index", expand = c(0, 0)) +
    ylab("log GDP from year 2000") +
    theme_classic()+
    theme(text = element_text(family = "Times"),
          legend.position = "none") +
    facet_wrap(~cont_africa)
  
  #adding a dummy variable
    #without dummy:
    b7.3 <-
      update(b7.1,
             newdata = dd)
    #with dummy:
    b7.4 <-
      update(b7.3,
             newdata = dd,
             formula = log_gdp ~ 1 + rugged + cont_africa) 
    
  #compare models using loo and waic:
    b7.3 <- add_criterion(b7.3, c("loo", "waic"))
    b7.4 <- add_criterion(b7.4, c("loo", "waic"))
    
    loo_compare(b7.3, b7.4,
                criterion = "waic")
    loo_compare(b7.3, b7.4,
                criterion = "loo") # loo and waic are in agreement
    
    model_weights(b7.3, b7.4,
                  weights = "waic") %>% 
      round(digits = 3) #all weight goes to b7.4 (the one with the dummy variable) 
    
  #wrangling the model for plotting it:
    nd <- 
      tibble(rugged      = seq(from = 0, to = 6.3, length.out = 30) %>% 
               rep(., times = 2),
             cont_africa = rep(0:1, each = 30))
    
    f <-
      fitted(b7.4, newdata = nd) %>%
      as_tibble() %>%
      bind_cols(nd) %>%
      mutate(cont_africa = ifelse(cont_africa == 1, "Africa", "not Africa"))
    
  #create figure 7.3:
    dd %>%
      mutate(cont_africa = ifelse(cont_africa == 1, "Africa", "not Africa")) %>%
      
      ggplot(aes(x = rugged)) +
      geom_smooth(data = f,
                  aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                      fill = cont_africa, color = cont_africa),
                  stat = "identity", 
                  alpha = 1/4, size = 1/2) +
      geom_point(aes(y = log_gdp, color = cont_africa),
                 size = 2/3) +
      scale_x_continuous("Terrain Ruggedness Index", expand = c(0, 0)) +
      ylab("log GDP from year 2000") +
      theme_classic() + 
      theme(text = element_text(family = "Times"),
            legend.position  = c(.69, .94),
            legend.title     = element_blank(),
            legend.direction = "horizontal")
    
  #7.1.2 - adding a linear interaction does work
    # fit the model: 
    b7.5 <-
      update(b7.4,
             formula = log_gdp ~ 1 + rugged*cont_africa)
    
    #compare it to the previous models:
    b7.5 <- add_criterion(b7.5, c("loo", "waic"))
    
    l <- loo_compare(b7.3, b7.4, b7.5, b7.5b,
                     criterion = "loo")
    
      print(l, simplify = F)
      # a traditional way of showing loo's:
      cbind(loo_diff = l[, 1] * -2,
            se       = l[, 2] *  2)
      
      # weigh the models based on loo's:
      model_weights(b7.3, b7.4, b7.5, b7.5b,
                    weights = "loo") %>% 
        round(digits = 3)
    
    # a traditional notation of model b7.5:
    b7.5b <-    
      update(b7.5,
             formula = log_gdp ~ 1 + rugged + cont_africa + rugged:cont_africa)
      #confirm that it is actually the same model using loo:
        b7.5b <- add_criterion(b7.5b, c("loo", "waic"))
        l <- loo_compare(b7.5, b7.5b,
                         criterion = "loo")
        print(l, simplify = F) #as its a result from sampling, its not exactly the same
        
  #7.1.3 - plotting the interaction: 
    f <-
      fitted(b7.5, newdata = nd) %>%  # we can use the same `nd` data from last time
      as_tibble() %>%
      bind_cols(nd) %>%
      mutate(cont_africa = ifelse(cont_africa == 1, "Africa", "not Africa"))
    
    #interpreting an interaction just by numbers is tricky: plot the predictions!
    dd %>%
      mutate(cont_africa = ifelse(cont_africa == 1, "Africa", "not Africa")) %>%
      
      ggplot(aes(x = rugged, color = cont_africa)) +
      geom_smooth(data = f,
                  aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                      fill = cont_africa),
                  stat = "identity", 
                  alpha = 1/4, size = 1/2) +
      geom_point(aes(y = log_gdp),
                 size = 2/3) +
      scale_x_continuous("Terrain Ruggedness Index", expand = c(0, 0)) +
      ylab("log GDP from year 2000") +
      theme_classic() + 
      theme(text = element_text(family = "Times"),
            legend.position = "none") +
      facet_wrap(~cont_africa)
    
  #7.1.4 - parameters change meaning
    posterior_summary(b7.5)
    #gamma (look up in rethinking!) does not appear, we have to compute it
    
    #inside africa:  
    fixef(b7.5)[2,1] + fixef(b7.5)[4,1] * 1
      fixef(b7.5)      #beta coefficients for each predictor
      fixef(b7.5)[2,1] #slope for rugged terrain (outside africa)
      fixef(b7.5)[4,1] #change in slope when inside africa
   
    #outside africa:
    fixef(b7.5)[2,1] + fixef(b7.5)[4,1] * 0 
    
  #7.1.4.2 - incorporating uncertainty
    #gamma depends on beta coefficients, which have a distribution
    #so gamma has a distribution:
    post <- as_draws_matrix(b7.5) %>%
      as.data.frame() # a slightly different way, as posterior_samples is deprecated
  
    post %>%
      transmute(gamma_Africa    = b_rugged + `b_rugged:cont_africa`,
                gamma_notAfrica = b_rugged) %>%
      gather(key, value) %>%
      group_by(key) %>%
      summarise(mean = mean(value))
        
  #figure 7.5: (doesnt work)
    post %>%
      transmute(gamma_Africa    = b_rugged + `b_rugged:cont_africa`,
                gamma_notAfrica = b_rugged) %>%
      gather(key, value) %>%
      
      ggplot(aes(x = value, group = key, color = key, fill = key)) +
      geom_density(alpha = 1/4) +
      scale_x_continuous(expression(gamma), expand = c(0, 0)) +
      scale_y_continuous(NULL, breaks = NULL) +
      ggtitle("Terraine Ruggedness slopes",
              subtitle = "Blue = African nations, Green = others") +
      theme_classic + 
      theme(text = element_text(family = "Times"),
            legend.position = "none")
  
  #what proportion is below zero?    
  post %>%
    mutate(gamma_Africa    = b_rugged + `b_rugged:cont_africa`,
           gamma_notAfrica = b_rugged) %>% 
      mutate(diff            = gamma_Africa -gamma_notAfrica) %>%
      summarise(Proportion_of_the_difference_below_0 = sum(diff < 0) / length(diff))  
   
#7.2 symmetry of the linear interaction
  
  
  
  