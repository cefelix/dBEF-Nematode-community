load()

#model plotting:
#new data to create regression curve:
#thats only block1
nd <- tibble(sowndiv = seq(from = 1, 
                           to = 60, 
                           length.out = 30) %>% rep(., times = 3),
             treatment = rep(1:3, each = 30),
             block = rep("B4", 90))
#fitted values
f <-
  fitted(m.Fu.12, newdata = nd) %>%
  as_tibble() %>%
  bind_cols(nd) %>%
  mutate(treatment = as.factor(treatment)) # treatment has to a factor to be plotted

#plotting fitted stuff
ggplot(dBEF_nem21, aes(x=log(sowndiv)))+
  geom_jitter(aes(y = Fu_per100g, color = treatment), width = 0.125, alpha=0.5)+
  geom_smooth(data = f,
              aes(y = exp(Estimate),  #exp estimate to not plot logged data
                  color = treatment),
              stat = "identity")+
  scale_x_continuous()+
  ylab("Fungivores per 100g DW")+
  theme_classic()


