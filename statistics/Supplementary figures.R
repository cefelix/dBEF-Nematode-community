# 
####S1 : soil moisture per week #### 

p.S1.a <- ggplot(dat, aes(x = week, y = SWC_gravimetric, col = week))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width= 0.2, alpha=0.3)+
  ylab("gravimetric soil water \ncontent (percent)")+
  ggtitle( "differences between weeks" ,subtitle = "a")+
  theme_bw()+
  theme(legend.position = "none")
  
p.S1.b <- ggplot(dat, aes(x = week, y = total_nematodes, col = week))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width= 0.2, alpha=0.3)+
  ylab("counted nematodes")+
  theme_bw()+
  ggtitle( "", subtitle = "b")+
  theme(legend.position = "none")

p.S1 <- ggarrange(p.S1.a, p.S1.b)

ggsave(p.S1, 
       filename = "./plots/sup_240229_SWCdens_byweek.png",
       dpi=300, width = 15, height = 8, units = "cm")
