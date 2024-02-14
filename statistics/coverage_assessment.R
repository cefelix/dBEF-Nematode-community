#analyzing coverage between weeks:
dat <- subset(dBEF_nem, year == 2021)


one.way <- aov(Coverage ~ week, data = dat)
summary(one.way)
plot(one.way)

two.way <- aov(Coverage ~ week*treatment, data = dat)
summary(two.way)

library(ggplot2)
ggplot(dat, aes(x=Coverage, y=total_nematodes, col = treatment))+
  #geom_jitter(width = 0.2)+
  geom_point()+
  facet_wrap(~week)

ggplot(dat, aes(x=Coverage, y=Hill_q1, col = treatment))+
  #geom_jitter(width = 0.2)+
  geom_point()+
  facet_wrap(~week)

ggplot(dat, aes(x=Rhabditidae.dauer.larvae, y=Rhabditis, col = treatment))+
  #geom_jitter(width = 0.2)+
  geom_point()+
  facet_wrap(~week)


####delete ####
nemaplex_fam <- read.csv("./wrangling/nemaplex_fam.csv")
nemaplex_fam$family %>% unique()
nemaplex_fam$genus %>% unique()

sum(dat$Dolichodoridae != 0)
sum(dat$Hoplolaimidae != 0)
sum(dat$Tylenchidae != 0)
sum(dat$Cephalobidae != 0)
sum(dat$Rhabditidae.dauer.larvae != 0)
sum(dat$Qudsianematidae != 0)

#Ba
p.Ba1 <- ggplot(dat, aes(x=sowndivLog, y=Ba1, col=treatment))+
  geom_jitter(width = .2)

p.Ba2 <- ggplot(dat, aes(x=sowndivLog, y=Ba2, col=treatment))+
  geom_jitter(width = .2)

p.Ba3 <- ggplot(dat, aes(x=sowndivLog, y=Ba3, col=treatment))+
  geom_jitter(width = .2)

p.Ba4 <- ggplot(dat, aes(x=sowndivLog, y=Ba4, col=treatment))+
  geom_jitter(width = .2)

sum(dat$Ba5 != 0)

grid.arrange(p.Ba1, p.Ba2, p.Ba3, p.Ba4)

#fu
p.Fu2 <- ggplot(dat, aes(x=sowndivLog, y=Fu2, col=treatment))+
  geom_jitter(width = .2)

p.Fu3 <- ggplot(dat, aes(x=sowndivLog, y=Fu3, col=treatment))+
  geom_jitter(width = .2)

p.Fu4 <- ggplot(dat, aes(x=sowndivLog, y=Fu4, col=treatment))+
  geom_jitter(width = .2)

grid.arrange(p.Fu2, p.Fu3, p.Fu4)

#Pl
p.Pl2 <- ggplot(dat, aes(x=sowndivLog, y=Pl2, col=treatment))+
  geom_jitter(width = .2)

p.Pl3 <- ggplot(dat, aes(x=sowndivLog, y=Pl3, col=treatment))+
  geom_jitter(width = .2)

p.Pl4 <- ggplot(dat, aes(x=sowndivLog, y=Pl4, col=treatment))+
  geom_jitter(width = .2)

p.Pl45 <- ggplot(dat, aes(x=sowndivLog, y=(Pl4 + Pl5), col=treatment))+
  geom_jitter(width = .2)

grid.arrange(p.Pl2, p.Pl3, p.Pl4, p.Pl45)

#Pr
dat$Pr3 %>% max()

p.Pr4 <- ggplot(dat, aes(x=sowndivLog, y=Pr4, col=treatment))+
  geom_jitter(width = .2)

p.Pr5 <- ggplot(dat, aes(x=sowndivLog, y=Pr5, col=treatment))+
  geom_jitter(width = .2)

grid.arrange(p.Pr4, p.Pr5)

#Om
p.Om4 <- ggplot(dat, aes(x=sowndivLog, y=Om4, col=treatment))+
  geom_jitter(width = .2)

dat$Om5 %>% max()
p.Om5 <- ggplot(dat, aes(x=sowndivLog, y=Om5, col=treatment))+
  geom_jitter(width = .2)

grid.arrange(p.Om4, p.Om5)
