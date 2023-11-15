library(ggplot2)
library(dplyr)

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
dBEF_nem$year <- dBEF_nem$year %>% as.factor()
dBEF_nem %>%
  str()

dBEF_nem %>%
  ggplot(aes(y=ind_per100g, x=log(sowndiv), col = year))+
  #geom_jitter()
  geom_point()

dBEF_nem %>%
  ggplot(aes(y=ind_per100g, x=SH, col = sowndiv))+
  #geom_jitter()
  geom_point()

##per trophic group
dBEF_nem %>%
  ggplot(aes(y=Fu_per100g, x = sowndiv, col = year))+
  geom_point()
dBEF_nem %>%
  ggplot(aes(y=Ba_per100g, x = sowndiv, col = year))+
  geom_point()
dBEF_nem %>%
  ggplot(aes(y=Pl_per100g, x = sowndiv, col = year))+
  geom_point()
dBEF_nem %>%
  ggplot(aes(y=Om_per100g, x = sowndiv, col = year))+
  geom_point()
dBEF_nem %>%
  ggplot(aes(y=Pr_per100g, x = sowndiv, col = year))+
  geom_point()

##per cp group
dBEF_nem %>%
  ggplot(aes(y=cp1_per100g, x = sowndiv, col = year))+
  geom_point(alpha = 0.5)

dBEF_nem %>%
  ggplot(aes(y=cp2_per100g, x = sowndiv, col = year))+
  geom_point()


dBEF_nem %>%
  ggplot(aes(y=cp3_per100g, x = sowndiv, col = year))+
  geom_point()

dBEF_nem %>%
  ggplot(aes(y=cp4_per100g, x = sowndiv, col = year))+
  geom_point()

dBEF_nem %>%
  ggplot(aes(y=cp5_per100g, x = sowndiv, col = year))+
  geom_point()
