#more elaborate exploratory plots

#always run https://github.com/cefelix/dBEF-Nematode-community/blob/main/wrangling/amynang-JenaXP_SP6_2021-wrangling-nematodes.R
#and https://github.com/cefelix/dBEF-Nematode-community/blob/main/wrangling/community%20indices.R
#before, as this will put the data into the necessary format

####libraries####
library(tidyverse)
library(ggplot2)

d <- data.analysis #less typing
  d1 <- subset(d, treatment == 1) #-S-P (?)
  d2 <- subset(d, treatment == 2) #+S-P (?)
  d3 <- subset(d, treatment == 3) #+S+P (?)
  d <- subset(d, sowndiv == 16)
  
ggplot(data = d, aes(x= SI, y= EI, col= block))+
  geom_point(aes(size=abundance))+
  scale_x_continuous(limits = c(0,100))+
  scale_y_continuous(limits = c(0,100))+
  facet_wrap(~treatment)
  

