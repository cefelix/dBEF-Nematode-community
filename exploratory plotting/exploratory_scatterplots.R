#creating some exploratory scatterplots.

#always run https://github.com/cefelix/dBEF-Nematode-community/blob/main/wrangling/amynang-JenaXP_SP6_2021-wrangling-nematodes.R
#and https://github.com/cefelix/dBEF-Nematode-community/blob/main/wrangling/community%20indices.R
#before, as this will put the data into the necessary format

####libraries####
library(tidyverse)

####Enrichment index####
#scatterplot
ggplot(data.indices, aes(x = sowndiv, y = EI, color = treatment))+
  geom_point()+
  scale_x_continuous(trans = 'log2')

####Channel ratio####
#scatterplot
ggplot(data.indices, aes(x = sowndiv, y = CI, color = treatment))+
  geom_point()+
  scale_x_continuous(trans = 'log2')
  #here something with the maRcel::nemaplex function must have gone wrong

####Structure index####
ggplot(data.indices, aes(x = sowndiv, y = SI, color = treatment))+
  geom_point()+
  scale_x_continuous(trans = 'log2')

####Maturity index####
ggplot(data.indices, aes(x = sowndiv, y = MI, color = treatment))+
  geom_point()+
  scale_x_continuous(trans = 'log2')

