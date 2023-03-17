#creating some exploratory scatterplots.

#always run https://github.com/cefelix/dBEF-Nematode-community/blob/main/wrangling/amynang-JenaXP_SP6_2021-wrangling-nematodes.R
#and https://github.com/cefelix/dBEF-Nematode-community/blob/main/wrangling/community%20indices.R
#before, as this will put the data into the necessary format

####libraries####
library(tidyverse)

####Enrichment index####
#scatterplot
ggplot(data.analysis, aes(x = sowndiv, y = EI, color = treatment))+
  geom_point()+
  scale_x_continuous(trans = 'log2')

#one graph for each treatment using facet_wrap()
ggplot(data.analysis, aes(x = sowndiv, y = EI, color = block))+
  geom_point()+
  facet_wrap(~treatment)+
  scale_x_continuous(trans = 'log2')

#one graph for each block
ggplot(data.analysis, aes(x = sowndiv, y = EI, color = treatment))+
  geom_point()+
  facet_wrap(~block)+
  scale_x_continuous(trans = 'log2')



####Channel index####
#scatterplot
ggplot(data.analysis, aes(x = sowndiv, y = CI, color = treatment))+
  geom_point()+
  scale_x_continuous(trans = 'log2')
  #here something with the maRcel::nemaplex function must have gone wrong

#one graph for each treatment
ggplot(data.analysis, aes(x = sowndiv, y = CI, color = block))+
  geom_point()+
  facet_wrap(~treatment)+
  scale_x_continuous(trans = 'log2')

#one graph for each block
ggplot(data.analysis, aes(x = sowndiv, y = CI, color = treatment))+
  geom_point()+
  facet_wrap(~block)+
  scale_x_continuous(trans = 'log2')

####Channel ratio####
#scatterplot
ggplot(data.analysis, aes(x = sowndiv, y = CR, color = treatment))+
  geom_point()+
  scale_x_continuous(trans = 'log2')
#here something with the maRcel::nemaplex function must have gone wrong

#one graph for each treatment
ggplot(data.analysis, aes(x = sowndiv, y = CR, color = block))+
  geom_point()+
  facet_wrap(~treatment)+
  scale_x_continuous(trans = 'log2')

#one graph for each block
ggplot(data.analysis, aes(x = sowndiv, y = CR, color = treatment))+
  geom_point()+
  facet_wrap(~block)+
  scale_x_continuous(trans = 'log2')



####Structure index####
ggplot(data.analysis, aes(x = sowndiv, y = SI, color = treatment))+
  geom_point()+
  scale_x_continuous(trans = 'log2')

#one graph for each treatment
ggplot(data.analysis, aes(x = sowndiv, y = SI, color = block))+
  geom_point()+
  facet_wrap(~treatment)+
  scale_x_continuous(trans = 'log2')

#one graph for each block
ggplot(data.analysis, aes(x = sowndiv, y = SI, color = treatment))+
  geom_point()+
  facet_wrap(~block)+
  scale_x_continuous(trans = 'log2')


####Maturity index####
ggplot(data.analysis, aes(x = sowndiv, y = MI, color = treatment))+
  geom_point()+
  scale_x_continuous(trans = 'log2')

#one  graph for each treatment
ggplot(data.analysis, aes(x = sowndiv, y = MI, color = block))+
  geom_point()+
  facet_wrap(~treatment)+
  scale_x_continuous(trans = 'log2')

#one graph for each block
ggplot(data.analysis, aes(x = sowndiv, y = MI, color = treatment))+
  geom_point()+
  facet_wrap(~block)+
  scale_x_continuous(trans = 'log2')

####share of cp1 nematodes in sample####
ggplot(data.analysis, aes(x = sowndiv, y = cp1, color = treatment ))+
  geom_point()+
  scale_x_continuous(trans='log2')



####soil moisture####
ggplot(data.analysis, aes(x = block, y = percent_water))+
  geom_boxplot()
  #soil moisture should be included in the model, if we use block as a random effect
  



