#analyzing the 2021 data with simple frequentist models
#reminder on stats: https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf

####packages####
library(lme4)
library(tidyverse)

####generalized mixed effect models####
#for now i will go with a binomial distribution, even though i am not sure how i should treat ecological indices
#(they originate from count data though)

####CR####
data.indices$CR %>%
  hist()

CR_glmer <- glmer(CR ~ sowndiv + treatment + (1|block) , data= data.indices, family = "binomial")
summary(CR_glmer)

