#analyzing the 2021 data with simple frequentist models
#reminder on stats: https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf

####packages####
library(lme4)
library(tidyverse)

####generalized mixed effect models####
#for now i will go with a binomial distribution, even though i am not sure how i should treat ecological indices
#(they originate from count data though)

####CR####
data.analysis$CR %>%
  hist()
  #right skewed beta distribution might be suitable

CR_glmer <- glmer(CR ~ sowndiv * treatment + (1|block) , data= data.analysis, family = "binomial")
summary(CR_glmer)

####CI####
data.analysis$CI %>%
  hist() #can this be informative?

####EI####
data.analysis$EI %>%
  hist()



