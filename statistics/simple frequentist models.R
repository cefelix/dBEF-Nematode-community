#analyzing the 2021 data with simple frequentist models
#reminder on stats: https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf

####packages####
library(lme4)
library(tidyverse)

####generalized mixed effect models####
#for now i will go with a binomial distribution, even though i am not sure how i should treat ecological indices
#(they originate from count data though)

####CP-CLASSES and overall abundance####
####
#I start here, just because these distributions look like negative binomial, and this seems for now easy to me

#lets look at the distributions
data.analysis$abundance_anja %>%
  hist() #negative binomial

data.analysis$cp1 %>%
  hist() #negative binomial
data.analysis$cp2 %>%
  hist() #negative binomial
data.analysis$cp3 %>%
  hist() #negative binomial
data.analysis$cp4 %>%
  hist()
data.analysis$cp5 %>%
  hist() #negative binomial

####cp-1, cp-2####
cp1_glmer <- glmer(cp1 ~ sowndiv * treatment + (1|block), data = data.analysis, family = "poisson") #this gives warnings as we input non-integer values
warnings()
summary(cp1_glmer) #we also should try out random slopes depending on treatment, according to Reich et al. 2012 (if that makes sence...)

#but, as the warnings complain about non-integer values, lets not use half nematodes and round our data (this is probably highly inappropriate):
data.analysis.round <- data.analysis
data.analysis.round$cp1 <- data.analysis$cp1 %>%
  round()
data.analysis.round$cp2 <- data.analysis$cp2 %>%
  round()
data.analysis.round$cp3 <- data.analysis$cp3 %>%
  round()
data.analysis.round$cp4 <- data.analysis$cp4 %>%
  round()
data.analysis.round$cp5 <- data.analysis$cp5 %>%
  round()
data.analysis.round$abundance_anja <- data.analysis$abundance_anja %>%
  round()

#repeat the code from the start of this section:
cp1_glmer <- glmer(cp1 ~ sowndiv * treatment + (1|block), data = data.analysis.round, family = "poisson") 
summary(cp1_glmer) #we also should try out random slopes depending on treatment, according to Reich et al. 2012 (if that makes sence...)

cp2_glmer <- glmer(cp2 ~ sowndiv * treatment + (1|block), data = data.analysis.round, family = "poisson") #large eigenvalue warning
summary(cp2_glmer)

cp2_glmer.nb <- glmer.nb(cp2 ~ sowndiv*treatment + (1|block), data =data.analysis.round) #large eigenvalue warning, model failed to converge warning
summary(cp2_glmer.nb) 

#hmm. Internet gives some ideas on troubleshooting: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html


####overall abundance####
abundance_glmer.nb <- glmer.nb(abundance_anja ~ treatment * sowndiv + (1|block), data = data.analysis.round)
summary(abundance_glmer.nb)

ggplot(data=data.analysis.round, aes(x=sowndiv, y=abundance_anja, color= treatment))+
  geom_point()+
  scale_x_continuous(trans = "log2")






####CR####
data.analysis$CR %>%
  hist()
  #right skewed beta distribution might be suitable

CR_glmer <- glmer(CR ~ sowndiv * treatment + (1|block) , data= data.analysis, family = "binomial")
summary(CR_glmer)




