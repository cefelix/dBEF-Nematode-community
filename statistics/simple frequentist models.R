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



####
####histograms of the cp1-cp5 nematodes. Assessing the abundances of nematodes by cp-class####
#this section might be moved to ./exploratory plotting

#lets create a dataframe with only 2 columns: cp-class and respective abundance within that class.
#each row represents the abundance cp-x nematodes in a sample -> each sample has 5 rows in this df
  #this is just for the convenience of plotting the histograms for each cp-class next to each other, using facet_wrap()!

cp.hist <- data.analysis$cp1 %>%
  append(data.analysis$cp2) %>%
  append(data.analysis$cp3) %>%
  append(data.analysis$cp4) %>%
  append(data.analysis$cp5)
cp.hist <- cbind(cp.hist, c(rep("cp1", 240 ), rep("cp2", 240), rep("cp3", 240), rep("cp4", 240), rep("cp5", 240))) %>%
  data.frame()
colnames(cp.hist) <- c("abun", "class")
cp.hist$abun <- cp.hist$abun %>% as.numeric()

#lets plot the histograms for each cp-class next to each other:  
ggplot(cp.hist, aes(x=abun, color=class))+
  geom_histogram(binwidth = 15)+
scale_x_continuous(limits = c(0,400))+
  facet_wrap(~class) #this is probably not very informative concerning cp1 and cp5 nematodes...

#let's check what the overall abundances in each cp-class are:
  data.analysis$cp1 %>% sum() #244
  data.analysis$cp2 %>% sum() #19102 
  data.analysis$cp3 %>% sum() #12203
  data.analysis$cp4 %>% sum() #5287
  data.analysis$cp5 %>% sum() #113
  #okay, this might be a problem: the very low abundances in cp1 and cp5 might have a big influence in the data analysis, while probably not being very reliable...
  
  #They will clearly have a huge influence on CI (maRcel::Channel(), as it jumps to 1 if there are no Ba1 nematodes in the sample 
  #as all cp1-nematodes are bacterivores, "no cp1 nematodes in sample" is synonymous with "no Ba1 nematodes in sample"
  #this also explains this:
  summary(data.analysis$CI)





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




