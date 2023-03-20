#analyzing the 2021 data with simple frequentist models
#reminder on stats: https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf

####packages####
library(lme4)
library(tidyverse)

####INAPPROPRIATE glmm for cp2####
#Negative binomial and poisson distribution are not appropriate! Our data are densities, not counts! Still, I will keep this section as a reminder. 

cp2_glmer <- glmer(cp2 ~ sowndiv * treatment + (1|block), data = data.analysis, family = "poisson") #large eigenvalue warning
warnings() #NOT TO USE. poisson distribution for continuous variable doesn't make sense. 
summary(cp2_glmer) 
  #NOTE: we also should try out random slopes depending on treatment, according to Reich et al. 2012 (if that makes sense...)

cp2_glmer.nb <- glmer.nb(cp2 ~ sowndiv*treatment + (1|block), data =data.analysis) #large eigenvalue warning, model failed to converge warning
warnings() #NEITHER USE THIS. negative binomial distribution also needs discrete responses!
summary(cp2_glmer.nb) 

#For future negative binomial distribution related warnings, the internet gives some ideas on troubleshooting: 
#https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html



####
####histograms of the cp1-cp5 nematodes. Assessing the abundances of nematodes by cp-class####
#this section might be moved to ./exploratory plotting

#lets create a dataframe with only 2 columns: cp-class and respective abundance within that class.
#each row represents the abundance cp1/2/.../5 nematodes in a sample -> each sample has 5 rows in this df
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


####GLMM for bacterivores, fungivores, channel ratio####
  #trying to recreate analysis as Dietrich et al. 2021 reportet in fig.2b/c and table 2:
  #check out beta distribution (?)







