#analyzing the 2021 data with simple frequentist models
#reminder on stats: https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf

####01packages####
library(lme4)
library(tidyverse)
library(lmerTest)

####02 INAPPROPRIATE glmm for cp2####
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
####11 histograms of the cp1-cp5 nematodes. Assessing the abundances of nematodes by cp-class####
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


####20 GLMM's for bacterivores, fungivores, channel ratio####
  #trying to recreate analysis as Dietrich et al. 2021 reportet in fig.2b/c and table 2:
  #check out beta distribution for Channel ratio
  
#lets check out how many nematodes should live in our 24kg of soil:
data.analysis$bacterivores %>% sum() #a lot
data.analysis$fungivores %>% sum() #a lot * 3
  
  
#histograms first
ggplot(data.analysis, aes(x=bacterivores))+
  geom_histogram(binwidth = 25) #unimodal, right skewed (="left leaning")
  #gamma distribution

ggplot(data.analysis, aes(x=fungivores))+
  geom_histogram(binwidth = 25) #unimodal, right skewed
  #gamma distribution   

ggplot(data.analysis, aes(x=CR))+
  geom_histogram(binwidth = 0.05) #left skewed
  #CR = Fu / (Fu + Ba)

####21 fitting a glmm on ba/fu density using gamma distribution####
ba.glmm.gamma <- glmer(bacterivores ~ sowndiv * treatment + (1|block), data = data.analysis, family = Gamma) #error: no zeros allowed
 summary(data.analysis$bacterivores == 0) # 16 out of 240 :(
fu.glmm.gamma <- glmer(fungivores ~ sowndiv * treatment + (1|block), data = data.analysis, family = Gamma)
  summary(data.analysis$fungivores == 0) # 4 out of 240 :(
  
#just add 0.1**6 to every density to get rid of the problem:
data.FuBa <- data.analysis
  data.FuBa$bacterivores <- data.FuBa$bacterivores +  0.1**6
  data.FuBa$fungivores <- data.FuBa$fungivores +  0.1**6
  
#retry:   
  summary(data.FuBa$bacterivores == 0)
  ba.glmm.gamma <- glmer(bacterivores ~ sowndiv * treatment + (1|block), data = data.FuBa, family = Gamma(link=log)) 
    #specifying link function to log prevents error: "PIRLS loop resulted in NaN value"
    summary(ba.glmm.gamma)
  fu.glmm.gamma <- glmer(fungivores ~ sowndiv * treatment + (1|block), data = data.FuBa, family = Gamma(link=log)) 
    summary(fu.glmm.gamma)

####22 check assumptions on ba.glmm.gamma####
ba.glmm.gamma %>%
  residuals() %>%
  plot() #slightky curved in the right quarter(?), like bottom of a quadratic relationship
  
  #residuals vs fitted
  ba.pred <- predict(ba.glmm.gamma)
  ba.resid <- residuals(ba.glmm.gamma)
  plot(ba.pred, ba.resid) #maybe remove the outliers, to see the pattern of the remaining data better
  
  #histogram of residuals
  hist(ba.resid) # skewed... the 16 outliers on the left...
      
        
####23 check assumptions on fu.glmm.gamma####            
fu.glmm.gamma %>%
  residuals() %>%
  plot() #this looks okay

  #residuals vs fitted
  fu.pred <- predict(fu.glmm.gamma)
  fu.resid <- residuals(fu.glmm.gamma)
  plot(fu.pred, fu.resid)
  
  #histogram of residuals
  hist(fu.resid) #left skewed... 4 outliers left?
  
  
####31 fitting a glmmm on CR using beta distribution####
  #they told me that it's possible to use package::glmmTMB
  #something on beta distributions in plain language https://rpubs.com/nicoleknight/936037
    #there they also mention to read chapter 10 and 11 from McElreath 

library(glmmTMB)  
  
CR.glmm.beta <- glmmTMB(CR ~ sowndiv * treatment + (1|block), data = data.analysis, family = beta_family)  #error: y values must be 0<y<1
  summary(data.analysis$CR == 0) #3
  summary(data.analysis$CR == 1) #15
  #...and 1 NA (so one sample with neither a fungivore, nor a bacterivore)
  
  #lets add/substract 0.1**6 to the CR's of 0/1
  data.analysis.CR <- data.analysis
  data.analysis.CR <- data.analysis.CR[is.na(data.analysis.CR$CR) == FALSE,] #dropping the NA row
  data.analysis.CR[data.analysis.CR$CR==0,]$CR <- rep(0.1**6, 3) #replacing zeros with very low value
  data.analysis.CR[data.analysis.CR$CR==1,]$CR <- rep(1-0.1**6, 15) #replacing 1 with almost 1
  
#re-run the stuff from above
CR.glmm.beta <- glmmTMB(CR ~ sowndiv * treatment + (1|block), data = data.analysis.CR, family = beta_family) 
  summary(CR.glmm.beta)
  
  
####32 checking assumptions of the CR.glmm.beta#### 
  CR.glmm.beta %>%
    residuals() %>%
    plot()
  
  #fitted vs residuals
  CR.pred <- predict(CR.glmm.beta)
  CR.resid <- residuals(CR.glmm.beta) 
  plot(CR.pred, CR.resid)
  
  #histogram of residuals
  hist(CR.resid) #left skewed
  

####41 fitting a log-transformed lmer on ba density , as Dietrich 2021####
ba.lmer <- lmer(log(bacterivores) ~ sowndiv * treatment + (1|block), data = data.analysis) #error NA/NaN/Inf in y
  log(data.analysis$bacterivores) %>%
    summary() #minimum -Inf, due to log(0), 16 rows in total 
  #lets drop the lines where we have zero bacterivores:
  data.analysis.ba.lmer <- data.analysis[(data.analysis$bacterivores > 0) == TRUE,]

#re-run the model above
ba.lmer <- lmer(log(bacterivores) ~ sowndiv*treatment + (1|block), data = data.analysis.ba.lmer)
summary(ba.lmer) #Dietrich doesn't specify how they got their p-values, as lmer does not output them by default


  #fitted vs residuals
  ba.lmer.pred <- predict(ba.lmer)
  ba.lmer.resid <- residuals(ba.lmer)
  plot(ba.lmer.pred, ba.lmer.resid)
  
  #histogram of residuals and response
  hist(ba.lmer.resid) #thats normal 
  data.analysis.ba.lmer$bacterivores %>%
    log() %>%
    hist() #thats not centered at zero, but the shape could be almost gaussian bell like
  
  
  
  
####42 fitting a lmer on CR, as Dietrich 2021####
  CR.lmer <- lmer(CR ~ sowndiv * treatment + (1|block), data = data.analysis)
  summary(CR.lmer)
  
  #fitted vs residuals:
  plot(CR.lmer)
  
  #histogram of residuals and response:
  CR.lmer %>%
    residuals() %>%
    hist() #that looks normal
  data.analysis$CR %>%
    hist() #that certainly doesn't
  
    
####51 just a check on gamma distribution based models####
#stole this code from https://stats.stackexchange.com/questions/390063/what-are-the-assumptions-of-a-gamma-glm-or-glmm-for-hypothesis-testing
library(tidyverse)

#Simulate
set.seed(0)
N <- 250
x <- runif(N, -1, 1)
a <- 0.5
b <- 1.2
y_true <- exp(a + b * x)
shape <- 2
y <- rgamma(N, scale = y_true/shape, shape = shape)

#Model
plot(x,y)
model <- glm(y ~ x, family = Gamma(link = "log"))
summary(model)

#Deviance residuals
ypred = predict(model)
res = residuals(model, type = 'deviance')
plot(ypred,res)
hist(res)

#Deviance GOF

deviance = model$deviance
p.value = pchisq(deviance, df = model$df.residual, lower.tail = F)
# Fail to reject; looks like a good fit.

 
####61 glmm bacterivore/fungivore density using normal distribution####
ba.glmm.normal <- lmer(bacterivores ~ sowndiv * treatment + (1|block), data = data.FuBa)
summary(ba.glmm.normal)

  ba.glmm.normal %>% plot() # 3 clouds, problem?
  ba.glmm.normal %>% 
    residuals() %>%
    hist() #thats not normal, but rather gamma
  
  qqnorm(ba.glmm.normal, ~ranef(., level = 0)) #error

fu.glmm.normal <- lmer(fungivores ~ sowndiv * treatment + (1|block), data = data.FuBa)
summary(fu.glmm.normal)

  fu.glmm.normal
  
  fu.glmm.normal %>%
    residuals() %>%
    hist()
  
  qqnorm(fu.glmm.normal)
  
    
  
  
  

