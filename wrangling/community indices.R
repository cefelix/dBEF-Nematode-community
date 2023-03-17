#calculating the community indices using maRcel

####required packages####
  library(devtools)
  install_github("amynang/marcel")
  #install.packages("RSelenium")
  #install_github("BEXIS2/rBExIS", subdir = "rBExIS") 

####loading libraries#### 
  library(maRcel)
  library(tidyverse)
  library(RSelenium)
  lsf.str("package:maRcel")
  library(httr)
  library(jsonlite)
  library(XML)
  library(rBExIS)


####data checking and re-arranging####
#get data by re-running code from: https://github.com/amynang/JenaXP_SP6_2021/blob/main/wrangling/nematodes.R

#rename the first raw data df:  
data.1 <- data 
  
head(data.3) #nematode densities per soil sample
head(data.4) #nematode densities per 100g soil dry weight
#head(data.7) #nematode densities per 1 square meter, in 0-10 cm depth


bexis.options(base_url = "https://jexis.uni-jena.de")
#get data from bexis from dataset with id = xy
main.plot = bexis.get.dataset_by(id = 90)

data.4 <- data.4 %>% 
  # we change plot to be "blockplot" by splitting Sample at D and keeping the 1st half
  mutate(plot = str_split(.$Sample, "D", simplify = T)[,1]) %>% 
  # we add a sowndiv column whose elements are filled from main.plot
  # based on plot matching
  mutate(.after = plot,
         sowndiv = as.character(main.plot$sowndiv[match(.$plot, main.plot$plotcode)]))

head(data.4) #now sown diversity is the 4th column



####maRcel testing####
#test_taxa <- c(colnames(data.4[6:8]))
#query_nemaplex(test_taxa) #Error in if (file.access(phantompath, 1) < 0) { : argument is of length zero
#system("taskkill /im java.exe /f", intern=FALSE, ignore.stdout=FALSE)

#query_nemaplex("Actinolaimidae") #Error in if (file.access(phantompath, 1) < 0) { : argument is of length zero
#system("taskkill /im java.exe /f", intern=FALSE, ignore.stdout=FALSE)



####indices calculation ####
data.nplx <- read.csv("./wrangling/nemaplex.csv", row.names = 1)


#CP proportions
data.4 %>%
  C_P(nemaplex = data.nplx) %>%
  summary()

#channel index
data.4 %>%
  Channel(nemaplex = data.nplx) %>% 
  summary() #thats a Median CI of 100, seems sketchy, should look into it
            #after looking into it, a CI of 100 simply means that there are no Ba1 nematodes
            #in the samples

#maturity index
data.4 %>% 
  Maturity(nemaplex = data.nplx) %>%
  summary()

#enrichment index
data.4 %>% 
  Enrichment(nemaplex = data.nplx) %>%
  summary()

#structure index
data.4 %>%
  Structure(nemaplex = data.nplx) %>%
  summary()

####create table with all FUNCTIONAL INDICES included in maRcel####
#using the fancy function which does all indices at once:
data.4$sowndiv <- data.4$sowndiv %>%
  as.factor() #this prevents sowndiv being dropped by all.índices()
data.indices <- data.4 %>%
  all.indices(nemaplex = data.nplx) #creating a new df with all indices and bloc/plot/treatment info
data.indices$sowndiv <- data.indices$sowndiv %>%
  as.character() %>%
  as.numeric() #converting sowndiv back to numeric
head(data.indices)
str(data.indices)

####add CHANNEL RATIO as described by Dietrich et al. 2021####

#This function calculates the channel ratio as Fu/(Fu+Ba):
ChannelRatio <- function(df, nemaplex)
  {
    Ba = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 
                  3)
    Fu = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 
                  2)
    index = df %>% mutate(CR = (rowSums(.[Fu])/(rowSums(.[Ba]) + rowSums(.[Fu]))), .keep = "none")
    out = cbind(index)
    return(out)
  }

#we add a new column to our index data frame which shows the channel ratio (CR)
data.indices <- cbind(data.indices, ChannelRatio(data.4, data.nplx))

#this is a test comment, as i am struggling with merging errors
#this is another test comment


####Adding SOIL MOISTURE ####
#calculating soil moisture as proportion of water in the fresh soil used for extraction
moisture <- soil$water.content / soil$init.weight 
data.indices <- cbind(data.indices, moisture)


####changing data type of variables####

data.indices$block <- data.indices$block %>% 
  as.factor()
data.indices$treatment <- data.indices$treatment %>%
  as.factor()

str(data.indices)



####NEW approach for a table with units appropriate to respective response variable### 

#what variables do I want, are they response/explanatory/covariates, which unit:
  #EI, response, index 0-100
  #SI
  #CI
  #MI

  #cp1 - cp5, response, ind. per 100g DW
  #abundance, response, ind. per 100g DW

  #soil moisture, covariate, %water of FW
  

#lets start with CI (channel index)
  data.analysis <- data.4%>%
    Channel(data.nplx) 

#add CR (channel ratio)        
  data.analysis <- cbind(data.analysis, ChannelRatio(data.4, data.nplx))
  
#add EI (enrichment index)
  data.analysis <- cbind(data.analysis, Enrichment(data.4, data.nplx)[6])
  
#add SI (structure index)
  data.analysis <- cbind(data.analysis, Structure(data.4, data.nplx)[6])

#add MI (maturity index) 
  data.analysis <- cbind(data.analysis, Maturity(data.4, data.nplx))

#extract abundances per sample as counted by Anja, and add join by Sample ID:
  abun_anja <- abun %>%
    data.frame() 
  names(abun_anja) <- c("Date", "Sample", "Plot", "Subplot", "abundance_anja")
  abun_anja <- abun_anja[-c(1,3:4)]
  data.analysis <- full_join(data.analysis, abun_anja, by = "Sample")
  
#calculate cumulated identifications by Marcel, and join by corresponding Sample IDs:
  abun_marcel <- data.1[2:65] %>%
    rowSums() 
  abun_marcel <- cbind(abun_marcel, data.1$Sample) %>%
    data.frame() 
  colnames(abun_marcel) <- c("abundance_marcel","Sample")
  abun_marcel$abundance_marcel <- abun_marcel$abundance_marcel %>%
    as.numeric()
  data.analysis <- full_join(data.analysis, abun_marcel, by="Sample")
  
  data.analysis$abundance_anja %>%
    summary()
  data.analysis$abundance_marcel %>%
    summary()
  summary(data.analysis)  

#have a brief look into the discrepancy between Anja's counts and Marcel's cumulated identifications
  ggplot(data.analysis, aes(x=abundance_anja, y=abundance_marcel))+
    geom_point()+
    geom_point(aes(x=abundance_anja, y=abundance_anja, color="red"))+
    scale_y_continuous(limits = c(0,100))+
    scale_x_continuous(limits = c(0,100))
  
#add cp classes (as individuals per 100g dw)  
   
  