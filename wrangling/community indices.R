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

#also, read in the nemaplex data - this is important, KEEP THIS LINE! (message to myself)
data.nplx <- read.csv("./wrangling/nemaplex.csv", row.names = 1)
  
  #I added family level information on specimens identified till genus level:
  #genus identification only changes information in families: 
    #Aporcelimidae, Dorylaimidae, Qudsianematidae, Tylenchidae
  


####maRcel testing####
#test_taxa <- c(colnames(data.4[6:8]))
#query_nemaplex(test_taxa) #Error in if (file.access(phantompath, 1) < 0) { : argument is of length zero
#system("taskkill /im java.exe /f", intern=FALSE, ignore.stdout=FALSE)

#query_nemaplex("Actinolaimidae") #Error in if (file.access(phantompath, 1) < 0) { : argument is of length zero
#system("taskkill /im java.exe /f", intern=FALSE, ignore.stdout=FALSE)



####add a function to calculate CHANNEL RATIO as described by Dietrich et al. 2021####
####
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

####create functions to determine fungal and and bacterial feeder proportions####
####

#proportion of bacterial feeders
bacterial <- function(df, nemaplex) {
  Ba = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3)
  abun.Ba = df %>% mutate(bacterivores = rowSums(.[Ba]), .keep = "none")
  abun.all = df %>% 
    select_if(is.numeric,) %>%
    rowSums()
  density = abun.Ba / abun.all
  return(density)
}

bacterial(data.4, data.nplx)  

#proportion of fungal feeders
fungal <- function(df, nemaplex) {
  Fu = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2)
  abun.Fu = df %>% mutate(fungivores = rowSums(.[Fu]), .keep = "none")
  abun.all = df %>% 
    select_if(is.numeric,) %>%
    rowSums()
  density = abun.Fu / abun.all
  return(density)
}

#lets generalize these functions to one single function, with an additional argument for feeding type:
#UNDER CONSTRUCTION
feeding_type <- function(df, nemaplex, feeding) {

}
#feeding types are coded according to yeates et al. 1993 (non listed categories dont appear in data)
  #1 - herbivore (plant feeder)
  #2 - fungivore (hyphal feeder)
  #3 - bacterivore 
  #5 - predator
  #8 - omnivore
  


####create a table with appropriate units#### 

#what variables do I want, are they response/explanatory/covariates, which unit:
  
  #Indices (0-100 / 0-1), all response:
  #EI
  #SI
  #CI
  #MI
  
  #Densities (ind. per 100g DW), all response
  #cp1 - cp5
  #abundance
  #bacterial feeders
  #fungal feeders

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
  
#add relative cp abundances
  data.analysis <- cbind(data.analysis, C_P(data.4, data.nplx)[6:10])

#extract abundances per sample as counted by Anja, and add join by Sample ID:
  abun_anja <- abun %>%
    data.frame() 
  names(abun_anja) <- c("Date", "Sample", "Plot", "Subplot", "abundance")
  abun_anja <- abun_anja[-c(1,3:4)]
  data.analysis <- full_join(data.analysis, abun_anja, by = "Sample")
  
#calculate cumulated identifications by Marcel, and join by corresponding Sample IDs:
  abun_marcel <- data.1[2:65] %>%
    rowSums() 
  abun_marcel <- cbind(abun_marcel, data.1$Sample) %>%
    data.frame() 
  colnames(abun_marcel) <- c("counts.cummulated","Sample")
  abun_marcel$counts.cummulated <- abun_marcel$counts.cummulated %>%
    as.numeric()
  abun_marcel$Sample[240] <- "B1A12D3" #replacing the unlabeled sample, without this line data.analysis gets 241 rows
  data.analysis <- full_join(data.analysis, abun_marcel, by="Sample") 
  
  data.analysis$abundance %>%
    summary()
  data.analysis$counts.cummulated %>%
    summary()
  summary(data.analysis)  

#Convert to no. of individuals per 100g DW 
  #first check whether net.weight is the net weight of the dried soil:
  soil$net.weight - (soil$dry.weight - soil$pot.weight) < 0.01 #presumably yes
  #add DW per sample to data.analysis:
  DW_sample <- soil[-c(2:6,8,9)] %>%
    data.frame()
  names(DW_sample) <- c("Sample", "DW")
  data.analysis <- full_join(data.analysis, DW_sample, by="Sample")
  #perform the magic, drop the DW column afterwards:
  data.analysis$abundance <- data.analysis$abundance*100/data.analysis$DW
  data.analysis$counts.cummulated <- data.analysis$counts.cummulated*100/data.analysis$DW
  data.analysis <- data.analysis[-18]
    
  
  
  
  
#have a brief look into the discrepancy between Anja's counts and Marcel's cumulated identifications 
#(only up to 100 individuals as more have not been identified)
  ggplot(data.analysis, aes(x=abundance, y=counts.cummulated, color= block))+
    geom_point()+
    geom_line(aes(x=abundance, y=abundance, color= "y = abundance"))+
    scale_y_continuous(limits = c(0,100))+
    scale_x_continuous(limits = c(0,100))
  
#convert cp1 - cp5 from proportions to individuals per 100g DW
  data.analysis$cp1 <- data.analysis$cp1 * data.analysis$abundance
  data.analysis$cp2 <- data.analysis$cp2 * data.analysis$abundance
  data.analysis$cp3 <- data.analysis$cp3 * data.analysis$abundance
  data.analysis$cp4 <- data.analysis$cp4 * data.analysis$abundance
  data.analysis$cp5 <- data.analysis$cp5 * data.analysis$abundance
  #check whether there are big differences in the sums:
  rowSums(data.analysis[11:15]) - data.analysis$abundance < 0.1 #looks good
  
#calculate abundance of fungivores and bacterivores per 100g DW
  data.analysis <- cbind(data.analysis, bacterial(data.4, data.nplx) * data.analysis$abundance)
  data.analysis <- cbind(data.analysis, fungal(data.4, data.nplx) * data.analysis$abundance)
  
  
#calculating soil moisture as %water in relation to fresh weight and join by corresponding sample number
  soil$percent_water <- soil$water.content / soil$init.weight * 100
  soil_moisture <- soil[-c(2:8)]
  data.analysis <- full_join(data.analysis, soil_moisture, by="Sample")
  

  
  
####adjust classes of columns to appropriate ones####
####
#block and treatment to factor, sown diversity to numeric:
data.analysis$block <- data.analysis$block %>% 
  as.factor()
data.analysis$treatment <- data.analysis$treatment %>%
  as.factor()
data.analysis$sowndiv <- data.analysis$sowndiv %>%
  as.numeric()
str(data.analysis)
  
  