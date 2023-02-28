#calculating the community indices (maybe using maRcel?)


####loading libraries#### 
#library(devtools)
#install_github("amynang/marcel")
#install.packages("RSelenium")
library(maRcel)
library(tidyverse)
library(RSelenium)
lsf.str("package:maRcel")

####data checking####
#get data by re-running code from: https://github.com/amynang/JenaXP_SP6_2021/blob/main/wrangling/nematodes.R

head(data.3) #nematode densities per soil sample
head(data.4) #nematode densities per 100g soil dry weight
head(data.7) #nematode densities per 1 square meter, in 0-10 cm depth

####maRcel testing####
test_taxa <- c(colnames(data.4[5:7]))
query_nemaplex("Actinolaimidae") #Error in if (file.access(phantompath, 1) < 0) { : argument is of length zero

system("taskkill /im java.exe /f", intern=FALSE, ignore.stdout=FALSE)

####maturity index - manual calculation####
#as maRcel produces an error atm, i will create a dataframe with c-p-values manually:

#extract all observed taxa from our data frame
taxa <- data.4[5:68] %>% 
  colnames() %>%
  as.data.frame() 

#create vectors for different indices of the same length
cp <- rep(NA, nrow(taxa))
feeding <- rep(NA, nrow(taxa))

#combine the df 'taxa' with the newly created vectors
taxa <- cbind(taxa, cp, feeding)



