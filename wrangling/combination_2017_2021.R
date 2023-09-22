#combining Vogel's 2017 and Amyntas' 2021 data
#to be safe, re-run the wrangling files for 2017/2021's data previous to this!

#libraries
library(dplyr)
library(tidyr)
library(maRcel)

amyntas2021 <- read.csv("./wrangling/Amyntas2021.csv", row.names = 1)
vogel2017 <- read.csv("./wrangling/Vogel2017.csv", row.names = 1)

###
####merge into one table####
###

#add the taxons which are found only in amyntas' data to vogel's data 
taxons <- setdiff(colnames(amyntas2021), colnames(vogel2017))
add_vo <- data.frame(matrix(0,
                            nrow = nrow(vogel2017),
                            ncol = length(taxons)) ) 
colnames(add_vo) = taxons
vogel2017 <- cbind(vogel2017, add_vo)

#vice versa for vogels data
taxons  <- setdiff(colnames(vogel2017), colnames(amyntas2021))
add_am <- data.frame(matrix(0,
                            nrow = nrow(amyntas2021),
                            ncol = length(taxons)) ) 
colnames(add_am) = taxons
amyntas2021 <- cbind(amyntas2021, add_am)

#clean the mess:
rm(add_am, add_vo, taxons)

#merge both data sets
dBEF_nem <- full_join(amyntas2021, vogel2017)


#### independent variables####

#creating a function to combine "year" and "treatment" into "years of Soil history":
calc_history <- function(data ,treatment, year, type){
  #this calculates years of soil history based on treatment and year of sampling
  if (type == "soil"){
    SH <- rep(NA, nrow(data))
    for (i in 1:nrow(data)) {
      if(treatment[i] == "1") {         # -SH-PH
        SH[i] = year[i]-2016
      } else if (treatment[i] == "2") { # +SH-PH
        SH[i] = year[i]-2002
      } else if (treatment[i] == "3") { # +SH+PH
        SH[i] = year[i]-2002
      } 
      else{
        return("treatment must be 1, 2, or 3")
      }
    } 
    return(SH)
  }
  #this calculates years of plant history based on treatment and year of sampling
  else if (type == "plant") {
    PH <- rep(NA, nrow(data))
    for (i in 1:nrow(data)) {
      if(treatment[i] == "1") {         #treatment 1: -SH-PH
        PH[i] = year[i]-2016
      } else if (treatment[i] == "2") { #treatment 2: +SH-PH
        PH[i] = year[i]-2016
      } else if (treatment[i] == "3") { #treatment 3: +SH+PH
        PH[i] = year[i]-2002
      } 
      else{
        return("treatment must be 1, 2, or 3")
      }
    } 
    return(PH)
    
    
  }
}

#adding columns with years of PH / SH:
dBEF_nem <- dBEF_nem %>%
  mutate(., SH = calc_history(., .$treatment,. $year, type="soil"), .after= sowndiv) %>%
  mutate(., PH = calc_history(., .$treatment,. $year, type="plant"), .after= sowndiv) 



#### response variables####

#total abundance (individuals per 100g DW)
#column indices of taxa:
taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))

dBEF_nem <- dBEF_nem %>%
  mutate(ind_per100g = rowSums(dBEF_nem[taxa]), .before = Acrobeles)


  
  grep("Acrobeles", colnames(.))

rowSums()
which(colnames(dBEF_nem)== "Acrobeles")


#densities of trophic guilds:

  # Ba per 100g
  
  # Fu per 100g
  
  # Pl per 100g 
  
  # Om per 100g
  
  # Ca per 100g


#densities by life strategy:
  
  # cp-1 per 100g

  # cp-2 per 100g

  # cp-3 per 100g

  # cp-4 per 100g

  # cp-5 per 100g










####how to deal with dauer larvae?####
#dauer larvae problem: they could be from each of the three occurring genera of rhabditidae:
dBEF_nem$Rhabditidae.dauer.larvae %>%
  sum() #9118
dBEF_nem$Mesorhabditis %>%
  sum() #31 (Ba-cp1)
dBEF_nem$Protorhabditis %>%
  sum() #67 (Ba-cp1)
dBEF_nem$Rhabditis %>%
  sum() #9170 (Ba-cp1) -> but most likely Rhabditis

####save the whole file as .csv####
write.csv(dBEF_nem, "./dBEF_nem.csv")



####trash####

nemaplex <- read.csv("./wrangling/nemaplex.csv")
nemaplex$X

setdiff(nemaplex$X, colnames(dBEF_nem))
setdiff(colnames(dBEF_nem), nemaplex$X)


