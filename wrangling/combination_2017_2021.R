#combining Vogel's 2017 and Amyntas' 2021 data
#to be safe, re-run the wrangling files for 2017/2021's data previous to this!

#libraries
library(dplyr)
library(tidyr)

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

# a function to combine "year" and "treatment" into "years of Soil history":
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
# treatment codes: 1 = -SH-PH, 2 = +SH-PH, 3 = +SH+PH
# years of setup: 2016, 2002

#testing it:
a <- calc_history(dBEF_nem, treatment = dBEF_nem$treatment, year = dBEF_nem$year, type = "plant")
hist(a) #looks good





#### response variables####



