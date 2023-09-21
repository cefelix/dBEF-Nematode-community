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




#### response variables####



