#re-sampling the nematode assemblages as n (=number of nematodes in sample) draws with replacement,
#with the probability pi of drawing species i being the relative abundance of species i in the sample

#maRcel::extrapolate.robust by Angelos Amnytas

#e.g.: if a sample of 100 individuals has 20 identified individuals from 3 different species; 
  #of which 15 individuals are species 1, 4 are species 2 and 1 is species 3;
  #then n=100, and pi is {p1 = 0.75, p2=0.20, p3=0.05}



#libraries:
library(dplyr)
library(tidyr)
library(tidyverse)
library(maRcel)
library(vegan)

#data:
dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)


####standardize taxon abundances:####
taxa <- grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem)
dBEF_nem[taxa] <- decostand(dBEF_nem[taxa], method = "total", MARGIN = 1)
  #check it:
  rowSums(dBEF_nem[taxa])
  


####bootstrap function doing almost the same as maRcel::extrapolate.robust####
#creating data:
taxa = matrix(c(.1,0,.2,.7,0,0,.4,.6,.8,.1,.1,0),
                  byrow = T, nrow = 3, ncol = 4,
                  dimnames = list(c("site1","site2","site3"),
                                  c("sp1","sp2","sp3","sp4"))) %>% as.data.frame()
counts <- c(10,14,32)
n = 3

#the function's parameters:
  #taxa = a data frame of standardized taxon frequencies (sum=1), 
    #with colnames being taxons,
    #with rownames being sampling sites
  #counts = a vector of counted individuals per site
  #n = number of iterations

bootstrap_abuns <- function(taxa, counts, n) {
  #an array to store the bootstrapped abundances for each of the n iterations:
  bs_array <- array(dim = c( nrow(taxa), ncol(taxa), n),
                    dimnames = list( c(rownames(taxa)),
                                     c(colnames(taxa)),
                                     seq(1:n)) )#taxa has to be a df
  #this loop is for the n iterations
  for (i in 1:n){
    #this loop re-samples a community for each site:
    for (j in 1:nrow(taxa)){
      
      draws <- sample(colnames(taxa), size = counts[j], replace = TRUE, prob = taxa[j,]) %>%
        as.factor() %>%
        summary() # a vector with the number of draws for each species
      site_freq <- draws[colnames(taxa)] # bootsprapped composition at site j
      site_freq <- site_freq %>% replace(is.na(.), 0) #the NA's mean that we didn't sample any individual 
      
      bs_array[j, , i] <- site_freq
      
    }
    print(i) #to see progress
  }
  #return the filled array
  return(bs_array)
}


#testing it on dBEF_nem
taxa <- dBEF_nem[c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))]
counts <- dBEF_nem$total_nematodes

dBEF_bs <- bootstrap_abuns(taxa, counts, n=100)

dBEF_bs[1:3,2:6,1:5]

#save the result:
save(dBEF_bs, file = "./dBEF_bootstrapped.RData")

#this is a test comment to resolve error: git did not send all necessary objects




