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
  #a list to store the df's of bootstrapped abundances for each of the n iterations:
  bs_list <- list(as.data.frame(matrix(nrow = nrow(taxa),
                                       ncol = ncol(taxa),
                                       dimnames = list( c(rownames(taxa)),
                                                        c(colnames(taxa)) )
                                       )) )
  #this loop is for the n iterations
  for (i in 1:n){
    #initialize a data frame:
    df <- as.data.frame(matrix(nrow = nrow(taxa),
                               ncol = ncol(taxa)))
    #this loop re-samples a community for each site:
    for (j in 1:nrow(taxa)){
      
      draws <- sample(colnames(taxa), size = counts[j], replace = TRUE, prob = taxa[j,]) %>%
        as.factor() %>%
        summary() # a named vector with the number of draws for each species
      site_freq <- draws[colnames(taxa)] # bootsprapped composition at site j
      site_freq <- site_freq %>% replace(is.na(.), 0) #the NA's mean that we didn't sample any individual 
      
      df[j,] <- site_freq
      
    }
    rownames(df) <- rownames(taxa)
    colnames(df) <- colnames(taxa)
    bs_list[[i]] <- df
    print(i) #to see progress
  }
  #return the list of data frames
  return(bs_list)
}


#### standardizing dBEF_nem ####
  dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)
  #standardize taxon abundances:#
  stand <- grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem)
    dBEF_nem[stand] <- decostand(dBEF_nem[stand], method = "total", MARGIN = 1)
    #check it:
    rowSums(dBEF_nem[stand])
  
#### bootstrapping dBEF_nem ####  
  #subset dBEF_nem to feed into bootstrap_abuns():
    taxa <- dBEF_nem[c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))]
    counts <- dBEF_nem$total_nematodes
  #re-sample:    
  resampled <- bootstrap_abuns(taxa, counts, n=100)
  
  #Insert independent variables to each bootstrapped df:
  independent <- dBEF_nem[,1:15] #select independent vars
  resampled <- lapply(resampled, cbind, independent) #add them
  resampled <- lapply(resampled, relocate, c(81:95), .before = 1) #move to front
  
  #save the result:
  save(resampled, file = "./wrangling/taxons_bootstrapped_list.RData")



####repeat calculation of functional indices for each iteration####

dBEF_bs <- array(rep(dBEF_nem, 100))





#this is a test comment to resolve error: git did not send all necessary objects




