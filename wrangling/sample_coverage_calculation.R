#calculating sample coverage according to chao and jost (2012) eq.4a (in roswell 2021)
  #coverage depends on number of singletons, doubletons and individuals in a sample
#the calculated coverages are appended to dBEF_nem.csv and stored in dBEF_nem_Coverage.csv

library(readxl)
library(dplyr)

abundances2021 <- read_xlsx("./wrangling/abundances2021b.xlsx")
abundances2021[,2:ncol(abundances2021)] %>% rowSums(na.rm = TRUE) #all identified nematodes
  abundances2021$Sample <- gsub("_", "D" ,abundances2021$Sample) #correct Sample IDs

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)


#### a function to calculate coverage according to Chao and Jost (2012), eq. 4a ####

coverage <- function(f1, f2, n) {
  C = 1 - (f1/n) * ( ((n-1)*f1) / ((n-1)*f1 + 2*f2) ) #This is from Chao and Jost (2012), eq. 4a (Or: Roswell et al. (2021), eq. B1)
  return(C)
}  

#### for 2021 ####
  #presence/absence matrix of singletons
  singletons21 <- abundances2021 %>%
    mutate(across(Acrobeles:Xiphinema, function(x) {ifelse(x !=1 , 0, 1)} )) 
  #the number of singletons
  singletons21 <- singletons21 %>% 
    rowwise() %>%
    mutate(no.singletons = sum(c_across(Acrobeles:Xiphinema), na.rm = TRUE), .before = Acrobeles)
  #select only sample id and no. of singletons
  singletons21 <- select(singletons21, Sample, no.singletons)

  
  #presence/absence doubletons:
  doubletons21 <- abundances2021 %>%
    mutate(across(Acrobeles:Xiphinema, function(x) {ifelse(x !=2 , 0, 1)} )) #presence/absence matrix of doubletons
  #the number of doubletons
  doubletons21 <- doubletons21 %>% 
    rowwise() %>%
    mutate(no.doubletons = sum(c_across(Acrobeles:Xiphinema), na.rm = TRUE), .before = Acrobeles)
  #select only sample id and no. ofdoubletons
  doubletons21 <- select(doubletons21, Sample, no.doubletons)

  
  #the number of individuals in each sample: 
  no.ind <- select(dBEF_nem, Sample, total_nematodes, year)
  no.ind21 <- subset(no.ind, year == 2021)
  
  #a matrix with a column of singletons, doubletons, total_individuals:
  mat21 <- merge(singletons21, doubletons21, by="Sample")
  mat21 <- merge(mat21, no.ind21, by="Sample")
  
  #calculating coverage: 
  mat21$Cov <- rep(NA, 240)
  for(i in 1:nrow(mat21)){
    mat21[i,]$Cov <- coverage(mat21[i,]$no.singletons, mat21[i,]$no.doubletons, mat21[i,]$total_nematodes  )
  }
  
  plot(density(mat21$Cov))
    


#### for 2017 ####
abundances2017 <- read.csv("./wrangling/abundances2017.csv")
  #presence/absence matrix of singletons
  singletons17 <- abundances2017 %>%
    mutate(across(Acrobeles:Xiphinema, function(x) {ifelse(x !=1 , 0, 1)} )) 
  #the number of singletons
  singletons17 <- singletons17 %>% 
    rowwise() %>%
    mutate(no.singletons = sum(c_across(Acrobeles:Xiphinema), na.rm = TRUE), .before = Acrobeles)
  #select only sample id and no. of singletons
  singletons17 <- select(singletons17, Sample, no.singletons)
  
  
  #presence/absence doubletons:
  doubletons17 <- abundances2017 %>%
    mutate(across(Acrobeles:Xiphinema, function(x) {ifelse(x !=2 , 0, 1)} )) #presence/absence matrix of doubletons
  #the number of doubletons
  doubletons17 <- doubletons17 %>% 
    rowwise() %>%
    mutate(no.doubletons = sum(c_across(Acrobeles:Xiphinema), na.rm = TRUE), .before = Acrobeles)
  #select only sample id and no. ofdoubletons
  doubletons17 <- select(doubletons17, Sample, no.doubletons)
  
  
  #the number of individuals in each sample: 
  no.ind <- select(dBEF_nem, Sample, total_nematodes, year)
  no.ind17 <- subset(no.ind, year == 2017)
  
  #a matrix with a column of singletons, doubletons, total_individuals:
  mat17 <- merge(singletons17, doubletons17, by="Sample")
  mat17 <- merge(mat17, no.ind17, by="Sample")
  
  #calculating coverage: 
  mat17$Cov <- rep(NA, nrow(mat17))
  for(i in 1:nrow(mat17)){
    mat17[i,]$Cov <- coverage(mat17[i,]$no.singletons, mat17[i,]$no.doubletons, mat17[i,]$total_nematodes  )
  }

  plot(density(mat17$Cov))
  
  
#### attach mat21 and mat17 to dBEF_nem ####
  coverage.merged <- rbind(mat21, mat17)
  coverage.merged <- select(coverage.merged, Sample, no.singletons, no.doubletons, year, Cov)
  
  dBEF_nem <- merge(coverage.merged, dBEF_nem, by=c("Sample", "year"))
  dBEF_nem <- dBEF_nem %>%
    relocate(c(no.singletons, no.doubletons, Cov), .after = total_nematodes) %>%
    rename(Coverage = Cov)
  
  write.csv(dBEF_nem, file = "./dBEF_nem_Coverage.csv")
  
  
  