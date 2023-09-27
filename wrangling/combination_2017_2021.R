#combining Vogel's 2017 and Amyntas' 2021 data
#to be safe, re-run the wrangling files for 2017/2021's data previous to this!
#
#
# in dBEF_nem: taxa densities always have to be the last columns of dBEF_nem to 
# make this code work!

#libraries
library(dplyr)
library(tidyr)
library(maRcel)
library(vegan)

amyntas2021 <- read.csv("./wrangling/Amyntas2021.csv", row.names = 1)
vogel2017 <- read.csv("./wrangling/Vogel2017.csv", row.names = 1)
nemaplex <- read.csv("wrangling/nemaplex.csv", row.names = 1)

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

####fixing NA's in nemaplex####
#
#check on NA's in nemaplex:
nemaplex %>% filter(if_any(everything(), is.na))

  #merge Macroposthonia and Criconemoides as one genus
  dBEF_nem$Criconemoides <- dBEF_nem$Criconemoides+dBEF_nem$Macroposthonia
  #and drop $Macroposhonia
  dBEF_nem <- dBEF_nem %>%
    subset(., select = -c(Macroposthonia))
  
  #Nordiidae can either be herbivores or Omnivores (Bongers 1999: The Maturity index, table 2)
  #but are always cp-4:
  nemaplex[rownames(nemaplex)=="Nordiidae",]$cp_value <- 4
  #assigning them to Omnivores:
  nemaplex[rownames(nemaplex)=="Nordiidae",]$feeding <- 8
  #to see if results differ, repeat analysis with Nordiidae as herbivores
  #nemaplex[rownames(nemaplex)=="Nordiidae",]$feeding <- 1
  
  #how to treat Rhabditidae dauer larvae?
  #keeping them as own taxon for now, as the ecological implication presence of dauerlarvae differ
  #from those of non-dauerlarvae rhabditidae being present: 
    #if dauerlarvae present, a microbial bloom might have occured previuosly;
    #if rhabditidae in normal form present, this might indicate an occuring microbial bloom
    dBEF_nem$Rhabditidae.dauer.larvae %>%
      sum() #9118
    dBEF_nem$Mesorhabditis %>%
      sum() #31 (Ba-cp1)
    dBEF_nem$Protorhabditis %>%
      sum() #67 (Ba-cp1)
    dBEF_nem$Rhabditis %>%
      sum() #9170 (Ba-cp1) -> but most likely Rhabditis
    (dBEF_nem$Rhabditidae.dauer.larvae - dBEF_nem$Rhabditis) %>% 
      summary() #--> maybe check whether in dry week there are more dauerlarvae
    
  
  #This will include Rhabditidae in the calculation of all trophic/cp guild density calculations:
  nemaplex[rownames(nemaplex)=="Rhabditidae.dauer.larvae",]$cp_value <- 1
  nemaplex[rownames(nemaplex)=="Rhabditidae.dauer.larvae",]$feeding <- 3


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

#functional groups numbers:
main.plot <- read.csv("./main_plot_information.csv")
dBEF_nem <- dBEF_nem %>% 
  mutate(.after = sowndiv,
         func.group = as.numeric(main.plot$func.group[match(.$plot, main.plot$plotcode)]))

#presence/absence of functional groups:
main.plot$numgrass[main.plot$numgrass > 0] <- 1
  main.plot$numleg[main.plot$numleg > 0]     <- 1
  main.plot$numsherb[main.plot$numsherb > 0] <- 1
  main.plot$numtherb[main.plot$numtherb > 0] <- 1
dBEF_nem <- dBEF_nem %>%
  mutate(.after = SH,
         grasses = as.character(main.plot$numgrass[match(.$plot, main.plot$plotcode)]),
         legumes = as.character(main.plot$numleg[match(.$plot, main.plot$plotcode)]),
         s.herbs = as.character(main.plot$numsherb[match(.$plot, main.plot$plotcode)]),
         t.herbs = as.character(main.plot$numtherb[match(.$plot, main.plot$plotcode)]),
    
  )




#### response variables####

#total abundance (individuals per 100g DW)
#column indices of taxa:
taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))

dBEF_nem <- dBEF_nem %>%
  mutate(ind_per100g = rowSums(dBEF_nem[taxa]), .before = "Acrobeles")


#densities of trophic guilds:

  #a function to filter out trophic guilds
  trophic_guilds <- function(df, nemaplex) {
    #plant feeders:
    Pl = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 1)
    abun.Pl = df %>% mutate(herbivores = rowSums(.[Pl]), .keep = "none")
    
    #fungal feeders:
    Fu = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2)
    abun.Fu = df %>% mutate(fungivores = rowSums(.[Fu]), .keep = "none")
    
    #bacterial feeders:
    Ba = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3)
    abun.Ba = df %>% mutate(bacterivores = rowSums(.[Ba]), .keep = "none")
    
    #Predators:
    Pr = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 5)
    abun.Pr = df %>% mutate(predators = rowSums(.[Pr]), .keep = "none")
    
    
    #Omnivores:
    Om = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 8)
    abun.Om = df %>% mutate(omnivores = rowSums(.[Om]), .keep = "none")
    
    densities_tr = cbind(abun.Pl, abun.Fu, abun.Ba, abun.Pr, abun.Om)
    return(densities_tr)
  }
  
  #apply trophic_guilds on dBEF_nem
  taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))
  dBEF_nem <- dBEF_nem %>% 
    mutate(trophic_guilds(.[taxa], nemaplex), .before = "Acrobeles") %>%  
    rename("Pl_per100g" = "herbivores",
           "Fu_per100g" = "fungivores" ,
           "Ba_per100g" = "bacterivores",
           "Pr_per100g" = "predators",
           "Om_per100g" = "omnivores")

  
#densities by cp-value:
  
  #determine columns where taxonomic data is:
  taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))
  
  #adding 5 columns, density of cp-groups in ind/100g, added before the taxa
  dBEF_nem <- dBEF_nem %>%
    mutate(maRcel::C_P(.[taxa],nemaplex)*.$ind_per100g, .before = "Acrobeles") %>%
    rename("cp1_per100g" = "cp1",
           "cp2_per100g" = "cp2",
           "cp3_per100g" = "cp3",
           "cp4_per100g" = "cp4",
           "cp5_per100g" = "cp5")
  
  
#combination of trophic guild and cp-value:
  
  # a function filtering the density of each trophic guild and the cp-value:
  trophic_guildsCP <- function(df, nemaplex) {
    
    #Plant feeders: Pl-2
    Pl2 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 1 &
                  nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 2)
    abun.Pl2 = df %>% mutate(Pl2 = rowSums(.[Pl2]), .keep = "none")
      #Pl-3
      Pl3 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 1 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 3)
      abun.Pl3 = df %>% mutate(Pl3 = rowSums(.[Pl3]), .keep = "none")
      #Pl-4
      Pl4 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 1 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 4)
      abun.Pl4 = df %>% mutate(Pl4 = rowSums(.[Pl4]), .keep = "none")
      #Pl-5
      Pl5 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 1 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 5)
      abun.Pl5 = df %>% mutate(Pl5 = rowSums(.[Pl5]), .keep = "none")
      
    #fungal feeders: Fu-2
    Fu2 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2 &
                  nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 2)
    abun.Fu2 = df %>% mutate(Fu2 = rowSums(.[Fu2]), .keep = "none")
      #Fu-3
      Fu3 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 3)
      abun.Fu3 = df %>% mutate(Fu3 = rowSums(.[Fu3]), .keep = "none")
      #Fu-4
      Fu4 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 4)
      abun.Fu4 = df %>% mutate(Fu4 = rowSums(.[Fu4]), .keep = "none")
      #Fu-5
      Fu5 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 5)
      abun.Fu5 = df %>% mutate(Fu5 = rowSums(.[Fu5]), .keep = "none")

    #bacterial feeders: Ba-1
    Ba1 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3 &
                  nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 1)
    abun.Ba1 = df %>% mutate(Ba1 = rowSums(.[Ba1]), .keep = "none")
      #Ba-2
      Ba2 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3 &
                      nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 2)
        abun.Ba2 = df %>% mutate(Ba2 = rowSums(.[Ba2]), .keep = "none")
      #Ba-3
      Ba3 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 3)
      abun.Ba3 = df %>% mutate(Ba3 = rowSums(.[Ba3]), .keep = "none")
      #Ba-4
      Ba4 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 4)
      abun.Ba4 = df %>% mutate(Ba4 = rowSums(.[Ba4]), .keep = "none")
      #Ba-5
      Ba5 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 5)
      abun.Ba5 = df %>% mutate(Ba5 = rowSums(.[Ba5]), .keep = "none")
      
    #Predators: Pr-3
    Pr3 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 5 &
                  nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 3)
    abun.Pr3 = df %>% mutate(Pr3 = rowSums(.[Pr3]), .keep = "none")
      #Pr-4
      Pr4 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 5 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 4)
      abun.Pr4 = df %>% mutate(Pr4 = rowSums(.[Pr4]), .keep = "none")
      #Pr-5
      Pr5 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 5 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 5)
      abun.Pr5 = df %>% mutate(Pr5 = rowSums(.[Pr5]), .keep = "none")
    
    #Omnivores: Om-4
    Om4 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 8 &
                 nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 4)
    abun.Om4 = df %>% mutate(Om4 = rowSums(.[Om4]), .keep = "none")
      #Om-5
      Om5 = which(nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 8 &
                    nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 5)
      abun.Om5 = df %>% mutate(Om5 = rowSums(.[Om5]), .keep = "none")
      
    
    densities_tr_cp = cbind(abun.Pl2, abun.Pl3, abun.Pl4, abun.Pl5,
                            abun.Fu2, abun.Fu3, abun.Fu4, abun.Fu5,
                            abun.Ba1, abun.Ba2, abun.Ba3, abun.Ba4, abun.Ba5,
                            abun.Pr3, abun.Pr4, abun.Pr5,
                            abun.Om4, abun.Om5)
    return(densities_tr_cp)
  }
  
  #apply trophic_guildsCP
  taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))
  dBEF_nem <- dBEF_nem %>%
    mutate(trophic_guildsCP(.[taxa], nemaplex), .before = "Acrobeles")

    
####indices####
  #shannon's diversity H':
  taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))
  dBEF_nem <- dBEF_nem %>%
    mutate(Shannon_H = diversity(.[taxa], index = "shannon"), .before = "Acrobeles" ) 
    
  #species number:
  taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))
  presence <- dBEF_nem[taxa] > 0
  dBEF_nem <- dBEF_nem %>%
    mutate(Hill_q0 = rowSums(presence == TRUE), .before = "Acrobeles")
    rm(presence)

  #effective species number (exp(H'), Hill number with q=1):
  taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))
  dBEF_nem <- dBEF_nem %>%
    mutate(Hill_q1 = exp(.$Shannon_H), .before = "Acrobeles")

  #nematode specific indices from maRcel
  taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))
  indices_nematodes <- cbind(Enrichment(dBEF_nem[taxa], nemaplex),
                            Structure(dBEF_nem[taxa], nemaplex),
                            Maturity(dBEF_nem[taxa], nemaplex),
                            Channel(dBEF_nem[taxa], nemaplex) ) 
  dBEF_nem <- dBEF_nem %>% 
    mutate(., indices_nematodes, .before = "Acrobeles")
    rm(indices_nematodes)
  
  #nematode specific indices as in Eisenhauer 2011
  #Fu/(Ba+Fu):
    taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))
    dBEF_nem <- dBEF_nem %>%
      mutate("Fu/(Fu+Ba)" = .$Fu_per100g / (.$Fu_per100g + .$Ba_per100g), .before = "Acrobeles")
  #(Fu+Ba)/Pl:
    taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))
    dBEF_nem <- dBEF_nem %>%
      mutate("(Fu+Ba)/Pl" = (.$Fu_per100g + .$Ba_per100g) / .$Pl_per100g, .before = "Acrobeles")
  #Pr/Pl:
    taxa <- c(grep("Acrobeles", colnames(dBEF_nem)):ncol(dBEF_nem))
    dBEF_nem <- dBEF_nem %>%
    mutate("Pr/Pl" = .$Pr_per100g / .$Pl_per100g, .before = "Acrobeles")  

  

####save the whole file as .csv####
write.csv(dBEF_nem, "./dBEF_nem.csv")
    
    





