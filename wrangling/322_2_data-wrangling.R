#structuring the vogel data

#libraries:
library(tidyr)
library(dplyr)
library(readxl)
library(stringr)
library(rBExIS)
library(vegan)
#functions:
`%not_in%` <- purrr::negate(`%in%`)

#load plot info from jexis:
bexis.options(base_url = "https://jexis.uni-jena.de")
#main.plot = bexis.get.dataset_by(id = 90)
#as the previous line throws an error, just load the file manually:
main.plot <- read.csv("./main_plot_information.csv")

####community composition data as in dataset 322 in JEXIS####
raw <- read.csv("./322_2_data.csv") #as in JEXIS
  str(raw)

#change from long format to wide format:  
composition_vog <- raw %>% 
  group_by(genus) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = genus, values_from = number) %>%
  select(-row) %>%
  arrange(plotcode)

#change misleading column name, as this is not the total abundance of the Sample:
colnames(composition_vog)[3] <- "total_identified"

#create Jexis conform Sample names:
composition_vog <- mutate(composition_vog, Sample = paste(plotcode, treatment, sep = "")) %>%
  arrange(Sample) %>%
  relocate(Sample, .before = plotcode)



####density data from nematodes_Jun2017.csv ####
# ---
#note: The nematodes_Jun2017.csv file has some flaws as:
#   1   Two columns with the supposedly same measure (nematodes per g DW), called
#       `indiv.per.mass..gdw.1` and `tot.abundance.per.soildw`
#   2   intransparent calculation of these the above mentioned columns
#   3   deviations of the two reaching more than 1 individual per g DW in 14 Samples
#
#THUS: Use the density data we calculated in the section below:

raw <- read.csv("./nematodes_Jun2017.csv")

#select the density columns:
densities_vog <- subset(raw, select = c(plotcode,
                                        treatment,
                                       `indiv.per.mass..gdw.1.`,
                                       `tot.abundance.per.soildw`))
#create Jexis conform Sample names:
densities_vog <- mutate(densities_vog, Sample = paste(plotcode, treatment, sep = "D")) %>% 
  arrange(Sample) %>%
  relocate(Sample, .before = plotcode)

densities_vog %>%
  filter(abs(indiv.per.mass..gdw.1. - tot.abundance.per.soildw) > 1) #14 Samples -.-
    #`indiv.per.mass..gdw.1.` is matching my calculations from below



####DW, abundances and densities as in `complete dataset nematodes.xlsx`####
# ---
#note: use the densities calculated in this section

raw <- read_xlsx("./complete dataset nematodes.xlsx", sheet = "Tabelle1")

#keep only plotcode, DW, total abundance and density; get correct column names:
DW_vog <- raw[4:nrow(raw), 1:5] 
  colnames(DW_vog) <- c(DW_vog[1, 1:3], "total_nematodes", "nem_gDW")
  DW_vog <- DW_vog[-1,] %>% #droping first row, it's the rownames
    filter(total_nematodes != -9999) #probably some weird excel error
  

#Jexis conform Sample names:  
DW_vog <- mutate(DW_vog, Sample =paste(plotcode, treatment, sep="D")) %>%
  arrange(Sample) %>%
  relocate(Sample, .before = plotcode) 

#recalculate nematodes per 100g DW:
DW_vog$total_nematodes <- as.numeric(DW_vog$total_nematodes)
DW_vog$`soil (gdw)` <- as.numeric(DW_vog$`soil (gdw)`)
DW_vog$nem_gDW <- as.numeric(DW_vog$nem_gDW)
DW_vog$nem_100gDW <- DW_vog$total_nematodes / (DW_vog$`soil (gdw)`/100)

#add gravimetric soil water content (WARNING: as we lack the original measuring sheet, 
  #we assume that each sample consisted of 25g fresh soil initially!)
#Thus: ("mass of moist soil" − "mass of oven-dried soil")/"mass of oven-dried soil" * 100:
DW_vog <- DW_vog %>%
  mutate(soil.water.gravimetric = (25 - DW_vog$`soil (gdw)`)/DW_vog$`soil (gdw)` *100, 
         .before = total_nematodes)


#check the differences between my calculated densities and Vogel's provided ones:
str(DW_vog)
DW_vog %>%
  filter(abs(nem_gDW*100 - nem_100gDW) > 0.001) #B4A0D2 has a -9999 error



####vogel2017a - merge density and composition data####

#define rows which have to be dropped due to inconsistent data:
  #following Samples are missing (Vogel script, line 465 to 470, additionally B2A03D2):
  composition_missing <- c("B2A04D2", "B1A04D2", "B2A01D2","B2A16D2", "B3A03D2", "B2A03D2")
  #following Samples occurred twice (472-474):
  composition_duplicate <- c("B3A23D2", "B1A22D2", "B2A22D2")
  #following Samples are from plots which dont have the target plant richness (476-478):
  wrong_plot <-  c("B2A23D1", "B2A23D2", "B2A23D3", "B2A20D1", "B2A20D2", "B2A20D3")

#merge into one table:
vogel2017a <- DW_vog %>% full_join(composition_vog, by = join_by(Sample)) %>%
  filter(Sample %not_in% composition_missing) %>%
  filter(Sample %not_in% composition_duplicate) %>%
  filter(Sample %not_in% wrong_plot)

#remove/rename duplicate plotcode and treatment columns:
vogel2017a$plotcode.x == vogel2017a$plotcode.y
vogel2017a <- vogel2017a  %>% rename(c(plot = plotcode.x, treatment = treatment.x))

#rename Rhabditidae.dauer.larva into (...)larvae (congruent with data from 2021)
vogel2017a <- vogel2017a %>%
  mutate(Rhabditidae.dauer.larvae = Rhabditidae.dauer.larva, .after=Rhabditidae.dauer.larva) %>%
  subset(select = -c(Rhabditidae.dauer.larva))


#Split blockplot into block and plot:
vogel2017a <- vogel2017a %>% 
  mutate(block = str_remove(.$plot, pattern = "A.*"), .after = 1) %>%
  #add year:
  mutate(year = c("2017"), .after = treatment)
  
vogel2017a[1,1:7]

#remove columns plotcode.y, treatment.y, blockplot, nem_gDW:
vogel2017a <- vogel2017a %>% subset(select = -c(plotcode.y, treatment.y,  nem_gDW))

head(vogel2017a)

#add sowndiv, whose elements are filled from main.plot based on plot matching
vogel2017a <- vogel2017a %>% 
  mutate(.after = plot,
         sowndiv = as.character(main.plot$sowndiv[match(.$plot, main.plot$plotcode)]))

write.csv(vogel2017a, file = "./wrangling/abundances2017.csv")


####vogel2017d - calculate taxon abundances per 100g soil DW####

#transform abundances to frequencies:
vogel2017b <- vogel2017a
#select taxon's column indices:
taxons <- grep("Acrobeles", colnames(vogel2017a)):ncol(vogel2017a)
#transform abundances to frequencies
vogel2017b[,taxons] <- vegan::decostand(vogel2017a[,taxons], "total", 1)

  #check whether row sums are equal to 1:
  rowSums(vogel2017b[,taxons]) == 1 

#multiply frequencies with absolute numbers of extracted nematodes:
vogel2017c <- vogel2017b
taxons <- grep("Acrobeles", colnames(vogel2017b)):ncol(vogel2017b)
vogel2017c[,taxons] <- vogel2017b[,taxons]*vogel2017b$total_nematodes

  #check whether rowsums are equal to total_nematodes
  abs(rowSums(vogel2017c[,taxons]) - vogel2017c$total_nematodes) <0.0001
  
#multiply frequencies with nematode densities per 100g DW
vogel2017d <- vogel2017b
taxons <- grep("Acrobeles", colnames(vogel2017b)):ncol(vogel2017b)
vogel2017d[,taxons] <- vogel2017b[,taxons]*vogel2017b$nem_100gDW

  #check whether rowsums are equal to nem_100gDW
  abs(rowSums(vogel2017d[,taxons]) - vogel2017d$nem_100gDW) <0.0001
  
  #drop unnecessary columns
  vogel2017d <- vogel2017d %>%
    subset(select = -c(`soil (gdw)`, total_nematodes, total_identified, nem_100gDW))
  
  #this is the same format as data.4 from our 2021's data wrangling




####saving as .csv####
#better to re-run each time
getwd()
write.csv(vogel2017d, "./wrangling/Vogel2017.csv")



