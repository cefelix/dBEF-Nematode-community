#structuring the vogel data

#libraries:
library(tidyr)
library(dplyr)
library(readxl)
#functions:
`%not_in%` <- purrr::negate(`%in%`)


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

#change misleading column name, as this is not the total abundance of the sample:
colnames(composition_vog)[3] <- "total_identified"

#create Jexis conform sample names:
composition_vog <- mutate(composition_vog, sample = paste(plotcode, treatment, sep = "")) %>%
  arrange(sample) %>%
  relocate(sample, .before = plotcode)



####density data from nematodes_Jun2017.csv ####
# ---
#note: The nematodes_Jun2017.csv file has some flaws as:
#   1   Two columns with the supposedly same measure (nematodes per g DW), called
#       `indiv.per.mass..gdw.1` and `tot.abundance.per.soildw`
#   2   intransparent calculation of these the above mentioned columns
#   3   deviations of the two reaching more than 1 individual per g DW in 14 samples
#
#THUS: Use the density data we calculated in the section below:

raw <- read.csv("./nematodes_Jun2017.csv")

#select the density columns:
densities_vog <- subset(raw, select = c(plotcode,
                                        treatment,
                                       `indiv.per.mass..gdw.1.`,
                                       `tot.abundance.per.soildw`))
#create Jexis conform sample names:
densities_vog <- mutate(densities_vog, sample = paste(plotcode, treatment, sep = "D")) %>% 
  arrange(sample) %>%
  relocate(sample, .before = plotcode)

densities_vog %>%
  filter(abs(indiv.per.mass..gdw.1. - tot.abundance.per.soildw) > 1) #14 samples -.-
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
  

#Jexis conform sample names:  
DW_vog <- mutate(DW_vog, sample =paste(plotcode, treatment, sep="D")) %>%
  arrange(sample) %>%
  relocate(sample, .before = plotcode) 

#recalculate nematodes per 100g DW:
DW_vog$total_nematodes <- as.numeric(DW_vog$total_nematodes)
DW_vog$`soil (gdw)` <- as.numeric(DW_vog$`soil (gdw)`)
DW_vog$nem_gDW <- as.numeric(DW_vog$nem_gDW)
DW_vog$nem_100gDW <- DW_vog$total_nematodes / (DW_vog$`soil (gdw)`/100)

#check the differences between my calculated densities and Vogel's provided ones:
str(DW_vog)
DW_vog %>%
  filter(abs(nem_gDW*100 - nem_100gDW) > 0.001) #B4A0D2 has a -9999 error








####merge density and composition data####

#define rows which have to be dropped due to inconsistent data:
  #following samples are missing (Vogel script, line 465 to 470):
  composition_missing <- c("B2A04D2", "B1A04D2", "B2A01D2","B2A16D2", "B3A03D2")
  #following samples occurred twice (472-474):
  composition_duplicate <- c("B3A23D2", "B1A22D2", "B2A22D2")
  #following samples are from plots which dont have the target plant richness (476-478):
  wrong_plot <-  c("B2A23D1", "B2A23D2", "B2A23D3", "B2A20D1", "B2A20D2", "B2A20D3")

#merge into one table
genusDW_vog <- DW_vog %>% full_join(composition_vog, by = join_by(sample)) %>%
  #filter(sample %not_in% composition_missing) %>%
  #filter(sample %not_in% composition_duplicate) %>%
  filter(sample %not_in% wrong_plot)



####Vogel excluded some plots where data was missing or occured twice####

    #following treatment D2 samples are missing (Vogel script line 465 to 470):
    a <- c("B2A04D2", "B1A04D2", "B2A01D2","B2A16D2", "B3A03D2")
      #B2A04 D2: composition missing, density available
      #B1A04 D2: composition missing, density available
      #B2A01 D2: composition missing, density available
    
    #the following plots occured twice and thus were removed (472 line to 474):
    c(B3A23, B1A22, B2A22)
      #B3A23 D2: twice in community composition with differing values -> DROP
      #B1A22 D2: twice in community composition with differing values -> DROP
      #B2A22 D2: twice in community composition, twice in densities -> DROP
        #(-> could be kept for some analysis)
    
    #the following plots should not have been sampled (line 476 to 478):
    c(B2A23, B2A20)










