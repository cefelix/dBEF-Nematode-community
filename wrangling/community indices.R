#calculating the community indices using maRcel

####required packages####
  library(devtools)
  install_github("amynang/marcel")
  #install.packages("RSelenium")
  #install_github("BEXIS2/rBExIS", subdir = "rBExIS") 

####loading libraries#### 
  library(maRcel)
  library(tidyverse)
  library(RSelenium)
  lsf.str("package:maRcel")
  library(httr)
  library(jsonlite)
  library(XML)
  library(rBExIS)


####data checking and re-arranging####
#get data by re-running code from: https://github.com/amynang/JenaXP_SP6_2021/blob/main/wrangling/nematodes.R

head(data.3) #nematode densities per soil sample
head(data.4) #nematode densities per 100g soil dry weight
#head(data.7) #nematode densities per 1 square meter, in 0-10 cm depth


bexis.options(base_url = "https://jexis.uni-jena.de")
#get data from bexis from dataset with id = xy
main.plot = bexis.get.dataset_by(id = 90)

data.4 <- data.4 %>% 
  # we change plot to be "blockplot" by splitting Sample at D and keeping the 1st half
  mutate(plot = str_split(.$Sample, "D", simplify = T)[,1]) %>% 
  # we add a sowndiv column whose elements are filled from main.plot
  # based on plot matching
  mutate(.after = plot,
         sowndiv = main.plot$sowndiv[match(.$plot, main.plot$plotcode)])

head(data.4) #now sown diversity is the 4th column



####maRcel testing####
#test_taxa <- c(colnames(data.4[6:8]))
#query_nemaplex(test_taxa) #Error in if (file.access(phantompath, 1) < 0) { : argument is of length zero
#system("taskkill /im java.exe /f", intern=FALSE, ignore.stdout=FALSE)

#query_nemaplex("Actinolaimidae") #Error in if (file.access(phantompath, 1) < 0) { : argument is of length zero
#system("taskkill /im java.exe /f", intern=FALSE, ignore.stdout=FALSE)



####indices calculation ####
#as maRcel produces an error atm, i will use the data which maRcel::query_nemaplex() should produce:
data.nplx <- read.csv("./wrangling/nemaplex.csv")
rownames(data.nplx) <- data.nplx$X  #this fixes an NaN error, when applying the index functions of maRcel
                                    #for details: https://github.com/cefelix/dBEF-Nematode-community/blob/main/wrangling/nemaplexCSV_error_solution.R

#CP proportions
data.4 %>%
  C_P(nemaplex = data.nplx) %>%
  summary()

#channel ratio
data.4 %>%
  Channel(nemaplex = data.nplx) %>% 
  summary() #thats a Median CI of 100, seems sketchy, should look into it

#maturity index
data.4 %>% 
  Maturity(nemaplex = data.nplx) %>%
  summary()

#enrichment index
data.4 %>% 
  Enrichment(nemaplex = data.nplx) %>%
  summary()

#structure index
data.4 %>%
  Structure(nemaplex = data.nplx) %>%
  summary()

####create table with all response variables####
#using the fancy function which does all indices at once:
data.4$sowndiv <- data.4$sowndiv %>%
  as.factor() #this prevents sowndiv being dropped by all.índices()
data.indices <- data.4 %>%
  all.indices(nemaplex = data.nplx) #creating a new df with all indices and bloc/plot/treatment info
data.indices$sowndiv <- data.indices$sowndiv %>%
  as.character() %>%
  as.numeric() #converting sowndiv back to numeric
head(data.indices)
str(data.indices)


