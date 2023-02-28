# library(devtools)
# install_github("BEXIS2/rBExIS", subdir = "rBExIS") 
library(httr)
library(jsonlite)
library(XML)
library(rBExIS)
bexis.options(base_url = "https://jexis.uni-jena.de")
#get data from bexis from dataset with id = xy
main.plot = bexis.get.dataset_by(id = 90)

data.4 %>% 
  # we change plot to be "blockplot" by splitting Sample at D and keeping the 1st half
  mutate(plot = str_split(.$Sample, "D", simplify = T)[,1]) %>% 
  # we add a sowndiv column whose elements are filled from main.plot
  # based on plot matching
  mutate(.after = plot,
         sowndiv = main.plot$sowndiv[match(.$plot, main.plot$plotcode)])
