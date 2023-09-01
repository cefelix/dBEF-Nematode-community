#structuring the vogel data

library(tidyr)

raw <- read.csv("./322_2_data.csv") #as in JEXIS
  str(raw)

data <- raw %>% 
  group_by(genus) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = genus, values_from = number) %>%
  select(-row)

#now we need soil FW and DW in order to work with the densities



