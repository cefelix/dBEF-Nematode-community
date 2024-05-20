#bringing the data into the long format
library(dplyr)
library(tidyr)

dBEF_nem <- read.csv("./dBEF_nem.csv", row.names = 1)

dat21 <- dBEF_nem %>% filter(year == 2021)

taxa <- c(grep("Acrobeles", colnames(dat21)):ncol(dat21))
a <- dat21[,taxa]

counts <- read.csv("./wrangling/abundances2021.csv")
ident <- read.csv("./wrangling/Amyntas2021.csv")  
head(ident)
head(counts)

ident.long <- pivot_longer(ident, cols = all_of(taxa))
