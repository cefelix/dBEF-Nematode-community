#disentangling some errors which might occur when using any index calculation function on the maRcel::query_nemaplex output
  #run https://github.com/cefelix/dBEF-Nematode-community/blob/main/wrangling/amynang-JenaXP_SP6_2021-wrangling-nematodes.R before you continue!
  #(this provides you the necessary data!)

library(devtools)
install_github("amynang/marcel")

####Error: only NaN's in output####
data.nplx <- read.csv("./wrangling/nemaplex.csv")
Enrichment(data.4 ,nemaplex = data.nplx) #produces NaN


####Error solution: rename nemaplex rows####
#maRcel::Enrichment() relies on match() function to match rownames of nemaplex with colnames of df

rownames(data.nplx) #these are just numbers
data.nplx$X #these are the taxons we want to have as rownames

#Thus, set the appropriate rownames
rownames(data.nplx) <- data.nplx$X
rownames(data.nplx) #now we got it




####Error solution details####
#lets see what maRcel::Enrichment() does:
  #function (df, nemaplex) 
  #  {
  #    Ba1 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 
  #                  1 & nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 
  #                  3)
  #   [...]
  #   }

#initialize df and nemaplex
df <- data.4
nemaplex <- read.csv("./wrangling/nemaplex.csv")

#start from the very inside of the function:
colnames(df) #that's legitimate
rownames(nemaplex) #that's meh
match(colnames(df), rownames(nemaplex)) #that returns NA's





