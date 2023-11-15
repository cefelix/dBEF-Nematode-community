####Hyptheses finding####
library(ggplot2)

lgnd_default <- legend

####abundances####
Pl_div <- c(1,2,4,8,16)

Sh_Ph <- log(Pl_div)*0.9+0.3 
Sh_nP <- log(Pl_div)*0.7+0.2
nS_nP <- log(Pl_div)*0.5+0.1

Sh_Ph2 <- log(Pl_div)*0.8+0.3 
Sh_nP2 <- log(Pl_div)*0.6+0.2
nS_nP2 <- log(Pl_div)*0.4+0.1

plot(log(Pl_div), Sh_Ph, col="blue", pch=16, type="l",
     ylim = c(0,max(Sh_Ph)+0.5),
     xlim =c(0, log(17)),
     ylab = "Nematode abundance",
     xlab = "Plant diversity",
     xaxt= 'n')
axis(1, at=log(Pl_div), labels = c(1,2,4,8,16))
lines(log(Pl_div), Sh_nP, col="red", pch=16)
lines(log(Pl_div), nS_nP, col="green", pch=16)
legend("topleft", legend=c("+SH +PH","+SH -PH", "-SH -PH"),
       col=c("blue", "red", "green"), pch=16
       )
#hypotheses/abun_01.jpg

#now the higher trophic levels
lines(log(Pl_div), Sh_Ph2, lty="dashed", col="blue")
lines(log(Pl_div), Sh_nP2, lty="dashed", col="red")
lines(log(Pl_div), nS_nP2, lty="dashed", col="green")
legend("topright", legend=c("lower trophic level",
                               "higher trophic level"), 
       lty=c(1,2))
#hypotheses/abun_02.jpg

#Expectation:
#higher plant diversity --> more abundance of nematodes overall, 
  #as well as in different trophic group
#mechanistic pathways: 


####taxon richness####

#hypothesis 2: The mature community consists of 2 parts:
  #1) a non-random subset of the initial (naive) community
  #2) a turnover part, distinct from the initial community

#hypothesis 2a: migrated taxa < extinct naive part

Pl_div <- c(1,2,4,8,16)

Sh_Ph <- log(Pl_div)*0.6+0.6 
Sh_nP <- log(Pl_div)*0.5+0.5
nS_nP <- log(Pl_div)*0.8+0.7

plot(log(Pl_div), Sh_Ph, col="blue", pch=16, type="l",
     ylim = c(0, 3.3),
     xlim =c(0, log(17)),
     ylab = "Number of nematode taxa",
     xlab = "Plant diversity",
     xaxt= 'n')
axis(1, at=log(Pl_div), labels = c(1,2,4,8,16))
lines(log(Pl_div), Sh_nP, col="red", pch=16)
lines(log(Pl_div), nS_nP, col="green", pch=16)
legend("topleft", legend=c("+SH +PH","+SH -PH", "-SH -PH"),
       col=c("blue", "red", "green"), pch=16
)

#Explanation:
#green line (-SH -PH): initial community in the soil from arable field is very diverse, 
#but through 

#hypothesis 2b: migrated taxa > extinct naive part
Sh_Ph <- log(Pl_div)*0.6+0.6 
Sh_nP <- log(Pl_div)*0.5+0.5
nS_nP <- log(Pl_div)*0.4+0.4

plot(log(Pl_div), Sh_Ph, col="blue", pch=16, type="l",
     ylim = c(0, 3.3),
     xlim =c(0, log(17)),
     ylab = "Number of nematode taxa",
     xlab = "Plant diversity",
     xaxt= 'n')
axis(1, at=log(Pl_div), labels = c(1,2,4,8,16))
lines(log(Pl_div), Sh_nP, col="red", pch=16)
lines(log(Pl_div), nS_nP, col="green", pch=16)
legend("topleft", legend=c("+SH +PH","+SH -PH", "-SH -PH"),
       col=c("blue", "red", "green"), pch=16
)


#Explanaition: naive community is a frequently disturbed community 
#thus we expecta higher number of migrated new taxons than extinct naive taxons

####taxon richness in different trophic groups along the cp scale####

#primary/secondary consumers:
Sh_Ph <- log(Pl_div)*0.6+0.5 
Sh_nP <- log(Pl_div)*0.5+0.4
nS_nP <- log(Pl_div)*0.8+0.8

plot(log(Pl_div), Sh_Ph, col="blue", pch=16, type="l",
     ylim = c(0,max(nS_nP)+0.5),
     xlim =c(0, log(17)),
     ylab = "Number of nematode taxa",
     xlab = "Plant diversity",
     xaxt= 'n')
axis(1, at=log(Pl_div), labels = c(1,2,4,8,16))
lines(log(Pl_div), Sh_nP, col="red", pch=16)
lines(log(Pl_div), nS_nP, col="green", pch=16)
legend("topleft", legend=c("+SH +PH","+SH -PH", "-SH -PH"),
       col=c("blue", "red", "green"), pch=16
)



#### composition####

nS_nP <- c(40, 25, 20, 7.5, 7.5)
Sh_nP <- c(25, 37, 23, 7.5, 7.5)
Sh_Ph <- c(19, 40, 26, 7.5, 7.5) 



composition <- cbind(nS_nP, Sh_nP, Sh_Ph)
  rownames(composition) <- c("Pl", "Ba", "Fu", "Ca", "Om")
  colnames(composition) <- c("-SH -PH", "+SH -PH", "+SH +PH")
  colSums(composition)
  
barplot(composition,
        col = c("lightgreen", "orange", "navajowhite3", "violet","pink"),
        ylab = "relative abundance",
        legend = rownames(composition))

  
  

####Shannon 

#H4a Shannon for entire communtiy

Sh_Ph <- log(Pl_div)*0.6+0.6 
Sh_nP <- log(Pl_div)*0.5+0.55
nS_nP <- log(Pl_div)*0.3+0.4

plot(log(Pl_div), Sh_Ph, col="blue", pch=16, type="l",
     ylim = c(0,2.5),
     xlim =c(0, log(17)),
     ylab = "Shannon index H'",
     xlab = "Plant diversity",
     xaxt= 'n')
axis(1, at=log(Pl_div), labels = c(1,2,4,8,16))
lines(log(Pl_div), Sh_nP, col="red", pch=16)
lines(log(Pl_div), nS_nP, col="green", pch=16)
title("Entire Nematode community")
legend("topleft", legend=c("+SH +PH","+SH -PH", "-SH -PH"),
       col=c("blue", "red", "green"), pch=16
)

#H4a2 Shannon for Fu/Ba

Sh_Ph <- log(Pl_div)*0.6+0.6 
Sh_nP <- log(Pl_div)*0.5+0.55
nS_nP <- log(Pl_div)*0.3+0.4

plot(log(Pl_div), Sh_Ph, col="blue", pch=16, type="l",
     ylim = c(0,2.5),
     xlim =c(0, log(17)),
     ylab = "Shannon index H'",
     xlab = "Plant diversity",
     xaxt= 'n')
axis(1, at=log(Pl_div), labels = c(1,2,4,8,16))
lines(log(Pl_div), Sh_nP, col="red", pch=16)
lines(log(Pl_div), nS_nP, col="green", pch=16)
title("Fungi- / Bacterivores")
legend("topleft", legend=c("+SH +PH","+SH -PH", "-SH -PH"),
       col=c("blue", "red", "green"), pch=16
)


#H4b Shannon for only Plant feeders

Sh_Ph <- log(Pl_div)*0.5+0.7 
Sh_nP <- log(Pl_div)*0.4+0.65
nS_nP <- log(Pl_div)*0.8+0.4

plot(log(Pl_div), Sh_Ph, col="blue", pch=16, type="l",
     ylim = c(0,2.5),
     xlim =c(0, log(17)),
     ylab = "Shannon index H'",
     xlab = "Plant diversity",
     xaxt= 'n')
axis(1, at=log(Pl_div), labels = c(1,2,4,8,16))
lines(log(Pl_div), Sh_nP, col="red", pch=16)
lines(log(Pl_div), nS_nP, col="green", pch=16)
title("Only Plant feeders")
legend("topleft", legend=c("+SH +PH","+SH -PH", "-SH -PH"),
       col=c("blue", "red", "green"), pch=16
)


