setwd("/workdir/gmv23/peppers/pheno/asreml")

library(asreml)
library(reshape2)

#Read and code data
pep <- read.csv("data/Ratings_filt.csv")
pep$Rep <- as.factor(pep$Rep)
pep$Block <- as.factor(pep$Block)
pep$Tray <- as.factor(pep$Tray)
pep$Column <- as.factor(pep$Column)

#Make dummy variables
check_rows <- which(pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"))
pep$New <- rep(2, nrow(pep))
pep$New[check_rows] <- 1
pep$New <- as.factor(pep$New)

#For now make pep.sub Aristotle
pep.sub <- pep[pep$Isolate == "Aristotle"]
