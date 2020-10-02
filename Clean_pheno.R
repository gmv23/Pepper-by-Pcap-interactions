setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/pheno/")
source("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/scripts/Functions_for_disease_analysis.R")

library(reshape2)

#Import data
pep <- read.csv("data/Ratings_all.csv", na.strings = "NA") 
pep$Rep <- as.factor(pep$Rep)
pep$Block <- as.factor(pep$Block)
pep$Tray <- as.factor(pep$Tray)
pep$Column <- as.factor(pep$Column)
colnames(pep)[colnames(pep) == "Isolate.Name"] <- "Isolate"
#Get rid of rows with NAs (trays that weren't inoculated bc of inoculum problems)
na_rows <- which(apply(pep[,8:13], 1, function(x) any(is.na(x))))
pep <- pep[-na_rows,]
pep$Isolate <- droplevels(pep$Isolate)

#Look at number of times each pepper and isolate appears
table(table(pep$Pepper))
table(pep$Pepper)
table(table(pep$Isolate))
table(pep$Isolate)

#Make some variables used in functions to check for errors 
pep$Plot_no <- 6
dpi <- as.integer(gsub("Dpi", "", colnames(pep[,8:13])))

#Check errors
check_errors(pep, 8:13, 14, F, F)
pep <- check_errors(pep, 8:13, 14, T, T) #Fix handful of observations where disease decreases over time

#Calculate AUDPC
pep$audpc <- apply(pep[,8:13], 1, my_audpc, time_points = dpi)

#Get rid of water blank and nonused columns
any(pep$audpc[pep$Isolate == "H2O"] != 0)
pep <- pep[pep$Isolate != "H2O",]
pep$Isolate <- droplevels(pep$Isolate)
pep$Isolate.CC <- NULL

#What percentage of isolates have >0 virulence for each pepper?
interact_means.long <- aggregate(audpc~Isolate*Pepper,data=pep, FUN=mean)
interact_means.wide <- dcast(interact_means.long, Isolate ~ Pepper, value.var = 'audpc')
percent_virulent <- apply(interact_means.wide[,-1],2,function(x) sum(x>0, na.rm=T)/sum(!is.na(x)))

#Filter out non-informative peppers
#Where to draw the line? 10% pathogenic? 20%? 30%
#Lets get rid of peppers where fewer than 15% of the isolates can cause disease
informative_peppers <- names(which(percent_virulent > 0.15))
length(informative_peppers)
pep.filt <- pep[pep$Pepper %in% informative_peppers,]

#Pull out non-pathogenic isolates --- remove them from data for further analysis
pathogenic <- apply(interact_means.wide[,-1], 1, function(x) all(x==0))
pathogenic_isolates <- interact_means.wide$Isolate[pathogenic]
length(pathogenic_isolates)
as.character(pathogenic_isolates)
pep.filt <- pep.filt[!pep.filt$Isolate %in% pathogenic_isolates,]

#Write clean pheno data
write.csv(pep, "data/Ratings_clean.csv", quote=F, row.names = F)
write.csv(pep.filt, "data/Ratings_filt.csv", quote=F, row.names = F)



