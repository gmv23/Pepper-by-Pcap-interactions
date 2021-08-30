setwd("~/Documents/Cornell/Pepper_Interactions/paper/pheno/")
source("~/Documents/Cornell/Pepper_Interactions/paper/Pepper-by-Pcap-interactions/Functions_for_disease_analysis.R")

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
pep$Isolate <- as.factor(pep$Isolate)
pep$Pepper <- as.factor(pep$Pepper)
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
checks <- c("CHECK1", "CHECK2", "CHECK3")

#How many observations excluding checks?
total_obs <- nrow(pep); total_obs
exp_obs <- sum(!pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3")); exp_obs
total_obs_all_time_points <- total_obs*6; total_obs_all_time_points
exp_obs_all_time_points <- exp_obs*6; exp_obs_all_time_points
unique_combinations <- length(unique((pep$Isolate:pep$Pepper)[!is.na(pep$audpc) & 
                                                              !as.character(pep$Isolate) %in% checks]))
unique_combinations

#Turn phenotype into cause disease or not cause disease
#Get sum of number reps an isolate caused any disease on a given pepper
disease_reps <- aggregate(audpc~Isolate*Pepper,data=pep, FUN=function(x) sum(x>0)/length(x))
colnames(disease_reps) <- c("Isolate", "Pepper", "rep_incidence")
disease_reps <- dcast(disease_reps, Isolate ~ Pepper, value.var = 'rep_incidence')
disease_reps <- disease_reps[-grep("CHECK", disease_reps$Isolate),]

#Pull out non-pathogenic isolates --- remove them from data for further analysis
non_pathogenic <- apply(disease_reps[,-1], 1, function(x) all(x==0, na.rm = T))
non_pathogenic_isolates <- disease_reps$Isolate[non_pathogenic]
length(non_pathogenic_isolates)
as.character(non_pathogenic_isolates)
pep.filt <- pep[!pep$Isolate %in% non_pathogenic_isolates,]
disease_reps <- disease_reps[!disease_reps$Isolate %in% non_pathogenic_isolates,]

#What percentage of isolates consistently cause disease (ie in both reps) on each pepper?
percent_virulent <- apply(disease_reps[,-1],2,function(x) sum(x==1, na.rm=T)/sum(!is.na(x)))
number_virulent <- apply(disease_reps[,-1],2,function(x) sum(x==1, na.rm=T))

#Plot percent virulent on each pepper
pdf("plots/percent_virulent_per_pepper.pdf")
old.par <- par(no.readonly = T)
par(mar=c(7,4,4,2))
percent_virulent <- sort(percent_virulent)
number_virulent <- sort(number_virulent)
barplot.data <- barplot(percent_virulent*100, plot=F)
barplot(percent_virulent*100,las=2,ylab="Isolates Pathogenic in All Reps (%)",ylim=c(0,100))
text(barplot.data, percent_virulent*100 + 5, labels=number_virulent)
abline(h=20, lty=2, col = 'black', lwd=2)
par(old.par)
dev.off()

#Lets get rid of peppers where fewer than 20% of the isolates can cause disease in both reps
cutoff <- 0.20
informative_peppers <- names(which(percent_virulent > cutoff))
length(informative_peppers)
pep.filt <- pep.filt[pep.filt$Pepper %in% informative_peppers,]
non_informative_peppers <- names(percent_virulent)[!(names(percent_virulent) %in% informative_peppers)]
disease_reps <- disease_reps[,!colnames(disease_reps) %in% non_informative_peppers]

#Write pheno data -- clean (i.e. errors fixed, no water blank, no NAs)
# and filtered --- no nonpathogenic isolates, no peppers where fewer than 20% of isolates consistently cause disease
write.csv(pep, "data/Ratings_clean.csv", quote=F, row.names = F)
write.csv(pep.filt, "data/Ratings_filt.csv", quote=F, row.names = F)
write.csv(disease_reps, "data/Disease_reps.csv", quote=F, row.names = F)



