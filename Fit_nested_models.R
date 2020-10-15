setwd("/workdir/gmv23/peppers/models")
library(asreml)

###################################       Read and clean up data           ########################################

#Read data and code factors as factors
pep <- read.csv("data/Ratings_filt.csv")
pep$Rep <- as.factor(pep$Rep)
pep$Block <- as.factor(pep$Block)
pep$Tray <- as.factor(pep$Tray)
pep$Column <- as.factor(pep$Column)

#Read blues from full dataset
blues <- read.csv("tables/virulence_blues.csv")
blues$main <- NULL
rownames(blues) <- blues$Isolate
blues$Isolate <- NULL

#Create dummy variable for check rows
check_rows <- which(pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"))
pep$New <- rep(2, nrow(pep))
pep$New[check_rows] <- 1
pep$New <- as.factor(pep$New)

#####################    Loop through peppers and fit separate model for each pepper    #########################

peppers <- levels(pep$Pepper)

#Set up matrices to store heritabilities
H2s <- matrix(NA,ncol=2, nrow=length(peppers), dimnames=list(peppers, c("H2", "SE")))

for(i in 1:length(peppers)){

	#Subset data for specific pepper
	pepper <- peppers[i]
	pep.sub <- pep[pep$Pepper == pepper,]
	pep.sub$Isolate <- droplevels(pep.sub$Isolate)

	#Fit model
	mod <- asreml(fixed = audpc ~ 1 + Rep + at(New, 1):Isolate,
			 random = ~ at(New, 2):Isolate + Rep:Block,
			 residual = ~ idv(units),
			 data = pep.sub, na.action = na.method(x="include",y="include"), maxit=500)

	#Get harmonic mean of number reps
	r <- table(pep.sub$Isolate)
	n_harm <- length(r)/sum(1/r)
	
	#Get H2
	if(round(n_harm) > 1){ #only calculate H2 for peppers with more than one rep
		H2s[i,] <- unlist(vpredict(mod, H2 ~ V2/(V2+V3/n_harm)))
	}	
}

#Add median BLUE to data frame and write table
medians <- apply(blues, 2, median, na.rm=T)
medians <- round(medians,2)
H2s_round <- apply(H2s,2,round,digits=2)
summary_table <- data.frame("Isolate" = levels(pep$Pepper))
summary_table$median <- medians[match(summary_table$Isolate, names(medians))]
summary_table <- cbind(summary_table, H2s_round[match(summary_table$Isolate, rownames(H2s_round)),])
write.csv(summary_table, "tables/medians_and_H2s.csv", quote=F, row.names=F) 
