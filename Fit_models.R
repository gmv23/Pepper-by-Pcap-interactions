setwd("/workdir/gmv23/peppers/models")
library(asreml)

###################################       Read and clean up data           ########################################

#Read data and code factors as factors
pep <- read.csv("data/Ratings_filt.csv")
pep$Rep <- as.factor(pep$Rep)
pep$Block <- as.factor(pep$Block)
pep$Tray <- as.factor(pep$Tray)
pep$Column <- as.factor(pep$Column)

#Read relationship matrix and set row and column names
K <- read.csv("data/K_full.csv")
rownames(K) <- K$X
K$X <- NULL
colnames(K) <- rownames(K)
K <- as.matrix(K)

#### Put in alphanumerical order
K <- K[order(rownames(K)), order(colnames(K))]
###

check_rows <- which(pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"))
pep$New <- rep(2, nrow(pep))
pep$New[check_rows] <- 1
pep$New <- as.factor(pep$New)

#####################    Loop through peppers and fit separate model for each pepper    #########################

peppers <- levels(pep$Pepper)
#peppers <- "RedKnight"

#Set up matrices to store heritabilities
H2s <- matrix(NA,ncol=2, nrow=length(peppers), dimnames=list(peppers, c("H2", "SE")))
h2s <- matrix(NA,ncol=2, nrow=length(peppers), dimnames=list(peppers, c("h2", "SE")))

#Set up vector to store correlations of predictions from both models
blup_correlations <- rep(NA, length(peppers))
names(blup_correlations) <- peppers

#Set up vector to store x-validation results
xval_accuracies <- rep(NA, length(peppers))
names(xval_accuracies) <- peppers

for(i in 1:length(peppers)){

	pepper <- peppers[i]
	pep.sub <- pep[pep$Pepper == pepper,]
	pep.sub$Isolate <- droplevels(pep.sub$Isolate)

	#Need to get rid of isolates in K that are not in data
	K_sub <- K
	K_sub <- K[rownames(K) %in% levels(pep.sub$Isolate),colnames(K) %in% levels(pep.sub$Isolate)]

	#Fit model WITH and WITHOUT relationship matrix

	mod <- asreml(fixed = audpc ~ at(New, 1):Isolate + Rep,
			 random = ~ at(New, 2):Isolate + Rep:Block,
			 data = pep.sub, na.action = na.method(x="include",y="include"))

	mod.GBLUP <- asreml(fixed = audpc ~ at(New, 1):Isolate + Rep,
			 random = ~ Rep:Block +
			 at(New, 2):vm(Isolate, K_sub),
			 data = pep.sub, na.action = na.method(x="include",y="include"))

	#Calculate heritabilities
	
	H2s[i,] <- unlist(vpredict(mod, H2 ~ V2/(V2+V3)))
	h2s[i,] <- unlist(vpredict(mod.GBLUP, h2 ~ V2/(V2+V3)))
	
	#Look at predictions

	mod.predict <- as.data.frame(predict(mod, classify="Isolate:New", levels=list("New"=2))$pvals)
	mod.GBLUP.predict <- as.data.frame(predict(mod.GBLUP, classify="Isolate:New", levels=list("New"=2))$pvals)
	blup_correlations[i] <- cor(mod.predict$predicted.value, mod.GBLUP.predict$predicted.value)

	#Cross_validation
	n_vals <- 20
	accuracies <- rep(NA, n_vals)
	for(j in 1:n_vals){
		pep.miss <- pep.sub
		experimental_isolates <- levels(pep.miss$Isolate)
		experimental_isolates <- experimental_isolates[-which(experimental_isolates %in% c("CHECK1", "CHECK2", "CHECK3"))]
		missing_isolates <- sample(experimental_isolates, round(length(experimental_isolates)*0.20), replace=F)
		pep.miss$audpc[pep.miss$Isolate %in% missing_isolates] <- NA
		mod.miss <- asreml(fixed = audpc ~ at(New, 1):Isolate + Rep,
				 random = ~ Rep:Block +
				 at(New, 2):vm(Isolate, K_sub),
				 data = pep.miss)
		mod.miss.predict <- as.data.frame(predict(mod.miss, classify="Isolate:New", levels=list("New"=2))$pvals)
		accuracies[j] <- cor(mod.miss.predict$predicted.value[mod.miss.predict$Isolate %in% missing_isolates],
		         	 mod.predict$predicted.value[mod.predict$Isolate %in% missing_isolates])
		}
	xval_accuracies[i] <- median(accuracies)

}
