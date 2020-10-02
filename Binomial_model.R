setwd("/workdir/gmv23/peppers/pheno/asreml")

library(asreml)

###################################       Read and clean up data           ########################################

#Read and code data as factors
pep <- read.csv("data/Ratings_filt.csv")
pep$Rep <- as.factor(pep$Rep)
pep$Block <- as.factor(pep$Block)
pep$Tray <- as.factor(pep$Tray)
pep$Column <- as.factor(pep$Column)

#Read in data frame with lower triangle of inverse relationshp matrix
Ginv <- readRDS("data/K_ginv.rds")

#Make dummy variables
check_rows <- which(pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"))
pep$New <- rep(2, nrow(pep))
pep$New[check_rows] <- 1
pep$New <- as.factor(pep$New)

#####################    Loop through peppers and fit separate model for each pepper    #########################

#Scramble K for the fuck of it
#G_rowNames <- attr(Ginv, "rowNames")[1:(length(attr(Ginv, "rowNames"))-3)]
#G_rowNames_scramble <- sample(G_rowNames, size=length(G_rowNames), replace=F)
#attr(Ginv, "rowNames") <- c(G_rowNames_scramble, "CHECK1", "CHECK2", "CHECK3")

peppers <- "RedKnight"
#peppers <- levels(pep$Pepper)

#Set up lists to store residuals and fitted values
res <- vector("list", length(peppers))
names(res) <- peppers
fit <- vector("list", length(peppers))
names(fit) <- peppers

#Set up vector to store heritabilities
H2s <- rep(NA, length(peppers))
h2s <- rep(NA, length(peppers))

for(i in 1:length(peppers)){

	pepper <- peppers[i]
	pep.sub <- pep[pep$Pepper == pepper,]

	#Fit model with non-transformed and log-transformed audpc
	mod.bin <- asreml(fixed = Dpi15 ~ at(New, 1):Isolate + Rep,
	       random = ~ at(New, 2):Isolate + Rep:Block,
	       data = pep.sub,
	       family = asr_binomial(link="logit", total=Plot_no))

	#Append residuals and fitted values to lists
#	res[[i]] <- mod$residuals
#	fit[[i]] <- mod$fitted.values

	#Get harmonic mean of number reps

	#Calculate broad-sense heritability

	#Only save H2 for traits that are replicated

	#Now run model with genomic relationship matrix

	#Get narrow sense heritability
}
