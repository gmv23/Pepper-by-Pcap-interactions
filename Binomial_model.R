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
<<<<<<< HEAD

=======
K <- readRDS("data/K.rds")
>>>>>>> 7385d7814540331318232030da43136a8aac8a81
#Make dummy variables
check_rows <- which(pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"))
pep$New <- rep(2, nrow(pep))
pep$New[check_rows] <- 1
pep$New <- as.factor(pep$New)

<<<<<<< HEAD
#####################    Loop through peppers and fit separate model for each pepper    #########################

#Scramble K for the fuck of it
=======
######################Test
pep$Isolate[pep$Isolate=="CHECK1"] <- "SJV_CAA"
pep$Isolate[pep$Isolate=="CHECK2"] <- "17EH01C"
pep$Isolate[pep$Isolate=="CHECK3"] <- "G_3A_4A_C5" 
pep$Isolate <- droplevels(pep$Isolate)

#####################    Loop through peppers and fit separate model for each pepper    #########################

#Scramble K to see what happens
>>>>>>> 7385d7814540331318232030da43136a8aac8a81
#G_rowNames <- attr(Ginv, "rowNames")[1:(length(attr(Ginv, "rowNames"))-3)]
#G_rowNames_scramble <- sample(G_rowNames, size=length(G_rowNames), replace=F)
#attr(Ginv, "rowNames") <- c(G_rowNames_scramble, "CHECK1", "CHECK2", "CHECK3")

<<<<<<< HEAD
peppers <- "RedKnight"
=======
peppers <- "Aristotle"
>>>>>>> 7385d7814540331318232030da43136a8aac8a81
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

<<<<<<< HEAD
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
=======
	#Fit model without relationship matrix
	mod <- asreml(fixed = Dpi15 ~ at(New, 1):Isolate + Rep,
	       random = ~ at(New, 2):Isolate + Rep:Block,
	       data = pep.sub, maxit=200,
	       family = asr_binomial(link="logit", total=Plot_no))

	mod2 <- asreml(fixed = Dpi15 ~ Rep,
	       random = ~ Isolate + Rep:Block,
	       data = pep.sub, maxit=200,
	       family = asr_binomial(link="logit", total=Plot_no))

	if(!mod$converge){
		stop(paste("Broad-sense model for", pepper, "failed to converge"))
	}

	#Append residuals and fitted values to lists
	res[[i]] <- mod$residuals
	fit[[i]] <- mod$fitted.values

	#Get harmonic mean of number reps
	non_check_isolates <- pep.sub$Isolate[pep.sub$New==2]
	non_check_isolates <- droplevels(non_check_isolates)
	r <- table(non_check_isolates)
	n_harm <- length(r)/sum(1/r)

	#Get  error variance in logistic regression
	error_variance <- (pi^2)/3

	#Calculate broad-sense heritability
	H2 <- mod$vparameters["at(New, 2):Isolate"]/
	      (mod$vparameters["at(New, 2):Isolate"] + error_variance/n_harm)

	#Only save H2 for traits that are replicated
	if(round(n_harm) > 1){
		H2s[i] <- H2
	}

	#Now run model with genomic relationship matrix
	mod.GBLUP <- asreml(fixed = Dpi15 ~ at(New, 1):Isolate + Rep,
		 random = ~ at(New, 2):vm(Isolate, Ginv)  + Rep:Block,
		 data = pep.sub, maxit=200,
		 family = asr_binomial(link="logit", total=Plot_no))	

	mod.GBLUP2 <- asreml(fixed = Dpi15 ~ Rep,
		 random = ~ vm(Isolate, Ginv) + Isolate + Rep:Block,
		 data = pep.sub, maxit=200,
		 family = asr_binomial(link="logit", total=Plot_no, dispersion=NA))	

	mod.GBLUP3 <- asreml(fixed = Dpi15 ~ Rep,
		 random = ~ vm(Isolate, K)  + Rep:Block,
		 data = pep.sub, maxit=200,
		 family = asr_binomial(link="logit", total=Plot_no))	

	if(!mod.GBLUP$converge){
		stop(paste("Narrow-sense model for", pepper, "failed to converge"))
	}

	#Get narrow sense heritability
	h2 <- mod.GBLUP$vparameters["at(New, 2):vm(Isolate, Ginv)"]/
	      (mod.GBLUP$vparameters["at(New, 2):vm(Isolate, Ginv)"] + error_variance)
	h2s[i] <- h2

>>>>>>> 7385d7814540331318232030da43136a8aac8a81
}
