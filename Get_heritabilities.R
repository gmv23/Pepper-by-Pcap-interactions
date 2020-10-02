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
rownames(Ginv) <- NULL
#Make dummy variables
check_rows <- which(pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"))
pep$New <- rep(2, nrow(pep))
pep$New[check_rows] <- 1
pep$New <- as.factor(pep$New)

############################################# test ###########################
#Give the checks their actual names
pep$Isolate[pep$Isolate == "CHECK1"] <- "SJV_CAA"
pep$Isolate[pep$Isolate == "CHECK2"] <- "17EH01C"
pep$Isolate[pep$Isolate == "CHECK3"] <- "G_3A_4A_C5"
pep$Isolate <- droplevels(pep$Isolate)

#####################    Loop through peppers and fit separate model for each pepper    #########################

#Scramble K for the fuck of it
#G_rowNames <- attr(Ginv, "rowNames")[1:(length(attr(Ginv, "rowNames"))-3)]
#G_rowNames_scramble <- sample(G_rowNames, size=length(G_rowNames), replace=F)
#attr(Ginv, "rowNames") <- c(G_rowNames_scramble, "CHECK1", "CHECK2", "CHECK3")

#Function to get residual error as percentage of all variance components
get_residual_error_percent <- function(model.as){
	res_error <- model.as$vparameters["units!R"]
	total_error <- sum(model.as$vparameters)
	return(res_error/total_error)
}

<<<<<<< HEAD
peppers <- "RedKnight"
=======
peppers <- "Aristotle"
>>>>>>> 7385d7814540331318232030da43136a8aac8a81
#peppers <- levels(pep$Pepper)

#Set up data frame to store residual errors of different models
res_errors <- data.frame("Pepper" = peppers, "Raw" = NA, "Log" = NA)

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

	mod.raw <- asreml(fixed = audpc ~ at(New, 1):Isolate + Rep,
			 random = ~ at(New, 2):Isolate + Rep:Block,
			 data = pep.sub)

	mod.log <- asreml(fixed = log(audpc + 1) ~ at(New, 1):Isolate + Rep,
		random = ~ at(New, 2):Isolate + Rep:Block,
		data = pep.sub)

	mod2 <- asreml(fixed = audpc ~ Rep,
		       random = ~ Isolate + Rep:Block,
		       data=pep.sub)

	#Pull out residual error variance component as a percentage
	res_error.raw <- get_residual_error_percent(mod.raw)
	res_error.log <- get_residual_error_percent(mod.log)
	res_errors[i,c("Raw", "Log")] <- c(res_error.raw, res_error.log)

	#Choose model with smaller residual error
	if(res_error.raw < res_error.log){
		mod <- mod.raw
	}else{
		mod <- mod.log
#		pep.sub$audpc <- log(pep.sub$audpc + 1)
	}

	##############
	mod <- mod.raw #also uncomment above in choose model
	###############

	#Append residuals and fitted values to lists
	res[[i]] <- mod$residuals
	fit[[i]] <- mod$fitted.values

	#Get harmonic mean of number reps
	non_check_isolates <- pep.sub$Isolate[pep.sub$New==2]
	non_check_isolates <- droplevels(non_check_isolates)
	r <- table(non_check_isolates)
	n_harm <- length(r)/sum(1/r)

	#Calculate broad-sense heritability
	H2 <- mod$vparameters["at(New, 2):Isolate"]/
	      (mod$vparameters["at(New, 2):Isolate"] + mod$vparameters["units!R"])

	#Only save H2 for traits that are replicated
	if(round(n_harm) > 1){
		H2s[i] <- H2
	}

	#Now run model with genomic relationship matrix
	mod.GBLUP <- asreml(fixed = audpc ~ at(New, 1):Isolate + Rep,
			 random = ~ Rep:Block +
			 at(New, 2):vm(Isolate, Ginv),
			 data = pep.sub)

	mod.GBLUP2 <- asreml(fixed = audpc ~ Rep,
		       random = ~ vm(Isolate, Ginv) + Rep:Block,
		       data=pep.sub)

	#Get narrow sense heritability
	h2 <- mod.GBLUP$vparameters["at(New, 2):vm(Isolate, Ginv)"]/
	      (mod.GBLUP$vparameters["at(New, 2):vm(Isolate, Ginv)"] +
	       mod.GBLUP$vparameters["units!R"])
	h2s[i] <- h2
}
