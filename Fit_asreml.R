#usr/bin/R env

setwd("/workdir/gmv23/peppers/pheno/asreml")

library(asreml)

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

#Look at full model
pep.full <- asreml(fixed = audpc ~ at(New, 1):Isolate + at(New, 1):Isolate:Pepper,
		   random = ~ at(New, 2):Isolate + Rep/Block/Tray + Pepper + at(New, 2):Isolate:Pepper,
		   data = pep)

pep.full.predict <- predict(pep.full,
		    classify = "Isolate:Pepper:New",
		    levels = list(New=2),
		    present = c("Isolate", "Pepper", "New"))
pep.full.blups <- pep.full.predict$predictions$pvals


pep.fixed <- asreml(fixed = audpc ~ Isolate + Pepper + Isolate:Pepper,
	     random = ~ Rep/Block/Tray,
	     data = pep)


#Get residual error as percentage of all variance components
get_residual_error_percent <- function(model.as){
	res_error <- model.as$gammas["R!variance"]
	total_error <- sum(model.as$gammas)
	return(res_error/total_error)
}

res_errors <- matrix(NA, ncol=2, nrow=nlevels(pep$Pepper))
rownames(res_errors) <- levels(pep$Pepper)

#Loop through peppers
for(i in 1:nlevels(pep$Pepper)){
	pepper <- levels(pep$Pepper)[i]

	pep.sub <- pep[pep$Pepper == pepper,]

	mod <- asreml(fixed = audpc ~ at(New, 1):Isolate,
			 random = ~ at(New, 2):Isolate + Rep:Block,
			 data = pep.sub)

	mod.log <- asreml(fixed = log(audpc + 1) ~ at(New, 1):Isolate,
		random = ~ at(New, 2):Isolate + Rep:Block,
		data = pep.sub)
	res_errors[i,1] <- get_residual_error_percent(mod)
	res_errors[i,2] <- get_residual_error_percent(mod.log)
}

