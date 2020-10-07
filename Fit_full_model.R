#usr/bin/R env

setwd("/workdir/gmv23/peppers/models/")

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

#pep <- pep[order(pep$Pepper),]
#Fit full model

mod <- asreml(fixed = audpc ~ Isolate + Pepper + Isolate:Pepper + Rep,
		   random = ~ Rep:Block/Tray,
		   data = pep, na.action=na.method(x="include",y="include"))

#mod2 <- asreml(fixed = audpc ~ Isolate + Pepper + Isolate:Pepper + Rep,
#		   random = ~ Rep:Block/Tray,
#		   residual = ~ dsum(~id(units) | Pepper),
#		   data = pep, maxit=500)


#Get pepper-specific BLUPs and main-effect BLUPs
pep.full.predict <- predict(mod,
		    classify = "Isolate:Pepper",
		    present = c("Isolate", "Pepper"))
pep.full.blues <- pep.full.predict$pvals
pep.full.blues <- pep.full.blues[!pep.full.blues$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"),]
pep.full.blues <- dcast(Isolate~Pepper, value.var="predicted.value", data=pep.full.blues)

pep.main.predict <- predict(mod, classify="Isolate")
pep.main.blues <- pep.main.predict$pvals
pep.main.blues <- pep.main.blues[!pep.main.blues$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"),]

pep.full.blues$main <- pep.main.blues$predicted.value[match(pep.full.blues$Isolate, pep.main.blues$Isolate)]

write.csv(pep.full.blups, "data/virulence_blues.csv", quote=F, row.names=F)


skip <- function(){

#Function to supply nested asreml models and get LRT p-value for random term
LRT <- function(mod.full, mod.red, df){
	LR <- 2*(mod.full$loglik - mod.red$loglik)
	pval <- pchisq(q=LR, df = df, lower.tail=F)
	return(c(LR,pval))
}

#Reduced model without Tray
pep.redTray <- asreml(fixed = audpc ~ at(New, 1):Isolate + at(New, 1):Isolate:Pepper + Rep,
		   random = ~ at(New, 2):Isolate + Rep:Block + Pepper + at(New, 2):Isolate:Pepper,
		   data = pep)

#Reduced model without Block
pep.redBlock <- asreml(fixed = audpc ~ at(New, 1):Isolate + at(New, 1):Isolate:Pepper + Rep,
		   random = ~ at(New, 2):Isolate + Rep:Block:Tray + Pepper + at(New, 2):Isolate:Pepper,
		   data = pep)

#Reduced model without Pepper
pep.redPepper <- asreml(fixed = audpc ~ at(New, 1):Isolate + at(New, 1):Isolate:Pepper + Rep,
		   random = ~ at(New, 2):Isolate + Rep:Block + Rep:Block:Tray + at(New, 2):Isolate:Pepper,
		   data = pep)

#Reduced model without Isolate
pep.redIsolate <- asreml(fixed = audpc ~ at(New, 1):Isolate + at(New, 1):Isolate:Pepper + Rep,
		   random = ~ Rep:Block + Rep:Block:Tray + Pepper + at(New, 2):Isolate:Pepper,
		   data = pep)

#Reduced model without  Pepper:Isolate interaction
pep.redInteraction <- asreml(fixed = audpc ~ at(New, 1):Isolate + at(New, 1):Isolate:Pepper + Rep,
		   random = ~ at(New, 2):Isolate + Rep:Block + Rep:Block:Tray + Pepper,
		   data = pep)

#Get variance components
pep.full.vcs <- as.data.frame(pep.full$vparameters*pep.full$sigma2)
colnames(pep.full.vcs) <- "Variance"
pep.full.vcs$Variance_percent <- pep.full.vcs$Variance/sum(pep.full.vcs$Variance)*100

#Add LRTs to variance components data frame
pep.full.vcs$LRT <- NA
pep.full.vcs$pval <- NA
pep.full.vcs[c("Rep:Block:Tray",
	       "Rep:Block",
	       "Pepper",
	       "at(New, 2):Isolate",
	       "at(New, 2):Isolate:Pepper"),
	      c("LRT","pval")] <- rbind(LRT(pep.full, pep.redTray, 1),
					LRT(pep.full, pep.redBlock, 1),
					LRT(pep.full, pep.redPepper, 1),
					LRT(pep.full, pep.redIsolate, 1),
					LRT(pep.full, pep.redInteraction, 1))

#Write variance components
write.csv(pep.full.vcs, "out/variance_components.csv", quote=F, row.names=T)

#Conduct Wald test for fixed effects
fixed_tests <- as.data.frame(wald(pep.full))
write.csv(fixed_tests, "out/Fixed_effects_tests.csv", quote=F, row.names=T)

#Get BLUPs
pep.full.predict <- predict(pep.full,
		    classify = "Isolate:Pepper:New",
		    levels = list("New"=2),
		    present = c("Isolate", "Pepper"))
pep.full.blups <- pep.full.predict$pvals
pep.full.blups <- pep.full.blups[!pep.full.blups$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"),]
pep.full.blups.wide <- dcast(Isolate~Pepper, value.var="predicted.value", data=pep.full.blups)

write.csv(pep.full.blups.wide, "out/full_blups.csv", quote=F, row.names=F)

pep.full.predict <- predict(mod,
		    classify = "Isolate:Pepper")
pep.full.blups <- pep.full.predict$pvals
pep.full.blups <- pep.full.blups[!pep.full.blups$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"),]
pep.full.blups <- dcast(Isolate~Pepper, value.var="predicted.value", data=pep.full.blups)

pep.sub.blups <- pep.full.blups
pep.sub.blups[,2:ncol(pep.sub.blups)] <- NA

peppers <- levels(pep$Pepper)
for(i in 1:length(peppers)){
	
	pepper <- peppers[i]
	pep.sub <- pep[pep$Pepper == pepper,]
	pep.sub$Isolate <- droplevels(pep.sub$Isolate)
	
	mod.sub <- asreml(fixed = audpc ~ Isolate + Rep,
			  random = ~ Rep:Block,
			  data= pep.sub)

	mod.sub.predict <- predict(mod.sub, classify="Isolate")
	mod.sub.blups <- mod.sub.predict$pvals
	
	pep.sub.blups[,pepper] <- mod.sub.blups$predicted.value[match(pep.sub.blups$Isolate, mod.sub.blups$Isolate)]
}



}
