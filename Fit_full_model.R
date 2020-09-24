#usr/bin/R env

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

#Fit full model
pep.full <- asreml(fixed = audpc ~ at(New, 1):Isolate + at(New, 1):Isolate:Pepper + Rep,
		   random = ~ at(New, 2):Isolate + Rep:Block + Rep:Block:Tray + Pepper + at(New, 2):Isolate:Pepper,
		   data = pep)

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
		   random = ~ at(New, 2):Isolate + Block:Tray + Pepper + at(New, 2):Isolate:Pepper,
		   data = pep)

#Reduced model without Pepper
pep.redPepper <- asreml(fixed = audpc ~ at(New, 1):Isolate + at(New, 1):Isolate:Pepper + Rep,
		   random = ~ at(New, 2):Isolate + Rep:Block + Block:Tray + at(New, 2):Isolate:Pepper,
		   data = pep)

#Reduced model without Isolate
pep.redIsolate <- asreml(fixed = audpc ~ at(New, 1):Isolate + at(New, 1):Isolate:Pepper + Rep,
		   random = ~ Rep:Block + Block:Tray + Pepper + at(New, 2):Isolate:Pepper,
		   data = pep)

#Reduced model without  Pepper:Isolate interaction
pep.redInteraction <- asreml(fixed = audpc ~ at(New, 1):Isolate + at(New, 1):Isolate:Pepper + Rep,
		   random = ~ at(New, 2):Isolate + Rep:Block + Block:Tray + Pepper,
		   data = pep)
interaction.pval <- LRT(pep.full.unconstrained, pep.redInteraction, 1)

#Get variance components
pep.full.vcs <- as.data.frame(pep.full$gammas*pep.full$sigma2)
colnames(pep.full.vcs) <- "Variance"
pep.full.vcs$Variance_percent <- pep.full.vcs$Variance/sum(pep.full.vcs$Variance)*100

#Add LRTs to variance components data frame
pep.full.vcs$LRT <- NA
pep.full.vcs$pval <- NA
pep.full.vcs[c("Rep:Block:Tray!Rep.var",
	       "Rep:Block!Rep.var",
	       "Pepper!Pepper.var",
	       "at(New, 2):Isolate!Isolate.var",
	       "at(New, 2):Isolate:Pepper!Isolate.var"),
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
		    levels = list(New=2),
		    present = c("Isolate", "Pepper", "New"))
pep.full.blups <- pep.full.predict$predictions$pvals
pep.full.blups.wide <- dcast(Isolate~Pepper, value.var="predicted.value", data=pep.full.blups)

write.csv(pep.full.blups.wide, "out/full_blups.csv", quote=F, row.names=F)



