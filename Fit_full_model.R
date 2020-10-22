#usr/bin/R env

setwd("/workdir/gmv23/peppers/models/")

library(asreml)
library(reshape2)

##################################    Read and clean    #############################

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

#Make balanced data set excluding peppers only included in 1 rep
pep.bal <- pep[!pep$Pepper %in% c("Archimedes", "Revolution", "Vanguard"),]
pep.bal$Pepper <- droplevels(pep.bal$Pepper)
#################################  Fit  models   #############################

mod <- asreml(fixed = audpc ~ 1 + Isolate + Pepper + Isolate:Pepper + Rep,
		   random = ~ Rep:Block/Tray,
		   data = pep, na.action=na.method(x="include",y="include"))

mod.bal <- asreml(fixed = audpc ~ 1 + Isolate + Pepper + Isolate:Pepper + Rep,
		   random = ~ Rep:Block/Tray,
		   data = pep.bal, na.action=na.method(x="include",y="include"))

#################################  Get BLUEs   ################################

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

write.csv(pep.full.blues, "tables/virulence_blues.csv", quote=F, row.names=F)

#Get LS means for peppers as well
pepper.means <- predict(mod,
		classify = "Pepper",
		present = c("Isolate", "Pepper"))$pvals
write.csv(pepper.means, "tables/pepper_blues.csv", quote=F, row.names=F)

#################################  Get model terms and p-values  ################################

#Function to supply nested asreml models and get LRT p-value for random term
LRT <- function(mod.full, mod.red, df){
	LR <- 2*(mod.full$loglik - mod.red$loglik)
	pval <- pchisq(q=LR, df = df, lower.tail=F)
	return(c(LR,pval))
}

#Reduced model without Tray
mod.redTray <- asreml(fixed = audpc ~ Isolate + Pepper + Isolate:Pepper + Rep,
		   random = ~ Rep:Block,
		   data = pep.bal, na.action=na.method(x="include",y="include"))

#Reduced model without Block
mod.redBlock <- asreml(fixed = audpc ~ Isolate + Pepper + Isolate:Pepper + Rep,
		   random = ~ Rep:Block:Tray,
		   data = pep.bal, na.action=na.method(x="include",y="include"))

#Get variance components
pep.full.vcs <- as.data.frame(mod.bal$vparameters*mod.bal$sigma2)
colnames(pep.full.vcs) <- "Variance"
pep.full.vcs$Variance_percent <- pep.full.vcs$Variance/sum(pep.full.vcs$Variance)*100

#Add LRTs to variance components data frame
pep.full.vcs$LRT <- NA
pep.full.vcs$pval <- NA
pep.full.vcs[c("Rep:Block:Tray",
	       "Rep:Block"),
	      c("LRT","pval")] <- rbind(LRT(mod.bal, mod.redTray, 1),
					LRT(mod.bal, mod.redBlock, 1))
pep.full.vcs_round <- as.data.frame(apply(pep.full.vcs[1:3],2,round,digits=2))
pep.full.vcs_round$pval <- signif(pep.full.vcs$pval,2)
write.csv(pep.full.vcs_round, "tables/variance_components.csv", quote=F, row.names=T)

#Test for fixed effects
fixed_tests <- as.data.frame(wald(mod.bal, denDF="default")$Wald)
fixed_tests_round <- as.data.frame(apply(fixed_tests[,1:3],2,round,digits=2))
fixed_tests_round$Pr <- signif(fixed_tests$Pr, 2)	
write.csv(fixed_tests_round, "tables/Fixed_effects_tests.csv", quote=F, row.names=T)





