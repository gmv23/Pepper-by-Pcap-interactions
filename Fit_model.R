setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/pheno/")

library(reshape2)

#Import data
pep <- read.csv("data/Ratings_filt.csv", na.strings = "NA") 
pep$Rep <- as.factor(pep$Rep)
pep$Block <- as.factor(pep$Block)
pep$Tray <- as.factor(pep$Tray)
pep$Column <- as.factor(pep$Column)

#Make some dummy variables
pep$checks <- rep(99, nrow(pep))
pep$checks[pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3")] <- as.character(pep$Isolate[pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3")])
pep$new <- rep(0, nrow(pep))
pep$new[pep$checks==99] <- 1
pep$checks <- as.factor(pep$checks)
pep$block <- pep$Rep:pep$Block

#Random effects augmented model
pep.mm <- lmer(audpc ~ checks + checks:Pepper + (1|Isolate:new) + (1|block) + (1|block:Isolate:new)
               + (1|Pepper) + (1|Pepper:Isolate:new), data=pep,
               control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

jpeg("plots/residuals.jpeg")
hist(residuals(pep.mm), freq=F)
dev.off()
jpeg("plots/fitted_vs_residuals.jpeg")
plot(fitted(pep.mm), residuals(pep.mm))
dev.off()

#Pull out main and main+interaction effects
pep.ranef <- ranef(pep.mm)
main_effects <- pep.ranef$`Isolate:new`
main_effects <- data.frame("Isolate" = sapply(rownames(main_effects), function(x) unlist(strsplit(x, ":"))[1]),
                           "effect" = main_effects$`(Intercept)`)
rownames(main_effects) <- NULL
main_effects <- main_effects[-grep("CHECK", main_effects$Isolate),]

interactions <- pep.ranef$`Pepper:Isolate:new`
interactions <- data.frame("Isolate" = sapply(rownames(interactions), function(x) unlist(strsplit(x, ":"))[2]),
                           "Pepper" = sapply(rownames(interactions), function(x) unlist(strsplit(x, ":"))[1]),
                           "Check" = sapply(rownames(interactions), function(x) unlist(strsplit(x, ":"))[3]),
                           "effect" = interactions$`(Intercept)`)
interactions.mat <- dcast(interactions[interactions$Check==1,], Isolate ~ Pepper, value.var = 'effect')

#Check in same order
all(interactions.mat$Isolate == main_effects$Isolate)
#Make main+interactions data frame
blups <- interactions.mat
for(i in 2:ncol(blups)){
  blups[,i] <- blups[,i] + main_effects$effect
}

#Save variance components
vcs <- as.data.frame(VarCorr(pep.mm))
vcs <- vcs[,c(1,4)]
colnames(vcs) <- c("Term", "Variance")
vcs$Term <- gsub(":new", "", vcs$Term)
vcs$Percent <- round(vcs$Variance/sum(vcs$Variance) * 100,2)
write.csv(vcs, "tables/Variance_components.csv", quote=F, row.names = F)

#Save BLUPs
write.csv(blups, "tables/isolate_blups.csv", quote=F, row.names=F)
