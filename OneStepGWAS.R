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


#Read in data frame with lower triangle of inverse relationshp matrix
Ginv <- readRDS("data/K_ginv.rds")
rownames(Ginv) <- NULL

#Pcs
pcs <- read.csv("data/pcs.csv")
pcs$Isolate <- pcs$X

#Geno
geno <- read.csv("data/geno_filt.csv")
rownames(geno) <- geno$X
geno$X <- NULL

#Make dummy variables
check_rows <- which(pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"))
pep$New <- rep(2, nrow(pep))
pep$New[check_rows] <- 1
pep$New <- as.factor(pep$New)

#Give the checks their actual names
pep$Isolate[pep$Isolate == "CHECK1"] <- "SJV_CAA"
pep$Isolate[pep$Isolate == "CHECK2"] <- "17EH01C"
pep$Isolate[pep$Isolate == "CHECK3"] <- "G_3A_4A_C5"
pep$Isolate <- droplevels(pep$Isolate)

#For now make pep.sub Aristotle
pep.sub <- pep[pep$Pepper == "RedKnight",]
pep.sub <- cbind(pep.sub, pcs[match(pep.sub$Isolate, pcs$Isolate),-c(1,6)])

#Get the AICs of a bunch of different models
skip <- function(x){
aics <- rep(NA, 10)

mod1 <- asreml(fixed=Dpi15 ~ Rep, random = ~ Isolate + Rep:Block,
data = pep.sub, maxit=200, family=asr_binomial(link="logit", total=Plot_no))
aics[1] <- summary(mod1)$aic
mod2 <- asreml(fixed=Dpi15 ~ Rep + PC1, random = ~ Isolate + Rep:Block,
data = pep.sub, maxit=200, family=asr_binomial(link="logit", total=Plot_no))
aics[2] <- summary(mod2)$aic
mod3 <- asreml(fixed=Dpi15 ~ Rep + PC1 + PC2, random = ~ Isolate + Rep:Block,
data = pep.sub, maxit=200, family=asr_binomial(link="logit", total=Plot_no))
aics[3] <- summary(mod3)$aic
mod4 <- asreml(fixed=Dpi15 ~ Rep + PC1 + PC2 + PC3, random = ~ Isolate + Rep:Block,
data = pep.sub, maxit=200, family=asr_binomial(link="logit", total=Plot_no))
aics[4] <- summary(mod4)$aic
mod5 <- asreml(fixed=Dpi15 ~ Rep + PC1 + PC2 + PC3 + PC4, random = ~ Isolate + Rep:Block,
data = pep.sub, maxit=200, family=asr_binomial(link="logit", total=Plot_no))
aics[5] <- summary(mod5)$aic
mod6 <- asreml(fixed=Dpi15 ~ Rep, random = ~ vm(Isolate, Ginv) + Rep:Block,
data = pep.sub, maxit=200, family=asr_binomial(link="logit", total=Plot_no))
aics[6] <- summary(mod6)$aic
mod7 <- asreml(fixed=Dpi15 ~ Rep + PC1, random = ~ vm(Isolate, Ginv) + Rep:Block,
data = pep.sub, maxit=200, family=asr_binomial(link="logit", total=Plot_no))
aics[7] <- summary(mod7)$aic
mod8 <- asreml(fixed=Dpi15 ~ Rep + PC1 + PC2, random = ~ vm(Isolate, Ginv) + Rep:Block,
data = pep.sub, maxit=200, family=asr_binomial(link="logit", total=Plot_no))
aics[8] <- summary(mod8)$aic
mod9 <- asreml(fixed=Dpi15 ~ Rep + PC1 + PC2 + PC3, random = ~ vm(Isolate, Ginv) + Rep:Block,
data = pep.sub, maxit=200, family=asr_binomial(link="logit", total=Plot_no))
aics[9] <- summary(mod9)$aic
mod10 <- asreml(fixed=Dpi15 ~ Rep + PC1 + PC2 + PC3 + PC4, random = ~ vm(Isolate, Ginv) + Rep:Block,
data = pep.sub, maxit=200, family=asr_binomial(link="logit", total=Plot_no))
aics[10] <- summary(mod10)$aic
}

null_mod <- asreml(fixed=Dpi15 ~ 1 + PC1 + PC4, random = ~ Isolate + Rep:Block,
data = pep.sub, maxit=200, family=asr_binomial(link="logit", total=Plot_no), na.action = na.method(x = "include"))
null_G_param <- null_mod$G.param
null_R_param <- null_mod$R.param

pvals <- rep(NA, ncol(geno))

for(i in 1:length(pvals)){
	if(i %% 500 == 0){
		print(round(paste(i/ncol(geno)*100,2), "% done", sep=''))
	}
	pep.gwas <- cbind(pep.sub, "SNP" = geno[match(pep.sub$Isolate, rownames(geno)),i])
	pep.gwas$SNP <- as.factor(pep.gwas$SNP)
	mod <- asreml(fixed=Dpi15 ~ 1 + PC1 + PC4 + SNP, random = ~ Isolate + Rep:Block,
	data = pep.gwas, maxit=200, family=asr_binomial(link="logit", total=Plot_no), na.action = na.method(x = "include"),
	G.param = null_G_param, R.param = null_R_param, trace=F)
	fix_test <- as.data.frame(wald(mod))
	pvals[i] <- fix_test["SNP",4]
}

write.csv(pvals, "out/RK_PC14.csv")
