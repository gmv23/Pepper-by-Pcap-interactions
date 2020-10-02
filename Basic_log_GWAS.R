setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/pheno/")

library(reshape2)
library(lme4)

#Import data
pep <- read.csv("data/Ratings_filt.csv", na.strings = "NA") 
pep$Rep <- as.factor(pep$Rep)
pep$Block <- as.factor(pep$Block)
pep$Tray <- as.factor(pep$Tray)
pep$Column <- as.factor(pep$Column)

#Reformat data
pep.counts <- aggregate(Dpi15 ~ Isolate*Pepper, data=pep, FUN=sum)
counts.wide <- dcast(Isolate ~ Pepper, value.var = "Dpi15", data=pep.counts)
counts.wide <- counts.wide[!counts.wide$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"),]

#Read in genotypes
geno <- read.csv("data/geno_filt.csv")
rownames(geno) <- geno$X
geno$X <- NULL
geno <- geno[match(counts.wide$Isolate, rownames(geno)),]

nonpatho <- apply(counts.wide[,-1], 1, function(x) all(x==0))
counts.wide <- counts.wide[-which(nonpatho),]
geno <- geno[-which(nonpatho),]

pvals <- rep(NA, ncol(geno))
for(i in 1:ncol(geno)){
  if(i %% 100 == 0){
    print(paste(round(i/ncol(geno), 2)*100, "%", sep=""))
  }
  snp_data <- data.frame("Dead" = counts.wide$Aristotle,
                    "Total" = 12,
                    "SNP" = geno[,i])
  snp_data$SNP <- as.factor(snp_data$SNP)
  snp.model <- glm(cbind(Dead,Total-Dead) ~ 1 + SNP, data=snp_data, family=quasibinomial)
  snp.anova <- Anova(snp.model)
  pvals[i] <- snp.anova$`Pr(>Chisq)`
}

test_plot <- data.frame("BP" = snps$V2, "CHR" = snps$V1, "P" = pvals)
qq(pvals)
manhattan(test_plot, genomewideline = F, suggestiveline = -log10(.05/length(pvals)))

