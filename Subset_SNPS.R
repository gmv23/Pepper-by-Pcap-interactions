setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/effectors/")

snps <- read.table("../geno/data/capsici_pepper_subset.012.pos")
colnames(snps) <- c("Chrom", "Bp")

eff <- read.table("data/Pc_effectors_blastout.txt")
colnames(eff) <- c("qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore")
eff$sseqid <- sapply(as.character(eff$sseqid), function(x) unlist(strsplit(x, "_"))[2])

snps$eff <- FALSE
distance <- 12000

for(i in 1:nrow(snps)){
  print(i)
  chr <- snps$Chrom[i]
  bp <- snps$Bp[i]
  if(chr %in% eff$sseqid){
    eff.chr <- eff[eff$sseqid == chr,]
    for(j in 1:nrow(eff.chr)){
      coords <- eff.chr[j,c("sstart","send")]
      start <- min(coords) - distance
      stop <- max(coords) + distance
      if(bp > start & bp < stop){
        snps$eff[i] <- TRUE
      }
    }
  }
}
write.csv(snps, "data/snps_annotated.csv")

pvals <- read.csv("../gwas/data/gwas_pvalues.csv")
phenos <- read.csv("../gwas/tables/transformed_blues.csv")
rownames(phenos) <- phenos$X
phenos$X <- NULL
geno <- read.table("../geno/data/capsici_pepper_subset.012")
indvs <- read.table("../geno/data/capsici_pepper_subset.012.indv")
geno$V1 <- NULL
rownames(geno) <- indvs$V1
geno <- as.matrix(geno)
geno[geno==-1] <- NA
geno <- geno - 1

geno.sub <- geno[,snps$eff]
library(rrBLUP)

K <- A.mat(geno)
K.sub <- A.mat(geno.sub)
n <- nrow(geno)
traits <- colnames(phenos)
accuracies.full <- matrix(NA, nrow=50,ncol=length(traits))
colnames(accuracies.full) <- traits
accuracies.sub <- matrix(NA, nrow=50,ncol=length(traits))
colnames(accuracies.sub) <- traits

for(j in 1:length(traits)){
  for(i in 1:50){
    trait <- traits[j]
    which.train <- sample(n, round(0.8*n), replace=F)
    which.test <- (1:n)[!(1:n) %in% which.train]
    pheno_col <- phenos[,trait]
    pheno.fold <- data.frame("Sample" = rownames(phenos),
                             "audpc" = pheno_col)
    pheno.fold[which.test,'audpc'] <- NA
    pheno.fold$intercept <- 1
    rrblup.full <- kin.blup(data = pheno.fold, geno = "Sample", pheno = "audpc", K=K, fixed="intercept")
    rrblup.sub <- kin.blup(data = pheno.fold, geno = "Sample", pheno = "audpc", K=K.sub, fixed="intercept")
    gebvs.full <- rrblup.full$g
    gebvs.sub <- rrblup.sub$g
    accuracies.full[i,j] <- cor(pheno_col[which.test], gebvs.full[which.test], use="pairwise.complete.obs")
    accuracies.sub[i,j] <- cor(pheno_col[which.test], gebvs.sub[which.test], use="pairwise.complete.obs")
  }
}
