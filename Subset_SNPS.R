setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/effectors/")

snps <- read.table("../geno/data/capsici_pepper_subset.012.pos")
colnames(snps) <- c("Chrom", "Bp")

eff <- read.table("data/Pc_effectors_blastout.txt")
colnames(eff) <- c("qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore")
eff$sseqid <- sapply(as.character(eff$sseqid), function(x) unlist(strsplit(x, "_"))[2])

snps$eff <- FALSE
distance <- 0

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

genes <- read.delim("../../../isolate_collection/analysis/genome/annotations.txt", sep='\t')
genes$chrom <- sapply(as.character(genes$chrom), function(x) unlist(strsplit(x, "_"))[2])
genes.ran <- genes[sample(1:nrow(genes), size=nrow(eff), replace=F),1:4]

snps$ran <- FALSE
distance <- 0

for(i in 1:nrow(snps)){
  print(i)
  chr <- snps$Chrom[i]
  bp <- snps$Bp[i]
  if(chr %in% genes.ran$chrom){
    genes.chr <- genes.ran[genes.ran$chrom == chr,]
    for(j in 1:nrow(genes.chr)){
      coords <- genes.chr[j,c("start","end")]
      start <- min(coords) - distance
      stop <- max(coords) + distance
      if(bp > start & bp < stop){
        snps$ran[i] <- TRUE
      }
    }
  }
}


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
geno.ran <- geno[,snps$ran]
library(rrBLUP)

K <- A.mat(geno)
K.sub <- A.mat(geno.sub)
K.ran <- A.mat(geno.ran)
n <- nrow(geno)
traits <- colnames(phenos)
accuracies.full <- matrix(NA, nrow=50,ncol=length(traits))
colnames(accuracies.full) <- traits
accuracies.sub <- matrix(NA, nrow=50,ncol=length(traits))
colnames(accuracies.sub) <- traits
accuracies.ran <- matrix(NA, nrow=50,ncol=length(traits))
colnames(accuracies.ran) <- traits

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
    rrblup.ran <- kin.blup(data = pheno.fold, geno = "Sample", pheno = "audpc", K=K.ran, fixed="intercept")
    gebvs.full <- rrblup.full$g
    gebvs.sub <- rrblup.sub$g
    gebvs.ran <- rrblup.ran$g
    accuracies.full[i,j] <- cor(pheno_col[which.test], gebvs.full[which.test], use="pairwise.complete.obs")
    accuracies.sub[i,j] <- cor(pheno_col[which.test], gebvs.sub[which.test], use="pairwise.complete.obs")
    accuracies.ran[i,j] <- cor(pheno_col[which.test], gebvs.ran[which.test], use="pairwise.complete.obs")
  }
}

library(RColorBrewer)
all_accuracies <- matrix(NA, nrow=50, ncol=ncol(accuracies.full)*3)
all_accuracies[,seq(1,25,by=3)] <- accuracies.full
all_accuracies[,seq(2,26,by=3)] <- accuracies.sub
all_accuracies[,seq(3,27,by=3)] <- accuracies.ran
jpeg("plots/prediction_accuracies.jpeg")
old.par <- par(no.readonly = T)
par(mar=c(7,3,1,8.5), xpd=NA)
draw_at <- (1:36)[!(1:36 %% 4 == 0)]
colrs <- rep(brewer.pal(3, "Dark2"),9)
boxplot(all_accuracies, at =draw_at, xaxt='n', col=colrs, ylim=c(min(all_accuracies)*1.1,1))
label_at <- seq(2,34,by=4)
axis(1, at=label_at,labels=traits, las=2)
legend(37,0.5,fill=colrs[1:3],legend=c("All SNPs", "Effectors", "Random genes"), bty='n')
par(old.par)
dev.off()

