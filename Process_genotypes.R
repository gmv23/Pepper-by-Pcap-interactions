##Filter genotype data for isolates in experiment and MAF > 0.05
##Create relationship matrix and calculate principal components

setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/geno/")
library(rrBLUP)
library(matrixcalc)

#Load geno data
geno <- read.table("data/capsici_diversity_PG.012")
indvs <- read.table("data/capsici_diversity_PG.012.indv", stringsAsFactors = F)
snps <- read.table("data/capsici_diversity_PG.012.pos")

#Load pheno data
pep <- read.csv("../pheno/data/Ratings_filt.csv")

#Clean up data
indvs <-unlist(indvs$V1)
indvs <- sapply(indvs, function(x) unlist(strsplit(x,":"))[1])

geno <- geno[,-1]
geno <- as.matrix(geno)
geno[geno==-1] <- NA

#There are two isolates that were inoculated where their clones were included in the clone-corrected data set
#Therefore the genotype for 17PZ18A needs to be relabeled 17PZ21A
#And the genotype for 14_55C relabled 14_55
#They are clones so the genotypes are essentially identical
indvs[indvs=="17PZ18A"] <- "17PZ21A"
indvs[indvs=="14_55C"] <- "14_55"
rownames(geno) <- indvs

#Subset geno for those that were phenotyped and pathogenic
geno <- geno[indvs %in% pep$Isolate,]
indvs <- indvs[indvs %in% pep$Isolate]

#Filter for MAF > minMAF
minMAF <- 0.05

maf <- function(x){
  x <- x[!is.na(x)]
  af <- (sum(x==0)*2 + sum(x==1)) / (length(x)*2)
  afs <- c(af, 1-af)
  return(min(afs))
}
mafs <- apply(geno, 2, maf)

sum(mafs<minMAF)
geno <- geno[,-which(mafs<minMAF)]
snps <- snps[-which(mafs<minMAF),]

#Calculate genomic relationship matrix
K <- A.mat(geno-1)

#For one-step models, need to add arbitrary elements for Checks 1-3 so that all isolates are represented
K.full <- cbind(K, matrix(0, nrow=nrow(K), ncol=3))
K.full <- rbind(K.full, matrix(0, nrow=3, ncol=ncol(K.full)))
pad_indices <- (ncol(K.full)-2):ncol(K.full)
diag(K.full)[pad_indices] <- rep(1,3)
rownames(K.full) <- c(indvs, "CHECK1", "CHECK2", "CHECK3")
colnames(K.full) <- c(indvs, "CHECK1", "CHECK2", "CHECK3")

#If matrix is not positive semi definite, then add 10^-6 to the diagonal
if(!is.positive.definite(K.full)){
  diag(K.full) <- diag(K.full) + 10^-6
}

#Get Pcs
library(pcaMethods)
geno.pca <- pca(geno, nPcs = 4, method = "nipals", scale="uv")
write.csv(geno.pca@scores, "data/pcs.csv", quote=F, row.names=T)

write.table(snps, "data/GWAS_positions.txt", col.names = F, row.names = F, quote=F)
gwas_indvs <- rownames(geno)
gwas_indvs[gwas_indvs == "14_55"] <- "14_55C"
gwas_indvs[gwas_indvs == "17PZ21A"] <- "17PZ18A"
write.table(gwas_indvs, "data/GWAS_individuals.txt", col.names=F, row.names=F, quote=F)

write.csv(K, "data/K.csv", quote=F, row.names=F)
write.csv(K.full, "data/K_full.csv", quote=F, row.names=T)

#Now call vcftools to filter vcf file and get new filtered vcf file
system(paste("vcftools --vcf data/capsici_diversity_PG.recode.vcf", 
             "--keep data/GWAS_individuals.txt",
             "--positions data/GWAS_positions.txt",
             "--recode --out data/capsici_pepper_subset_rename"))

system(paste("sed 's/17PZ18A/17PZ21A/' data/capsici_pepper_subset_rename.recode.vcf |",
             "sed 's/14_55C/14_55/' > data/capsici_pepper_subset.vcf"))

system("rm data/capsici_pepper_subset_rename.recode.vcf")

system("vcftools --vcf data/capsici_pepper_subset.vcf --012 --out")


