setwd("/workdir/gmv23/peppers/pheno")
library(rrBLUP)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(matrixcalc)

#################################   Load and clean data   ################################
set.seed(385832)
#Reformat VCF file as GDS file
snpgdsVCF2GDS(vcf.fn="data/capsici_pepper_subset.recode.vcf", out.fn="capsici_pepper_subset.gds")

#Load data
geno <- read.csv("data/geno_filt.csv")
rownames(geno) <- geno$X
geno$X <- NULL

snps <- read.csv("data/snps_filt.csv")
colnames(snps) <- c("CHROM", "BP")
snps$MARKER <- paste(snps$CHROM, snps$BP, sep="_")

phenos <- read.csv("data/full_blups.csv")

geno.gds <- GdsGenotypeReader("capsici_pepper_subset.gds")

geno.pca <- read.csv("data/pcs.csv")
rownames(geno.pca) <- geno.pca$X
geno.pca$X <- NULL

#Clean up snps data frame
colnames(snps) <- c("CHROM", "BP", "MARKER")

#For now need to rename 2 samples that are the same but differently labeled between datasets
rownames(geno)[rownames(geno) == "14_55"] <- "14_55C"
rownames(geno)[rownames(geno) == "17PZ21A"] <- "17PZ18A"
rownames(geno.pca)[rownames(geno.pca) == "14_55"] <- "14_55C"
rownames(geno.pca)[rownames(geno.pca) == "17PZ21A"] <- "17PZ18A"
phenos$Isolate <- as.character(phenos$Isolate)
phenos$Isolate[phenos$Isolate == "14_55"] <- "14_55C"
phenos$Isolate[phenos$Isolate == "17PZ21A"] <- "17PZ18A"

#Put phenos in VCF sample order
sample_order <- getScanID(geno.gds)
phenos <- phenos[match(sample_order, phenos$Isolate),]
geno <- geno[match(sample_order, rownames(geno)),]
geno.pca <- geno.pca[match(sample_order, rownames(geno.pca)),]

#################################   Get phenotypes and covariates   ################################

#Make ScanAnnotationDataFrame object with phenotypes and covariates
scanAnnot.df <- data.frame("scanID" = phenos$Isolate,
                           "PC1" = geno.pca[,1],
                           "PC2" = geno.pca[,2],
                           "PC3" = geno.pca[,3],
                           "PC4" = geno.pca[,4],
                           phenos[,-1])

scanAnnot <- ScanAnnotationDataFrame(scanAnnot.df)

#Make GenotypeData object
genoData <- GenotypeData(geno.gds, scanAnnot = scanAnnot)

#Get genomic relationship matrix
K <- A.mat(geno-1)
rownames(K) <- scanAnnot.df$scanID
colnames(K) <- scanAnnot.df$scanID

#Fill in heritability estimates
h2s <- data.frame("Pepper" = c("Archimedes", "Aristotle", "EarlyJalapeno", 
                               "NMRIL-N", "Paladin","Perennial", "RedKnight", "Revolution", "Vanguard"),
                  "h2" = NA,
                  "Lower" = NA,
                  "Upper" = NA)

#Fill in empirical p-value thresholds
pval_thresholds <- data.frame("Pepper" = c("Archimedes", "Aristotle", "EarlyJalapeno", 
                               "NMRIL-N", "Paladin","Perennial", "RedKnight", "Revolution", "Vanguard"),
                              "Threshold5" = NA,
                              "Threshold1" = NA)

#Fill in pvalues
pvalues <- matrix(NA, nrow=ncol(geno), ncol=9)

n_permutations <- 1000
###########################   Run Archimedes   ##########################
print("Starting Archimedes")
#Model Testing
mod_simple <- fitNullModel(scanAnnot, outcome = "Archimedes",
                           family = gaussian)
mod_K <- fitNullModel(scanAnnot, outcome = "Archimedes",
                      cov.mat = K, family = gaussian)
mod_KQ1 <- fitNullModel(scanAnnot, outcome = "Archimedes", covars = c("PC1"),
                        cov.mat = K, family = gaussian)
mod_Q1 <- fitNullModel(scanAnnot, outcome = "Archimedes", covars = c("PC1"),
                       family = gaussian)
mod_KQ2 <- fitNullModel(scanAnnot, outcome = "Archimedes", covars = c("PC1","PC2"),
                        cov.mat = K, family = gaussian)
mod_Q2 <- fitNullModel(scanAnnot, outcome = "Archimedes", covars = c("PC1", "PC2"),
                       family = gaussian)
mod_KQ3 <- fitNullModel(scanAnnot, outcome = "Archimedes", covars = c("PC1","PC2", "PC3"),
                        cov.mat = K, family = gaussian)
mod_Q3 <- fitNullModel(scanAnnot, outcome = "Archimedes", covars = c("PC1", "PC2", "PC3"),
                       family = gaussian)
mod_KQ4 <- fitNullModel(scanAnnot, outcome = "Archimedes", covars = c("PC1","PC2", "PC3", "PC4"),
                        cov.mat = K, family = gaussian)
mod_Q4 <- fitNullModel(scanAnnot, outcome = "Archimedes", covars = c("PC1","PC2", "PC3", "PC4"),
                       family = gaussian)

plot(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC),
     ylab = "AIC", xlab = "", xaxt="n")
axis(1, at=1:10, las="2",
     labels = c("intercept", "K", "K+PC1", "PC1", "K+PC1-2", "PC1-2",
                "K+PC1-3", "PC1-3", "K+PC1-4", "PC1-4"))
text(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC)+2,
     labels = round(c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
                      mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC), 1))

varcomp <- as.data.frame(varCompCI(mod_K))
h2s[1,2:4] <- unlist(varcomp[1,])

#Get genoIterator object dividing SNPs into groups of 10k
genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)

#Run association tests
assoc <- assocTestSingle(gdsobj=genoIterator, null.model=mod_Q4)
pvalues[,1] <- assoc$Score.pval

#Perms
scrambled_pvals <- matrix(NA, nrow=n_permutations, ncol=nrow(snps))
for(i in 1:n_permutations){
  scanAnnot.scramble.df <- scanAnnot.df
  scramble_order <- sample(1:nrow(scanAnnot.scramble.df), nrow(scanAnnot.scramble.df), replace=F)
  scanAnnot.scramble.df[,2:ncol(scanAnnot.scramble.df)] <- scanAnnot.scramble.df[scramble_order,2:ncol(scanAnnot.scramble.df)]
  K_scramble <- K[scramble_order, scramble_order]
  scanAnnot.scramble <- ScanAnnotationDataFrame(scanAnnot.scramble.df)
  genoScramble <- GenotypeData(geno.gds, scanAnnot = scanAnnot.scramble)
  null_scramble <- fitNullModel(scanAnnot.scramble, outcome = "Archimedes", covars = c("PC1", "PC2", "PC3", "PC4"),
                                family = gaussian)
  genoIterator <- GenotypeBlockIterator(genoScramble,snpBlock = 10000)
  assocSCRAMBLE <- assocTestSingle(gdsobj=genoIterator, null.model=null_scramble)
  scrambled_pvals[i,] <- assocSCRAMBLE$Score.pval
}

permute_mins <- apply(scrambled_pvals,1,min)
pval_thresholds[1,2:3] <- c(quantile(permute_mins, .05), quantile(permute_mins, .01))

###########################   Run Aristotle   ##########################
print("Starting Aristotle")
#Model Testing
mod_simple <- fitNullModel(scanAnnot, outcome = "Aristotle",
                           family = gaussian)
mod_K <- fitNullModel(scanAnnot, outcome = "Aristotle",
                      cov.mat = K, family = gaussian)
mod_KQ1 <- fitNullModel(scanAnnot, outcome = "Aristotle", covars = c("PC1"),
                        cov.mat = K, family = gaussian)
mod_Q1 <- fitNullModel(scanAnnot, outcome = "Aristotle", covars = c("PC1"),
                       family = gaussian)
mod_KQ2 <- fitNullModel(scanAnnot, outcome = "Aristotle", covars = c("PC1","PC2"),
                        cov.mat = K, family = gaussian)
mod_Q2 <- fitNullModel(scanAnnot, outcome = "Aristotle", covars = c("PC1", "PC2"),
                       family = gaussian)
mod_KQ3 <- fitNullModel(scanAnnot, outcome = "Aristotle", covars = c("PC1","PC2", "PC3"),
                        cov.mat = K, family = gaussian)
mod_Q3 <- fitNullModel(scanAnnot, outcome = "Aristotle", covars = c("PC1", "PC2", "PC3"),
                       family = gaussian)
mod_KQ4 <- fitNullModel(scanAnnot, outcome = "Aristotle", covars = c("PC1","PC2", "PC3", "PC4"),
                        cov.mat = K, family = gaussian)
mod_Q4 <- fitNullModel(scanAnnot, outcome = "Aristotle", covars = c("PC1","PC2", "PC3", "PC4"),
                       family = gaussian)

plot(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC),
     ylab = "AIC", xlab = "", xaxt="n")
axis(1, at=1:10, las="2",
     labels = c("intercept", "K", "K+PC1", "PC1", "K+PC1-2", "PC1-2",
                "K+PC1-3", "PC1-3", "K+PC1-4", "PC1-4"))
text(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC)+2,
     labels = round(c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
                      mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC), 1))

varcomp <- as.data.frame(varCompCI(mod_K))
h2s[2,2:4] <- unlist(varcomp[1,])

#Get genoIterator object dividing SNPs into groups of 10k
genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)

#Run association tests
assoc <- assocTestSingle(gdsobj=genoIterator, null.model=mod_Q3)
assoc$chr <- snps$CHROM #Provide correct chromosome names
pvalues[,2] <- assoc$Score.pval

#Perms
scrambled_pvals <- matrix(NA, nrow=n_permutations, ncol=nrow(snps))
for(i in 1:n_permutations){
  scanAnnot.scramble.df <- scanAnnot.df
  scramble_order <- sample(1:nrow(scanAnnot.scramble.df), nrow(scanAnnot.scramble.df), replace=F)
  scanAnnot.scramble.df[,2:ncol(scanAnnot.scramble.df)] <- scanAnnot.scramble.df[scramble_order,2:ncol(scanAnnot.scramble.df)]
  K_scramble <- K[scramble_order, scramble_order]
  scanAnnot.scramble <- ScanAnnotationDataFrame(scanAnnot.scramble.df)
  genoScramble <- GenotypeData(geno.gds, scanAnnot = scanAnnot.scramble)
  null_scramble <- fitNullModel(scanAnnot.scramble, outcome = "Aristotle", covars = c("PC1", "PC2", "PC3"),
                                family = gaussian)
  genoIterator <- GenotypeBlockIterator(genoScramble,snpBlock = 10000)
  assocSCRAMBLE <- assocTestSingle(gdsobj=genoIterator, null.model=null_scramble)
  scrambled_pvals[i,] <- assocSCRAMBLE$Score.pval
}

permute_mins <- apply(scrambled_pvals,1,min)
pval_thresholds[2,2:3] <- c(quantile(permute_mins, .05), quantile(permute_mins, .01))


###########################   Run EJ   ##########################
print("Starting Early Jalapeno")
#Model Testing
mod_simple <- fitNullModel(scanAnnot, outcome = "EarlyJalapeno",
                           family = gaussian)
mod_K <- fitNullModel(scanAnnot, outcome = "EarlyJalapeno",
                      cov.mat = K, family = gaussian)
mod_KQ1 <- fitNullModel(scanAnnot, outcome = "EarlyJalapeno", covars = c("PC1"),
                        cov.mat = K, family = gaussian)
mod_Q1 <- fitNullModel(scanAnnot, outcome = "EarlyJalapeno", covars = c("PC1"),
                       family = gaussian)
mod_KQ2 <- fitNullModel(scanAnnot, outcome = "EarlyJalapeno", covars = c("PC1","PC2"),
                        cov.mat = K, family = gaussian)
mod_Q2 <- fitNullModel(scanAnnot, outcome = "EarlyJalapeno", covars = c("PC1", "PC2"),
                       family = gaussian)
mod_KQ3 <- fitNullModel(scanAnnot, outcome = "EarlyJalapeno", covars = c("PC1","PC2", "PC3"),
                        cov.mat = K, family = gaussian)
mod_Q3 <- fitNullModel(scanAnnot, outcome = "EarlyJalapeno", covars = c("PC1", "PC2", "PC3"),
                       family = gaussian)
mod_KQ4 <- fitNullModel(scanAnnot, outcome = "EarlyJalapeno", covars = c("PC1","PC2", "PC3", "PC4"),
                        cov.mat = K, family = gaussian)
mod_Q4 <- fitNullModel(scanAnnot, outcome = "EarlyJalapeno", covars = c("PC1","PC2", "PC3", "PC4"),
                       family = gaussian)

plot(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC),
     ylab = "AIC", xlab = "", xaxt="n")
axis(1, at=1:10, las="2",
     labels = c("intercept", "K", "K+PC1", "PC1", "K+PC1-2", "PC1-2",
                "K+PC1-3", "PC1-3", "K+PC1-4", "PC1-4"))
text(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC)+2,
     labels = round(c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
                      mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC), 1))

varcomp <- as.data.frame(varCompCI(mod_K))
h2s[3,2:4] <- unlist(varcomp[1,])

#Get genoIterator object dividing SNPs into groups of 10k
genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)

#Run association tests
assoc <- assocTestSingle(gdsobj=genoIterator, null.model=mod_K)
assoc$chr <- snps$CHROM #Provide correct chromosome names
pvalues[,3] <- assoc$Score.pval

#Perms
scrambled_pvals <- matrix(NA, nrow=n_permutations, ncol=nrow(snps))
for(i in 1:n_permutations){
  scanAnnot.scramble.df <- scanAnnot.df
  find_order <- F
  while(find_order == F){
    scramble_order <- sample(1:nrow(scanAnnot.scramble.df), nrow(scanAnnot.scramble.df), replace=F)
    K_scramble <- K[scramble_order, scramble_order]
    if(is.positive.semi.definite(K_scramble)){
      find_order <- T
    }
  }
  scanAnnot.scramble.df[,2:ncol(scanAnnot.scramble.df)] <- scanAnnot.scramble.df[scramble_order,2:ncol(scanAnnot.scramble.df)]
  scanAnnot.scramble <- ScanAnnotationDataFrame(scanAnnot.scramble.df)
  genoScramble <- GenotypeData(geno.gds, scanAnnot = scanAnnot.scramble)
  null_scramble <- fitNullModel(scanAnnot.scramble, outcome = "EarlyJalapeno", cov.mat = K_scramble,
                                family = gaussian)
  write.csv(K_scramble, "K_scramble.csv", row.names=T, quote=F)
  genoIterator <- GenotypeBlockIterator(genoScramble,snpBlock = 10000)
  assocSCRAMBLE <- assocTestSingle(gdsobj=genoIterator, null.model=null_scramble)
  scrambled_pvals[i,] <- assocSCRAMBLE$Score.pval
}

permute_mins <- apply(scrambled_pvals,1,min)
pval_thresholds[3,2:3] <- c(quantile(permute_mins, .05), quantile(permute_mins, .01))

###########################   Run NMRILN   ##########################
print("Starting NMRILN")
#Model Testing
mod_simple <- fitNullModel(scanAnnot, outcome = "NMRIL.N",
                           family = gaussian)
mod_K <- fitNullModel(scanAnnot, outcome = "NMRIL.N",
                      cov.mat = K, family = gaussian)
mod_KQ1 <- fitNullModel(scanAnnot, outcome = "NMRIL.N", covars = c("PC1"),
                        cov.mat = K, family = gaussian)
mod_Q1 <- fitNullModel(scanAnnot, outcome = "NMRIL.N", covars = c("PC1"),
                       family = gaussian)
mod_KQ2 <- fitNullModel(scanAnnot, outcome = "NMRIL.N", covars = c("PC1","PC2"),
                        cov.mat = K, family = gaussian)
mod_Q2 <- fitNullModel(scanAnnot, outcome = "NMRIL.N", covars = c("PC1", "PC2"),
                       family = gaussian)
mod_KQ3 <- fitNullModel(scanAnnot, outcome = "NMRIL.N", covars = c("PC1","PC2", "PC3"),
                        cov.mat = K, family = gaussian)
mod_Q3 <- fitNullModel(scanAnnot, outcome = "NMRIL.N", covars = c("PC1", "PC2", "PC3"),
                       family = gaussian)
mod_KQ4 <- fitNullModel(scanAnnot, outcome = "NMRIL.N", covars = c("PC1","PC2", "PC3", "PC4"),
                        cov.mat = K, family = gaussian)
mod_Q4 <- fitNullModel(scanAnnot, outcome = "NMRIL.N", covars = c("PC1","PC2", "PC3", "PC4"),
                       family = gaussian)

plot(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC),
     ylab = "AIC", xlab = "", xaxt="n")
axis(1, at=1:10, las="2",
     labels = c("intercept", "K", "K+PC1", "PC1", "K+PC1-2", "PC1-2",
                "K+PC1-3", "PC1-3", "K+PC1-4", "PC1-4"))
text(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC)+2,
     labels = round(c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
                      mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC), 1))

varcomp <- as.data.frame(varCompCI(mod_K))
h2s[4,2:4] <- unlist(varcomp[1,])

#Get genoIterator object dividing SNPs into groups of 10k
genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)

#Run association tests
assoc <- assocTestSingle(gdsobj=genoIterator, null.model=mod_Q3)
assoc$chr <- snps$CHROM #Provide correct chromosome names
pvalues[,4] <- assoc$Score.pval

#Perms
scrambled_pvals <- matrix(NA, nrow=n_permutations, ncol=nrow(snps))
for(i in 1:n_permutations){
  scanAnnot.scramble.df <- scanAnnot.df
  scramble_order <- sample(1:nrow(scanAnnot.scramble.df), nrow(scanAnnot.scramble.df), replace=F)
  scanAnnot.scramble.df[,2:ncol(scanAnnot.scramble.df)] <- scanAnnot.scramble.df[scramble_order,2:ncol(scanAnnot.scramble.df)]
  K_scramble <- K[scramble_order, scramble_order]
  scanAnnot.scramble <- ScanAnnotationDataFrame(scanAnnot.scramble.df)
  genoScramble <- GenotypeData(geno.gds, scanAnnot = scanAnnot.scramble)
  null_scramble <- fitNullModel(scanAnnot.scramble, outcome = "NMRIL.N", covars = c("PC1", "PC2", "PC3"),
                                family = gaussian)
  genoIterator <- GenotypeBlockIterator(genoScramble,snpBlock = 10000)
  assocSCRAMBLE <- assocTestSingle(gdsobj=genoIterator, null.model=null_scramble)
  scrambled_pvals[i,] <- assocSCRAMBLE$Score.pval
}

permute_mins <- apply(scrambled_pvals,1,min)
pval_thresholds[4,2:3] <- c(quantile(permute_mins, .05), quantile(permute_mins, .01))

###########################   Run Paladin   ##########################
print("Starting Paladin")
#Model Testing
mod_simple <- fitNullModel(scanAnnot, outcome = "Paladin",
                           family = gaussian)
mod_K <- fitNullModel(scanAnnot, outcome = "Paladin",
                      cov.mat = K, family = gaussian)
mod_KQ1 <- fitNullModel(scanAnnot, outcome = "Paladin", covars = c("PC1"),
                        cov.mat = K, family = gaussian)
mod_Q1 <- fitNullModel(scanAnnot, outcome = "Paladin", covars = c("PC1"),
                       family = gaussian)
mod_KQ2 <- fitNullModel(scanAnnot, outcome = "Paladin", covars = c("PC1","PC2"),
                        cov.mat = K, family = gaussian)
mod_Q2 <- fitNullModel(scanAnnot, outcome = "Paladin", covars = c("PC1", "PC2"),
                       family = gaussian)
mod_KQ3 <- fitNullModel(scanAnnot, outcome = "Paladin", covars = c("PC1","PC2", "PC3"),
                        cov.mat = K, family = gaussian)
mod_Q3 <- fitNullModel(scanAnnot, outcome = "Paladin", covars = c("PC1", "PC2", "PC3"),
                       family = gaussian)
mod_KQ4 <- fitNullModel(scanAnnot, outcome = "Paladin", covars = c("PC1","PC2", "PC3", "PC4"),
                        cov.mat = K, family = gaussian)
mod_Q4 <- fitNullModel(scanAnnot, outcome = "Paladin", covars = c("PC1","PC2", "PC3", "PC4"),
                       family = gaussian)

plot(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC),
     ylab = "AIC", xlab = "", xaxt="n")
axis(1, at=1:10, las="2",
     labels = c("intercept", "K", "K+PC1", "PC1", "K+PC1-2", "PC1-2",
                "K+PC1-3", "PC1-3", "K+PC1-4", "PC1-4"))
text(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC)+2,
     labels = round(c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
                      mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC), 1))

varcomp <- as.data.frame(varCompCI(mod_K))
h2s[5,2:4] <- unlist(varcomp[1,])

#Get genoIterator object dividing SNPs into groups of 10k
genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)

#Run association tests
assoc <- assocTestSingle(gdsobj=genoIterator, null.model=mod_K)
assoc$chr <- snps$CHROM #Provide correct chromosome names
pvalues[,5] <- assoc$Score.pval

#Perms
scrambled_pvals <- matrix(NA, nrow=n_permutations, ncol=nrow(snps))
for(i in 1:n_permutations){
  scanAnnot.scramble.df <- scanAnnot.df
  find_order <- F
  while(find_order == F){
    scramble_order <- sample(1:nrow(scanAnnot.scramble.df), nrow(scanAnnot.scramble.df), replace=F)
    K_scramble <- K[scramble_order, scramble_order]
    if(is.positive.semi.definite(K_scramble)){
      find_order <- T
    }
  }
  scanAnnot.scramble.df[,2:ncol(scanAnnot.scramble.df)] <- scanAnnot.scramble.df[scramble_order,2:ncol(scanAnnot.scramble.df)]
  scanAnnot.scramble <- ScanAnnotationDataFrame(scanAnnot.scramble.df)
  genoScramble <- GenotypeData(geno.gds, scanAnnot = scanAnnot.scramble)
  null_scramble <- fitNullModel(scanAnnot.scramble, outcome = "Paladin", cov.mat=K_scramble,
                                family = gaussian)
  genoIterator <- GenotypeBlockIterator(genoScramble,snpBlock = 10000)
  assocSCRAMBLE <- assocTestSingle(gdsobj=genoIterator, null.model=null_scramble)
  scrambled_pvals[i,] <- assocSCRAMBLE$Score.pval
}

permute_mins <- apply(scrambled_pvals,1,min)
pval_thresholds[5,2:3] <- c(quantile(permute_mins, .05), quantile(permute_mins, .01))


###########################   Run Perennial   ##########################
print("Starting Pernennial")
#Model Testing
mod_simple <- fitNullModel(scanAnnot, outcome = "Perennial",
                           family = gaussian)
mod_K <- fitNullModel(scanAnnot, outcome = "Perennial",
                      cov.mat = K, family = gaussian)
mod_KQ1 <- fitNullModel(scanAnnot, outcome = "Perennial", covars = c("PC1"),
                        cov.mat = K, family = gaussian)
mod_Q1 <- fitNullModel(scanAnnot, outcome = "Perennial", covars = c("PC1"),
                       family = gaussian)
mod_KQ2 <- fitNullModel(scanAnnot, outcome = "Perennial", covars = c("PC1","PC2"),
                        cov.mat = K, family = gaussian)
mod_Q2 <- fitNullModel(scanAnnot, outcome = "Perennial", covars = c("PC1", "PC2"),
                       family = gaussian)
mod_KQ3 <- fitNullModel(scanAnnot, outcome = "Perennial", covars = c("PC1","PC2", "PC3"),
                        cov.mat = K, family = gaussian)
mod_Q3 <- fitNullModel(scanAnnot, outcome = "Perennial", covars = c("PC1", "PC2", "PC3"),
                       family = gaussian)
mod_KQ4 <- fitNullModel(scanAnnot, outcome = "Perennial", covars = c("PC1","PC2", "PC3", "PC4"),
                        cov.mat = K, family = gaussian)
mod_Q4 <- fitNullModel(scanAnnot, outcome = "Perennial", covars = c("PC1","PC2", "PC3", "PC4"),
                       family = gaussian)

plot(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC),
     ylab = "AIC", xlab = "", xaxt="n")
axis(1, at=1:10, las="2",
     labels = c("intercept", "K", "K+PC1", "PC1", "K+PC1-2", "PC1-2",
                "K+PC1-3", "PC1-3", "K+PC1-4", "PC1-4"))
text(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC)+2,
     labels = round(c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
                      mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC), 1))

varcomp <- as.data.frame(varCompCI(mod_K))
h2s[6,2:4] <- unlist(varcomp[1,])

#Get genoIterator object dividing SNPs into groups of 10k
genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)

#Run association tests
assoc <- assocTestSingle(gdsobj=genoIterator, null.model=mod_Q1)
assoc$chr <- snps$CHROM #Provide correct chromosome names
pvalues[,6] <- assoc$Score.pval

#Perms
scrambled_pvals <- matrix(NA, nrow=n_permutations, ncol=nrow(snps))
for(i in 1:n_permutations){
  scanAnnot.scramble.df <- scanAnnot.df
  scramble_order <- sample(1:nrow(scanAnnot.scramble.df), nrow(scanAnnot.scramble.df), replace=F)
  scanAnnot.scramble.df[,2:ncol(scanAnnot.scramble.df)] <- scanAnnot.scramble.df[scramble_order,2:ncol(scanAnnot.scramble.df)]
  K_scramble <- K[scramble_order, scramble_order]
  scanAnnot.scramble <- ScanAnnotationDataFrame(scanAnnot.scramble.df)
  genoScramble <- GenotypeData(geno.gds, scanAnnot = scanAnnot.scramble)
  null_scramble <- fitNullModel(scanAnnot.scramble, outcome = "Perennial", covars = c("PC1"),
                                family = gaussian)
  genoIterator <- GenotypeBlockIterator(genoScramble,snpBlock = 10000)
  assocSCRAMBLE <- assocTestSingle(gdsobj=genoIterator, null.model=null_scramble)
  scrambled_pvals[i,] <- assocSCRAMBLE$Score.pval
}

permute_mins <- apply(scrambled_pvals,1,min)
pval_thresholds[6,2:3] <- c(quantile(permute_mins, .05), quantile(permute_mins, .01))


###########################   Run Red Knight   ##########################
print("Starting RK")
#Model Testing
mod_simple <- fitNullModel(scanAnnot, outcome = "RedKnight",
                           family = gaussian)
mod_K <- fitNullModel(scanAnnot, outcome = "RedKnight",
                      cov.mat = K, family = gaussian)
mod_KQ1 <- fitNullModel(scanAnnot, outcome = "RedKnight", covars = c("PC1"),
                        cov.mat = K, family = gaussian)
mod_Q1 <- fitNullModel(scanAnnot, outcome = "RedKnight", covars = c("PC1"),
                       family = gaussian)
mod_KQ2 <- fitNullModel(scanAnnot, outcome = "RedKnight", covars = c("PC1","PC2"),
                        cov.mat = K, family = gaussian)
mod_Q2 <- fitNullModel(scanAnnot, outcome = "RedKnight", covars = c("PC1", "PC2"),
                       family = gaussian)
mod_KQ3 <- fitNullModel(scanAnnot, outcome = "RedKnight", covars = c("PC1","PC2", "PC3"),
                        cov.mat = K, family = gaussian)
mod_Q3 <- fitNullModel(scanAnnot, outcome = "RedKnight", covars = c("PC1", "PC2", "PC3"),
                       family = gaussian)
mod_KQ4 <- fitNullModel(scanAnnot, outcome = "RedKnight", covars = c("PC1","PC2", "PC3", "PC4"),
                        cov.mat = K, family = gaussian)
mod_Q4 <- fitNullModel(scanAnnot, outcome = "RedKnight", covars = c("PC1","PC2", "PC3", "PC4"),
                       family = gaussian)

plot(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC),
     ylab = "AIC", xlab = "", xaxt="n")
axis(1, at=1:10, las="2",
     labels = c("intercept", "K", "K+PC1", "PC1", "K+PC1-2", "PC1-2",
                "K+PC1-3", "PC1-3", "K+PC1-4", "PC1-4"))
text(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC) + 1,
     labels = round(c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
                      mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC), 1))

varcomp <- as.data.frame(varCompCI(mod_K))
h2s[7,2:4] <- unlist(varcomp[1,])

#Get genoIterator object dividing SNPs into groups of 10k
genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)

#Run association tests
assoc <- assocTestSingle(gdsobj=genoIterator, null.model=mod_Q2)
assoc$chr <- snps$CHROM #Provide correct chromosome names
pvalues[,7] <- assoc$Score.pval

#Perms
scrambled_pvals <- matrix(NA, nrow=n_permutations, ncol=nrow(snps))
for(i in 1:n_permutations){
  scanAnnot.scramble.df <- scanAnnot.df
  scramble_order <- sample(1:nrow(scanAnnot.scramble.df), nrow(scanAnnot.scramble.df), replace=F)
  scanAnnot.scramble.df[,2:ncol(scanAnnot.scramble.df)] <- scanAnnot.scramble.df[scramble_order,2:ncol(scanAnnot.scramble.df)]
  K_scramble <- K[scramble_order, scramble_order]
  scanAnnot.scramble <- ScanAnnotationDataFrame(scanAnnot.scramble.df)
  genoScramble <- GenotypeData(geno.gds, scanAnnot = scanAnnot.scramble)
  null_scramble <- fitNullModel(scanAnnot.scramble, outcome = "RedKnight", covars = c("PC1", "PC2"),
                                family = gaussian)
  genoIterator <- GenotypeBlockIterator(genoScramble,snpBlock = 10000)
  assocSCRAMBLE <- assocTestSingle(gdsobj=genoIterator, null.model=null_scramble)
  scrambled_pvals[i,] <- assocSCRAMBLE$Score.pval
}

permute_mins <- apply(scrambled_pvals,1,min)
pval_thresholds[7,2:3] <- c(quantile(permute_mins, .05), quantile(permute_mins, .01))

###########################   Run Revolution  ##########################
print("Starting Revolution")
#Model Testing
mod_simple <- fitNullModel(scanAnnot, outcome = "Revolution",
                           family = gaussian)
mod_K <- fitNullModel(scanAnnot, outcome = "Revolution",
                      cov.mat = K, family = gaussian)
mod_KQ1 <- fitNullModel(scanAnnot, outcome = "Revolution", covars = c("PC1"),
                        cov.mat = K, family = gaussian)
mod_Q1 <- fitNullModel(scanAnnot, outcome = "Revolution", covars = c("PC1"),
                       family = gaussian)
mod_KQ2 <- fitNullModel(scanAnnot, outcome = "Revolution", covars = c("PC1","PC2"),
                        cov.mat = K, family = gaussian)
mod_Q2 <- fitNullModel(scanAnnot, outcome = "Revolution", covars = c("PC1", "PC2"),
                       family = gaussian)
mod_KQ3 <- fitNullModel(scanAnnot, outcome = "Revolution", covars = c("PC1","PC2", "PC3"),
                        cov.mat = K, family = gaussian)
mod_Q3 <- fitNullModel(scanAnnot, outcome = "Revolution", covars = c("PC1", "PC2", "PC3"),
                       family = gaussian)
mod_KQ4 <- fitNullModel(scanAnnot, outcome = "Revolution", covars = c("PC1","PC2", "PC3", "PC4"),
                        cov.mat = K, family = gaussian)
mod_Q4 <- fitNullModel(scanAnnot, outcome = "Revolution", covars = c("PC1","PC2", "PC3", "PC4"),
                       family = gaussian)

plot(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC),
     ylab = "AIC", xlab = "", xaxt="n")
axis(1, at=1:10, las="2",
     labels = c("intercept", "K", "K+PC1", "PC1", "K+PC1-2", "PC1-2",
                "K+PC1-3", "PC1-3", "K+PC1-4", "PC1-4"))
text(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC) + 2,
     labels = round(c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
                      mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC), 3))

varcomp <- as.data.frame(varCompCI(mod_K))
h2s[8,2:4] <- unlist(varcomp[1,])

#Get genoIterator object dividing SNPs into groups of 10k
genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)

#Run association tests
assoc <- assocTestSingle(gdsobj=genoIterator, null.model=mod_simple)
assoc$chr <- snps$CHROM #Provide correct chromosome names
pvalues[,8] <- assoc$Score.pval

#Perms
scrambled_pvals <- matrix(NA, nrow=n_permutations, ncol=nrow(snps))
for(i in 1:n_permutations){
  scanAnnot.scramble.df <- scanAnnot.df
  scramble_order <- sample(1:nrow(scanAnnot.scramble.df), nrow(scanAnnot.scramble.df), replace=F)
  scanAnnot.scramble.df[,2:ncol(scanAnnot.scramble.df)] <- scanAnnot.scramble.df[scramble_order,2:ncol(scanAnnot.scramble.df)]
  K_scramble <- K[scramble_order, scramble_order]
  scanAnnot.scramble <- ScanAnnotationDataFrame(scanAnnot.scramble.df)
  genoScramble <- GenotypeData(geno.gds, scanAnnot = scanAnnot.scramble)
  null_scramble <- fitNullModel(scanAnnot.scramble, outcome = "Revolution",
                                family = gaussian)
  genoIterator <- GenotypeBlockIterator(genoScramble,snpBlock = 10000)
  assocSCRAMBLE <- assocTestSingle(gdsobj=genoIterator, null.model=null_scramble)
  scrambled_pvals[i,] <- assocSCRAMBLE$Score.pval
}

permute_mins <- apply(scrambled_pvals,1,min)
pval_thresholds[8,2:3] <- c(quantile(permute_mins, .05), quantile(permute_mins, .01))


###########################   Run Vanguard  ##########################
print("Starting Vanguard")
#Model Testing
mod_simple <- fitNullModel(scanAnnot, outcome = "Vanguard",
                           family = gaussian)
mod_K <- fitNullModel(scanAnnot, outcome = "Vanguard",
                      cov.mat = K, family = gaussian)
mod_KQ1 <- fitNullModel(scanAnnot, outcome = "Vanguard", covars = c("PC1"),
                        cov.mat = K, family = gaussian)
mod_Q1 <- fitNullModel(scanAnnot, outcome = "Vanguard", covars = c("PC1"),
                       family = gaussian)
mod_KQ2 <- fitNullModel(scanAnnot, outcome = "Vanguard", covars = c("PC1","PC2"),
                        cov.mat = K, family = gaussian)
mod_Q2 <- fitNullModel(scanAnnot, outcome = "Vanguard", covars = c("PC1", "PC2"),
                       family = gaussian)
mod_KQ3 <- fitNullModel(scanAnnot, outcome = "Vanguard", covars = c("PC1","PC2", "PC3"),
                        cov.mat = K, family = gaussian)
mod_Q3 <- fitNullModel(scanAnnot, outcome = "Vanguard", covars = c("PC1", "PC2", "PC3"),
                       family = gaussian)
mod_KQ4 <- fitNullModel(scanAnnot, outcome = "Vanguard", covars = c("PC1","PC2", "PC3", "PC4"),
                        cov.mat = K, family = gaussian)
mod_Q4 <- fitNullModel(scanAnnot, outcome = "Vanguard", covars = c("PC1","PC2", "PC3", "PC4"),
                       family = gaussian)

plot(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC),
     ylab = "AIC", xlab = "", xaxt="n")
axis(1, at=1:10, las="2",
     labels = c("intercept", "K", "K+PC1", "PC1", "K+PC1-2", "PC1-2",
                "K+PC1-3", "PC1-3", "K+PC1-4", "PC1-4"))
text(1:10, c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
             mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC) + 2,
     labels = round(c(mod_simple$AIC, mod_K$AIC, mod_KQ1$AIC, mod_Q1$AIC, mod_KQ2$AIC, mod_Q2$AIC,
                      mod_KQ3$AIC, mod_Q3$AIC, mod_KQ4$AIC, mod_Q4$AIC), 3))

varcomp <- as.data.frame(varCompCI(mod_K))
h2s[9,2:4] <- unlist(varcomp[1,])

#Get genoIterator object dividing SNPs into groups of 10k
genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)

#Run association tests
assoc <- assocTestSingle(gdsobj=genoIterator, null.model=mod_Q3)
assoc$chr <- snps$CHROM #Provide correct chromosome names
pvalues[,9] <- assoc$Score.pval

#Perms
scrambled_pvals <- matrix(NA, nrow=n_permutations, ncol=nrow(snps))
for(i in 1:n_permutations){
  scanAnnot.scramble.df <- scanAnnot.df
  scramble_order <- sample(1:nrow(scanAnnot.scramble.df), nrow(scanAnnot.scramble.df), replace=F)
  scanAnnot.scramble.df[,2:ncol(scanAnnot.scramble.df)] <- scanAnnot.scramble.df[scramble_order,2:ncol(scanAnnot.scramble.df)]
  K_scramble <- K[scramble_order, scramble_order]
  scanAnnot.scramble <- ScanAnnotationDataFrame(scanAnnot.scramble.df)
  genoScramble <- GenotypeData(geno.gds, scanAnnot = scanAnnot.scramble)
  null_scramble <- fitNullModel(scanAnnot.scramble, outcome = "Vanguard", covars = c("PC1", "PC2", "PC3"),
                                family = gaussian)
  genoIterator <- GenotypeBlockIterator(genoScramble,snpBlock = 10000)
  assocSCRAMBLE <- assocTestSingle(gdsobj=genoIterator, null.model=null_scramble)
  scrambled_pvals[i,] <- assocSCRAMBLE$Score.pval
}

permute_mins <- apply(scrambled_pvals,1,min)
pval_thresholds[9,2:3] <- c(quantile(permute_mins, .05), quantile(permute_mins, .01))

########################################    Save data   #########################

write.csv(pval_thresholds, "data/empirical_pval_thresholds.csv", quote=F, row.names=F)
write.csv(h2s, "data/h2_estimates.csv", quote=F, row.names=F)
pvalues <- data.frame(cbind(snps), pvalues)
colnames(pvalues)[4:ncol(pvalues)] <- as.character(pval_thresholds$Pepper)
write.csv(pvalues, "data/pvalues.csv", quote=F, row.names=F)
