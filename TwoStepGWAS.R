#setwd("/workdir/gmv23/peppers/pheno")
setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/gwas/")
library(rrBLUP)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(matrixcalc)

<<<<<<< HEAD
#################################   Load and clean data   ################################
set.seed(385832)
#Reformat VCF file as GDS file
snpgdsVCF2GDS(vcf.fn="data/capsici_pepper_subset.recode.vcf", out.fn="capsici_pepper_subset.gds")
=======
#######################################          Load and clean data       ######################################

#Reformat VCF file as GDS file and import
snpgdsVCF2GDS(vcf.fn="../geno/data/capsici_pepper_subset.vcf", out.fn="data/capsici_pepper_subset.gds")
>>>>>>> 493a77abfcc787459c36342ec063f0e9399f1fa1

#Read geno as 012 file
geno <- read.table("../geno/data/capsici_pepper_subset.012")
geno$V1 <- NULL #Get rid of row names
geno <- as.matrix(geno)
geno[geno==-1] <- NA

#Read marker positions
snps <- read.table("../geno/data/capsici_pepper_subset.012.pos")
colnames(snps) <- c("CHROM", "BP")
snps$MARKER <- paste(snps$CHROM, snps$BP, sep="_")

#Read sample names
indvs <- read.table("../geno/data/capsici_pepper_subset.012.indv")
indvs <- as.character(unlist(indvs$V1))
rownames(geno) <- indvs

#Read virulence BLUEs
phenos <- read.csv("data/virulence_blues.csv")
rownames(phenos) <- phenos$Isolate
phenos$Isolate <- NULL

#Read marker PCs
geno.pca <- read.csv("../geno/data/pcs.csv")
rownames(geno.pca) <- geno.pca$X
geno.pca$X <- NULL

#Put all datasets in VCF sample order
sample_order <- getScanID(geno.gds)
phenos <- phenos[match(sample_order, rownames(phenos)),]
geno <- geno[match(sample_order, rownames(geno)),]
geno.pca <- geno.pca[match(sample_order, rownames(geno.pca)),]

###################################        Functions for GWAS        ########################################

#######  Test if log transform improves normality #######

try_log <- function(x){
  
  #Shift distribution so minimum value = 1
  x.min <- min(x, na.rm=T)
  if(x.min < 1){
    raw <- x + (1-x.min)
  }else{
    raw <- x
  }
  
    transform <- log(raw)
    
    #Shapiro-Wilk test
    raw.W <- shapiro.test(raw)$statistic
    transform.W <- shapiro.test(transform)$statistic
   
    #Return data and whether or not it was transformed
    if(raw.W < transform.W){
      return(list("data" = transform, "transform" = TRUE))
    }else{
      return(list("data" = x, "transform" = FALSE))
    }
}

#######  Perform model testing in GENESIS with K and PCs #######

test_models <- function(scanAnnot, outcome, plot=T){
  
  models <- list()
  
  models$simple <- fitNullModel(x=scanAnnot, outcome = outcome,
                             family = gaussian)
  models$K <- fitNullModel(x=scanAnnot, outcome = outcome,
                        cov.mat = K, family = gaussian)
  models$KQ1 <- fitNullModel(x=scanAnnot, outcome = outcome, covars = c("PC1"),
                          cov.mat = K, family = gaussian)
  models$Q1 <- fitNullModel(x=scanAnnot, outcome = outcome, covars = c("PC1"),
                         family = gaussian)
  models$KQ2 <- fitNullModel(x=scanAnnot, outcome = outcome, covars = c("PC1","PC2"),
                          cov.mat = K, family = gaussian)
  models$Q2 <- fitNullModel(x=scanAnnot, outcome = outcome, covars = c("PC1", "PC2"),
                         family = gaussian)
  models$KQ3 <- fitNullModel(x=scanAnnot, outcome = outcome, covars = c("PC1","PC2", "PC3"),
                          cov.mat = K, family = gaussian)
  models$Q3 <- fitNullModel(x=scanAnnot, outcome = outcome, covars = c("PC1", "PC2", "PC3"),
                         family = gaussian)
  models$KQ4 <- fitNullModel(x=scanAnnot, outcome = outcome, covars = c("PC1","PC2", "PC3", "PC4"),
                          cov.mat = K, family = gaussian)
  models$Q4 <- fitNullModel(x=scanAnnot, outcome = outcome, covars = c("PC1","PC2", "PC3", "PC4"),
                         family = gaussian)
  
  model_names <- names(models)
  model_AICs <- sapply(models, function(x) x$AIC)
  
  if(plot){
    plot(1:10, model_AICs,
         ylab = "AIC", xlab = "", xaxt="n", col='white', 
         main = outcome)
    axis(1, at=1:10, las="2", labels = model_names)
    text(1:10, model_AICs, labels=round(model_AICs,2))
  }
  
  model_choose <- which(model_AICs == min(model_AICs))
  return(list(name=model_names[model_choose],
              model=models[[model_choose]]))
  
}

#######  Find FDR Threshold #######

fdr_cutoff <- function(x, alpha){
  x.sort <- sort(x)
  x.adjust <- p.adjust(x.sort, method="fdr")
  significant <- any(x.adjust <= alpha)
  if(significant){
    threshold <- x.sort[max(which(x.adjust < alpha))]
  }else{
    threshold <- NA
  }
}

###################################        Put together phenotypes and covariates        ########################################

#Make new phenos data frame with transformed phenos
transformations <- rep(NA, ncol(phenos))
names(transformations) <- colnames(phenos)
phenos_transformed <- phenos
for(i in 1:ncol(phenos_transformed)){
  transform_test <- try_log(phenos[,i])
  phenos_transformed[,i] <- transform_test$data
  transformations[i] <- transform_test$transform
  hist(transform_test$data, 
       main=paste(colnames(phenos_transformed)[i],transform_test$transform, sep=": "))
}

#Make ScanAnnotationDataFrame object with phenotypes and covariates
scanAnnot.df <- data.frame("scanID" = rownames(phenos),
                           "PC1" = geno.pca[,1],
                           "PC2" = geno.pca[,2],
                           "PC3" = geno.pca[,3],
                           "PC4" = geno.pca[,4],
                           phenos_transformed)

scanAnnot <- ScanAnnotationDataFrame(scanAnnot.df)

#Make GenotypeData object
genoData <- GenotypeData(geno.gds, scanAnnot = scanAnnot)

#Get genomic relationship matrix
K <- A.mat(geno-1)
rownames(K) <- scanAnnot.df$scanID
colnames(K) <- scanAnnot.df$scanID

<<<<<<< HEAD
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
=======
###########################################        RUN GWAS        ########################################
>>>>>>> 493a77abfcc787459c36342ec063f0e9399f1fa1

models <- rep(NA, ncol(phenos_transformed))
names(models) <- colnames(phenos_transformed)

<<<<<<< HEAD
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
=======
pvals <- matrix(NA, nrow=nrow(snps), ncol=ncol(phenos_transformed))
colnames(pvals) <- colnames(phenos_transformed)
pvals <- as.data.frame(pvals)

for(i in 1:ncol(phenos_transformed)){
  trait <- colnames(phenos_transformed)[i]
  best_model <- test_models(scanAnnot=scanAnnot, outcome=trait, plot=T)
  null_model <- best_model$model
  models[i] <- best_model$name
  
  genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)
  assoc <- assocTestSingle(gdsobj=genoIterator, null.model=null_model)
  pvals[,i] <- assoc$Score.pval
}

close(geno.gds)

###########################################        Save results        ########################################

gwas_models <- data.frame("Trait" = colnames(phenos_transformed),
                          "Log_transform" = transformations,
                          "model" = models)
rownames(gwas_models) <- NULL
write.csv(gwas_models, "tables/gwas_models.csv", quote=F, row.names=F)

pvals_write <- cbind(snps, pvals)
write.csv(pvals_write, "data/gwas_pvalues.csv", quote=F, row.names=F)

###########################################        Plot results        ########################################

#Function to draw stacked manhattan and qq plots
draw_plots <- function(pvals.plot, name){
  
  require(qqman)  
  
  n.plots <- ncol(pvals.plot)
  jpeg(name, width=7, height=(1.5*n.plots), units="in",res=100)
  
  old.par <- par(no.readonly = T)
  par(mar=c(6.5,4.5,1,1.5), oma=c(0,1,1.5,1))

  #Get layout
  m <- c(rep(seq(1,((n.plots*2)-1),by=2),4),
         seq(2,(n.plots*2),by=2))
  m <- sort(m)
  m <- matrix(m, nrow=n.plots, ncol=5, byrow=T)
  layout(m)
  
  #Plot
  for(i in 1:n.plots){
    trait <- colnames(pvals.plot)[i]
    man.df <- data.frame("BP"=snps$BP, "CHR"=snps$CHROM,"P"=pvals.plot[,i])
    threshold <- -log10(fdr_cutoff(man.df$P,alpha=0.10))
    if(is.na(threshold)){
      genomewideline <- FALSE
    }else{
      genomewideline <- threshold
    }
    if(i < n.plots){
      xlab.man <- ""
      xlab.qq <- ""
    }else{
      xlab.man <- "Scaffold"
      xlab.qq <- expression(paste("Observed -log"[10], "(",italic("p"), ")", sep=""))
    }
    manhattan(man.df,main=trait,suggestiveline=FALSE,
              genomewideline=genomewideline,
              xlab = xlab.man,
              ylim=c(0,-log10(min(pvals.plot))*1.1))
    qq(man.df$P, main=trait, xlab=xlab.qq, ylab="")
    mtext("Expected", side=2, line=3.5, cex=0.65)
    mtext(expression(paste("-log"[10], "(",italic("p"),")")), side=2, line=2.25, cex=0.65)
  }
  
  par(old.par)
  
  dev.off()
  
>>>>>>> 493a77abfcc787459c36342ec063f0e9399f1fa1
}

#Move 'main' to beginning of pvals
pvals.order <- cbind("Across-pepper" = pvals$main, pvals)
pvals.order$main <- NULL

#Separate into peppers with significant hits and those without
thresholds <- apply(pvals.order, 2, fdr_cutoff, alpha = 0.10)
pvals.sig <- pvals.order[,!is.na(thresholds)]
pvals.insig <- pvals.order[,is.na(thresholds)]

<<<<<<< HEAD
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
=======
#Make plots

draw_plots(pvals.sig, "plots/GWAS_sig.jpeg")
draw_plots(pvals.insig, "plots/GWAS_insig.jpeg")
>>>>>>> 493a77abfcc787459c36342ec063f0e9399f1fa1



