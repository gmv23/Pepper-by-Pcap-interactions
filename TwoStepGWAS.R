#setwd("/workdir/gmv23/peppers/pheno")
setwd("~/Documents/Cornell/Pepper_Interactions/paper/gwas/")
library(rrBLUP)
library(GENESIS)
library(GWASTools)
library(SNPRelate)

#######################################          Load and clean data       ######################################

#Reformat VCF file as GDS file and import
snpgdsVCF2GDS(vcf.fn="../geno/data/capsici_pepper_subset.vcf", out.fn="data/capsici_pepper_subset.gds")
geno.gds <- GdsGenotypeReader("data/capsici_pepper_subset.gds")

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
  n.models <- length(model_names)
  if(plot){
    plot(1:n.models, model_AICs,
         ylab = "AIC", xlab = "", xaxt="n", col='white', 
         main = outcome)
    axis(1, at=1:n.models, las="2", labels = model_names)
    text(1:n.models, model_AICs, labels=round(model_AICs,2))
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

#Save transformed phenos
write.csv(phenos_transformed, "tables/transformed_blues.csv", quote=F, row.names=T)

###########################################        RUN GWAS        ########################################

models <- rep(NA, ncol(phenos_transformed))
names(models) <- colnames(phenos_transformed)

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

#### Get and write SNP p-values, R2s, and additive allelic effects ###

#Which is the significant SNP we identified
sig_snp <- which(pvals$RedKnight == min(pvals$RedKnight))
sig_geno <- geno[,sig_snp]

p.summary <- data.frame("Pepper" = colnames(pvals),
                        "P" = unlist(pvals[sig_snp,]),
                        "R2" = NA,
                        "Allelic_effect" = NA)

for(i in 1:nrow(p.summary)){
  mod <- lm(phenos_transformed[,i] ~ as.integer(as.character(sig_geno)))
  p.summary$R2[i] <- summary(mod)$r.squared
  effect <- mod$coefficients[2]
  p.summary$Allelic_effect[i] <- effect
}

p.summary$P <- sapply(p.summary$P, signif, digits=2)
p.summary$R2 <- sapply(p.summary$R2, round, digits=2)
p.summary$Allelic_effect <- sapply(p.summary$Allelic_effect, round, digits=2)

write.csv(p.summary, "tables/pval_summary.csv", quote=F, row.names = F)

###########################################        Random calculations       ########################################

#This SNP is top X % in each of the traits
sig_ps <- unlist(pvals[sig_snp,])
percentiles <- rep(NA, length(sig_ps))
names(percentiles) = colnames(pvals)
for(i in 1:ncol(pvals)){
  percentiles[i] <- sum(pvals[,i] <= sig_ps[i])/nrow(pvals)
}
percentiles*100

#MAF
sig_calls <- sig_geno[!is.na(sig_geno)]
(sum(sig_calls==1)*2 + sum(sig_calls=2)) / (2*length(sig_calls))

