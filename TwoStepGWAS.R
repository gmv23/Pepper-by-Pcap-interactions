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

#Read phenotypic PCs and add to phenos
pheno.pca <- read.csv("../pop_structure/tables/phenotypic_pcs.csv")
rownames(pheno.pca) <- pheno.pca$X
pheno.pca$X <- NULL
pheno.pca <- pheno.pca[match(rownames(phenos), rownames(pheno.pca)),]
colnames(pheno.pca) <- paste("Pheno", colnames(pheno.pca), sep="")
phenos <- cbind(phenos, pheno.pca)

#Put all datasets in VCF sample order
sample_order <- getScanID(geno.gds)
phenos <- phenos[match(sample_order, rownames(phenos)),]
geno <- geno[match(sample_order, rownames(geno)),]
geno.pca <- geno.pca[match(sample_order, rownames(geno.pca)),]

###################################        Functions for GWAS        ########################################

#######  Perform model testing in GENESIS with K and PCs #######

test_models <- function(scanAnnot, outcome, plot=T, K){
  
  covars.list <- list(NULL,
                c("PC1"),
                c("PC1", "PC2"),
                c("PC1", "PC2", "PC3"),
                c("PC1", "PC2", "PC3", "PC4"))
  
  cov.mat.list <- list(NULL,K)
  
  results <- list()
  counter <- 0
  for(i in 1:length(covars.list)){
    for(j in 1:length(cov.mat.list)){
      counter <- counter + 1
      covars <- covars.list[[i]]
      cov.mat <- cov.mat.list[[j]]
      model <- fitNullModel(x=scanAnnot, outcome = outcome, family = gaussian,
                            covars = covars, cov.mat = cov.mat)
      results[counter] <- list(list(covars = covars, cov.mat = cov.mat, model = model))                    
    }
  }
  
  aics <- data.frame("Model" = 1:counter,
                     "Covars" = sapply(results, function(x) paste(x$covars, collapse = " ")),
                     "Covar_matrix" = sapply(results, function(x) !is.null(x$cov.mat)),
                     "AIC" = sapply(results, function(x) x$model$AIC))
  if(plot){
    pdf(paste("plots/model_testing_", outcome, ".pdf", sep=""))
    old.par <- par(no.readonly = T)
    par(mar=c(12,4,2,2))
    plot(1:counter, aics$AIC,
         ylab = "AIC", xlab = "", col='white', xaxt='n',
         main = outcome)
    axis(1, at=1:counter, las="2", labels = paste(aics$Covars, aics$Covar_matrix), cex=0.75)
    text(1:counter, aics$AIC, labels=round(aics$AIC, 3))
    par(old.par)
    dev.off()
  }
  
  model_choose <- which(aics$AIC == min(aics$AIC))
  if(is.null(results[[model_choose]]$cov.mat)){
    res <- results[[model_choose]]$model$resid.marginal
  }else{
    res <- results[[model_choose]]$model$resid.conditional
  }
  return(list(covars = results[[model_choose]]$covars,
              cov.mat = results[[model_choose]]$cov.mat,
              res = res))
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
  return(threshold)
}

#######  Shift distribution and log transform #######

transform <- function(x){
  x.min <- min(x, na.rm=T)
  if(x.min < 1){
    x.shift <- x + (1-x.min)
  }else{
    x.shift <- x
  }
  return(log(x.shift))
}

###################################        Model test and Run GWAS        ########################################

#Get genomic relationship matrix
K <- A.mat(geno-1)
rownames(K) <- rownames(phenos)
colnames(K) <- rownames(phenos)

#Make starting scanAnnot object
scanAnnot.df <- data.frame("scanID" = rownames(phenos),
                           "PC1" = geno.pca[,1],
                           "PC2" = geno.pca[,2],
                           "PC3" = geno.pca[,3],
                           "PC4" = geno.pca[,4],
                           phenos)
scanAnnot <- ScanAnnotationDataFrame(scanAnnot.df)

#Make variables used in loop or to store results from looping through traits
p <- ncol(phenos)
model_summaries <- data.frame("Trait" = colnames(phenos),
                              "Covariates" = NA,
                              "K" = NA,
                              "Transformed" = NA)
pvals <- matrix(NA, nrow=nrow(snps), ncol=p)
colnames(pvals) <- colnames(phenos)
pvals <- as.data.frame(pvals)

for(i in 1:p){
  trait <- colnames(phenos)[i]
  
  #Find best covariates
  test_out <- test_models(scanAnnot = scanAnnot, outcome = trait, plot = T, K=K)
  hist(test_out$res, main=trait)
  
  #Test if residuals are normal
  shap.test <- shapiro.test(test_out$res)
  
  transformed <- FALSE
  if(shap.test$p.value < 0.05){
    scanAnnot.df.test <- scanAnnot.df
    
    #If residuals non normal, transform phenotype and compare normality of residuals 
    scanAnnot.df.test[,trait] <- transform(scanAnnot.df.test[,trait])
    scanAnnot.test <- ScanAnnotationDataFrame(scanAnnot.df.test)
    log_out <- fitNullModel(x=scanAnnot.test, outcome = trait, family = gaussian,
                 covars = test_out$covars, cov.mat = test_out$cov.mat)
   
     if(is.null(test_out$cov.mat)){
      log_out_res <- log_out$resid.marginal
    }else{
      log_out_res <- log_out$resid.conditional
    }
    shap.log <- shapiro.test(log_out_res)
    if(shap.log$statistic > shap.test$statistic){
      transformed <- TRUE
      scanAnnot.df[,trait] <- scanAnnot.df.test[,trait]
    }else{
      transformed <- FALSE
    }
    scanAnnot <- ScanAnnotationDataFrame(scanAnnot.df)
  }
  
  #Find best model again, using 'correct' phenotype (transformed or non transformed)
  best_model <- test_models(scanAnnot = scanAnnot, outcome = trait, plot = T, K=K)
  model_summaries[i,2:4] <- c(paste(best_model$covars, collapse=" "),
                      !is.null(best_model$cov.mat),
                      transformed)
  
  #Run GWAS
  null_model <- fitNullModel(x=scanAnnot, outcome = trait, family = gaussian,
                             covars = best_model$covars, cov.mat = best_model$cov.mat)
  genoData <- GenotypeData(geno.gds, scanAnnot = scanAnnot)
  genoIterator <- GenotypeBlockIterator(genoData,snpBlock = 10000)
  assoc <- assocTestSingle(gdsobj=genoIterator, null.model=null_model)
  pvals[,i] <- assoc$Score.pval
}

close(geno.gds)

###########################################        Save results        ########################################

write.csv(model_summaries, "tables/gwas_models.csv", quote=F, row.names=F)

pvals_write <- cbind(snps, pvals)
write.csv(pvals_write, "data/gwas_pvalues.csv", quote=F, row.names=F)



