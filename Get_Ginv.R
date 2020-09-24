##Filter genotype data for isolates in experiment and MAF > 0.05
##Create relationship matrix and inverse relationship matrix in asreml-readable format

setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/pheno/")
library(rrBLUP)

#Load geno data
geno <- read.table("../../../isolate_collection/paper/pop_structure/capsici_diversity_PG.012")
indvs <- read.table("../../../isolate_collection/paper/pop_structure/capsici_diversity_PG.012.indv", stringsAsFactors = F)
snps <- read.table("../../../isolate_collection/paper/pop_structure/capsici_diversity_PG.012.pos")

#Load pheno data
pep <- read.csv("data/Ratings_filt.csv")

#Clean up data
indvs <-unlist(indvs$V1)
indvs <- sapply(indvs, function(x) unlist(strsplit(x,":"))[1])

geno <- geno[,-1]
geno <- as.matrix(geno)
geno[geno==-1] <- NA

#Genotype for 17PZ18A in the CC set, but inoculated with 17PZ21A
#For 14-55, genotype is for 14_55C but inoculated with 14_55
#Theyre clones so rename
indvs[indvs=="17PZ18A"] <- "17PZ21A"
indvs[indvs=="14_55C"] <- "14_55"

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

K <- A.mat(geno-1)
rownames(K) <- indvs
colnames(K) <- indvs
colnames(geno_gwas)[4:ncol(geno_gwas)] <- indvs

Get_ginv <- function(mat, digits=10){
  require(MASS)
  ################################################################
  # This function takes a relationship matrix                    #
  # and creates 3-column data frame with non-zero elements of    #
  # lower triangle of inverse relationhsip matrix                #
  # Code is based off Chris Hernandez script write.asremlMat     #
  ################################################################
  #Get inverse
  mat.inv <- solve(mat)
  # Round matrix to a certain number of digits to decrease file size
  mat.inv <- round(mat.inv, digits) 
  # Reformat matrix into three column format representing
  # the lower triangular portion of the matrix and diagnal 
  # elements.
  lg.mat <- data.frame(Row = rep(1:nrow(mat.inv), nrow(mat.inv)),
                       Column = rep(1:nrow(mat.inv), each=nrow(mat.inv)),
                       coeff = as.numeric(mat.inv),
                       lower = as.logical(lower.tri(mat.inv, diag=T)))
  # Only need lower triangular
  lg.mat <- lg.mat[lg.mat$lower, c("Row", "Column", "coeff")]
  # Sort by row and column within row. This is an 
  # ASReml requirement.
  lg.mat <- lg.mat[order(lg.mat$Row, lg.mat$Column),]
  # Remove entries that are zero on the off diag.
  # ASReml assumes that non diagonal entries ommitted are
  # zero, so removing them is a form of compression.
  lg.mat <- lg.mat[!lg.mat$coeff == 0,]
  #Create row name attribute
  K.rownames <- row.names(mat)
  names(K.rownames) <- NULL
  attr(lg.mat, "rowNames") <- K.rownames
  # Return data frame
  return(lg.mat)
}

K.ginv <- Get_ginv(K)

saveRDS(K.ginv, "data/K_ginv.rds")
