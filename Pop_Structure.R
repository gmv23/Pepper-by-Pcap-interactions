setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/pop_structure/")
library(pcaMethods)
library(RColorBrewer)

#Read virulence BLUEs
phenos_full <- read.csv("../gwas/data/virulence_blues.csv")
rownames(phenos_full) <- phenos_full$Isolate
phenos_full$Isolate <- NULL
#Get rid of 'main' phenotype
phenos <- phenos_full
phenos$main <- NULL

#Read assignment of isolates to fields/subpopulations
pops <- read.csv("../pheno/tables/pop_assignments.csv")

#Read geno
geno <- read.table("../geno/data/capsici_pepper_subset.012")
geno$V1 <- NULL #Get rid of row names
geno <- as.matrix(geno)
geno[geno==-1] <- NA
indvs <- read.table("../geno/data/capsici_pepper_subset.012.indv")
indvs <- unlist(indvs$V1)
rownames(geno) <- indvs

#Put phenos in same order as genos
phenos <- phenos[match(rownames(geno), rownames(phenos)),]

###################              geno PCA vs PCA-biplot       ######################################

#Calculate PCs
#Impute missing with mean
impute <- function(x){
  x[is.na(x)] <- mean(x,na.rm=T)
  return(x)
}
geno.imp <- apply(geno,2,impute)
geno.pc <-prcomp(geno.imp, center = T, scale=T)
screeplot(geno.pc)
pcs.g <- geno.pc$x[,1:4]

max.clusters <- 10
wss <- rep(NA,max.clusters)
wss[1] <- (nrow(pcs.g)-1)*sum(apply(pcs.g,2,var))
for (i in 2:15){
  km.i <- kmeans(pcs.g,centers=i,nstart=20,iter.max=1000)
  wss[i] <- sum(km.i$withinss)
}
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

n.clusters <- 5
pcs.km <- kmeans(pcs.g, n.clusters, nstart=20, iter.max=1000)
color_choices <- brewer.pal(n.clusters, "Set1")
color_assignments <- color_choices[match(pcs.km$cluster, 1:n.clusters)]

pheno.imp <- apply(phenos, 2, impute)
pheno.pc <- prcomp(pheno.imp, center = T, scale=T)
pcs.p <- pheno.pc$x[,1:2]

pdf("plots/geno_and_phenos_pca.pdf")
old.par <- par(no.readonly = T)
m <- cbind(c(1,2),c(3,3))
layout(m)
plot(pcs.g[,1], pcs.g[,2],
     col = 'black', bg=color_assignments, pch=21, cex=1.5,
     xlab = paste("PC1: ", round(geno.pc$sdev[1]^2/sum(geno.pc$sdev^2)*100,1), "%", sep=""),
     ylab = paste("PC2: ", round(geno.pc$sdev[2]^2/sum(geno.pc$sdev^2)*100,1), "%", sep=""))
plot(pcs.g[,3], pcs.g[,4],
     col = 'black', bg=color_assignments, pch=21, cex=1.5,
     xlab = paste("PC3: ", round(geno.pc$sdev[3]^2/sum(geno.pc$sdev^2)*100,1), "%", sep=""),
     ylab = paste("PC4: ", round(geno.pc$sdev[4]^2/sum(geno.pc$sdev^2)*100,1), "%", sep=""))

plot(pcs.p[,1],
     pcs.p[,2],
     pch=21,cex=1.5,col='black',asp=1,
     bg=color_assignments,
     xlab = paste("PC1: ", round(pheno.pc$sdev[1]^2/sum(pheno.pc$sdev^2)*100,1), "%", sep=""),
     ylab = paste("PC2: ", round(pheno.pc$sdev[2]^2/sum(pheno.pc$sdev^2)*100,1), "%", sep=""))
segments(x0=rep(0,ncol(phenos)),y0=rep(0,ncol(phenos)), 
       x1=pheno.pc$rotation[,1]*pheno.pc$sdev[1]*4, 
       y1=pheno.pc$rotation[,2]*pheno.pc$sdev[2]*4)
text(pheno.pc$rotation[,1]*3*66.9/14.8, pheno.pc$rotation[,2]*5, rownames(pheno.pc$rotation), cex=0.8)
par(old.par)
dev.off()



