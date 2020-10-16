setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/pop_structure/")
library(pcaMethods)
library(RColorBrewer)
library(grid)

#####################################          Read data          ####################################

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

#Function to impute missing with mean
impute <- function(x){
  x[is.na(x)] <- mean(x,na.rm=T)
  return(x)
}
geno.imp <- apply(geno,2,impute)

#Calculate PCs and visualize screeplot
geno.pc <- prcomp(geno.imp, center = T, scale=T)
screeplot(geno.pc)
pcs.g <- geno.pc$x[,1:4]

#K-means clusters ---- first look at how within sum of square varies with increasing K 
#in order to choose number clusters
set.seed(5834785)
max.clusters <- 10
wss <- rep(NA,max.clusters)
wss[1] <- (nrow(pcs.g)-1)*sum(apply(pcs.g,2,var))
for (i in 2:15){
  km.i <- kmeans(pcs.g,centers=i,nstart=20,iter.max=1000)
  wss[i] <- sum(km.i$withinss)
}
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#Now identify clusters
n.clusters <- 5
pcs.km <- kmeans(pcs.g, n.clusters, nstart=20, iter.max=1000)
color_choices <- brewer.pal(n.clusters, "Set1")
color_assignments <- color_choices[match(pcs.km$cluster, 1:n.clusters)]

#Impute missing data in phenotypes and calculate phenotype PCs
pheno.imp <- apply(phenos, 2, impute)
pheno.pc <- prcomp(pheno.imp, center = T, scale=T)
pcs.p <- pheno.pc$x[,1:4]

#function to turn usr plot coordinates into relative plot coordinates to plot pane letters
get_coords <- function(x=-0.1,y=1.1){
  x.coord <- grconvertX(x, "npc", "user")
  y.coord <- grconvertY(y, "npc", "user")
  return(c(x.coord, y.coord))
}

pdf("plots/geno_and_phenos_pca.pdf")
old.par <- par(no.readonly = T)
m <- rbind(c(1,3),
           c(1,3),
           c(1,3),
           c(2,4),
           c(2,4),
           c(2,4),
           c(5,5))
layout(m)
par(mar=c(4,3.5,2,2.5), oma=c(0,2,1,0), xpd=NA)

plot(pcs.g[,1], pcs.g[,2],
     col = 'black', bg=color_assignments, pch=21, cex=1.5,
     xlab = paste("Genotype PC 1: ", round(geno.pc$sdev[1]^2/sum(geno.pc$sdev^2)*100,1), "%", sep=""),
     ylab = paste("Genotype PC 2: ", round(geno.pc$sdev[2]^2/sum(geno.pc$sdev^2)*100,1), "%", sep=""))
text(get_coords()[1], get_coords()[2], "A", cex=1.5)

plot(pcs.g[,3], pcs.g[,4],
     col = 'black', bg=color_assignments, pch=21, cex=1.5,
     xlab = paste("Genotype PC 3: ", round(geno.pc$sdev[3]^2/sum(geno.pc$sdev^2)*100,1), "%", sep=""),
     ylab = paste("Genotype PC 4: ", round(geno.pc$sdev[4]^2/sum(geno.pc$sdev^2)*100,1), "%", sep=""))
text(get_coords()[1], get_coords()[2], "B", cex=1.5)

plot(pcs.p[,1],
     pcs.p[,2],
     pch=21,cex=1.5,col='black',asp=1,
     bg=color_assignments,
     xlab = paste("Phenotype PC 1: ", round(pheno.pc$sdev[1]^2/sum(pheno.pc$sdev^2)*100,1), "%", sep=""),
     ylab = paste("Phenotype PC 2: ", round(pheno.pc$sdev[2]^2/sum(pheno.pc$sdev^2)*100,1), "%", sep=""))
text(get_coords()[1], get_coords()[2], "C", cex=1.5)

plot(pcs.p[,3],
     pcs.p[,4],
     pch=21,cex=1.5,col='black',asp=1,
     bg=color_assignments,
     xlab = paste("Phenotype PC 3: ", round(pheno.pc$sdev[3]^2/sum(pheno.pc$sdev^2)*100,1), "%", sep=""),
     ylab = paste("Phenotype PC 4: ", round(pheno.pc$sdev[4]^2/sum(pheno.pc$sdev^2)*100,1), "%", sep=""))
text(get_coords()[1], get_coords()[2], "D", cex=1.5)

#Add legend
par(mar=c(2,3.5,3.5,2.5))
plot(0,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='', xlim=c(1,10),ylim=c(1,10))
legend("center",
       fill=color_choices,
       legend=1:5,
       title = "Cluster",
       bty = 'n', cex=1.25,
       ncol=5, x.intersp=1, text.width=1)

par(old.par)
dev.off()

###################              ANOVA --- cluster vs phenotype       ######################################

clust <- as.factor(pcs.km$cluster)
p <- ncol(phenos)
pvals <- rep(NA,p)
for(i in 1:p){
  clust.lm <- lm(phenos[,i] ~ clust)
  clust.anova <- anova(clust.lm)
  pvals[i] <- clust.anova$`Pr(>F)`[1]
}
names(pvals) <- colnames(phenos)
pvals <- p.adjust(pvals, method="bonferroni")

#Write clusters
cluster_assignments <- data.frame("Isolate" = names(clust),
                                  "Cluster" = clust)
write.csv(cluster_assignments, "tables/cluster_assignments.csv", quote=F, row.names = F)

