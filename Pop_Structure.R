setwd("~/Documents/Cornell/Pepper_Interactions/paper/pop_structure/")
library(pcaMethods)
library(RColorBrewer)
library(grid)

#####################################          Read data          ####################################

#Read virulence BLUEs
phenos <- read.csv("../gwas/data/virulence_blues.csv")
rownames(phenos) <- phenos$Isolate
phenos$Isolate <- NULL

#Read assignment of isolates to fields/subpopulations
sites <- read.csv("../pheno/tables/renamed_field_metadata.csv")
sites$SampleSZ <- as.character(sites$SampleSZ)
sites$SampleSZ[sites$SampleSZ =="14_55C"] <- "14_55"

#Read geno
geno <- read.table("../geno/data/capsici_pepper_subset.012")
geno$V1 <- NULL #Get rid of row names
geno <- as.matrix(geno)
geno[geno==-1] <- NA
indvs <- read.table("../geno/data/capsici_pepper_subset.012.indv")
indvs <- unlist(indvs$V1)
rownames(geno) <- indvs

#Read finlay-wilkinson regression parameters
fw <- read.csv("../pheno/tables/fw_params.csv")
rownames(fw) <- fw$X
fw$X <- NULL
fw <- fw[match(rownames(geno), rownames(fw)),]

#Put phenos in same order as genos
phenos <- phenos[match(rownames(geno), rownames(phenos)),]

#Put sites in same order as genos
sites <- sites[match(rownames(geno), sites$SampleSZ),]

#Get rid of 'main' phenotype
main_effects <- phenos$main
phenos$main <- NULL

###################              geno PCA vs pheno PCA       ######################################

#Function to impute missing with mean
impute <- function(x){
  x[is.na(x)] <- mean(x,na.rm=T)
  return(x)
}
geno.imp <- apply(geno,2,impute)

#Calculate genotypic  PCs and visualize screeplot
geno.pc <- prcomp(geno.imp, center = T, scale=T)
screeplot(geno.pc)
pcs.g <- geno.pc$x[,1:4]

#K-means clusters ---- first look at how within sum of square varies with increasing K 
#in order to choose number clusters
set.seed(5834785)
max.clusters <- 15
wss <- rep(NA,max.clusters)
wss[1] <- (nrow(pcs.g)-1)*sum(apply(pcs.g,2,var))
for (i in 2:max.clusters){
  km.i <- kmeans(pcs.g,centers=i,nstart=20,iter.max=1000)
  wss[i] <- sum(km.i$withinss)
}
plot(1:max.clusters, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#Put screeplot and K-means cluster plot together in one supplemental plot
pdf("plots/geno_clustering.pdf", width=7, height=4)
get_coords <- function(x=-0.15,y=1.15){
  x.coord <- grconvertX(x, "npc", "user")
  y.coord <- grconvertY(y, "npc", "user")
  return(c(x.coord, y.coord))
}
old.par <- par(no.readonly = T)
par(mfrow=c(1,2), xpd=NA)
screeplot(geno.pc, type="lines", main="")
mtext("Principal component", side = 1, line=3)
text(get_coords()[1], get_coords()[2], "A", cex=1.5)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
text(get_coords()[1], get_coords()[2], "B", cex=1.5)
par(old.par)
dev.off()

#Now identify clusters
n.clusters <- 5
set.seed(1515)
pcs.km <- kmeans(pcs.g, n.clusters, nstart=20, iter.max=1000)
clust <- as.factor(pcs.km$cluster)
color_choices <- brewer.pal(n.clusters, "Set1")
color_assignments <- color_choices[match(clust, 1:n.clusters)]

#How do clusters relate to field sites
sites$cluster <- clust[match(sites$SampleSZ, names(clust))]
sites <- sites[order(sites$cluster),]
clusters_by_sites <- table(sites$Field, sites$cluster)
write.csv(clusters_by_sites, "tables/clusters_by_sites.csv", quote=F)

#Impute missing data in phenotypes and calculate phenotype PCs
pheno.imp <- apply(phenos, 2, impute)
pheno.pc <- prcomp(pheno.imp, center = T, scale=T)
pcs.p <- pheno.pc$x[,1:4]

#Save these PCs
write.csv(pcs.p, "tables/phenotypic_pcs.csv", quote=F, row.names = T)

#Plot loadings
pdf("plots/pepper_loadings.pdf", height=7, width=4)
old.par <- par(no.readonly=T)
par(mfrow=c(4,1), mar=c(2,5,2,2), oma=c(6,0,0,0), xpd=NA)
for(i in 1:4){
  barplot(pheno.pc$rotation[,i],
          names.arg=rep("",nrow(pheno.pc$rotation)),
          ylab="Loading")
  text(get_coords()[1], get_coords()[2], LETTERS[i], cex=1.5)
  if(i == 4){
    barplot.data <- barplot(pheno.pc$rotation[,i], plot=F)
    axis(1, at=barplot.data, labels=rownames(pheno.pc$rotation),las=2)
  }
}
par(old.par)
dev.off()

#What are the phenotypic PCS?
#PC 1 relates to general virulence
cor.test(pcs.p[,1], main_effects)

#PC 2 relates to FW slope
cor.test(pcs.p[,2], fw$Slope)

#PC 3 harder to define
pheno.pc$rotation

#######################      Make plot of genotypic and phenotypic PCs      ###################

#function to turn usr plot coordinates into relative plot coordinates to plot pane letters
get_coords <- function(x=-0.1,y=1.1){
  x.coord <- grconvertX(x, "npc", "user")
  y.coord <- grconvertY(y, "npc", "user")
  return(c(x.coord, y.coord))
}

pdf("plots/geno_and_phenos_pca.pdf", width=7, height=7)
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

phenos.all <- cbind(phenos, pcs.p)
p <- ncol(phenos.all)
cluster_results <- data.frame("Trait" = colnames(phenos.all),
                              "R2" = NA,
                              "p"  = NA,
                              "p.adjust" = NA)
pvals <- rep(NA,p)
for(i in 1:p){
  clust.lm <- lm(phenos.all[,i] ~ clust)
  clust.anova <- anova(clust.lm)
  cluster_results$R2[i] <- summary(clust.lm)$"r.squared"*100
  cluster_results$p[i] <- clust.anova$"Pr(>F)"[1]
}

cluster_results$p.adjust <- p.adjust(cluster_results$p, method="bonferroni")


#Write clusters
cluster_assignments <- data.frame("Isolate" = names(clust),
                                  "Cluster" = clust)
write.csv(cluster_assignments, "tables/cluster_assignments.csv", quote=F, row.names = F)

#Write results of regressing phenotypes on cluster
cluster_results.clean <- cluster_results
cluster_results.clean[,3:4] <- apply(cluster_results[,3:4], 2, signif, digits=2)
cluster_results.clean$R2 <- round(cluster_results$R2, digits = 2)
write.csv(cluster_results.clean, "tables/cluster_tests.csv", quote=F, row.names=F)
