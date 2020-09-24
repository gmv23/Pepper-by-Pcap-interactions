setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/pheno/")

library(reshape2)
library(RColorBrewer)

#Import data
blups <- read.csv("tables/isolate_blups.csv")
colnames(blups) <- gsub(".", "-", colnames(blups), fixed=T)
pep <- read.csv("data/Ratings_filt.csv")
pep$Rep <- as.factor(pep$Rep)
pep$Block <- as.factor(pep$Block)
pep$Tray <- as.factor(pep$Tray)
pep$Column <- as.factor(pep$Column)

###### Make interaction plot

#First rank peppers by their mean disease ratings
interact_means.long <- aggregate(audpc~Isolate*Pepper,data=pep, FUN=mean)
interact_means.wide <- dcast(interact_means.long, Isolate ~ Pepper, value.var = 'audpc')
mean_disease <- apply(interact_means.wide[,-1], 2, mean, na.rm=T)
mean_disease <- sort(mean_disease, decreasing=T)

blups.sort <- blups[,names(mean_disease)]
pdf("plots/interaction_matrix.pdf")
old.par <- par(no.readonly = T)
par(mfrow=c(11,11), oma=c(1,1,0.5,0.5), mar=c(0,0,0,0))
for(i in 1:nrow(blups.sort)){
  isolate_highlight <- blups$Isolate[i]
  plot(0, type='n', 
       xlim=c(1,length(mean_disease)), 
       ylim=range(blups.sort, na.rm=T),
       xaxt='n', yaxt='n', ylab='', xlab='')
  pep.background <- blups.sort[blups$Isolate!=isolate_highlight,]
  for(j in 1:nrow(blups.sort)){
    lines(1:length(mean_disease), 
          pep.background[j,], 
          main = isolate_highlight, lwd=0.5,
          col=adjustcolor("gray", alpha.f=0.6))
  }
  lines(1:length(mean_disease),
        apply(blups.sort,2,mean,na.rm=T),
        col='black', lwd=1, lty=2)
  lines(1:length(mean_disease),
        blups.sort[i,],
        col='red', lwd=1.5)
  mtext(isolate_highlight,side=3,line=-1.5, cex=0.7)
}
dev.off()

library(pcaMethods)
blups.pc <- pca(blups,nPcs = 4, method = "nipals")

rownames(blups.sort) <- blups$Isolate
heatmap(as.matrix(blups.sort))


isolates.d <- dist(as.matrix(interact_means.wide[,-1]))
peppers.d <- dist(t(as.matrix(interact_means.wide[,-1])))
isolates.clust <- hclust(isolates.d)
plot(isolates.clust)

peppers.clust <- hclust(peppers.d)

pdf("~/Downloads/isolate_dendrogram.pdf")
old.par <- par(no.readonly = T)
par(cex=0.5)
plot(isolates.clust, hang=-1)
par(old.par)
dev.off()

pdf("~/Downloads/pepper_dendrogram.pdf")
old.par <- par(no.readonly = T)
par(cex=0.5)
plot(peppers.clust, hang=-1)
par(old.par)
dev.off()

library(vioplot)

pdf("means_vioplot.pdf")
old.par <- par(no.readonly = T)
vioplot(interact_means.wide[,colnames(blups.sort)])
par(old.par)
dev.off()

vioplot(blups.sort)


isolate.assignments <- rect.hclust(isolates.clust, h=50)
isolates.clust <- hclust(isolates.d)

pdf("plots/dendrogram_cut.pdf")
plot(isolates.clust)
rect.hclust(isolates.clust, h=50, border=brewer.pal(length(isolate.assignments), "Set3"))
dev.off()

pdf("plots/pepper_dendrogram.pdf")
plot(peppers.clust)
dev.off()

#Turn isolate assignments into data.frame
assignments.df <- data.frame("Isolate" = interact_means.wide$Isolate,
                             "Assignment" = NA)
for(i in 1:length(isolate.assignments)){
  isolate.cluster <- unlist(isolate.assignments[i])
  assignments.df$Assignment[assignments.df$Isolate %in%names(isolate.cluster)] <- i
}
assignments.df$Isolate <- as.character(assignments.df$Isolate)
assignments.df$Assignment <- as.numeric(as.character(assignments.df$Assignment))


blups.sort <- blups[,names(mean_disease)]
pdf("plots/interaction_matrix_colored.pdf")
old.par <- par(no.readonly = T)
par(mfrow=c(11,11), oma=c(1,1,0.5,0.5), mar=c(0,0,0,0))
for(i in 1:nrow(blups.sort)){
  isolate_highlight <- blups$Isolate[i]
  plot(0, type='n', 
       xlim=c(1,length(mean_disease)), 
       ylim=range(blups.sort, na.rm=T),
       xaxt='n', yaxt='n', ylab='', xlab='')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],
       col = brewer.pal(length(isolate.assignments), "Set3")[assignments.df$Assignment[assignments.df$Isolate==isolate_highlight]])
  pep.background <- blups.sort[blups$Isolate!=isolate_highlight,]
  for(j in 1:nrow(blups.sort)){
    lines(1:length(mean_disease), 
          pep.background[j,], 
          main = isolate_highlight, lwd=0.5,
          col=adjustcolor("gray", alpha.f=0.6))
  }
  lines(1:length(mean_disease),
        apply(blups.sort,2,mean,na.rm=T),
        col='black', lwd=1, lty=2)
  lines(1:length(mean_disease),
        blups.sort[i,],
        col='red', lwd=1.5)
  mtext(isolate_highlight,side=3,line=-1.5, cex=0.7)
}
dev.off()


