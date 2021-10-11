setwd("~/Documents/Cornell/Pepper_Interactions/paper/gwas/")
library(grid)
library(qqman)

#######################################          Load and clean data       ######################################

#Read geno as 012 file
geno <- read.table("../geno/data/capsici_pepper_subset.012")
geno$V1 <- NULL #Get rid of row names
geno <- as.matrix(geno)
geno[geno==-1] <- NA

#Read sample names
indvs <- read.table("../geno/data/capsici_pepper_subset.012.indv")
indvs <- as.character(unlist(indvs$V1))
rownames(geno) <- indvs

#Read virulence BLUEs
phenos <- read.csv("data/virulence_blues.csv")
rownames(phenos) <- phenos$Isolate
phenos$Isolate <- NULL

#Read phenotypic PCs and add to phenos
pheno.pca <- read.csv("../pop_structure/tables/phenotypic_pcs.csv")
rownames(pheno.pca) <- pheno.pca$X
pheno.pca$X <- NULL
pheno.pca <- pheno.pca[match(rownames(phenos), rownames(pheno.pca)),]
colnames(pheno.pca) <- paste("Pheno", colnames(pheno.pca), sep="")
phenos <- cbind(phenos, pheno.pca)

#Read GWAS pvalues
pvals <- read.csv("data/gwas_pvalues.csv")

#Put genos and phenos in same order
phenos <- phenos[match(rownames(geno), rownames(phenos)),]

#Subset pvals into snps and pvals and rename "main" trait
snps <- pvals[,c("CHROM","BP","MARKER")]
pvals <- pvals[,!colnames(pvals) %in% c("CHROM","BP","MARKER")]
colnames(pvals)[colnames(pvals) == "main"] <- "Across-pepper"

###########################################        Plot manhattan and qq plots        ########################################

#Function to get FDR threshold

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

#Function to get plot coordinates for pane letters
get_coords <- function(x=-0.04,y=1.4){
  require(grid)
  x.coord <- grconvertX(x, "npc", "user")
  y.coord <- grconvertY(y, "npc", "user")
  return(c(x.coord, y.coord))
}

#Function to draw stacked manhattan and qq plots
draw_plots <- function(pvals.plot, name){
  
  n.plots <- ncol(pvals.plot)
  pdf(name, width=7, height=(1.5*n.plots))
  
  old.par <- par(no.readonly = T)
  par(mar=c(5.5,4.5,0.25,1.5), oma=c(0,1,2.5,1), xpd=NA)
  
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
    man.df$SNP <- 1:nrow(man.df)
    par(xpd=F)
    manhattan(man.df,suggestiveline=FALSE,
              genomewideline=genomewideline,
              xlab = xlab.man,
              ylim=c(0,-log10(min(pvals.plot))*1.1))
    par(xpd=NA)
    mtext(trait,side=3,line=0.75)
    text(get_coords()[1], get_coords()[2], LETTERS[i], cex=1.7)
    par(xpd=F)
    qq(man.df$P, xlab=xlab.qq, ylab="")
    par(xpd=NA)
    mtext(trait,side=3,line=0.75)
    mtext("Expected", side=2, line=3.5, cex=0.65)
    mtext(expression(paste("-log"[10], "(",italic("p"),")")), side=2, line=2.25, cex=0.65)
    
  }
  
  par(old.par)
  
  dev.off()
  
}

#Separate into peppers with significant hits and those without
thresholds <- apply(pvals, 2, fdr_cutoff, alpha = 0.10)
pvals.sig <- pvals[,!is.na(thresholds)]
pvals.insig <- pvals[,is.na(thresholds)]

#Make plots
draw_plots(pvals.sig, "plots/GWAS_sig.pdf")
draw_plots(pvals.insig, "plots/GWAS_insig.pdf")

###########################################        Plot SNP effects       ########################################

#Which is the significant SNP we identified
sig_snp <- which(pvals$RedKnight == min(pvals$RedKnight))
sig_geno <- as.factor(geno[,sig_snp])

get_coords <- function(x=-0.15,y=1.3){
  require(grid)
  x.coord <- grconvertX(x, "npc", "user")
  y.coord <- grconvertY(y, "npc", "user")
  return(c(x.coord, y.coord))
}

pdf("plots/snp_effects.pdf")

old.par <- par(no.readonly = T)
par(mfrow=c(3,3))

for(i in 1:9){
  boxplot.data <- boxplot(phenos[,i] ~ sig_geno, plot=F)
  boxplot(phenos[,i] ~ sig_geno, main=colnames(pvals)[i],
          xlab = "Genotype",
          ylab = "AUDPC", 
          names = c("CC", "CA", "AA"), outline=F)
  stripchart(phenos[,i] ~ sig_geno, vertical=T, pch=1, add=T, method="jitter", cex=1)
  par(xpd=NA)
  text(get_coords()[1], get_coords()[2], LETTERS[i], cex=1.5)
  par(xpd=F)
}

par(old.par)
dev.off()

####################################    Show LD to peak SNPs in significant regions    ##################################

peak_bp <- snps[sig_snp,"BP"]

#Show peak +/- this distance, and for this trait
distance <- 400000
trait <- "RedKnight"

filter <- snps$CHROM == 39 & snps$BP > (peak_bp-distance) & snps$BP <= (peak_bp+distance)
n <- sum(filter)

#Subset snps, genotypes, and pvals
snps.sub <- snps[filter,]
pvals.sub <- pvals[filter,trait]
geno.sub <- geno[,filter]

#get -log10 pvalues
logp <- -log10(pvals.sub)

#Get LD to peak SNP
peak <- which(logp == max(logp))
geno.peak <- geno.sub[,peak]
lds <- rep(NA, n)
for(i in 1:n){
  lds[i] <- cor(geno.sub[,i], geno.peak, use = "complete.obs")^2
}

#rescale r2 distribution to -log10 distribution for plotting
lds.rescale <- (max(logp)-min(logp))/(max(lds)-min(lds)) * (lds-max(lds)) + max(logp)

pdf("plots/LD_with_SNP.pdf")

old.par <- par(no.readonly=T)
par(mar=c(5,5,3,5))

#Plot 1: p-values and LD
plot(0, type="n", 
     xlim = range(snps.sub$BP), 
     ylim=c(0, max(logp)), xaxt="n",
     ylab = expression(paste("-log"[10], "(", italic("p"), ")")),
     xlab = "Physical position (Kb) on scaffold 39")
for(r in 1:n){
  lines(x=rep(snps.sub$BP[r],2), y = c(0,logp[r]), col='gray')
}
points(snps.sub$BP, lds.rescale, pch=2)
points(snps.sub$BP[peak], lds.rescale[peak], pch=17, cex=0.9, col='orange')
axis(side=4, at = seq(0,max(logp), length.out=6), labels = seq(0,max(lds), length.out=6))
axis(side=1, at = seq(100000,1000000,by=100000), labels = seq(100000,1000000,by=100000)/1000)
mtext(expression(italic("r")^2), side=4, line=3)

par(old.par)
dev.off()


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


