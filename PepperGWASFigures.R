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
#Get rid of pheno PC1 because it is basically across-pepper virulence with a different distribution
phenos$PhenoPC1 <- NULL

#Read GWAS pvalues
pvals <- read.csv("data/gwas_pvalues.csv")

#Put genos and phenos in same order
phenos <- phenos[match(rownames(geno), rownames(phenos)),]

#Subset pvals into snps and pvals and rename "main" trait
snps <- pvals[,c("CHROM","BP","MARKER")]
pvals <- pvals[,!colnames(pvals) %in% c("CHROM","BP","MARKER")]
colnames(pvals)[colnames(pvals) == "main"] <- "Across-pepper"

###########################################        Functions        ########################################

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


###########################################    Identify significant hits     ####################################3

#FDR threshold for each trait
thresholds <- apply(pvals, 2, fdr_cutoff, alpha = 0.10)

#Which markers are significant for each tray
sig_hits <- matrix(NA, nrow=nrow(pvals), ncol=ncol(pvals), dimnames = dimnames(pvals))
p <- length(thresholds)
for(i in 1:p){
  sig_hits[,i] <- pvals[,i] <= thresholds[i]
  print(which(sig_hits[,i] <- pvals[,i] <= thresholds[i]))
}

#Significant markers for one or more traits
sig_markers <- snps[apply(sig_hits, 1, any, na.rm=T),]

#Peak markers for each scaffold
sig_scaffolds <- unique(sig_markers$CHROM)
n.sig_scaffolds <- length(sig_scaffolds)
peak_markers <- matrix(NA, ncol=ncol(sig_markers), 
                       nrow=n.sig_scaffolds)
peak_markers <- as.data.frame(peak_markers)
colnames(peak_markers) <- colnames(sig_markers)
for(i in 1:n.sig_scaffolds){
  scaffold <- sig_scaffolds[i]
  sig_markers.scaff <- sig_markers[sig_markers$CHROM == scaffold,]
  peak_marker <- sig_markers.scaff[sig_markers.scaff$BP == 
                                          min(sig_markers.scaff$BP),]
  peak_markers[i,] <- peak_marker
}

##################################    Make table of allelic effects, R2, pvals, ranks     ####################################3

p.summary <- data.frame("Trait" = colnames(phenos))

for(i in 1:nrow(peak_markers)){ #Loop through peak markers
  
  #Pull marker genotype
  marker <- peak_markers$MARKER[i]
  marker.index <- which(snps$CHROM == peak_markers$CHROM[i] &
                          snps$BP == peak_markers$BP[i])
  marker.geno <- geno[,marker.index]

  #Store results for that marker
  marker.summary <- data.frame("Trait" = colnames(phenos),
                               "R2" = NA,
                               "Effect" = NA,
                               "Hit_no" = NA,
                               "Hit_percentile" = NA)

  for(j in 1:nrow(p.summary)){ #Loop through traits
    
    #Fit model -- fit marker as numerical to get allelic effect
    mod <- lm(phenos[,j] ~ marker.geno)
    
    #Populate data frame with data from simple linear model
    marker.summary$R2[j] <- summary(mod)$r.squared
    marker.summary$Effect[j] <- mod$coefficients[2]
    
    #Add info on GWAS pvalue rank
    marker.pval <- pvals[marker.index,j]
    pvals.order <- sort(pvals[,j], decreasing=F)
    marker.rank <- which(pvals.order == marker.pval)
    marker.summary$Hit_no[j] <- marker.rank
    marker.summary$Hit_percentile[j] <- marker.rank/length(pvals.order) * 100
  }
  
  #Round values and append info for that marker to dataframe
  marker.summary$R2 <- sapply(marker.summary$R2, round, digits=2)
  marker.summary$Effect <- sapply(marker.summary$Effect, round, digits=2)
  marker.summary$Hit_percentile <- sapply(marker.summary$Hit_percentile, round, digits=2)
  colnames(marker.summary)[2:ncol(marker.summary)] <- 
    paste(marker, "_", colnames(marker.summary)[2:ncol(marker.summary)], sep="")
  p.summary <- merge(p.summary,marker.summary, sort=F)
}

write.csv(p.summary, "tables/pval_summary.csv", quote=F, row.names = F)

###########################################    Make Manhattan and QQ plots     ####################################3

#Separate into traits with significant hits and those without
traits.sig <- pvals[,!is.na(thresholds)]
traits.insig <- pvals[,is.na(thresholds)]
#Make plots
draw_plots(traits.sig, "plots/GWAS_sig.pdf")
draw_plots(traits.insig, "plots/GWAS_insig.pdf")

####################################    Pull distance to nearest effector    ##################################

eff <- read.table("../effectors/data/Pc_effectors_blastout.txt")
colnames(eff) <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
                   "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
eff$sseqid <- as.integer(gsub("PHYCAscaffold_", "", eff$sseqid))

#For every SNP, calculate distance to nearest effector to get "null" distribution, then look at hits
nearest_effectors <- snps
nearest_effectors$distance <- NA
nearest_effectors$effector <- NA

for(i in 1:nrow(snps)){
  chrom <- snps$CHROM[i]
  bp <- snps$BP[i]
  if(chrom %in% eff$sseqid){
    eff.chrom <- eff[eff$sseqid == chrom,]
    #Find distance to nearest start or stop codon
    distances <- abs(bp - cbind(eff.chrom$sstart, eff.chrom$send))
    match.index <- which(distances == min(distances), arr.ind=T)
    eff.match <- eff.chrom[match.index[1],]
    #If SNP is actually INSIDE the gene, the distance with 0
    if( (bp >= eff.match$sstart & bp <= eff.match$send) | 
        (bp >= eff.match$send & bp <= eff.match$sstart)){
      distance <- 0
    }else{
      distance <- distances[match.index]
    }
    #Populate data frame with information
    nearest_effectors$distance[i] <- distance
    nearest_effectors$effector[i] <- eff.match$qseqid
  }
}

#Pull data just for peak SNPs
peak_effectors <- merge(peak_markers, nearest_effectors)
for(i in 1:nrow(peak_effectors)){
  distance <- peak_effectors$distance[i]
  peak_effectors$percentile[i] <- 
    round((sum(nearest_effectors$distance >= distance, na.rm=T) / 
             sum(!is.na(nearest_effectors$distance)))*100, 2)
}

write.csv(peak_effectors, "tables/peak_snp_contexts.csv", quote=F, row.names = F)

####################################    Pull lists of SNPs to get pairwise LDs    ##################################

# Get SNP sets to feed to vcftools in order to calculate pairwise LD between all sig SNPs
# and region +- 400 kb of peak SNPs

#First all significant positions
write.table(sig_markers[,1:2], "data/ld_snpsets/All_sig_markers.txt", sep="\t",
            quote=F, row.names=F, col.names=F)

#Now each region
for(i in 1:nrow(peak_markers)){
  chrom <- peak_markers$CHROM[i]
  bp <- peak_markers$BP[i]
  snps.peak <- snps[snps$CHROM == chrom &
                      (snps$BP > (bp - 400000) & snps$BP < (bp + 400000)),1:2]
  write.table(snps.peak,
              paste("data/ld_snpsets/chrom", chrom, "_peak.table", sep=""),
              quote=F, row.names = F, col.names = F, sep = "\t")
}





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


