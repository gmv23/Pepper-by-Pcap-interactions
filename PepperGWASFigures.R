setwd("~/Documents/Cornell/Pepper_Interactions/paper/gwas/")
library(grid)
library(qqman)
library(fields)

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
sig_markers.indices <- apply(sig_hits, 1, any, na.rm=T)
sig_markers <- snps[sig_markers.indices,]
sig_markers <- cbind(sig_markers, pvals[sig_markers.indices,])

#Peak markers for each scaffold
sig_scaffolds <- unique(sig_markers$CHROM)
n.sig_scaffolds <- length(sig_scaffolds)
peak_markers <- matrix(NA, ncol=3, 
                       nrow=n.sig_scaffolds)
peak_markers <- as.data.frame(peak_markers)
colnames(peak_markers) <- colnames(sig_markers)[1:3]
for(i in 1:n.sig_scaffolds){
  scaffold <- sig_scaffolds[i]
  sig_markers.scaff <- sig_markers[sig_markers$CHROM == scaffold,]
  peak_marker <- sig_markers.scaff[sig_markers.scaff$"Across-pepper" == 
                                     min(sig_markers.scaff$"Across-pepper"),1:3]
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

summary(nearest_effectors$distance)

#Pull data just for peak SNPs
peak_effectors <- merge(peak_markers, nearest_effectors)
for(i in 1:nrow(peak_effectors)){
  distance <- peak_effectors$distance[i]
  peak_effectors$percentile[i] <- 
    round((sum(nearest_effectors$distance >= distance, na.rm=T) / 
             sum(!is.na(nearest_effectors$distance)))*100, 2)
}

#Add MAFs
get_maf <- function(x){
  x <- x[!is.na(x)]
  p <- (2*sum(x==2) + sum(x==1)) / (2*length(x))
  pq <- c(p, 1-p)
  return(pq[which(pq == min(pq))][1])
}
mafs <- apply(geno,2,get_maf)
peak_effectors$MAF <- round(mafs[which(snps$MARKER %in% peak_effectors$MARKER)],2)

#Add traits that are significant for that marker
peak_effectors$trait <- NA
for(i in 1:nrow(peak_effectors)){
  x <- sig_hits[snps$MARKER == peak_effectors$MARKER[i],]
  peak_effectors$trait[i] <- paste(names(x)[which(x==TRUE)], collapse = ";")
}

#save table
peak_effectors.print <- peak_effectors[,c("CHROM", "BP", "trait", "MAF", "distance", "effector")]

write.csv(peak_effectors.print, "tables/peak_snp_contexts.csv", quote=F, row.names = F)

####################################    Look at pairwise LDs    ##################################

# Get SNP sets to feed to vcftools in order to calculate 
#pairwise LD between all sig SNPs
# and region +- 400 kb of peak SNPs

####### First look at LD between all significant markers
write.table(sig_markers[,1:2], "data/ld_snpsets/All_sig_markers.txt", sep="\t",
            quote=F, row.names=F, col.names=F)
#LD between sig hits on same scaffold
system_call <- paste("vcftools --vcf ../geno/data/capsici_pepper_subset.vcf --positions ",
                     "data/ld_snpsets/All_sig_markers.txt",
                     " --geno-r2 --out data/ld/sig_markers", sep="")
system(system_call)
#LD between sig hits on different scaffolds
system_call <- paste("vcftools --vcf ../geno/data/capsici_pepper_subset.vcf --positions ",
                                        "data/ld_snpsets/All_sig_markers.txt",
                                        " --interchrom-geno-r2 --out data/ld/sig_markers", sep="")
system(system_call)
#Put them together
ld.intra <- read.table("data/ld/sig_markers.geno.ld", header=T)
ld.inter <- read.table("data/ld/sig_markers.interchrom.geno.ld", header=T)
colnames(ld.intra)[colnames(ld.intra) == "CHR"] <- "CHR1"
ld.intra$CHR2 <- ld.intra$CHR1
ld <- rbind(ld.intra, ld.inter)
ld$snp1 <- paste(ld$CHR1, ld$POS1, sep="_")
ld$snp2 <- paste(ld$CHR2, ld$POS2, sep="_")
#Turn into matrix
markers <- sig_markers$MARKER
ld.pairwise <- matrix(NA,
                      nrow = length(markers),
                      ncol = length(markers), 
                      dimnames = list(markers,markers))
for(i in 1:nrow(ld)){
  ld.pairwise[ld$snp1[i], ld$snp2[i]] <- ld$R.2[i]
}

#Make heatplot
pdf("plots/Pairwise_ld_heatplot.pdf")
old.par <- par(no.readonly=T)
par(mar=c(5,4,1,4))
imagePlot(ld.pairwise, xaxt="n", yaxt="n", col=heat.colors(12)[12:1])
axis(1,at=seq(0,1,length.out=length(markers)),markers, las=2)
axis(2,at=seq(0,1,length.out=length(markers)),markers, las=2)
par(old.par)
dev.off()

#########  Now look at LD with peak marker in significant regions
region <- 400000 #how much distance on either side of peak SNP
trait <- "Across-pepper" #which trait to show p-values for

#Put all together in one plot
pdf("plots/LD_plots.pdf", height=6, width=7)

old.par <- par(no.readonly=T)
par(mar=c(5,5,1,5), mfrow=c(4,1), xpd=F, oma=c(0.5,0.5,0.5,0.5))

for(i in 1:nrow(peak_markers)){
  
  #Create table of markers in region
  chrom <- peak_markers$CHROM[i]
  bp <- peak_markers$BP[i]
  snps.peak <- snps[snps$CHROM == chrom &
                      (snps$BP > (bp - region) & snps$BP < (bp + region)),1:2]
  file_name <- paste("data/ld_snpsets/chrom", chrom, "_peak.txt", sep="")
  write.table(snps.peak,file_name,
              quote=F, row.names = F, col.names = F, sep = "\t")
  
  #Get pairwise LD values from VCFtools
  system_call <- paste("vcftools --vcf ../geno/data/capsici_pepper_subset.vcf --positions ",
                      file_name,
                      " --geno-r2 --out data/ld/chrom",
                      chrom, "_peak", sep="") 
 system(system_call)
 ld <- read.table(paste("data/ld/chrom", chrom, "_peak.geno.ld", sep=""), 
                  header=T)
 
 #Reformat to pull R2s for SNPs with peak SNP and p-values
 ld <- ld[ld$POS1 == bp | ld$POS2 == bp,]
 ld$POS <- apply(ld[,c("POS1", "POS2")], 1, function(x)
   x[which(x != bp)])
 for(j in 1:nrow(ld)){
   ld$pval[j] <- pvals[which(snps$CHROM == chrom & snps$BP == ld$POS[j]),
                    trait]
 }
 
 #Add peak SNP LD with itself
 ld <- rbind(ld,
             c(chrom, bp, bp, NA, 1, bp, 
               pvals[snps$CHROM == chrom & snps$BP == bp,trait]))
 
 #rescale r2 distribution to -log10 distribution for plotting
 ld$logp <- -log10(ld$pval)
 ld$ld.scale <- (max(ld$logp)-min(ld$logp))/(max(ld$R.2)-min(ld$R.2)) * (ld$R.2-max(ld$R.2)) + max(ld$logp)
 
 ### Make plot
 
 plot(0, type="n", 
      xlim = range(ld$POS), 
      ylim=c(0, max(ld$logp)*1.05), xaxt="n",
      ylab = expression(paste("-log"[10], "(", italic("p"), ")")),
      xlab = paste("Physical position (Kb) on scaffold", chrom))
 for(r in 1:nrow(ld)){
   lines(x=rep(ld$POS[r],2), y = c(0,ld$logp[r]), col='gray')
 }
 points(ld$POS, ld$ld.scale, pch=2)
 points(ld$POS[nrow(ld)], ld$ld.scale[nrow(ld)], pch=17, cex=0.9, col='orange')
 axis(side=4, at = seq(0,max(ld$logp), length.out=6), labels = seq(0,max(ld$R.2), length.out=6))

 mtext(expression(italic("r")^2), side=4, line=3)
 
 #Draw rectangle
 par(xpd=NA)
 xleft <- grconvertX(0, from = "npc", to = "user")
 xright <-grconvertX(1, from = "npc", to = "user")
 ybot <- grconvertY(-0.2, from = "npc", to = "user")
 ytop <- grconvertY(-0.05, from = "npc", to = "user")
 rect(xleft, ybot, xright, ytop)
 
 #####Add effector locations
 eff.chrom <- eff[eff$sseqid == chrom & eff$send <= (bp+region) & eff$send >= (bp-region),]
 for(j in 1:nrow(eff.chrom)){
   rect(eff.chrom$sstart[j], ybot, eff.chrom$send[j], ytop, col="red")
 }
 
 #Add X axis
 par(xpd=FALSE)
 axis(side=1, at = seq(0,2000000,by=100000)[seq(0,2000000,by=100000) >= (bp-region) &
                                              seq(0,2000000,by=100000) <= (bp+region)], 
      labels = (seq(0,2000000,by=100000)/1000)[seq(0,2000000,by=100000) >= (bp-region) &
                                                 seq(0,2000000,by=100000) <= (bp+region)],
      pos=ybot)
 
}

par(old.par)
dev.off()