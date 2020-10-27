setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/effectors/")

#########################################  Read data #####################################################

snps <- read.table("../geno/data/capsici_pepper_subset.012.pos")
colnames(snps) <- c("Chrom", "Bp")

eff <- read.table("data/Pc_effectors_blastout.txt")
colnames(eff) <- c("qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore")
eff$sseqid <- sapply(as.character(eff$sseqid), function(x) as.integer(unlist(strsplit(x, "_"))[2]))

genes <- read.delim("../../../isolate_collection/analysis/genome/annotations.txt", sep='\t')
genes$chrom <- sapply(as.character(genes$chrom), function(x) as.integer(unlist(strsplit(x, "_"))[2]))

pvals <- read.csv("../gwas/data/gwas_pvalues.csv")

####################################### Add effector predictions to genes ###########################

#If no overlap between a gene prediction and an effector prediction, make a new gene
#If overlap, make a consensus gene that is the union

new_eff <- data.frame("eff" = eff$qseqid,
                      "gene" = NA,
                      "chrom" = eff$sseqid,
                      "start" = NA,
                      "stop" = NA)

n <- nrow(new_eff)

for(i in 1:n){
  chrom <- eff$sseqid[i]
  start.eff <- eff$sstart[i]
  stop.eff <- eff$send[i]
  genes.chrom <- genes[genes$chrom==chrom,]
  #If no genes annotated on that chromosome
  if(nrow(genes.chrom) == 0){
    new_eff$start[i] <- start.eff
    new_eff$stop[i] <- stop.eff
    #Otherwise, look to see if any overlap
  }else{
    #Loop through genes on that chrom
    for(j in 1:nrow(genes.chrom)){
      start.gene <- genes.chrom$start[j]
      stop.gene <- genes.chrom$end[j]
      #check orientation
      if(stop.gene > start.gene){
        if((start.eff >= start.gene & start.eff <= stop.gene) |
           (stop.eff >= start.gene & stop.eff <= stop.gene)){ #If overlap, extend
          new_eff$start[i] <- min(c(start.eff, start.gene))
          new_eff$stop[i] <- max(c(stop.eff, stop.gene))
          new_eff$gene[i] <- genes.chrom$proteinID[j]
        }else{ #If no overlap with any genes
          new_eff$start[i] <- start.eff
          new_eff$stop[i] <- stop.eff
        }
      }else if(start.gene > stop.gene){
        if((start.eff <= start.gene & start.eff >= stop.gene) |
           (stop.eff <= start.gene & stop.eff >= stop.gene)){ #If overlap, extend
          new_eff$start[i] <- max(c(start.eff, start.gene))
          new_eff$stop[i] <- min(c(stop.eff, stop.gene))
          new_eff$gene[i] <- genes.chrom$proteinID[j]
        }else{ #If no overlap with any genes
          new_eff$start[i] <- start.eff
          new_eff$stop[i] <- stop.eff
        }
      }
    } 
  }
}

#Now append list of genes with no effector match, then sort to get new list 
genes_unmatch <- genes[!genes$proteinID %in% new_eff$gene,1:4]
genes_add <- data.frame("eff" = NA,
                        "gene" = genes_unmatch$proteinID,
                        "chrom" = genes_unmatch$chrom,
                        "start" = genes_unmatch$start,
                        "stop" = genes_unmatch$end)
all_genes <- rbind(new_eff, genes_add)
all_genes <- all_genes[order(all_genes$chrom, all_genes$start),]

####################################### Find closest genes to GWAS hits ###########################

#p data frame with chrom, bp, pval

p <- pvals[,c("CHROM", "BP", "RedKnight")] 
percent <- .01
return_gene_list <- function(p, genes, percent){
  n.top <- round(nrow(p)*percent)
  p <- p[order(p[,3]),]
  p <- p[1:n.top,]
  
  #Make list to fill in with genes
  closest_genes <- matrix(NA, nrow=n.top*2, ncol=ncol(genes)) #Make it as big as possible if every single match is a tie bw 2 genes
  colnames(closest_genes) <- colnames(genes)
  closest_genes <- as.data.frame(closest_genes)
  match_number <- 1
  for(i in 1:n.top){
    chrom <- p$CHROM[i]
    bp <- p$BP[i]
    #Find closest gene to SNP
    genes.chrom <- genes[genes$chrom==chrom,]
    distance_start <- abs(bp-genes.chrom$start)
    distance_stop <- abs(bp-genes.chrom$stop)
    gene_dist <- apply(cbind(distance_start,distance_stop),1,min)
    closest <- which(gene_dist == min(gene_dist))
    if(length(closest)==1){
      closest_genes[match_number,] <- unlist(genes.chrom[closest,])
      match_number <- match_number + 1
      added <- 1
    }else if(length(closest)==2){
      matches <- as.matrix(genes.chrom[closest,])
      if(length(unique(paste(genes$eff, genes$gene)))==1){
        closest_genes[match_number,] <- matches[1,]
        match_number <- match_number + 1
        added <- 1
      }else if(length(unique(paste(genes$eff, genes$gene)))==2){
        closest_genes[match_number:(match_number+1),] <- matches[]
        match_number <- match_number + 2
        added <- 2
      }else{
        print("hmmmm...?")
      }
    }
  }
  #Get rid of extra rows
  closest_genes <- closest_genes[1:(match_number-added),]
  #Get rid of duplicated entries
  rows_to_remove <- c()
  gene_ID <- paste(closest_genes$eff, closest_genes$gene)
  gene_ID.table <- table(gene_ID)
  gene_ID.dups <- names(gene_ID.table[gene_ID.table>1])
  for(dup in gene_ID.dups){
    matches <- which(gene_ID %in% dup)
    rows_to_remove <- c(rows_to_remove,matches[-1])
  }
  closest_genes <- closest_genes[-rows_to_remove,]
  return(closest_genes)
}

prop <- rep(NA,ncol(pvals)-3)
names(prop) <- colnames(pvals)[4:ncol(pvals)]
for(i in 1:length(prop)){
  p <- pvals[,c(1:2,(i+3))]
  closest_genes <- return_gene_list(p, all_genes, 0.01)
  x1 <- sum(!is.na(closest_genes$eff))
  x2 <- nrow(closest_genes)
  y1 <- sum(!is.na(all_genes$eff))
  y2 <- nrow(all_genes)
  con.table <- cbind(c(x1,x2),c(y1,y2))
  fish.out <- fisher.test(con.table)
  prop[i] <- fish.out$p.value
  #prop[i] <- sum(!is.na(closest_genes$eff))/nrow(closest_genes)
}
