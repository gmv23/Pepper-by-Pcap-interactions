setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/pheno/")
library(vioplot)
library(RColorBrewer)
library(stringr)
#####################################          Read data          ####################################

#Read virulence BLUEs
phenos_full <- read.csv("../gwas/data/virulence_blues.csv")
rownames(phenos_full) <- phenos_full$Isolate
phenos_full$Isolate <- NULL

#Get rid of 'main' phenotype
phenos <- phenos_full
phenos$main <- NULL
colnames(phenos)[colnames(phenos) == "RedKnight"] <- "Red\nKnight"
colnames(phenos)[colnames(phenos) == "EarlyJalapeno"] <- "Early\nJalapeno"

#Read assignment of isolates to fields/subpopulations
pops <- read.csv("../../../isolate_collection/paper/phenotypes_and_clones/isolate_plotting_metadata.csv")
pops$SampleSZ <- as.character(pops$SampleSZ)
pops$SampleSZ[pops$SampleSZ=="14_55C"] <- "14_55" #Rename
pops$SampleSZ <- as.factor(pops$SampleSZ)

#####################################          Violin plots          ####################################

pdf("plots/violinplot.pdf", width=7,height=5)
old.par <- par(no.readonly = T)
par(mar=c(6.5,4.5,4,2))
#Draw violin plots and fill in with color
vioplot(phenos,
        h=5,
        rectCol=NA, lineCol=NA,
        colMed=adjustcolor('lightgray',alpha.f=0),
        col=adjustcolor('lightgray',alpha.f=0.75),
        ylab = 'Virulence estimate', las=2)

#Add points on top
stripchart(phenos, 
           vertical=T,
           method='jitter',
           add=T, pch=1, cex=1.25,
           col='black')

#Add violin plots again so we can add median line on top
vioplot(phenos, 
        h=5,
        rectCol=NA, lineCol=NA,
        pchMed='_',cex=4, colMed="black",
        lwd=2, add=T, col='NA')

par(old.par)
dev.off()

#####################################          Interaction plots          ####################################

#Sort peppers by mean virulence
pep_means <- apply(phenos,2,mean, na.rm=T)
pep_means <- sort(pep_means, decreasing=F)
pep_order <- names(pep_means)
phenos <- phenos[,pep_order]

#Now lets put phenos in field order --- so isolates of the same field are next to each other
pops <- pops[pops$SampleSZ %in% rownames(phenos),]
pops <- pops[order(pops$Field),]
phenos <- phenos[match(pops$SampleSZ,rownames(phenos)),]


#Abbreviate county names
field_names <- as.character(pops$Field)
rename <- function(x, remove_year=T){
  #Rename fields that have years - those are individual fields and need to be renamed
  is.field <- regexpr("\\(.*\\)", x)
  if(is.field != -1){
    if(remove_year){
      x <- gsub("\\(.*\\)", "", x)
    }
    x.split <- unlist(strsplit(x, " "))
    county <- x.split[1]
    initials <- str_to_title(substr(county,1,2))
    x.new <- paste(c(initials, x.split[2:length(x.split)]), collapse=" ")
  }else{
    x.new <- x
  }
  return(x.new)
}
new_field_names <- sapply(field_names, rename)
names(new_field_names) <- NULL
pops$Field <- new_field_names
pops$Field <- as.factor(pops$Field)

n <- nrow(phenos) #Number of isolates
p <- ncol(phenos) #Number of peppers
isolates <- rownames(phenos)

#Make vector of alternating colors for fields
#field_no <- as.integer(pops$Field)
#col_choices <- brewer.pal(3,"Pastel2")[1:2]
#field_col <- rep(NA, n)
#field_col[field_no %% 2 == 0] <- col_choices[1]
#field_col[field_no %% 2 == 1] <- col_choices[2]

#Make vector of alternating colors for fields
field_no <- as.integer(pops$Field)
col_choices <- c(brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2"), brewer.pal(6, "Set3"))
field_col <- col_choices[match(field_no, 1:23)]

pdf("plots/interaction_matrix.pdf")
old.par <- par(no.readonly = T)
par(mfrow=c(10,11), oma=c(1,1,0.5,0.5), mar=c(0,0,0,0))
for(i in 1:n){
  isolate_highlight <- isolates[i]
  plot(0, type='n', 
       xlim=c(1,p), 
       ylim=range(phenos, na.rm=T),
       xaxt='n', yaxt='n', ylab='', xlab='')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = field_col[i])
  phenos.background <- phenos[-i,]
  for(j in 1:(n-1)){
    lines(1:p, 
          phenos.background[j,], 
          lwd=0.5,
          col=adjustcolor("gray", alpha.f=0.6))
  }
  lines(1:p,
        pep_means,
        col='black', lwd=1, lty=2)
  lines(1:p,
        phenos[i,],
        col='red', lwd=1.5)
  mtext(new_field_names[i],side=3,line=-1.5, cex=0.7)
}
par(old.par)
dev.off()

#####################################          Finlay-Wilkinson Regression         ####################################

pdf("plots/FW_regression.pdf")
plot(0,type='n',
     xlim=range(phenos,na.rm=T),
     ylim=range(phenos,na.rm=T),
     xlab = "Mean pepper susceptibility",
     ylab = "Isolate virulence")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = 'white')
for(i in 1:nrow(phenos)){
  abline(lm(unlist(phenos[i,]) ~ pep_means), 
         col=field_col[i],
         lwd=1.5)
}
abline(0,1,lwd=4,lty=2,col='red')
dev.off()



#####################################          Put it together         ####################################


pdf("plots/interactions.pdf", width=7, height=4)
old.par <- par(no.readonly = T)
par(oma=c(4,4,4,4), mar=c(0,0,0,0))
m1 <- matrix(c(1:105,107:110),nrow=10, ncol=11, byrow=T)
m2 <- matrix(106, nrow=10,ncol=7)
m <- cbind(m1,m2)
layout(m)

for(i in 1:n){
  isolate_highlight <- isolates[i]
  plot(0, type='n', 
       xlim=c(1,p), 
       ylim=range(phenos, na.rm=T),
       xaxt='n', yaxt='n', ylab='', xlab='')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = field_col[i])
  phenos.background <- phenos[-i,]
  for(j in 1:(n-1)){
    lines(1:p, 
          phenos.background[j,], 
          lwd=0.5,
          col=adjustcolor("gray", alpha.f=0.6))
  }
  lines(1:p,
        pep_means,
        col='black', lwd=1, lty=2)
  lines(1:p,
        phenos[i,],
        col='red', lwd=1.5)
  mtext(new_field_names[i],side=3,line=-1.5, cex=0.7)
}
par(mar=c(4,4,4,4), oma=c(0,0,0,0))
plot(0,type='n',
     xlim=range(phenos,na.rm=T),
     ylim=range(phenos,na.rm=T))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = 'white')
for(i in 1:nrow(phenos)){
  abline(lm(unlist(phenos[i,]) ~ pep_means), 
         col=field_col[i],
         lwd=1.5)
}
abline(0,1,lwd=4,lty=2,col='red')


dev.off()


#Write new pops table
write.csv(pops,"tables/pop_assignments.csv", quote=F, row.names = F)

