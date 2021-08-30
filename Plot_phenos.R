setwd("~/Documents/Cornell/Pepper_Interactions/paper/pheno/")
library(vioplot)
library(RColorBrewer)
library(stringr)
library(grid)
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
pops <- read.csv("data/isolate_plotting_metadata.csv")
pops$SampleSZ <- as.character(pops$SampleSZ)
pops$SampleSZ[pops$SampleSZ=="14_55C"] <- "14_55" #Rename
pops$SampleSZ <- as.factor(pops$SampleSZ)

#Read average pepper LS means
pepper_means <- read.csv("data/pepper_blues.csv")
pepper_means$Pepper <- as.character(pepper_means$Pepper)
pepper_means$Pepper[pepper_means$Pepper == "RedKnight"] <- "Red\nKnight"
pepper_means$Pepper[pepper_means$Pepper == "EarlyJalapeno"] <- "Early\nJalapeno"

rownames(phenos_full) <- phenos_full$Isolate
phenos_full$Isolate <- NULL

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

#################################          Violin with dendrogram          ####################################

#Make dendrogram
phenos.t <- t(phenos)
blues.d <- dist(t(phenos), method="euclidian")
blues.dendro <- as.dendrogram(hclust(blues.d))
plot(blues.dendro)

#Reorder phenos in dendro order
phenos <- phenos[,labels(blues.dendro)]

pdf("plots/violinplot_with_dendro.pdf", width=7,height=5)
old.par <- par(no.readonly = T)

#Create layout
m <- matrix(c(1,2,2),nrow=3)
layout(m)

#Par settings for first plot
par(mar=c(0,5,4,2))

#Plot dendorgram on top
plot(blues.dendro, leaflab='none', yaxt='n')

#Par settings for second plot
par(mar=c(8,5,0,2), cex.lab=1.5)

#Draw violin plots and fill in with color
vioplot(phenos,
        h=5,
        rectCol=NA, lineCol=NA,
        colMed=adjustcolor('lightgray',alpha.f=0),
        col=adjustcolor('lightgray',alpha.f=0.75),
        ylab = 'AUDPC', las=2, cex.names=1.25)

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
pep_means <- unlist(pepper_means$predicted.value)
names(pep_means) <- pepper_means$Pepper
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
write.csv(pops, "tables/renamed_field_metadata.csv", quote=F, row.names = F)

n <- nrow(phenos) #Number of isolates
p <- ncol(phenos) #Number of peppers
isolates <- rownames(phenos)

#Make vector of colors for fields
field_no <- as.integer(pops$Field)
col_choices <- c(brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2"), brewer.pal(6, "Set3"))
field_col <- col_choices[match(field_no, 1:23)]

pdf("plots/interaction_matrix.pdf", height=7, width=7)
old.par <- par(no.readonly = T)
m <- matrix(1:110,nrow=10,ncol=11,byrow=T)
m[10,7:11] <- 107
m[10,7] <- 106
m[10,11] <- 108
m <- cbind(m,109,110)
layout(m)
par(oma=c(6,4,1.5,1.5), mar=c(0,0,0,0))
for(i in 1:n){
  isolate_highlight <- isolates[i]
  
  #Draw emtpy plot
  plot(0, type='n', 
       xlim=c(1,p), 
       ylim=c(min(phenos,na.rm = T), max(phenos,na.rm=T)*1.25),
       xaxt='n', yaxt='n', ylab='', xlab='')
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = field_col[i])
  lines(1:p,
        pep_means,
        col='black', lwd=1, lty=2)
  lines(1:p,
        phenos[i,],
        col='red', lwd=1.5)
  points(1:p,
        phenos[i,],
        col='red', lwd=1)
  
  #Write isolate name and field
#  mtext(new_field_names[i],side=3,line=-2, cex=0.7)
  mtext(isolate_highlight,side=3,line=-1.25, cex=0.55)
  
  #Draw axis tick marks
  if(i %in% seq(1,100,by=11)){
    axis(2, at=seq(0,max(phenos,na.rm=T),by=25))
  }
  if(i %in% seq(95,105,by=1)){
    axis(1,at=1:p,labels=rep('',p), las=2, cex.axis=0.5)
  }

  #Add labels
  if(i == 45){
    mtext("AUDPC",side=2,line=2.5, cex=1)
  }
  if(i == 105){
    mtext("Pepper",side=1,line=3, cex=1)
  }
  
  #Get coordinates to draw dashed line to pepper 'legend'
  if(i == 97){
    y.end <- grconvertY(0, "npc", "ndc")
    x1.end <- grconvertX(0, "npc", "ndc")
    x2.end <- grconvertX(1, "npc", "ndc")
  }
}

#Add tick mark showing order of peppers
plot(0,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='', xlim=c(1,p),ylim=c(1,10))
plot(0,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='', xlim=c(1,p),ylim=c(1,10))
axis(side = 1, at = 1:p, labels=colnames(phenos), las=2)
x1.start <- grconvertX(0, "npc", "ndc")
x2.start <- grconvertX(1, "npc", "ndc")
y.start <- grconvertY(0, "npc", "ndc")
plot(0,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='', xlim=c(1,p),ylim=c(1,10))

#Add legend
par(xpd=NA)
plot(0,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='', xlim=c(0,10),ylim=c(0,10))
legend(0,8,fill=col_choices, legend=levels(pops$Field), bty='n', cex=1.25)

pushViewport(viewport())
grid.lines(x=c(x1.start, x1.end), y=c(y.start, y.end), gp=gpar(lty=2))
grid.lines(x=c(x2.start, x2.end), y=c(y.start, y.end), gp=gpar(lty=2))

par(old.par)
dev.off()

#####################################          Finlay-Wilkinson Regression         ####################################

#Get matrix of regression slopes, intercepts, and MSEs
fw_params <- matrix(NA, nrow=n, ncol=3)
colnames(fw_params) <- c("Intercept", "Slope", "MSE")
fw_params <- as.data.frame(fw_params)
rownames(fw_params) <- rownames(phenos)

r2s <- rep(NA, nrow=n)

for(i in 1:n){
  fw.lm <- lm(unlist(phenos[i,])~pep_means)
  fw_params[i,1:2] <- fw.lm$coefficients
  r2 <- summary(fw.lm)$r.squared
  fw_params[i,3] <- anova(fw.lm)$'Mean Sq'[2]
  r2s[i] <- r2
  print(paste(rownames(phenos)[i],r2, fw_params[i,3]))
}

#Save regressoin parameters
write.csv(fw_params,"tables/fw_params.csv",quote=F,row.names = T)

#function to turn usr plot coordinates into relative plot coordinates to plot pane letters
get_coords <- function(x=-0.1,y=1.075){
  x.coord <- grconvertX(x, "npc", "user")
  y.coord <- grconvertY(y, "npc", "user")
  return(c(x.coord, y.coord))
}

pdf("plots/FW_regression.pdf", width=7, height=4)

old.par <- par(no.readonly = T)
par(mar=c(5,4,3,1))
m <- rbind(c(1,1,1,1,2,2,2,2,3))
layout(m)

#Draw empty plot
plot(0,type='n',
     xlim=range(phenos,na.rm=T),
     ylim=range(phenos,na.rm=T),
     xlab = "Average Pepper AUDPC",
     ylab = "Individual Isolate AUDPC")

#Caclulate regression intercept and slope for all isolates and plot
for(i in 1:n){
  isolate <- rownames(phenos)[i]
  abline(fw_params$Intercept[i], fw_params$Slope[i], 
         col='gray',
         lwd=1)
  #Hard code in arrows to label four isolates: 
  if(isolate == "STK_5A"){
    fw.lm <- lm(unlist(phenos[i,])~pep_means)
    x.end <- -5
    y.end <- predict(fw.lm, list("pep_means"=x.end))
    x.start <- 0
    y.start <- y.end+3
    arrows(x.start, y.start,x.end,y.end, length = 0.15)
    text(x.start+2,y.start+2,isolate)
  }
  if(isolate == "LEWT3_2A"){
    fw.lm <- lm(unlist(phenos[i,])~pep_means)
    x.end <- 8
    y.end <- predict(fw.lm, list("pep_means"=x.end))
    x.start <- 15
    y.start <- y.end+10
    arrows(x.start, y.start,x.end,y.end, length = 0.15)
    text(x.start+2,y.start+2,isolate)
  }
  if(isolate == "RCZ_1B"){
    fw.lm <- lm(unlist(phenos[i,])~pep_means)
    x.end <- -9
    y.end <- predict(fw.lm, list("pep_means"=x.end))
    x.start <- x.end+5
    y.start <- y.end+8
    arrows(x.start, y.start,x.end,y.end, length = 0.15)
    text(x.start+2,y.start+2,isolate)
  }
  if(isolate == "4E_5A"){
    fw.lm <- lm(unlist(phenos[i,])~pep_means)
    x.end <- 0
    y.end <- predict(fw.lm, list("pep_means"=x.end))
    x.start <- x.end-1
    y.start <- y.end-7
    arrows(x.start, y.start,x.end,y.end, length = 0.15)
    text(x.start-2,y.start-2,isolate)
  }
}
abline(0,1,lwd=3,lty=2,col='red')
par(xpd=NA)
text(get_coords()[1], get_coords()[2], "A", cex=1.5)

#Plot 2: Intercept vs Slope, colored by MSE
#Turn MSES into color values ranging from blue to red
col_func <- colorRampPalette(c("blue", "red"))
col_func(100)
mses <- fw_params$MSE
mses_transform <- round(100/(max(mses)-min(mses)) * (mses-max(mses)) + 100)
mses_transform[mses_transform==0] <- 1
plot(fw_params$Intercept,fw_params$Slope,
     bg=col_func(100)[mses_transform],col='black',
     pch=21,cex=2,
     xlab = "Regression Intercept",
     ylab = "Regression Slope")

#Hard code in arrows to label five isolates: 

#14_55
isolate_vals <- fw_params["14_55",]
isolate_col <- col_func(100)[mses_transform[which(rownames(fw_params)=="14_55")]]
points(isolate_vals[1],isolate_vals[2],bg=isolate_col,col='black',pch=21,cex=2)
x.end <- unlist(isolate_vals)[1] + 1
y.end <- unlist(isolate_vals)[2] +.025
x.start <- 30
y.start <- 1.5
arrows(x.start, y.start, x.end, y.end, length = 0.15)
text(x.start+5,y.start+.08,"14_55")

#13EH04A
isolate_vals <- fw_params["13EH04A",]
isolate_col <- col_func(100)[mses_transform[which(rownames(fw_params)=="14_55")]]
points(isolate_vals[1],isolate_vals[2],bg=isolate_col,col='black',pch=21,cex=2)
x.end <- unlist(isolate_vals)[1] + 1
y.end <- unlist(isolate_vals)[2]
x.start <- 12
y.start <- 1.75
arrows(x.start, y.start, x.end, y.end, length = 0.15)
text(x.start+5,y.start+.08,"13EH04A")

#17EH15_2B
isolate_vals <- fw_params["17EH15_2B",]
isolate_col <- col_func(100)[mses_transform[which(rownames(fw_params)=="14_55")]]
points(isolate_vals[1],isolate_vals[2],bg=isolate_col,col='black',pch=21,cex=2)
x.end <- unlist(isolate_vals)[1] + 1
y.end <- unlist(isolate_vals)[2] + .03
x.start <- 20
y.start <- 1.6
arrows(x.start, y.start, x.end, y.end, length = 0.15)
text(x.start+5,y.start+.08,"17EH15_2B")

#17EH68A
isolate_vals <- fw_params["17EH68A",]
isolate_col <- col_func(100)[mses_transform[which(rownames(fw_params)=="14_55")]]
points(isolate_vals[1],isolate_vals[2],bg=isolate_col,col='black',pch=21,cex=2)
x.end <- unlist(isolate_vals)[1] +1
y.end <- unlist(isolate_vals)[2] - .015
x.start <- 30
y.start <- 1.1
arrows(x.start, y.start, x.end, y.end, length = 0.15)
text(x.start+10,y.start+.02,"17EH68A")


#3W_8A
isolate_vals <- fw_params["3W_8A",]
isolate_col <- col_func(100)[mses_transform[which(rownames(fw_params)=="14_55")]]
points(isolate_vals[1],isolate_vals[2],bg=isolate_col,col='black',pch=21,cex=2)
x.end <- unlist(isolate_vals)[1]
y.end <- unlist(isolate_vals)[2] - 0.03
x.start <- 10
y.start <- 0.75
arrows(x.start, y.start, x.end, y.end, length = 0.15)
text(x.start+5,y.start-.08,"3W_8A")
text(get_coords()[1], get_coords()[2], "B", cex=1.5)

#Plot 3: MSE legend
plot(0,type='n',xlim=c(0,10),ylim=c(0,10),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
par(xpd=NA)
legend(-50,9,
       fill=col_func(100)[seq(1,100,length.out=15)],
       legend = c(round(min(mses)), rep("",13), round(max(mses))),
       title="MSE",
       cex=1.25, bty='n', y.intersp=0.75)

par(old.par)
dev.off()


