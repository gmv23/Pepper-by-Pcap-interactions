setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/pheno/")

library(reshape2)

#Import data
pep <- read.csv("data/Ratings_clean.csv", na.strings = "NA") 
pep$Rep <- as.factor(pep$Rep)
pep$Block <- as.factor(pep$Block)
pep$Tray <- as.factor(pep$Tray)
pep$Column <- as.factor(pep$Column)

#Look at pairwise correlations and histograms, using means for peppers with two reps
pep.mat <- dcast(pep, Isolate ~ Pepper, value.var = 'audpc', fun.aggregate = mean)
pep.mat <- pep.mat[!pep.mat$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"),]

jpeg("plots/raw_distributions.jpeg", height=5, width=7, units="in", res=600)

old.par <- par(no.readonly = T)
par(mfrow=c(4,4), xpd=NA,
    mar=c(4,2,0.5,0.5), oma=c(1,1,1,1))

for(pepper in c("RedKnight", "Aristotle", "Paladin", "CM334",
                "Perennial", "EarlyJalapeno", "NMRIL-N", "NMRIL-G",
                "NMRIL-A", "NMRIL-Z", "NMRIL-I", "NMRIL-H",
                "Archimedes", "Intruder", "Revolution", "Vanguard")){
  audpc.pepper <- pep.mat[,pepper]
  hist(audpc.pepper,
       main = pepper,
       xlab = '', ylab='')
}
par(old.par)
dev.off()

#Visualize pairwise correlations
#Subset to remove peppers that basically have no virulent isolates
pep.mat.sub <- pep.mat[,-which(colnames(pep.mat) %in% c("CM334", "Intruder", "NMRIL-A", "NMRIL-G", "NMRIL-H", "NMRIL-I", "NMRIL-Z"))]

jpeg("plots/correlations.jpeg", height=5, width=8, units="in", res=600)

panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, use="complete.obs", method="spearman"), digits=2)
  txt <- paste0("r = ", r)
  cex.cor <- 1/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
text.plot <- function(x, y, labels, cex, font){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  txt <- labels
  text(0.5, 0.5, txt, cex = 0.8, font=2)
}

pairs(pep.mat.sub[,-1],
      lower.panel=panel.cor,
      text.panel=text.plot)

dev.off()

#Make interaction plot
pep_means <- apply(pep.mat.sub[,-1],2,mean, na.rm=T)
peppers_order <- names(sort(pep_means, decreasing=T))
pep_means <- sort(pep_means, decreasing=T)
pep.mat.sub.order <- pep.mat.sub[,peppers_order]
pdf("plots/interaction_matrix.pdf")
old.par <- par(no.readonly = T)
par(mfrow=c(11,11), oma=c(1,1,0.5,0.5), mar=c(0,0,0,0))
for(i in 1:nrow(pep.mat.sub.order)){
  isolate_highlight <- pep.mat.sub$Isolate[i]
  plot(0, type='n', 
       xlim=c(1,length(peppers_order)), 
       ylim=range(pep.mat.sub.order[,-1], na.rm=T),
       xaxt='n', yaxt='n', ylab='', xlab='')
  pep.background <- pep.mat.sub.order[pep.mat.sub$Isolate!=isolate_highlight,]
  for(j in 1:nrow(pep.mat.sub.order)){
    lines(1:length(peppers_order), 
          pep.background[j,], 
          main = isolate_highlight, lwd=0.5,
          col=adjustcolor("gray", alpha.f=0.6))
  }
  lines(1:length(peppers_order),
        apply(pep.mat.sub.order,2,mean,na.rm=T),
        col='black', lwd=1, lty=2)
  lines(1:length(peppers_order),
        pep.mat.sub.order[i,],
        col='red', lwd=1.5)
  mtext(isolate_highlight,side=3,line=-1.5, cex=0.7)
}
dev.off()






#Visualize check phenotypes

#First look at check phenotypes by pepper cultivar
checks <- pep[pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3"),]
pdf("plots/check_values_by_pepper.pdf")
old.par <- par(no.readonly = T)
par(mar=c(7,5,2,2))
peppers <- sort(unique(checks$Pepper))
check_names <- unique(checks$Isolate)
check_colors <- c('red', 'blue', 'gray')

plot(0, type='n', 
     xlim = c(0,16), 
     ylim = c(0, max(checks$audpc)),
     xaxt = 'n',
     ylab = "AUDPC", xlab = "")
for(i in 1:16){
  pepper <- peppers[i]
  for(j in 1:3){
    check <- check_names[j]
    check_pepper_values <- checks$audpc[checks$Isolate==check & 
                                          checks$Pepper==pepper]
    x_jitter <- rnorm(length(check_pepper_values),0,0.075)
    points(rep(i,length(check_pepper_values)) + x_jitter, 
           check_pepper_values, 
           pch=21, bg=check_colors[j], cex=1.5)
  }
}
axis(1, at=1:16, labels=peppers, las=2, cex=0.9)
legend("top", fill = check_colors, legend=check_names, bty='n', cex=0.9)
dev.off()

#Now separate by block
jpeg("plots/check_values_by_block.jpeg", width=5, height=5, units="in", res=600)
old.par <- par(no.readonly = T)
par(mar=c(3,4,0.5,4), oma=c(1,1,1,1), xpd=NA)
m <- matrix(rbind(1,1,2,2,3,3,4))
layout(m)
checks$rep_block <- checks$Rep:checks$Block
repblocks <- levels(checks$rep_block)
peppers <- c("Archimedes", "Aristotle", "EarlyJalapeno", "NMRIL-N", "Paladin", "Perennial", "RedKnight", "Revolution", "Vanguard")
check_names <- c("CHECK1", "CHECK2", "CHECK3")
pepper_colors <- brewer.pal(length(peppers), "Set1")
checks <- checks[checks$Pepper %in% peppers,]

for(i in 1:3){
  check <- check_names[i]
  checks.sub <- checks[checks$Isolate == check,]
  plot(0, type='n', 
       xlim = c(0,10), 
       ylim = c(0, max(checks.sub$audpc)),
       xaxt = 'n',
       ylab = "AUDPC", xlab = "")
  block_means <- rep(NA,10)
  for(j in 1:10){
    repblock <- repblocks[j]
    checks.sub.repblock <- checks.sub[checks.sub$rep_block == repblock,]
    check.values <- checks.sub.repblock$audpc
    x_jitter <- rnorm(length(check.values),0,0.075)
    points(rep(j,length(check.values)) + x_jitter, 
           check.values,
           pch=21, cex=1.5,
           bg= pepper_colors[match(checks.sub.repblock$Pepper, peppers)])
    block_means[j] <- mean(check.values)
  }
  lines(1:10, block_means, lty=2)
  points(1:10, block_means, pch=14)
  axis(1, at=1:10, labels=repblocks, las=1, cex=0.9)
  text(11.5, max(checks.sub$audpc)/2, check)
}

#legend
plot(0, type='n', xlab='',ylab='', xlim=c(0,10), ylim=c(0,10),
     bty='n',xaxt='n',yaxt='n')
legend(0,12,legend=peppers,fill=pepper_colors,
       ncol=3, bty='n', cex=1.2)
par(old.par)
dev.off()

#Fit check only model to look at block effects

#Just look at the 8 pepper accessions from Rep 1 and Rep 2
peppers_bal <- names(table(pep$Pepper)[table(pep$Pepper) == max(table(pep$Pepper))])
checks <- pep[pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3") & pep$Pepper %in% peppers_bal,
              c("Rep", "Block", "Isolate", "Pepper", "audpc")]
checks$Pepper <- droplevels(checks$Pepper)
#Let's 'impute' Block 5 Check 3 with the mean from the other blocks
check_means <- aggregate(checks$audpc~checks$Pepper*checks$Isolate, FUN=mean)
colnames(check_means) <- c("Pepper", "Isolate", "AUDPC")
check3_means <- check_means[check_means$Isolate == "CHECK3",]
add_to_df <- data.frame("Rep" = 1, "Block" = 5, "Isolate" = check3_means$Isolate, "Pepper" = check3_means$Pepper, "audpc" = check3_means$AUDPC)
checks <- rbind(checks, add_to_df)

checks$block <- checks$Rep:checks$Block

checks$Isolate <- droplevels(checks$Isolate)
checks.mm <- lmer(audpc ~ (1|block) + Pepper*Isolate + (1|block:Isolate), data=checks)
anova(checks.mm)
summary(checks.mm)
ranef(checks.mm)
