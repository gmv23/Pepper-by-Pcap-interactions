setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/pheno/")

old.par <- par(no.readonly = T)
par(mfrow=c(1,2))
for(pepper in levels(pep$Pepper)){
  hist(pep$audpc[pep$Pepper == pepper], main = pepper)
  hist(log(pep$audpc[pep$Pepper == pepper] + 1), main=pepper)
}
par(old.par)

blups <- read.csv("asreml_out/full_blups.csv")
blues <- read.csv("asreml_out/full_blues.csv")
colors <- rep('gray', nrow(blups))
colors[grep("CHECK", blups$Isolate)] <- 'red'
plot(blups$predicted.value, blues$predicted.value, col=colors)

old.par <- par(no.readonly = T)
par(mfrow=c(1,2))
for(pepper in levels(blues$Pepper)){
  hist(blues$predicted.value[blues$Pepper==pepper], main = pepper)
  hist(blups$predicted.value[blups$Pepper==pepper], main = pepper)
}
par(old.par)

residuals <- read.csv("asreml_out/residuals.csv")

library(reshape2)
blups.wide <- dcast(blups, Isolate ~ Pepper, value.var = 'predicted.value')
cv <- function(x){
  return(sqrt(var(x, na.rm=T))/mean(x, na.rm=T))
}
cvs <- apply(blups.wide[,-1],1,cv)

blues.wide <- dcast(blues, Isolate ~ Pepper, value.var = 'predicted.value')

res <- readRDS("asreml_out/model_residuals.rds")
fit <- readRDS("asreml_out/model_fittedvalues.rds")
for(pepper in names(res)){
  hist(res[[pepper]], main=pepper)
  plot(fit[[pepper]], res[[pepper]], main=pepper)
}

rk.gblup.res <- readRDS("data/gblup_res_rk.rds")
rk.gblup.fit <- readRDS("data/gblup_fit_rk.rds")

blups.plot <- as.matrix(blups[-grep("CHECK",blups$Isolate),-1])
rownames(blups.plot) <- blups$Isolate[-grep("CHECK", blups$Isolate)]
