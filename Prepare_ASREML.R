setwd("~/Documents/work/Smart_lab/P_capsici/Pepper_Interactions/paper/pheno/")

#Import data
pep <- read.csv("data/Ratings_filt.csv", na.strings = "NA") 
pep$Rep <- as.factor(pep$Rep)
pep$Block <- as.factor(pep$Block)
pep$Tray <- as.factor(pep$Tray)
pep$Column <- as.factor(pep$Column)

#Make some dummy variables
pep$checks <- rep(99, nrow(pep))
pep$checks[pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3")] <- as.character(pep$Isolate[pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3")])
pep$new <- rep(1, nrow(pep))
pep$new[pep$checks==99] <- 2
pep$checks <- as.factor(pep$checks)
pep$block <- pep$Rep:pep$Block

#Add 1 to all audpc values to make log transform possible
pep$audpc <- pep$audpc + 1

#Write files 

#Complete data
write.csv(pep, "data/asreml_in/pep_complete.csv", quote=F, row.names = F)

#Separate file for each pepper
for(pepper in levels(pep$Pepper)){
  pep.sub <- pep[pep$Pepper == pepper,]
  write.csv(pep.sub, paste("data/asreml_in/pep_", pepper, ".csv", sep=''), quote=F, row.names = F)
}


library(reshape2)
for(pepper in levels(pep$Pepper)){
  pep.sub <- pep[pep$Pepper == pepper & !(pep$Isolate %in% c("CHECK1", "CHECK2", "CHECK3")),]
  pep.long <- melt(pep.sub, measure.vars=c('Dpi4', 'Dpi6', 'Dpi8', 'Dpi11', 'Dpi13', 'Dpi15'))
  pep.long$variable <- as.numeric(gsub('Dpi','',pep.long$variable))
  pep.mean <- aggregate(value~variable*Isolate, data=pep.long, FUN=mean, na.rm=T)
  plot(0, type='n', xlim=c(0,max(pep.long$variable)), ylim=c(0,6),
       main = pepper)
  for(isolate in levels(pep.long$Isolate)){
    lines(c(0,pep.mean$variable[pep.mean$Isolate==isolate]),
          c(0,pep.mean$value[pep.mean$Isolate==isolate] + rnorm(sum(pep.mean$Isolate==isolate),0,0.25)),
          col=adjustcolor('gray', alpha.f=0.5))
    print(pep.mean$variable[pep.mean$Isolate==isolate])
    print(pep.mean$value[pep.mean$Isolate==isolate] + rnorm(sum(pep.mean$Isolate==isolate),0,0.25))
  }
}
