g <- t(geno)

ibs <- function(x,y){
  
  alleles_ibs <- 2 - abs(x-y)
  return(sum(alleles_ibs, na.rm = T)/(2*sum(!is.na(alleles_ibs))))
  
}


d <- ncol(g)

IBS_matrix <- matrix(nrow=d, ncol=d)

for(i in 1:(d-1)){
  for (j in (i +1):d){
    IBS_matrix[i,j] <- ibs(g[,i], g[,j])
  }
  print(i)
}

rownames(IBS_matrix) <- colnames(g)
colnames(IBS_matrix) <- colnames(g)

lows <- head(sort(IBS_matrix), 10)
for(i in 1:length(lows)){
  ind <- which(IBS_matrix==lows[i], arr.ind = T)
  print(paste(lows[i],colnames(IBS_matrix)[ind]))
}

highs <- head(sort(IBS_matrix, decreasing = T), 10)
for(i in 1:length(highs)){
  ind <- which(IBS_matrix==highs[i], arr.ind = T)
  print(paste(highs[i],colnames(IBS_matrix)[ind]))
}