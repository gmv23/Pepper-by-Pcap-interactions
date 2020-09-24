#Get residual error as percentage of all variance components
get_residual_error_percent <- function(model.as){
	res_error <- model.as$gammas["R!variance"]
	total_error <- sum(model.as$gammas)
	return(res_error/total_error)
}


#Loop through peppers
for(i in 1:nlevels(pep$Pepper)){
	pepper <- levels(pep$Pepper)[i]

	pep.sub <- pep[pep$Pepper == pepper,]

	mod <- asreml(fixed = audpc ~ at(New, 1):Isolate,
			 random = ~ at(New, 2):Isolate + Rep:Block,
			 data = pep.sub)

	mod.log <- asreml(fixed = log(audpc + 1) ~ at(New, 1):Isolate,
		random = ~ at(New, 2):Isolate + Rep:Block,
		data = pep.sub)
	res_errors[i,1] <- get_residual_error_percent(mod)
	res_errors[i,2] <- get_residual_error_percent(mod.log)
}
