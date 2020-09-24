### Functions for time series disease data either in GH or field 
# audpc function gives same results as agricolae audpc function, has some extra options
# like when to start calculating (time point 0 vs first disease rating)
# can do AUDPS or AUDPC, relative or absolute

check_errors <- function(scores, time_points_columns, germ_column, fix_decreases, fix_over){
  
  for(r in 1:nrow(scores)){
    germ_mistake <- F
    disease_mistake <- F
    rw <- as.numeric(scores[r, time_points_columns]) #Make current row numeric vector to work with
    germ <- scores[r,germ_column]
    
    #Check values that are greater than germ number
    over_values <- rw > germ
    if(any(over_values)){
      germ_mistake <- T
      print(scores[r,])
      rw[over_values] <- germ
    }
    
    #Check for decreases in disease over time (in my case user error) -> disease was called prematurely
    if (any(rw != sort(rw))){
      disease_mistake <- T
      if(!germ_mistake){ #only print line again if you havent already printed it from a mistake with the germ checker
        print(scores[r,])
      }
      if(fix_decreases){
        for(i in 1:length(rw)){ #Loop over values in row
          if(i == length(rw)){
            break
          } #Break at end of row
          else{ #If value is greater than minimum value in rest of row, replace with minimum value in rest of row
            if(rw[i] > min(rw[(i+1):length(rw)])){
              rw[i] <- min(rw[(i+1):length(rw)]) 
            } 
          } 
        }     
      } 
    }
    if( (fix_decreases | fix_over) & (germ_mistake | disease_mistake) ){ #If you want to fix mistakes, print fixed version
      scores[r,time_points_columns] <- rw
      print(scores[r,])
      print("..........................................................................................................")
    }
  } #End loop through rows
  if(fix_decreases | fix_over){
    return(scores)
  }
} #End function

#Option to exclude anything with germ less than a certain amount
convert_percentages <- function(scores, time_points_columns, germ_column, exclude=FALSE){
  
  scores[,time_points_columns] <- scores[,time_points_columns]/scores[,germ_column] * 100
  
  if(exclude != FALSE){
    for(r in 1:nrow(scores)){
      if(scores[r,germ_column] < exclude){
        scores[r,time_points_columns] <- rep(NA,length(time_points_columns))
      } 
    } 
  } 
  return(scores)
} 

# start_at = c("zero", order of column with first symptoms)
# method = c("curve", "steps")
# type = c("absolute", "standard", "relative")
my_audpc <- function(disease_values, time_points,
                     start_at = 1,
                     method = "curve",
                     type = "absolute"){
  
  if(length(disease_values) != length(time_points)){
    stop('Time points vector different length from disease ratings vector.')
  }
  
  #Where to start? Default is with first value. You can also specify a rating to start at or 
  #specify to start at day 0 (day of inoculation)
  if(start_at == "zero"){
    if(time_points[1] != 0){
      disease_values <- c(0, disease_values)
      time_points <- c(0, time_points)
    }
  }else if(is.numeric(start_at) & start_at >= 1 & start_at < length(disease_values)){
    disease_values <- disease_values[start_at:length(disease_values)]
    time_points <- time_points[start_at:length(time_points)]
  }else{
    stop("time_points must be zero or numeric value with first disease column")
  }
  
  #Calculate AUDPC for use in either AUDPC or AUDPS
  
  t <- length(time_points)
  audpc <- 0
  
  for(i in 1:t){
    if(i == t){
      break
    }
    else{
      trapezoid <- ((disease_values[i] + disease_values[i+1])/2) * (time_points[i+1] - time_points[i])
      audpc <- audpc + trapezoid
    } 
  }
  
  n <- length(time_points)
  D <- time_points[n] - time_points[1]
  #Add on weighted first and last observation for steps
  if(method == "curve"){
    disease_area <- audpc
  }else if(method == "steps"){
    disease_area <- audpc + (disease_values[1] + disease_values[n])/2*D/(n-1)
  }else{
    stop("Only curve or stairs are valid methods")
  }
  
  #Absolute or standard or relative?
  if(type == "absolute"){
    statistic <- disease_area
  }else if(type == "relative" & method == "steps"){
    statistic <- disease_area * (n-1) / (D * n * 100)  
  }else if(type == "relative" & method == "curve"){
    statistic <- disease_area/(D*100)
  }
  
  return(statistic)
  
} # End function



