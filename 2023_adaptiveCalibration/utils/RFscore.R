RFscore <- function(obs, series, index, custom_function = NULL, methods, weights = rep(1/length(index), length(index))){
  
  if (sum(weights) > 1) {
    
    stop("The sum of the weigths must be lower than 1")
    
  }
  
  if ((length(weights) >0) && (length(weights) != (length(index)+length(custom_function)))) {
    
    stop("Size of weights and index+custom functions are different")
  }
  

  
  index.obs <- c()
  
  for (i in c(1:length(index))) {
    
    index.obs <- c(index.obs, 
                   valueIndex(obs, index.code = index[i])$Index$Data)
    
    
  }
  aux.obs <- setTimeResolution(obs, 
                               time_resolution = "DD")
  aux.obs <- setGridDates.asPOSIXlt(aux.obs,
                                    tz = "UTC")
  
  
  index.list <- list()
  for (j in c(1:length(methods))) {
    index.cal <- c()
    for (i in c(1:length(index))) {
      
      index.cal <- c(index.cal,
                     valueIndex(grid = series[[j]], index.code = index[i])$Index$Data)
      
    }
    
    aux.obs$Dates$start <- as.POSIXct(aux.obs$Dates$start,
                                      tz = "CEST")
    aux.obs$Dates$end <- as.POSIXct(aux.obs$Dates$end,
                                    tz = "CEST")
    
    if (length(custom_function) > 0) {
      for (i in c(1:length(custom_function))) {
        index.obs <- c(index.obs, 
                       custom_function[[i]](aux.obs))
        index.cal <- c(index.cal,
                       custom_function[[i]](series[[j]]))
        
      }
      
    }
    
    index.list[[j]] <- index.cal
  }
  
  names(index.list) <- methods
  
  
  
  
  normalization <- function(measure){
    measure.norm <- c()
    for (i in c(1:length(measure))) {
      measure.norm <- c(measure.norm, 
                        1-((measure[i]-min(measure))/(max(measure)-min(measure))))
    }
    return(measure.norm)
  }
  
  measures <- list()
  
  for (i in c(1:length(index.cal))) {
    aux <- c()
    for (j in c(1:length(methods))) {
      aux <- c(aux, abs(index.list[[j]][i]-index.obs[i]))
    }
    measures[[length(measures)+1]] <- aux
  }
  
  aux <- NULL
  
  norm.vector <- list()
  for (i in c(1:length(measures))) {
    norm.vector[[length(norm.vector)+1]] <- normalization(measures[[i]]) 
  }
  
  
  
  scores <- c()
  for (j in c(1:(length(methods)))) {
    score <- c()
    for (i in c(1:length(measures))) {
      
      score <- c(score,norm.vector[[i]][j])
      
    }
    
    scores <- c(scores, score <- weighted.mean(score, w = weights))
  }
  
  names(scores) <- methods
  
  
  return(scores)
}