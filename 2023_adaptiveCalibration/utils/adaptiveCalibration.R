adaptiveCalibration <- function(obs, sim, clustering, methods,
                                  index, custom_function = NULL,
                                  weights = NULL, 
                                  scaling.type = "multiplicative",
                                  window = NULL, theta = .95){
  
  if (sum(weights) > 1) {
    
    stop("The sum of the weigths must be lower than 1")
    
  }
  
  if ((length(weights) >0) && (length(weights) != (length(index)+length(custom_function)))) {
    
    stop("Size of weights and index+custom functions are different")
  }
  
  if(length(which(methods == "gpqm")) != length(theta)){
    
    stop(paste0("Number of gpqm methods and theta size are different. Must be the same!", " Sizes: ",
                length(which(methods == "gpqm")), " != ", length(theta)))
    
  }
  
  final.list <- list()
  scores.list <- list()
  cal.method <- c()
  
  #loop on the WTs
  
  for (k in unique(clustering)) {
    
    gpqm.count <- 1
    
    aux.obs <- subsetDimension(obs, dimension = "time", 
                               indices = which(clustering == k))
    aux.sim <- subsetDimension(sim, dimension = "time", 
                               indices = which(clustering == k))
    
    n.folds <- trunc(length(aux.obs$Data)/275)
    
    if (n.folds <= 1) {
      n.folds <- 2
    }
    
    index.obs <- c()
    index.sim <- c()
    
    for (i in c(1:length(index))) {
      
      index.obs <- c(index.obs,
                     valueIndex(aux.obs, index.code = index[i])$Index$Data)
      index.sim <- c(index.sim, 
                     valueIndex(aux.sim, index.code = index[i])$Index$Data)
      
    }
    
    aux.obs <- setTimeResolution(aux.obs,
                                 time_resolution = "DD")
    aux.sim <- setTimeResolution(aux.sim, 
                                 time_resolution = "DD")
    
    aux.obs <- setGridDates.asPOSIXlt(aux.obs,
                                      tz = "UTC")
    aux.sim <- setGridDates.asPOSIXlt(aux.sim, 
                                      tz = "UTC")
    
    cal.list <- list()
    index.list <- list()
    
    #computation of the calibration methods for the WT-subseted data
    
    for (j in c(1:length(methods))) {
      
      if (methods[j] == "scaling") {
        
        cal <- biasCorrection(y = aux.obs, x = aux.sim, 
                              precipitation = TRUE, 
                              method = "scaling",
                              scaling.type = scaling.type,
                              window = window, 
                              cross.val = "kfold",
                              folds = n.folds)
        
      }else if (methods[j] == "eqm") {
        
        cal <- biasCorrection(y = aux.obs, x = aux.sim,
                              precipitation = TRUE,  
                              method = "eqm",
                              cross.val = "kfold", 
                              folds = n.folds)
        
      }else if (methods[j] == "pqm") {
        
        cal <- biasCorrection(y = aux.obs, x = aux.sim,
                              precipitation = TRUE, 
                              method = "pqm",
                              cross.val = "kfold",
                              folds = n.folds)
        
      }else if (methods[j] == "gpqm") {
        
        
        if(length(which(methods == "gpqm")) > 1){
          
          cal <- biasCorrection(y = aux.obs, x = aux.sim, 
                                precipitation = TRUE, 
                                method = "gpqm",
                                theta = theta[gpqm.count],
                                cross.val = "kfold", 
                                folds = n.folds)
          
          gpqm.count <- gpqm.count + 1
          
        }else{
          
          cal <- biasCorrection(y = aux.obs, x = aux.sim,
                                precipitation = TRUE, 
                                method = "gpqm",
                                theta = theta,
                                cross.val = "kfold",
                                folds = n.folds)
          
        }
        
      }
      
      
      
      cal <- subsetDimension(cal, dimension = "time",
                             indices = which(!is.na(cal$Data)))
      
      cal$Dates$start <- as.POSIXct(cal$Dates$start,tz = "GMT")
      cal$Dates$end <- as.POSIXct(cal$Dates$end,tz = "GMT")
      
      #computation of climate indices from VALUE
      
      index.cal <- c()
      for (i in c(1:length(index))) {
        
        index.cal <- c(index.cal, valueIndex(grid = cal, 
                                             index.code = index[i])$Index$Data)
        
      }
      
      aux.obs$Dates$start <- as.POSIXct(aux.obs$Dates$start, tz = "CEST")
      aux.obs$Dates$end <- as.POSIXct(aux.obs$Dates$end, tz = "CEST")
      
      #computation of custom functions as additional indices
      if (length(custom_function) > 0) {
        for (i in c(1:length(custom_function))) {
          index.obs <- c(index.obs, 
                         custom_function[[i]](aux.obs))
          index.cal <- c(index.cal,
                         custom_function[[i]](cal))
          
        }
        
      }
      
      cal.list[[j]] <- cal 
      index.list[[j]] <- index.cal
    }
    
    names(cal.list) <- methods
    names(index.list) <- methods
    
    #function which normalize the values for a certain index
    
    normalization <- function(measure){
      measure.norm <- c()
      for (i in c(1:length(measure))) {
        measure.norm <- c(measure.norm,
                          1-((measure[i]-min(measure))/(max(measure)-min(measure))))
      }
      return(measure.norm)
    }
    
    #computation of the absolute bias of the indices of each calibration method with respect to the observational values
    
    measures <- list()
    
    for (i in c(1:length(index.cal))) {
      aux <- c()
      for (j in c(1:length(methods))) {
        aux <- c(aux, abs(index.list[[j]][i]-index.obs[i]))
      }
      measures[[length(measures)+1]] <- aux
    }
    
    aux <- NULL
    
    #values normalization for the different indices
    
    norm.vector <- list()
    for (i in c(1:length(measures))) {
      norm.vector[[length(norm.vector)+1]] <- normalization(measures[[i]]) 
    }
    
    #calculation of the scores for the different methods
    
    scores <- c()
    for (j in c(1:(length(methods)))) {
      score <- c()
      for (i in c(1:length(measures))) {
        
        score <- c(score,norm.vector[[i]][j])
        
      }
      if (length(weights > 0)) {
        
        score <- weighted.mean(score, w = weights)
        
      }else{
        
        score <- mean(score)
        
      }
      
      scores <- c(scores, score)
    }
    
    names(scores) <- methods
    
    if (length(which(methods == "gpqm")) > 1) {
      
      idx <- which(names(scores) == "gpqm")
      for (i in c(1:length(which(methods == "gpqm")))) {
        
        names(scores)[idx[i]] <- paste0("gpqm","-",theta[i])
        
      }
    }
    
    
    scores.list[[length(scores.list)+1]] <- scores
    
    final.list[[length(final.list)+1]] <- cal.list[[order(scores, 
                                                          decreasing = T)[1]]]
    
    cal.method <- c(cal.method,names(scores)[order(scores, 
                                                   decreasing = T)][1])
  }
  
  #merging of the best calibrations per WT
  
  output <- bindGrid(final.list[[1]], 
                     final.list[[2]],
                     dimension = "time")
  
  #if the number of WTs is greater than 2 to merge the remainder is applied in a loop
  
  if (length(final.list) > 2) {
    
    for (i in c(3:length(final.list))) {
      output <- bindGrid(output,
                         final.list[[i]], 
                         dimension = "time")
    }
  }
  
  attr(output$Data, "dimensions") <- "time"
  attr(output$Data, "weather_types") <- unique(clustering)
  attr(output$Data, 'adaptive_calibration') <- cal.method
  attr(output$Data, "RF_scores") <- scores.list
  
  return(output)
}