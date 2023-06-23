adaptiveCalibration <- function(obs, sim, clustering, methods,
                                  index, custom_function = NULL,
                                  weights = rep(1/length(index), length(index)), 
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
    
    for (i in c(1:length(index))) {
      
      index.obs <- c(index.obs,
                     valueIndex(aux.obs, index.code = index[i])$Index$Data)

      
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
      
      cal.list[[j]] <- cal 
    }
    
    names(cal.list) <- methods

    scores <- RFscore(obs = obs,
                      series = cal.list,
                      index = index,
                      custom_function = custom_function,
                      methods = methods,
                      weights = weights)

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