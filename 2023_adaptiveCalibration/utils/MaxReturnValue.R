MaxReturnValue <- function(data, indices = NULL){
  require(evd)
  data$Dates$start <- as.POSIXct(data$Dates$start,
                                 tz = "GMT")
  data$Dates$end <- as.POSIXct(data$Dates$end,
                               tz = "GMT")
  if(!is.null(indices)){
    subsetDimension(grid = data, dimension = "time", 
                    indices = indices)
  }
  nyears <- 20
  data.maxYear <- aggregateGrid(redim(data), 
                                aggr.y = list(FUN = "max",na.rm = TRUE))
  
  data.maxrv <- climatology(data, 
                            clim.fun = list(FUN = "mean", na.rm = T))
  
  auxData <- data.maxYear$Data
  auxData[which(is.infinite(auxData))] <- NA
  
  if (any(!is.na(auxData))){
    auxGEV <- fgev(auxData)
    if ((auxGEV$estimate[3] - auxGEV$std.err[3] < 0) & (0  < auxGEV$estimate[3] + auxGEV$std.err[3])){
      auxGEV <- fgev(auxData, shape = 0)
      auxRV <- qgev(1-1/nyears, 
                    loc = auxGEV$estimate[1], 
                    scale = auxGEV$estimate[2], 
                    shape = 0)
    }else{
      auxRV <- qgev(1-1/nyears, 
                    loc = auxGEV$estimate[1], 
                    scale = auxGEV$estimate[2], 
                    shape = auxGEV$estimate[3])
    }
    data.maxrv <- as.numeric(auxRV)
    names(data.maxrv) <- paste("RV",nyears,"_max", sep = "")
  }
  return(data.maxrv)
  
}