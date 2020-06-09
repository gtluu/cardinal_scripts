#' Data Subset
#' 
#' Returns a new object from which certain roi(s) have been removed
#' 
#' @param dataObject is loaded data from a .imzml file
#' @param rm_roi is a vector of the rois that need to be removed
#' @return dataObject containin the new sub settes object without the specefied rois
#' @example 
#' 
#' subset <- dataSubset(omentum, c("DMEM","DMEMO"))
#' 
#' @export


dataSubset <- function(dataObject, rm_roi){
  for(i in rm_roi){
    dataObject <- dataObject[, which(pixelData(dataObject)$roi != i)]
  }
  return(dataObject)
}

