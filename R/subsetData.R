#' Data Subset
#' 
#' Returns a new object from which certain roi(s) have been removed
#' 
#' @param dataObject \code{MSImagingExperiment} loaded \code{.imzML} file
#' @param rois \code{vector} of the regions of interest (ROIs) to be removed
#' @return \code{MSImagingExperiment} with the specefied ROIs removed
#' @example 
#' 
#' subset <- dataSubset(dataset, c("roi1","roi2"))
#' 
#' @export

subsetData <- function(dataObject, rois){
  for(i in rois){
    dataObject <- dataObject[, which(pixelData(dataObject)$roi != i)]
  }
  return(dataObject)
}

