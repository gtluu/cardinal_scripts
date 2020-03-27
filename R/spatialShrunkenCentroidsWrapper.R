#' SpatialShrunkenCentroids Wrapper
#' 
#' Run \code{spatialShrunkenCentroids} to generate multiple segmentation models and
#' output a list of \code{spatialShrunkenCentroids2} objects containing one model per
#' object instead of one object with all segmentation models.
#' 
#' @param dataset \code{SpatialShrunkenCentroids2} or \code{list} of \code{SpatialShrunkenCentroids} object
#' @param sparam \code{integer} value or \code{vector} of values used for sparsity parameter in \code{spatialShrunkenCentroids}
#' @param rparam \code{integer} value or \code{vector} of values used for radius parameter in \code{spatialShrunkenCentroids}
#' @param kparam \code{integer} value or \code{vector} of values used for k parameter in \code{spatialShrunkenCentroids}
#' @return \code{list} of \code{spatialShrunkenCentroids2} objects
#' @examples
#' 
#' rparam <- c(1,2,3)
#' kparam <- c(2,4,6,8,10,12,14,16,18,20)
#' sparam <- c(0,3,6,9,12,15)
#' 
#' ssc <- spatialShrunkenCentroidsWrapper(data, r=rparam, k=kparam, s=sparam)
#' 
#' @export
spatialShrunkenCentroidsWrapper <- function(dataset, r, k, s) {
  sscList <- list()
  inc <- 1
  for (r in rparam) {
    for (k in kparam) {
      for (s in sparam) {
        ssc <- try(spatialShrunkenCentroids(dataset, r=r, k=k, s=s))
        if (class(ssc) != 'try-error') {
          sscList[[inc]] <- ssc
          inc <- inc + 1
        }
      }
    }
  }
  return(sscList)
}