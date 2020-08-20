#' SpatialShrunkenCentroids Wrapper
#' 
#' Run \code{spatialShrunkenCentroids} to generate multiple segmentation models and
#' output a list of \code{spatialShrunkenCentroids2} objects containing one model per
#' object instead of one object with all segmentation models.
#' 
#' @param dataset \code{SpatialShrunkenCentroids2} or \code{list} of \code{SpatialShrunkenCentroids} object
#' @param r The spatial neightborhood radius of nearby pixels to consider. This can be a vector
#' of multiple radii values.
#' @param k The maximum number of segments (clusters). This can be a vector to try initializing
#' the clustering wiht different numbers of maximum segments. The final number of segmenets may
#' differ.
#' @param s The sparsity thresholding parameter by which to shrink the t-statistics.
#' @param ... Parameters to be passed to \code{spatialShrunkenCentroids()}
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
spatialShrunkenCentroidsWrapper <- function(dataset, r=1, k=3, s=0, ...) {
  sscList <- list()
  inc <- 1
  for (r in rparam) {
    for (k in kparam) {
      for (s in sparam) {
        ssc <- try(spatialShrunkenCentroids(dataset, r=r, k=k, s=s, ...))
        if (class(ssc) != 'try-error') {
          sscList[[inc]] <- ssc
          inc <- inc + 1
        }
      }
    }
  }
  return(sscList)
}