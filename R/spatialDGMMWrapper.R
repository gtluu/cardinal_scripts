#' Spatial-DGMM Wrapper
#' 
#' Wrapper for \code{spatialDGMM} to apply univariate segmentation while ignoring problematic
#' features that cause errors.
#' 
#' @param dataset \code{MSImagingExperiment} to be analyzed.
#' @param r The spatial neightborhood radius of nearby pixels to consider. This can be a vector
#' of multiple radii values.
#' @param k The maximum number of segments (clusters). This can be a vector to try initializing
#' the clustering wiht different numbers of maximum segments. The final number of segmenets may
#' differ.
#' @param ... Parameters to be passed to \code{spatialDGMM()}.
#' @return \code{list} of \code{SpatialDGMM} objects.
#' @examples
#' 
#' sdgmmList <- spatialDGMMWrapper(data, r=1, k=4)
#' 
#' @export
spatialDGMMWrapper <- function(dataset, r=1, k=3, ...) {
  # Split MSImagingExperiment into list of MSImagingExperiments.
  # One feature per object.
  processedList <- list()
  for (i in 1:length(mz(dataset))) {
    processedList[[i]] <- dataset[features(dataset)[i],]
  }
  
  # Run spatial-DGMM for each individual object
  sdgmm <- function(processed, ...) {
    library(Cardinal)
    # Test if there are any bad ROIs for each feature first.
    for (j in 1:length(levels(pixelData(processed)$run))) {
      run <- as.character(levels(pixelData(processed)$run)[j])
      idata <- processed[,which(pixelData(processed)$run==run)]
      idata <- imageData(idata)[[1]][1,]
      if (length(unique(idata)) != 1) {
        badFeature <- FALSE
      } else {
        badFeature <- TRUE
        break
      }
    }
    print(paste('Feature:', as.character(i)))
    if (badFeature) {
      return(try(if (all(NA)) {}))
    } else {
      return(try(spatialDGMM(processed, r=r, k=k, ...)))
    }
  }
  
  if (class(bpparam()) == 'SerialParam') {
    # Linear
    sdgmmList <- list()
    for (i in 1:length(processedList)) {
      sdgmmList[[i]] <- sdgmm(processedList[[i]])
    }
  } else {
    # Parallel
    sdgmmList <- bplapply(processedList, FUN=sdgmm)
  }
  
  return(sdgmmList)
}
