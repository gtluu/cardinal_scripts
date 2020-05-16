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
#' @param method The method to use to calculate the spatial smoothing weights. The 'gaussian'
#' method refers to spatially-aware (SA) weights, and 'adaptive' refers to spatially-aware
#' structurally-adaptive (SASA) weights.
#' @param dist The type of distance metric to use when calculating neighboring pixels based on r.
#' The options are 'radial', 'manhattan', 'minkowski', and 'chebyshev' (the default).
#' @param annealing Should simulated annealing be used during the optimization process to
#' improve parameter estimates?
#' @param init Should the parameter estimates be initialized using k-means ('kmeans') or Gaussian
#' mixture model ('gmm')?
#' @param p0 A regularization parameter applied to estimated posterior class probabilities to
#' avoid singularities. Must be positive for successful gradient descent optimization. Changing
#' this value (within reason) shoudl have only minimal impact on values of parameter estimates,
#' but may greatly affect the algorithm's speed and stability.
#' @param iter.max The maximum number of EM-algorithm iterations.
#' @param tol The toelrance convergence criterion for the EM-algorithm. Corresponds to the change
#' in log-likelihood.
#' @param BPPARAM An optional instance of \code{BiocParallelParam}. See documentation for
#' \code{bpapply}.
#' @return \code{list} of \code{SpatialDGMM} objects.
#' @examples
#' 
#' sdgmmList <- spatialDGMMWrapper(data, r=1, k=4)
#' 
#' @export
spatialDGMMWrapper <- function(dataset, r=1, k=3, method='gaussian', dist='chebyshev',
                               annealing=TRUE, init='gmm', p0=0.05, iter.max=100,
                               tol=1e-9, BPPARAM=BiocParallel::bpparam()) {
  # Split MSImagingExperiment into list of MSImagingExperiments.
  # One feature per object.
  processedList <- list()
  for (i in 1:length(mz(dataset))) {
    processedList[[i]] <- dataset[features(dataset)[i],]
  }
  
  # Run spatial-DGMM for each individual object
  sdgmm <- function(processed) {
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
      return(try(spatialDGMM(processed, r=r, k=k, groups=pixelData(processed)$run,
                                       method=method, dist=dist, annealing=annealing,
                                       init=init, p0=p0, iter.max=iter.max, tol=tol,
                                       BPPARAM=BPPARAM)))
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
