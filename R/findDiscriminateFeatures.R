#' Find Discriminate Features
#' 
#' Find discriminate features using \code{spatialDGMM} to apply univariate segmentation followed
#' by \code{segmentationTest} to identify statistically significant values, returned in a
#' \code{dataframe}. Note that AdjP/FDR values using this wrapper are not reliable and should be
#' ignored.
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
#' @param fixedCondition A one-sided formula giving the fixed effects of the model on the RHS.
#' The response will added to the LHS, and the formula will be passed to the underlying modeling
#' function.
#' @param classControl Either the method used to match segmented classes to the fixed effects, or
#' a list where each element is a vector of name-value pairs giving the mapping between groups and
#' classes (e.g., c(group1=class1, group2=class2, ...)). For automated matching methods, 'Ymax'
#' means to use the classes with the highest mean response for each group, and 'Mscore' means to
#' select classes based on a match score quantifying the overlap between classes and fixed effects.
#' @param BPPARAM An optional instance of \code{BiocParallelParam}. See documentation for
#' \code{bpapply}.
#' @return \code{dataframe} containing features and their associated statistics.
#' @examples
#' 
#' featuresDf <- findDiscriminateFeatures(data, r=1, k=4, fixedCondition='treated',
#'                                        classControl='Ymax')
#' 
#' @export
findDiscriminateFeatures <- function(dataset, r=1, k=3, method='gaussian', dist='chebyshev',
                                     annealing=TRUE, init='gmm', p0=0.05, iter.max=100,
                                     tol=1e-9, fixedCondition, classControl='Ymax',
                                     BPPARAM=bpparam()) {
  # Split MSImagingExperiment into list of MSImagingExperiments.
  # One feature per object.
  processedList <- list()
  for (i in 1:length(mz(dataset))) {
    processedList[[i]] <- dataset[features(dataset)[i],]
  }
  
  # Run spatial-DGMM for each individual object
  sdgmmList <- list()
  for (i in 1:length(processedList)) {
    # Test if there are any bad ROIs for each feature first.
    for (j in 1:length(levels(run(processedList[[i]])))) {
      run <- as.character(levels(run(processedList[[i]]))[j])
      idata <- processedList[[i]][,which(run(processedList[[i]])==run)]
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
      sdgmmList[[i]] <- try(if (all(NA)) {})
    } else {
      sdgmmList[[i]] <- try(spatialDGMM(processedList[[i]], r=1, k=4,
                                        groups=run(processedList[[i]]), method=method, dist=dist,
                                        annealing=annealing, init=init, p0=p0, iter.max=iter.max,
                                        tol=tol, BPPARAM=BPPARAM))
    }
  }
  
  # Run segmentationTest for each successful spatial-DGMM run.
  df <- data.frame(mz=double(), r=double(), k=double(), feature=double(),
                   LR=double(), PValue=double(), AdjP=double())
  for (i in 1:length(sdgmmList)) {
    if (class(sdgmmList[[i]]) != 'try-error') {
      segTest <- segmentationTest(sdgmmList[[i]], as.formula(paste('~', fixedCondition)),
                                  classControl=classControl)
      segTestDf <- as.data.frame(topFeatures(segTest))
      segTestDf$feature <- c(i)
      df <- rbind(df, segTestDf)
    }
  }
  return(df)
}
