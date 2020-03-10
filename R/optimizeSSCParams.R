#' Optimize SpatialShrunkenCentroids Parameters
#' 
#' Optimize parameters used for Spatial Shrunken Centroids segmentation using
#' a heuristic algorithm.
#' 
#' @param ssc \code{SpatialShrunkenCentroids} object
#' @param sparam \code{vector} of values used for sparsity parameter in \code{spatialShrunkenCentroids}
#' @param rparam \code{vector} of values used for radius parameter in \code{spatialShrunkenCentroids}
#' @param kparam \code{vector} of values used for k parameter in \code{spatialShrunkenCentroids}
#' @return \code{list} with optimal \code{spatialShrunkenCentroids} parameters for r, k, and s
#' @examples
#' 
#' rparam <- c(1,2,3)
#' kparam <- c(2,4,6,8,10,12,14,16,18,20)
#' sparam <- c(0,3,6,9,12,15)
#' 
#' ssc <- spatialShrunkenCentroids(data, r=rparam, k=kparam, s=sparam)
#' 
#' sscParams <- optimizeSSCParams(ssc, r=rparam, k=kparam, s=sparam)
#' 
#' @export
optimizeSSCParams <- function(ssc, sparam, rparam, kparam) {
  # Optimize sparsity (s) parameter.
  # Get vectors for s values and number of unique classes per s value.
  ssc <- as.data.frame(summary(ssc))
  colnames (ssc) <- c('r', 'k', 's', 'classes', 'features_per_class', 'BIC')
  uniqueClasses <- c()
  for (s in sparam) {
    uniqueFreq <- as.data.frame(table(sort(ssc[,c('s', 'classes')][which(ssc$s==s),]$classes)))
    colnames(uniqueFreq) <- c('classes', 'freq')
    uniqueFreq$weight <- (as.numeric(uniqueFreq$freq) / as.numeric(sum(uniqueFreq$freq)))
    uniqueclass <- as.numeric(length(uniqueFreq$classes) / max(uniqueFreq$weight))
    uniqueClasses <- c(uniqueClasses, uniqueclass)
  }
  # Find the knee of the curve of uniqueClasses vs sparam.
  kneeParam <- 1
  while(TRUE) {
    knee <- kneed$KneeLocator(x=sparam, y=uniqueClasses, S=skneeparam, curve='convex', direction='decreasing')$knee
    if (!is.null(knee)) {
      optimalS <- knee
      break
    }
    kneeParam <- kneeParam - 0.1
    if (kneeParam <= 0) {
      # Error if no optimal s parameter is found.
      stop('No s parameter identified.')
    }
  }
  
  # Optimize r and k parameters.
  rk <- data.frame(r=double(), k=double(), slope=double(), knee=double())
  for (r in rparam) {
    # Find the knee for each combination of r and k.
    kneeParam <- 1
    while(TRUE) {
      knee <- kneed$KneeLocator(x=ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$s,
                                y=ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$classes,
                                S=kneeParam, curve='convex', direction='decreasing')$knee
      if (!is.null(knee)) {
        break
      }
      kneeParam <- kneeParam - 0.1
      if (kneeParam <= 0) {
        break
      }
      if (is.null(knee)) {
        knee <- NA
      }
      # Find slope * 100 of exponential curve for each combination of r and k.
      slope <- (lm(log(ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$classes) ~
                    ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$s)$coefficients[[2]]) * 100
      if (!is.na(slope)) {
        rk <- rbind(rk, c(r, k, slope, knee))
      }
    }
  }
  colnames(rk) <- c('r', 'k', 'slope', 'knee')
  # Find combinations of r and k with s==optimalS.
  # If no combination found, look for s + next increment of s.
  kneeValue <- optimalS
  while(TRUE) {
    srk <- rk[which(rk$knee==kneeValue),]
    if (nrow(srk) != 0) {
      break
    }
    if (which(sparam == kneeValue) < length(sparam)) {
      kneeValue <- sparam[which(sparam == kneeValue) + 1]
    } else {
      # Error if no kneeValue > max(sparam) parameter is found.
      stop('No s parameter identified.')
    }
  }
  
  # Find curvature * 100 (end derivative of exponential curve) at s == optimalS
  optimalParams <- data.frame(r=double(), k=double(), slope=double(), curvature=double())
  for (r in unique(srk$r)) {
    for (k in unique(srk$k)) {
      classes <- ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$classes
      sParamModel <- lm(log(classes) ~ sparam)
      sParamSeq <- seq(min(sparam), max(sparam), 0.1)
      classesModel <- exp(predict(sParamModel, sparam=list(sParamSeq)))
      dx <- numpy$gradient(sparam)
      dy <- numpy$grasdient(classes)
      d2x <- numpy$gradient(dx)
      d2y <- numpy$gradient(dy)
      curvature <- (numpy$abs(d2y) / (1 + (dy**2))**1.5) * 100
      optimalParams <- rbind(optimalParams, c(r, k, srk[which(srk$r==r & srk$k==k),]$slope,
                                              curvature[match(kneeValue, sparam)]))
    }
  }
  colnames(optimalParams) <- c('r', 'k', 'slope', 'curvature')
  
  # Remove combinations of r and k that have no slope or an increasing curve.
  optimalParams <- optimalParams[which(optimalParams$slope < 0),]
  # score = slope * curvature for each exponential curve.
  optimalParams$score <- numpy$abs(optimalParams$slope * optimalParams$curvature)
  # Choose the line where the combination of r and k has the best score.
  optimalParams <- optimalParams[which(optimalParams$score==max(optimalParams$score)),]
  optimalR <- optimalParams$r
  optimalK <- optimalParams$k
  return(list('r'=optimalR, 'k'=optimalK, 's'=optimalS))
}
