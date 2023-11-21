#' Get Summary DataFrame
#' 
#' Aggregate summary SSC information into single dataframe
#' 
#' @param sscList \code{list} of \code{SpatialShrunkenCentroids} object
#' @return \code{dataframe} with \code{summary()} information from \code{SpatialShrunkenCentroids2} objects
#' 
#' @export
getSSCDf <- function(sscList) {
  sscDf <- data.frame(r=double(), k=double(), s=double(), classes=double(), features_per_class=double())
  for (ssc in sscList) {
    tmp <- as.data.frame(Cardinal::summary(ssc))
    colnames(tmp) <- c('r', 'k', 's', 'classes', 'features_per_class')
    sscDf <- rbind(sscDf, tmp)
  }
  return(sscDf)
}

#'Makes a Sparsity plot
#'
#'a plot showing weighted cardinality score vs shrinkage parameter (s)
#'
#'@param sparam \code{vector} of values used for sparsity parameter in \code{spatialShrunkenCentroids}
#'@param uniqueClasses \code{dataframe} with plateaus removed resulting from k == classes obtained from \code{optimizeSParam}
#'@param optimalS value of optimized s
#'@return plot showing weighted cardinality score vs shrinkage parameter (s)
#'
#'@export
plotSparsity <- function(sparam, uniqueClasses, optimalS){
    figure <- ggplot() +
      xlab('sparsity (s)') +
      ylab('weighted cardinality score') +
      geom_point(aes(x=sparam, y=uniqueClasses)) +
      geom_line(aes(x=sparam, y=uniqueClasses)) +
      geom_vline(xintercept=optimalS)
    return(figure)
}

#' Optimize Sparsity (s) Parameter
#' 
#' Get the optimal parameter for s based on \code{spatialShrunkenCentroids} models
#' 
#' @param ssc \code{dataframe} with \code{summary()} information from \code{SpatialShrunkenCentroids2} objects
#' @param sparam \code{vector} of values used for sparsity parameter in \code{spatialShrunkenCentroids}
#' @param plot can be TREUE or FALSE to define if the plot showing weighted cardinality score vs s parameter is required or not
#' @return \code{integer} optimal value of s
#' 
#' @export
optimizeSParam <- function(ssc, sparam, plot=TRUE) {
  uniqueClasses <- c()
  for (s in sparam) {
    uniqueFreq <- as.data.frame(table(sort(ssc[,c('s', 'classes')][which(ssc$s==s),]$classes)))
    colnames(uniqueFreq) <- c('classes', 'freq')
    uniqueFreq$weight <- (as.numeric(uniqueFreq$freq) / as.numeric(sum(uniqueFreq$freq)))
    uniqueClass <- as.numeric(length(uniqueFreq$classes) / max(uniqueFreq$weight))
    uniqueClasses <- c(uniqueClasses, uniqueClass)
  }
  # Remove plateau resulting from k == classes.
  for (i in 2:length(sparam)) {
    if (length(sparam) == length(uniqueClasses)) {
      if (uniqueClasses[1] == uniqueClasses[2]) {
        uniqueClasses <- uniqueClasses[2:length(uniqueClasses)]
        sparam <- sparam[2:length(sparam)]
      }
    }
  }
  # Find the knee of the curve of uniqueClasses vs sparam.
  optimalS <- try(kneedle(x=sparam, y=uniqueClasses)[1], silent=TRUE)
  if (class(optimalS) != 'try-error') {
    if(plot) {
      sparsityPlot <- plotSparsity(sparam, uniqueClasses, optimalS)
    } else {
      sparsityPlot <- NULL
    }
    return(list('optimalS'=optimalS, 'sparsityPlot'=sparsityPlot))
  }
}

#' Optimize Sparsity (s) Parameter
#' 
#' Get the optimal parameter for s based on \code{spatialShrunkenCentroids} models
#' 
#' @param ssc \code{dataframe} with \code{summary()} information from \code{SpatialShrunkenCentroids2} objects
#' @param optimalS \code{integer} optimal value of s from \code{optimizeSParam}
#' @param sparam \code{vector} of values used for sparsity parameter in \code{spatialShrunkenCentroids}
#' @param rparam \code{vector} of values used for radius parameter in \code{spatialShrunkenCentroids}
#' @param kparam \code{vector} of values used for k parameter in \code{spatialShrunkenCentroids}
#' @return \code{dataframe} with row(s) of optimized parameter sets
#' 
#' @export
optimizeRKParams <- function(ssc, optimalS, rparam, kparam) {
  rk <- data.frame(r=double(), k=double(), slope=double(), knee=double())
  for (r in rparam) {
    for (k in kparam) {
      # Find the knee for each combination of r and k.
      kneeDf <- ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]
      kneeDfS <- kneeDf$s
      kneeDfClasses <- kneeDf$classes
      # Remove plateau resulting from k == classes.
      for (i in 2:length(kneeDfS)) {
        if (length(kneeDfS) == length(kneeDfClasses)) {
          if (kneeDfClasses[1] == kneeDfClasses[2]) {
            kneeDfClasses <- kneeDfClasses[2:length(kneeDfClasses)]
            kneeDfS <- kneeDfS[2:length(kneeDfS)]
          }
        }
      }
      # Only find knee if line is not horizontal.
      if (length(kneeDfS) > 1 & length(kneeDfClasses) > 1) {
        if (nrow(kneeDf) > 0) {
          knee <- try(kneedle(x=kneeDfS, y=kneeDfClasses)[1], silent=TRUE)
          if (class(knee) != 'try-error') {
            # Find slope * 100 of exponential curve for each combination of r and k.
            slope <- (lm(log(ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$classes) ~
                           ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$s)$coefficients[[2]]) * 100
            if (!is.na(slope) & !is.na(knee)) {
              rk <- rbind(rk, c(r, k, slope, knee))
            }
          }
        }
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
      sParamModel <- lm(log(classes) ~ ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$s)
      sParamSeq <- seq(min(sparam), max(sparam), 0.1)
      classesModel <- exp(predict(sParamModel, sparam=list(sParamSeq)))
      dx <- pracma::gradient(sparam)
      dy <- pracma::gradient(classes)
      d2x <- pracma::gradient(dx)
      d2y <- pracma::gradient(dy)
      curvature <- (abs(d2y) / (1 + (dy**2))**1.5) * 100
      optimalParams <- rbind(optimalParams, c(r, k, srk[which(srk$r==r & srk$k==k),]$slope,
                                              curvature[match(kneeValue, sparam)]))
    }
  }
  colnames(optimalParams) <- c('r', 'k', 'slope', 'curvature')
  
  # Remove combinations of r and k that have no slope or an increasing curve.
  optimalParams <- optimalParams[which(optimalParams$slope < 0),]
  # score = slope * curvature for each exponential curve.
  optimalParams$score <- abs(optimalParams$slope * optimalParams$curvature)
  # Choose the line where the combination of r and k has the best score.
  return(list('optimalParams'=optimalParams[which(optimalParams$score==max(optimalParams$score)),], 'rkScores'=optimalParams))
}

#'Plot of sparsity (s) parameter optimization
#'
#'Plot that shows predicted # of segments vs. shrinkage parameter s
#'
#'@param sscDf object output from \code{spatialShrunkenCentroids()} as a \code{dataframe}
#'@param plot is user defined either TRUE or FALSE
#'@return figure with plot of predicted # of segments for different s values
#'
#'@examples
#'
#'rparam <- c(1,2,3)
#'kparam <- c(2,4)
#'sparam <- c(0,3,6,9,12)
#'
#'ssc <- spatialShrunkenCentroids(data, r=rparam, k=kparam, s=sparam)
#'sscDf <- as.data.frame(Cardinal::summary(ssc))
#'SFig <- plotSSCLines(sscDf, plot=TRUE)
#'
#'@export
plotSSCLines <- function(sscDf, plot=TRUE){
  figure <- ggplot() +
    xlab('sparsity (s)') +
    ylab('predicted # of segments') +
    labs(colour='Line')
  
  #nested for loop that iterates through all combinations of r and k for each s value and adds a line to the plot
  i <- 1
  for (r in sort(unique(sscDf$r))) {
    for (k in sort(unique(sscDf$k))) {
      df <- sscDf[,c('s', 'classes')][which(sscDf$r==r & sscDf$k==k),]
      figure <- figure +
        geom_point(aes_(x=df$s, y=df$classes, colour = paste('r=', as.character(r), 'k=', as.character(k))), shape=i) +
        geom_line(aes_(x=df$s, y=df$classes,colour = paste('r=', as.character(r), 'k=', as.character(k))))
      i <- i + 1
    }
  }
  
  return(figure)
}

#' Optimize SpatialShrunkenCentroids Parameters
#' 
#' Optimize parameters used for Spatial Shrunken Centroids segmentation using
#' a heuristic algorithm.
#' 
#' @param x \code{SpatialShrunkenCentroids2} or \code{list} of \code{SpatialShrunkenCentroids} object
#' @param sparam \code{vector} of values used for sparsity parameter in \code{spatialShrunkenCentroids}
#' @param rparam \code{vector} of values used for radius parameter in \code{spatialShrunkenCentroids}
#' @param kparam \code{vector} of values used for k parameter in \code{spatialShrunkenCentroids}
#' @param plotLines can be TREUE or FALSE to define if the plot of s parameter optimization is required or not
#' @param plotS can be TREUE or FALSE to define if the plot showing weighted cardinality score vs s parameter is required or not
#' @param showRKScorescan be TREUE or FALSE to define if RK scores are required
#' @return \code{list} with optimal \code{spatialShrunkenCentroids} parameters for r, k, and s and figures from \code{sscLinesPlot} and
#' \code{plotSparcity} and rkScores
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
optimizeSSCParams <- function(x, sparam, rparam, kparam, plotLines=TRUE, plotS=TRUE, showRKScores=TRUE) {
  # Get vectors for s values and number of unique classes per s value.
  if (class(x) == "SpatialShrunkenCentroids2") {
    sscDf <- as.data.frame(Cardinal::summary(x))
  } else if (class(x) == "list") {
    sscDf <- getSSCDf(x)
  }
  colnames(sscDf) <- c('r', 'k', 's', 'classes', 'features_per_class')
  
  # Optimize sparsity (s) parameter.
  tmp <- optimizeSParam(sscDf, sparam, plotS)
  optimalS <- tmp$optimalS
  sparsityPlot <- tmp$sparsityPlot
  
  # Optimize r and k parameters.
  tmp2 <- optimizeRKParams(sscDf, optimalS, rparam, kparam)
  optimalR <- tmp2$optimalParams$r
  optimalK <- tmp2$optimalParams$k
  
  #plot of s parameter optimization
  if (plotLines) {
    sscLinesPlot <- plotSSCLines(sscDf, linesPlot)
  } else {
    sscLinesPlot <- NULL
  }
  
  if (showRKScores) {
    rkScoresDf <- tmp2$rkScores
  } else {
    rkScoresDf <- NULL
  }
  
  return(list('params'=list('r'=optimalR, 'k'=optimalK, 's'=optimalS), 'sparsityPlot'=sparsityPlot, 'sscLinesPlot'=sscLinesPlot, 'rkScoresDf'=rkScoresDf))
}