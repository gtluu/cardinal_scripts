#' Get Statistic Table 
#' 
#' Returns statistic table with optimized r,k and s parameters using 
#' \code{getStatisticTable}
#' 
#' @param sscObject \code{SpatialShrunkenCentroids2} object output from \code{spatialShrunkenCentroids()}
#' @param sscParams \code{named list} containing the optimized r, k, and s values
#' @return \code{dataframe} containing features and associated statistics for each segment
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
#' getStatisticTable(ssc, sscParams$params)
#' 
#' @export
getStatisticTable <- function(sscObject, sscParams) {
  # df is a dataframe that contains a complete list of the set of r, k, and s parameters tested
  df <- modelData(sscObject)
  rownames(df) <- c(1:nrow(df))
  
  # index stores the index of the model with the optimized r, k, and s values
  index <- as.numeric(rownames(df)[df$r==sscParams$r & df$k==sscParams$k & df$s==sscParams$s])
  
  # statTable stores the statistic table corresponding to the model with the optimiized r, k, and s values
  statTable <- as.data.frame(resultData(sscObject)[[index]]$statistic)
  return(statTable)
}
