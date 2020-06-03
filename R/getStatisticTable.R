#' Get Statistic Table 
#' 
#' Returns statistic table with optimized r,k and s parameters using 
#' \code{getStatisticTable}
#' 
#' @param sscObject should be the output from \code{spatialShrunkenCentroids}
#' @param sscParams should be a list of the three optimized r,k and s values. 
#' @param df is a data frame that contains a complete list of the set of r,k,s 
#' parameters tested
#' @param index stores the index of the model with the optimized r,k and s values
#' @param statTable stores the statistic table corresponding to the model with 
#' the optimized  r,k and s values
#' @return \code{getStatisticTable} containing features associated with statistics
#' @examples 
#' 
#' getStatisticTable(ssc, sscParams)
#' 
#' @export

getStatisticTable <- function(sscObject, sscParams) 
{
  df <- modelData(sscObject)
  rownames(df) <- c(1:nrow(df))
  index <- as.numeric(rownames(df)[df$r==sscParams[1] & df$k==sscParams[2] & df$s==sscParams[3]])
  statTable <- as.data.frame(resultData(sscObject)[[index]]$statistic)
  return(statTable)
}