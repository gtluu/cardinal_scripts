#' Generate Reference Peak List
#' 
#' Generate reference peak list from multiple \code{MSImagingExperiment} objects
#' to align and combine multiple datasets.
#' 
#' @param listOfDatasets List of \code{MSImagingExperiment} objects
#' @param tol m/z tolerance for aligning peaks
#' @return sorted \code{vector} of features
#' @examples 
#' refPeaks <- generateRefPeaksList(list(data1, data2), tol=0.2)
#' 
#' @export
generateRefPeaksList <- function(listOfDatasets, tol=0.2) {
  # Get feature dataframes from each dataset.
  peaks <- list()
  for (i in 1:length(listOfDatasets)) {
    peaks[[i]] <- data.frame('peaks'=mz(listOfDatasets[[i]]))
  }
  
  # Concatenate all feature lists into one dataframe (one long feature list) and sort.
  df <- data.frame('peaks'=c())
  for (i in 1:length(peaks)) {
    df <- rbind(df, peaks[[i]])
  }
  df$peaks <- sort(df$peaks)
  
  # Fuzzy inner join master feature list with itself.
  df2 <- difference_inner_join(df, df, by='peaks', max_dist=tol)
  colnames(df2) <- c('peaks', 'average')
  
  # Aggregate/average rows with similar features to generate list of averaged features.
  df3 <- aggregate(df2, by=list(df2$peaks), FUN=mean)
  
  # Return sorted list of unique features.
  return(sort(unique(df3$average)))
}
