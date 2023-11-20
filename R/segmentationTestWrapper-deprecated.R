#' Segmentation Test Wrapper
#' 
#' Wrapper for \code{segmentationTest} following \code{spatialDGMM} to find discriminate
#' features while ignoring problematic features causing error.
#' 
#' @param sdgmmList \code{list} of \code{SpatialDGMM} objects to be analyzed.
#' @param fixedCondition A one-sided formula giving the fixed effects of the model on the RHS.
#' The response will added to the LHS, and the formula will be passed to the underlying modeling
#' function.
#' @param ... Parameters to be passed to \code{segmentationTest()}
#' @return \code{dataframe} containing features and their associated statistics.
#' @examples
#' 
#' sdgmmList <- spatialDGMMWrapper(data, r=1, k=4)
#' 
#' featuresDf <- segmentationTestWrapper(sdgmmList, r=1, k=4, fixedCondition='treated',
#'                                       classControl='Ymax')
#' 
#' @export
segmentationTestWrapper <- function(sdgmmList, fixedCondition, ...) {
  # Initialize data.frame.
  df <- data.frame(mz=double(), r=double(), k=double(), feature=double(),
                   LR=double(), PValue=double(), AdjP=double())
  
  # Run segmentationTest for each successful spatial-DGMM run.
  for (i in 1:length(sdgmmList)) {
    if (class(sdgmmList[[i]]) != 'try-error') {
      segTest <- segmentationTest(sdgmmList[[i]], as.formula(paste('~', fixedCondition)), ...)
      segTestDf <- as.data.frame(topFeatures(segTest))
      segTestDf$feature <- c(i)
      df <- rbind(df, segTestDf)
    }
  }
  
  return(df)
}
