#' Segmentation Test Wrapper
#' 
#' Wrapper for \code{segmentationTest} following \code{spatialDGMM} to find discriminate
#' features while ignoring problematic features causing error.
#' 
#' @param sdgmmList \code{list} of \code{SpatialDGMM} objects to be analyzed.
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
#' sdgmmList <- spatialDGMMWrapper(data, r=1, k=4)
#' 
#' featuresDf <- segmentationTestWrapper(sdgmmList, r=1, k=4, fixedCondition='treated',
#'                                       classControl='Ymax')
#' 
#' @export
segmentationTestWrapper <- function(sdgmmList, fixedCondition, classControl='Ymax',
                                    BPPARAM=bpparam()) {
  # Initialize data.frame.
  df <- data.frame(mz=double(), r=double(), k=double(), feature=double(),
                   LR=double(), PValue=double(), AdjP=double())
  
  # Run segmentationTest for each successful spatial-DGMM run.
  for (i in 1:length(sdgmmList)) {
    if (class(sdgmmList[[i]]) != 'try-error') {
      segTest <- segmentationTest(sdgmmList[[i]], as.formula(paste('~', fixedCondition)),
                                  classControl=classControl, BPPARAM=BPPARAM)
      segTestDf <- as.data.frame(topFeatures(segTest))
      segTestDf$feature <- c(i)
      df <- rbind(df, segTestDf)
    }
  }
  
  return(df)
}
