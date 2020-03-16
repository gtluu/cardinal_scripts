spatialShrunkenCentroidsWrapper <- function(dataset, r, k, s) {
  sscList <- list()
  inc <- 1
  for (r in rparam) {
    for (k in kparam) {
      for (s in sparam) {
        ssc <- try(spatialShrunkenCentroids(dataset, r=r, k=k, s=s))
        if (class(ssc) != 'try-error') {
          sscList[[inc]] <- ssc
          inc <- inc + 1
        }
      }
    }
  }
  return(sscList)
}