#' Update Metadata
#' 
#' Update the metadata in an \code{MSImagingExperiment} to include
#' ROIs and user defined conditions for each ROI.
#' 
#' @param dataset \code{MSImagingExperiment} object
#' @param filedir Directory containing metadata files (default = current working directory)
#' @param spotFile Spot File containing XY coordinates and ROI names exported from Bruker flexImaging
#' @param roiFile CSV file created by the user containing ROIs and conditions for each ROI
#' @return \code{MSImagingExperiment} with updated metadata
#' @examples 
#' data <- as(readMSIData("data.imzML"), "MSContinuousImagingExperiment")
#' spots <- "spots.txt"
#' rois <- "rois.csv"
#' data <- updateMetadata(data, spotFile=spots, roiFile=rois)
#' 
#' @export
updateMetadata <- function(dataset, filedir=getwd(), spotFile, roiFile) {
  # Load in Spot File exported from flexImaging as dataframe.
  spots <- as.data.frame(read.table(paste0(filedir, '/', spotFile), sep=" ", col.names=c("x", "y", "spot", "roi")))
  spots$x <- as.integer(str_split(str_split(spots$spot, "X", simplify=TRUE)[,2], "Y", simplify=TRUE)[,1])
  spots$y <- as.integer(str_split(spots$spot, "Y", simplify=TRUE)[,2])
  spots$roi <- as.character(spots$roi)
  
  # Load in user created attribute file.
  rois <- read.csv(paste0(filedir, '/', roiFile), fileEncoding='UTF-8-BOM')
  
  # Exit with error if first column is not 'roi'.
  if (colnames(rois)[1] != 'roi') {
    stop("Attribute table should begin with column 'roi'.")
  }
  
  # Add columns for each condition to Spot File dataframe.
  for (cond in 2:length(rois)) {
    spots[names(rois)[cond]] <- spots$roi
  }
  for (i in 1:nrow(rois)) {
    roi <- as.list(rois[i,])
    for (cond in 2:length(roi)) {
      spots[[names(roi)[cond]]][spots[[names(roi)[cond]]] == as.character(roi[1][[1]])] <- as.character(roi[cond][[1]])
    }
  }
  
  # Use inner fuzzy join to reorder rows correctly according to xy coordinates.
  pd <- as.data.frame(pixelData(dataset))
  pd <- difference_inner_join(pd, spots, by=c('x','y'), max_dist=0)
  
  # Add roi and condition columns to pixelData.
  for (cond in 1:length(rois)) {
    pixelData(dataset)[[names(rois)[cond]]] <- as.factor(pd[[names(rois)[cond]]])
  }
  
  # Add column to metadata to reflect run + roi.
  pixelData(dataset)$run <- as.factor(paste(run(dataset), pixelData(dataset)$roi, sep='_'))
  
  # Return dataset with updated metadata.
  return(dataset)
}
