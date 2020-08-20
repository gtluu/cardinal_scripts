#' Spatial Shrunken Centroids Segmentation Image Report
#' 
#' Generate a \code{PDF} file containing all segmentaiton image(s) from a \code{SpatialShrunkenCentroids2} object.
#' 
#' @param ssc \code{SpatialShrunkenCentroids2} to generate segmentation images from.
#' @param filedir Output file directory.
#' @param filename Output PDF filename.
#' @param ... Parameters to be passed to \code{Cardinal::image()}
#' @return \code{NULL}; a PDF file will be in the specified location.
#' @example 
#' 
#' sscImageReport(processedData, filename='three_images')
#' 
#' sscImageReport(processedData, filename='three_images', layout=c(2, 2), xlim=c(0, 300), ylim=c(0, 150))
#'
#' @export
sscImageReport <- function(ssc, filedir=getwd(), filename, ...) {
  pdf(paste0(filename, '.pdf'))
  for (i in nrow(modelData(ssc))) {
    print(image(ssc, model=as.list(modelData(ssc)[i,]), ...))
  }
  dev.off()
}
