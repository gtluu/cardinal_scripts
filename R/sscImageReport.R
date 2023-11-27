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
#' sscImageReport(ssc, filename='three_images')
#' 
#' sscImageReport(ssc, filename='three_images', layout=c(2, 2), xlim=c(0, 300), ylim=c(0, 150))
#'
#' @export
sscImageReport <- function(ssc, filedir=getwd(), filename, ...) {
  pdf(paste0(filedir, '/', filename, '.pdf'))
  for (i in 1:nrow(modelData(ssc))) {
    print(paste0('Saving SSC image for',
                 ' r=', as.character(modelData(ssc)[i,]$r),
                 ' k=', as.character(modelData(ssc)[i,]$k),
                 ' s=', as.character(modelData(ssc)[i,]$s)))
    print(image(ssc, model=as.list(modelData(ssc)[i,]), ...))
  }
  dev.off()
}
