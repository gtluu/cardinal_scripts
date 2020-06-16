#' Ion Image Report
#' 
#' Generate a \code{PDF} file containing ion image(s) from a single or \code{vector} of m/z values.
#' 
#' @param dataset \code{MSImagingExperiment} to generate ion images from.
#' @param mz \code{numeric} or \code{vector} of \code{numerics} containing m/z values.
#' @param filedir Output file directory.
#' @param filename Output PDF filename.
#' @param ... Parameters to be passed to \code{Cardinal::image()}
#' @return \code{NULL}; a PDF file will be in the specified location.
#' @example 
#' 
#' ionImageReport(processedData, mz=c(211.2, 269.7, 609.8), filename='three_images')
#' 
#' ionImageReport(processedData, mz=c(211.2, 269.7, 609.8), filename='three_images', plusminus=0.2,
#'                normalize.image='linear', contrast.enhance='histogram')
#'
#' @export
ionImageReport <- function(dataset, mz, filedir=getwd(), filename, ...) {
  pdf(paste0(filename, '.pdf'))
  for (i in mz) {
    print(image(dataset, mz=i, ...))
  }
  dev.off()
}
