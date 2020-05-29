ionImageReport <- function(dataset, mz, filedir=getwd(), filename) {
  pdf(paste0(filename, '.pdf'))
  for (i in mz) {
    print(image(dataset, mz=i, plusminus=0.2))
  }
  dev.off()
}
