ionImageReport <- function(dataset, mz, filedir=getwd(), filename) {
  pdf(paste0(filename, '.pdf'))
  for (i in mz) {
    print(image(dataset, mz=iplusminus=0.2))
  }
  dev.off()
}
