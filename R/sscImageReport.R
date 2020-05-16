sscImageReport <- function(ssc, filedir=getwd(), filename) {
  pdf(paste0(filename, '.pdf'))
  for (i in nrow(modelData(ssc))) {
    print(image(ssc, model=as.list(modelData(ssc)[i,])))
  }
  dev.off()
}
