getStatisticTable <- function(sscObject, sscParams) 
{
  df <- modelData(sscObject)
  rownames(df) <- c(1:nrow(df))
  index <- as.numeric(rownames(df)[df$r==sscParams[1] & df$k==sscParams[2] & df$s==sscParams[3]])
  statTable <- as.data.frame(resultData(sscObject)[[index]]$statistic)
  return(statTable)
}