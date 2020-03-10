optimize_ssc_params <- function(ssc, sparam) {
  # optimize s param
  ssc <- as.data.frame(summary(ssc))
  colnames(ssc) <- c('r', 'k', 's', 'classes', 'features_per_class', 'BIC')
  unique_classes <- c()
  for (s in sparam) {
    unique_freq <- as.data.frame(table(sort(ssc[,c('s', 'classes')][which(ssc$s==s),]$classes)))
    colnames(unique_freq) <- c('classes', 'freq')
    unique_freq$weight <- (as.numeric(unique_freq$freq) / as.numeric(sum(unique_freq$freq)))
    unique_classes <- c(unique_classes, as.numeric(length(unique_freq$classes) / max(unique_freq$weight)))
  }
  skneeparam <- 1
  while(TRUE) {
    knee <- kneed$KneeLocator(x=sparam, y=unique_classes, S=skneeparam, curve='convex', direction='decreasing')$knee
    if (!is.null(knee)) {
      optimal_s <- knee
      break
    }
    skneeparam <- skneeparam - 0.1
    if (skneeparam <= 0) {
      return('No s parameter identified')
    }
  }
  # optimize r and k params
  rk_df <- data.frame(r=double(), k=double(), slope=double(), knee=double())
  for (r in rparam) {
    for (k in kparam) {
      # find knee for each r, k combination
      skneeparam <- 1
      while(TRUE) {
        knee <- kneed$KneeLocator(x=ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$s,
                                  y=ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$classes,
                                  S=skneeparam, curve='convex', direction='decreasing')$knee
        if (!is.null(knee)) {
          break
        }
        skneeparam <- skneeparam - 0.1
        if (skneeparam <= 0) {
          break
        }
      }
      if (is.null(knee)) {
        knee <- NA
      }
      # find slope*100 of exponential curve for each r, k combination
      slope <- (lm(log(ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$classes) ~
                     ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$s)$coefficients[[2]]) * 100
      if (!is.na(slope)) {
        rk_df <- rbind(rk_df, c(r, k, slope, knee))
      }
    }
  }
  colnames(rk_df) <- c('r', 'k', 'slope', 'knee')
  # find r, k combinations with s == optimal_s
  # if no combination found, look for s + next increment
  knee_value <- optimal_s
  while(TRUE) {
    optimal_rk_df <- rk_df[which(rk_df$knee==knee_value),]
    if (nrow(optimal_rk_df) != 0) {
      break
    }
    knee_value <- knee_value + 3
  }
  # find curvature*100 (2nd derivative of exponential curve) at s == optimal_s
  optimal_rk <- data.frame(r=double(), k=double(), slope=double(), curvature=double())
  for (r in unique(optimal_rk_df$r)) {
    for (k in unique(optimal_rk_df$k)) {
      classes <- ssc[,c('s', 'classes')][which(ssc$r==r & ssc$k==k),]$classes
      sparam_model <- lm(log(classes) ~ sparam)
      sparam_seq <- seq(min(sparam), max(sparam), 0.1)
      classes_model <- exp(predict(sparam_model, sparam=list(sparam_seq)))
      dx <- numpy$gradient(sparam)
      dy <- numpy$gradient(classes)
      d2x <- numpy$gradient(dx)
      d2y <- numpy$gradient(dy)
      curvature <- (numpy$abs(d2y) / (1 + (dy**2))**1.5) * 100
      optimal_rk <- rbind(optimal_rk, c(r, k, optimal_rk_df[which(optimal_rk_df$r==r & optimal_rk_df$k==k),]$slope,
                                        curvature[match(knee_value, sparam)]))
    }
  }
  colnames(optimal_rk) <- c('r', 'k', 'slope', 'curvature')
  optimal_rk <- optimal_rk[which(optimal_rk$slope < 0),]
  # score = slope * curvature for each exponential curve
  optimal_rk$score <- numpy$abs(optimal_rk$slope * optimal_rk$curvature)
  print(optimal_rk)
  # choose line where r, k combination has best score
  optimal_rk <- optimal_rk[which(optimal_rk$score==max(optimal_rk$score)),]
  optimal_r <- optimal_rk$r
  optimal_k <- optimal_rk$k
  return(list('r'=optimal_r, 'k'=optimal_k, 's'=optimal_s))
}
