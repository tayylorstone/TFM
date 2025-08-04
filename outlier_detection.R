### This function follows the method described in McCracken and Ng (2016) which not only detects outliers but
### replaces them. For each observation, we will determine if it deviates from the median more than ten times
### the interquartile range. If found outside the range, the function will replace it with the value that is
### exactly 10 IQRs above/below the median (depending on the original position of the observation).

handle_outliers = function(series) {
  M = median(series, na.rm = TRUE)
  IR = IQR(series, na.rm = TRUE)
  threshold = 10 * IR
  
  # Calculate deviations
  deviations = series - M
  
  # Find outlier indices
  outlier_idx = which(abs(deviations) > threshold)
  
  # Replace outliers
  series[outlier_idx] = M + threshold * sign(deviations[outlier_idx])
  
  return(series)
}
