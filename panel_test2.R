setwd("/Users/taylorstone/Desktop/TFM")
#pEPA Test

#saved files
forecasted_ar = read.csv("forecasted_ar", stringsAsFactors = FALSE)
forecasted_ar_gif = read.csv("forecasted_ar_gif", stringsAsFactors = FALSE)
forecasted_ucsv = read.csv("forecasted_ucsv", stringsAsFactors = FALSE)
forecasted_ucsv_gif = read.csv("forecasted_ucsv_gif", stringsAsFactors = FALSE)
inflation_arima = read.csv("Seasadj_log_inflation_rates.csv", stringsAsFactors = FALSE)

#drop date columns
forecasted_ar$X = as.Date(forecasted_ar$X)
forecasted_ar_gif$X = as.Date(forecasted_ar_gif$X)
forecasted_ucsv$X = as.Date(forecasted_ucsv$X)
forecasted_ucsv_gif$X = as.Date(forecasted_ucsv_gif$X)



force_numeric_matrix <- function(x) {
  result <- matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x))
  rownames(result) <- rownames(x)
  colnames(result) <- colnames(x)
  return(result)
}

panel_test_fixed = function(original_values,
                            forecasted_ar,
                            forecasted_ar_gif,
                            forecasted_ucsv,
                            forecasted_ucsv_gif) {
  
  # Force numeric matrices
  original_values <- force_numeric_matrix(as.matrix(original_values))
  forecasted_ar <- force_numeric_matrix(as.matrix(forecasted_ar))
  forecasted_ar_gif <- force_numeric_matrix(as.matrix(forecasted_ar_gif))
  forecasted_ucsv <- force_numeric_matrix(as.matrix(forecasted_ucsv))
  forecasted_ucsv_gif <- force_numeric_matrix(as.matrix(forecasted_ucsv_gif))
  
  # Verify they're numeric now
  cat("All numeric?", all(sapply(list(original_values, forecasted_ar, 
                                      forecasted_ar_gif, forecasted_ucsv, 
                                      forecasted_ucsv_gif), is.numeric)), "\n")
  
  # Rest of your function...
  models = list(
    ARIMA      = forecasted_ar,
    ARIMA_GIF  = forecasted_ar_gif,
    UCSV       = forecasted_ucsv,
    UCSV_GIF   = forecasted_ucsv_gif
  )
  
  pairs = combn(names(models), 2, simplify = FALSE)
  results = lapply(pairs, function(p) {
    loss1 = (original_values - models[[p[1]]])
    loss2 = (original_values - models[[p[2]]])
    realized   = original_values
    
    s1 = pool_av.S1.test(loss1, loss2, realized, loss.type = "SE")
    s3 = pool_av.S3.test(loss1, loss2, realized, loss.type = "SE")
    
    data.frame(
      Model1 = p[1],
      Model2 = p[2],
      S1_stat = s1$statistic,
      S1_pval = s1$p.value,
      S3_stat = s3$statistic,
      S3_pval = s3$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}


