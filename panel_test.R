#######
# Panel Test 
#######
# theory: https://www.sciencedirect.com/science/article/pii/S0169207023000092
# code adapted from: https://github.com/kdrachal/pEPA/blob/main/R/tests.R
#######

library(pEPA)

panel_test = function(original_values,
                      forecasted_ar,
                      forecasted_ar_gif,
                      forecasted_ucsv,
                      forecasted_ucsv_gif) {
  # Put forecasts into a list
  models = list(
    ARIMA      = forecasted_ar,
    ARIMA_GIF  = forecasted_ar_gif,
    UCSV       = forecasted_ucsv,
    UCSV_GIF   = forecasted_ucsv_gif
  )
  
  # Generate pairs of models
  pairs = combn(names(models), 2, simplify = FALSE)
  
  results = lapply(pairs, function(p) {
    # Loss differential (absolute error here, could be squared error too)
    loss1 = (original_values - models[[p[1]]])
    loss2 = (original_values - models[[p[2]]])
    
    # AE loss
    evaluated1 = loss1^2
    evaluated2 = loss2^2
    realized   = original_values
    
    # Panel tests
    s1 = pool_av.S1.test(evaluated1, evaluated2, realized, loss.type = "SE")
    s3 = pool_av.S3.test(evaluated1, evaluated2, realized, loss.type = "SE")
    
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