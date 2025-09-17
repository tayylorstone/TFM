# adapted from https://pkg.robjhyndman.com/forecast/reference/dm.test.html#references
################
# DM Test Loop #
################
dm_test2 = function(country, original_values,
                    forecasted_ar,
                    forecasted_ar_gif,
                    forecasted_ucsv,
                    forecasted_ucsv_gif
) {
  
  # Common index from original_values
  common_index = index(original_values)
  
  # Coerce everything to zoo with same index
  models = list(
    ARIMA      = zoo(forecasted_ar[, country], common_index),
    ARIMA_GIF  = zoo(forecasted_ar_gif[, country], common_index),
    UCSV       = zoo(forecasted_ucsv[, country], common_index),
    UCSV_GIF   = zoo(forecasted_ucsv_gif[, country], common_index)
  )
  
  orig = zoo(original_values[, country], common_index)
  
  # Align all series together and drop NAs
  all_data = na.omit(merge(orig, models$ARIMA, models$ARIMA_GIF, models$UCSV, models$UCSV_GIF))
  colnames(all_data) = c("orig", names(models))
  
  # regenerate models list with aligned data
  models = as.list(all_data[, -1])
  
  # generate pairs
  pairs = combn(names(models), 2, simplify = FALSE)
  
  # Run DM tests
  results = lapply(pairs, function(p) {
    e1 = all_data$orig - models[[p[1]]]
    e2 = all_data$orig - models[[p[2]]]
    
    dm_out = dm.test(as.numeric(e1), as.numeric(e2),
                      alternative = "two.sided", h = 1, power = 2)
    
    data.frame(
      Country = country,
      Model1 = p[1],
      Model2 = p[2],
      DM_stat = as.numeric(dm_out$statistic),
      p_value = as.numeric(dm_out$p.value),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, results)
}
