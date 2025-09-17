################
#### ARIMA ####
################
forecast_arima = function(y, frequency = 12) {
  if (!is.numeric(y)) stop("Input must be a numeric vector")
  
  # Parameters
  R = (2004 - 1991 + 1) * 12  # estimation window length (168 obs)
  OOS = (2023 - 2005 + 1) * 12  # out-of-sample periods (228 obs)
  
  # Storage vector for forecasts
  oos_forecasts = numeric(OOS)
  
  #progress bar
  pb = txtProgressBar(min = 0, max = OOS, style = 3)
  
  # Rolling window loop
  for (row in 1:OOS) {
    training_data = y[(1 + row - 1):(R + row - 1)]
    ts_training = ts(training_data, 
                      start = c(1991, 1), 
                      frequency = frequency)
    
    model_ARIMA = auto.arima(ts_training,
                              seasonal = TRUE,
                              stepwise = TRUE,
                              approximation = TRUE,
                              trace = FALSE)
    
    fc = forecast(model_ARIMA, h = 1)
    oos_forecasts[row] = fc$mean
    
    #for progress bar to update
    setTxtProgressBar(pb, row)
  }
  
  close(pb)
  return(oos_forecasts)
}


################
## ARIMA + GIF##
################
forecast_arima_gif = function(y, global_factor, frequency = 12) {
  library(forecast)
  
  if (!is.numeric(y)) stop("Input must be a numeric vector")
  if (!is.numeric(global_factor)) stop("Global factor must be a numeric vector")
  if (length(global_factor) != length(y)) stop("Global factor must be the same length as the series")
  
  # Parameters
  R = (2004 - 1991 + 1) * 12  # estimation window length (168 obs)
  OOS = (2023 - 2005 + 1) * 12  # out-of-sample periods (228 obs)
  
  # Storage vector for forecasts
  oos_forecasts = numeric(OOS)
  
  # Progress bar
  pb = txtProgressBar(min = 0, max = OOS, style = 3)
  
  # Rolling window loop
  for (row in 1:OOS) {
    # Training windows for y and global factor
    training_y   = y[(1 + row - 1):(R + row - 1)]
    training_gf  = global_factor[(1 + row - 1):(R + row - 1)]
    
    ts_training = ts(training_y, 
                      start = c(1991, 1), 
                      frequency = frequency)
    
    # Fit ARIMA with xreg = training global factor
    model_ARIMA = auto.arima(ts_training,
                              xreg = training_gf,
                              seasonal = TRUE,
                              stepwise = TRUE,
                              approximation = TRUE,
                              trace = FALSE)
    
    # Forecast 1-step ahead using *next period’s* global factor
    next_gf = global_factor[R + row]  # the regressor for the OOS point
    fc = forecast(model_ARIMA, h = 1, xreg = next_gf)
    
    oos_forecasts[row] = fc$mean
    
    # Update progress bar
    setTxtProgressBar(pb, row)
  }
  
  close(pb)
  return(oos_forecasts)
}

################
#### UCSVAR ####
################
forecast_ucsv = function(y, nloop = 1000, burnin = 200, frequency = 12) {
  source("UCSV.R")
  if (!is.numeric(y)) stop("Input must be a numeric vector")
  
  # Parameters
  R = (2004 - 1991 + 1) * 12  # initial estimation window
  OOS = (2023 - 2005 + 1) * 12  # out-of-sample periods
  forecast_vec = numeric(OOS)
  
  pb = txtProgressBar(min = 0, max = OOS, style = 3)
  
  # Rolling window
  for (row in 1:OOS) {
    training_data = y[(1 + row - 1):(R + row - 1)]
    
    # Fit UCSVAR on the training data
    fit = fit_ucsv(training_data, nloop = nloop, burnin = burnin, plot_results = FALSE)
    
    # Forecast next step as tauhat + chat for the next time point
    # Since UCSVAR is univariate, we approximate next point as last fitted value
    forecast_vec[row] = tail(fit$fitted_values, 1)
    
    setTxtProgressBar(pb, row)
  }
  
  close(pb)
  return(forecast_vec)
}

################
# UCSVAR + GIF #
################
forecast_ucsv_gif = function(y, global_factor, nloop = 1000, burnin = 200, frequency = 12) {
  if (!is.numeric(y)) stop("Input must be a numeric vector")
  if (!is.numeric(global_factor)) stop("Global factor must be numeric")
  if (length(global_factor) != length(y)) stop("y and global factor must be same length")
  
  # Parameters
  R = (2004 - 1991 + 1) * 12  # initial estimation window
  OOS = (2023 - 2005 + 1) * 12  # out-of-sample periods
  forecast_vec = numeric(OOS)
  
  pb = txtProgressBar(min = 0, max = OOS, style = 3)
  
  for (row in 1:OOS) {
    training_y = y[(1 + row - 1):(R + row - 1)]
    training_gf = global_factor[(1 + row - 1):(R + row - 1)]
    
    # Fit UCSVAR with global factor
    fit = fit_ucsv_gif(training_y, training_gf, nloop = nloop, burnin = burnin, plot_results = FALSE)
    
    # One-step-ahead forecast:
    next_gf = global_factor[R + row]  # global factor at t+1
    forecast_vec[row] = tail(fit$tauhat, 1) +
      tail(fit$chat, 1) +
      fit$betahat * next_gf
    
    setTxtProgressBar(pb, row)
  }
  
  close(pb)
  return(forecast_vec)
}

##################
# Plot Forecasts #
##################
plot_forecasts = function(forecast_matrix, method){
  for (i in 1:n_countries) {
    country <- g7_countries[i]
    
    plot_df <- data.frame(
      Date = as.Date(dates_oos),
      Actual = data_g7[dates_oos, country],
      Forecast = forecast_matrix[, country]
    )
    
    p <- ggplot(plot_df, aes(x = Date)) +
      geom_line(aes(y = Actual, color = "Actual"), lwd = 1) +
      geom_line(aes(y = Forecast, color = "Forecast"), lwd = 1) +
      theme_minimal() +
      labs(
        title = paste(method, "Forecast vs Actual:", country),
        y = "Inflation",
        x = "Date",
        color = ""
      )
    
    print(p)
  }
}

#######################
## Error Computation ##
#######################
compute_errors <- function(forecast_mat, actual_mat) {
  errors <- forecast_mat - actual_mat
  rmse <- sqrt(colMeans(errors^2, na.rm = TRUE))
  mae <- colMeans(abs(errors), na.rm = TRUE)
  list(rmse = rmse, mae = mae)
}