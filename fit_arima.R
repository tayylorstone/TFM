# Function to fit ARIMA to Time Series #

### No factors ###
fit_arima = function(y, frequency=12){
  #load required libraries
  library(forecast)
  
  if (!is.numeric(y)) stop("Input must be a numeric vector")
  
  #Create date sequence from 1991-01-01 to 2023-12-01
  start_date = as.Date("1991-01-01")
  end_date = as.Date("2023-12-01")
  
  ts_start = ts(y,
                start=c(1991,1),
                frequency = frequency)
  
  arima_model = auto.arima(ts_start,
                           seasonal = TRUE,
                           stepwise = FALSE,
                           approximation = FALSE,
                           trace = FALSE)
  
  #Results
  result = list(
    model = arima_model,
    time_series = ts_start,
    summary = summary(arima_model),
    fitted_values = fitted(arima_model),
    residuals = residuals(arima_model)
  )
  
  #Fitted Values
  new_ts = fitted(arima_model)
  #Residuals
  new_residuals = residuals(arima_model)
  
  #Model Summary
  cat("ARIMA Model Summary:\n")
  print(arima_model)
  cat("Date range: from", as.character(start_date), "to", as.character(end_date), "\n")
  
  #Plots
  par(mfrow = c(1, 2))
  
  #Fitted vs. Actual
  plot(ts_start, main = "Fitted vs Actual Values", 
       ylab = "Value", xlab = "Year", col = "cornflowerblue", lwd = 1.5)
  lines(result$fitted_values, col = "coral2", lwd = 2)
  legend("topleft", c("Actual", "Fitted"), col = c("cornflowerblue", "coral2"), 
         lty = 1, lwd = c(1.5, 2), bg = "white")
  grid(col = "lightgray", lty = "dotted")
  
  #ACF of Residuals
  acf(result$residuals, main = "ACF of Residuals", 
      col = "cornflowerblue", lwd = 2)
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
  
  
  return(result)
}

####################
### With Factors ###
####################

fit_arima_gif = function(y, global_factor, frequency = 12){
  # load required libraries
  library(forecast)
  
  if (!is.numeric(y)) stop("Input must be a numeric vector")
  if (!is.numeric(global_factor)) stop("Global factor must be a numeric vector")
  if (length(global_factor) != length(y)) stop("Global factor must be the same length as the series")
  
  # Create date sequence
  start_date = as.Date("1991-01-01")
  end_date   = as.Date("2023-12-01")
  
  ts_start = ts(y, start = c(1991,1), frequency = frequency)
  
  #Fit ARIMA with global_factor as exogenous regressor
  arima_model = auto.arima(ts_start,
                           xreg = global_factor,
                           seasonal = TRUE,
                           stepwise = FALSE,
                           approximation = FALSE,
                           trace = FALSE)
  
  #Results
  result = list(
    model = arima_model,
    time_series = ts_start,
    summary = summary(arima_model),
    fitted_values = fitted(arima_model),
    residuals = residuals(arima_model)
  )
  
  #Model Sumary
  cat("ARIMA + global_factor Model Summary:\n")
  print(arima_model)
  
  cat("Date range: from", as.character(start_date), "to", as.character(end_date), "\n")
  
  #Plots
  par(mfrow = c(1, 2))
  
  
  #Fitted vs Actual
  plot(ts_start, main = "Fitted vs Actual Values", 
       ylab = "Value", xlab = "Year", col = "cornflowerblue", lwd = 1.5)
  lines(result$fitted_values, col = "coral2", lwd = 2)
  legend("topleft", c("Actual", "Fitted"), col = c("cornflowerblue", "coral2"), 
         lty = 1, lwd = c(1.5, 2), bg = "white")
  grid(col = "lightgray", lty = "dotted")
  
  
  #ACF of residuals
  acf(result$residuals, main = "ACF of Residuals", 
      col = "cornflowerblue", lwd = 2)
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
  
  return(result)
}
  
