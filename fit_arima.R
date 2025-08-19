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
  
  # Create a comprehensive result list
  result <- list(
    model = arima_model,
    time_series = ts_start,
    summary = summary(arima_model),
    fitted_values = fitted(arima_model),
    residuals = residuals(arima_model),
    aic = AIC(arima_model),
    bic = BIC(arima_model),
    start_date = start_date,
    end_date = end_date,
    frequency = frequency,
    data_length = length(y)
  )
  
  # Store fitted values
  new_ts = fitted(arima_model)
  # Store residuals
  new_residuals = residuals(arima_model)
  
  # Print model summary
  cat("ARIMA Model Summary:\n")
  cat("===================\n")
  print(arima_model)
  cat("\nModel AIC:", result$aic, "\n")
  cat("Model BIC:", result$bic, "\n")
  cat("Data frequency:", frequency, "\n")
  cat("Time series length:", length(y), "\n")
  cat("Date range: from", as.character(start_date), "to", as.character(end_date), "\n")
  
  # Create diagnostic plots
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  # Plot 1: Original time series
  plot(ts_start, main = "Original Time Series", 
       ylab = "Value", xlab = "Year", col = "blue", lwd = 1.5)
  grid(col = "lightgray", lty = "dotted")
  
  # Plot 2: Fitted vs Actual
  plot(ts_start, main = "Fitted vs Actual Values", 
       ylab = "Value", xlab = "Year", col = "blue", lwd = 1.5)
  lines(result$fitted_values, col = "red", lwd = 2)
  legend("topleft", c("Actual", "Fitted"), col = c("blue", "red"), 
         lty = 1, lwd = c(1.5, 2), bg = "white")
  grid(col = "lightgray", lty = "dotted")
  
  # Plot 3: Residuals over time
  plot(result$residuals, main = "Residuals Over Time", 
       ylab = "Residuals", xlab = "Year", col = "darkred", pch = 16, cex = 0.7)
  abline(h = 0, col = "black", lty = 2, lwd = 2)
  grid(col = "lightgray", lty = "dotted")
  
  # Plot 4: ACF of residuals
  acf(result$residuals, main = "ACF of Residuals", 
      col = "darkblue", lwd = 2)
  
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
  
  # Fit ARIMA with global_factor as exogenous regressor
  arima_model = auto.arima(ts_start,
                           xreg = global_factor,
                           seasonal = TRUE,
                           stepwise = FALSE,
                           approximation = FALSE,
                           trace = FALSE)
  
  # Create a comprehensive result list
  result <- list(
    model = arima_model,
    time_series = ts_start,
    summary = summary(arima_model),
    fitted_values = fitted(arima_model),
    residuals = residuals(arima_model),
    aic = AIC(arima_model),
    bic = BIC(arima_model),
    start_date = start_date,
    end_date = end_date,
    frequency = frequency,
    data_length = length(y),
    global_factor = global_factor
  )
  
  # Print model summary
  cat("ARIMA + global_factor Model Summary:\n")
  cat("==========================\n")
  print(arima_model)
  cat("\nModel AIC:", result$aic, "\n")
  cat("Model BIC:", result$bic, "\n")
  cat("Data frequency:", frequency, "\n")
  cat("Time series length:", length(y), "\n")
  cat("Date range: from", as.character(start_date), "to", as.character(end_date), "\n")
  
  # Diagnostic plots
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  # Plot 1: Original time series
  plot(ts_start, main = "Original Time Series", 
       ylab = "Value", xlab = "Year", col = "blue", lwd = 1.5)
  grid(col = "lightgray", lty = "dotted")
  
  # Plot 2: Fitted vs Actual
  plot(ts_start, main = "Fitted vs Actual Values", 
       ylab = "Value", xlab = "Year", col = "blue", lwd = 1.5)
  lines(result$fitted_values, col = "red", lwd = 2)
  legend("topleft", c("Actual", "Fitted"), col = c("blue", "red"), 
         lty = 1, lwd = c(1.5, 2), bg = "white")
  grid(col = "lightgray", lty = "dotted")
  
  # Plot 3: Residuals over time
  plot(result$residuals, main = "Residuals Over Time", 
       ylab = "Residuals", xlab = "Year", col = "darkred", pch = 16, cex = 0.7)
  abline(h = 0, col = "black", lty = 2, lwd = 2)
  grid(col = "lightgray", lty = "dotted")
  
  # Plot 4: ACF of residuals
  acf(result$residuals, main = "ACF of Residuals", 
      col = "darkblue", lwd = 2)
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
  
  return(result)
}
  
