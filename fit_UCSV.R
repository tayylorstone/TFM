# --- Set workspace and source other functions --- #
setwd("/Users/taylorstone/Desktop/TFM")
source("SVRW.r")
source("gamrand.r")

####################
### No Factors ###
####################
fit_ucsv = function(y, nloop = 1000, burnin = 200, plot_results = TRUE) {
  # y: numeric vector
  # nloop, burnin: MCMC settings
  # plot_results: if TRUE, will show tauhat and hhat plots
  
  if (!is.numeric(y)) stop("Input must be a numeric vector")
  
  # Convert to column matrix
  y = matrix(y, ncol = 1)
  T = length(y)
  
  # Priors (as given in Stock and Watson 2007)
  Vtau = 9; Vh = 9
  atau = 10; ltau = .25^2 * (atau - 1)
  ah   = 10; lh   = .2^2  * (ah - 1)
  
  # Initialize Markov chain
  omega2tau = .25^2
  omega2h   = .2^2
  h = as.numeric(log(var(y) * .8)) * matrix(rep(1, T))
  
  # initialize AR component
  c = numeric(T)
  phi = 0.5
  
  # state space set up
  Sp = diag(T)
  d  = matrix(rep(0, T), 1)
  sparse = rbind(d, Sp)
  sparse = sparse[-(T+1), ]
  H = diag(T) - sparse
  
  # storage for MCMC draws
  store_omega2tau = matrix(0, nloop - burnin, 1)
  store_omega2h   = matrix(0, nloop - burnin, 1)
  store_tau       = matrix(0, (nloop - burnin), T)  # Fixed dimensions
  store_h         = matrix(0, (nloop - burnin), T)  # Fixed dimensions
  store_c         = matrix(0, (nloop - burnin), T)  # Store AR component
  store_phi       = matrix(0, (nloop - burnin), 1)  # Store AR coefficient
  
  # Precompute constants
  newatau = (T - 1) / 2 + atau
  newah   = (T - 1) / 2 + ah
  
  # --- Main MCMC loop ---
  for (loop in 1:nloop) {
    ## sample tau
    invOmegatau = diag(T) * c(1/Vtau, 1/omega2tau*rep(1, T-1))
    invSigy = diag(T) * c(exp(-h))
    Ktau = t(H) %*% invOmegatau %*% H + invSigy
    Ctau = t(chol(Ktau))
    tauhat_current = solve(Ktau, (invSigy %*% (y - c)))  # Subtract AR component
    tau = tauhat_current + solve(t(Ctau), matrix(rnorm(T),,1))
    
    ## Update AR component using current tau estimate
    c[1] = 0  # Initialize first observation
    for (t in 2:T) {
      c[t] = phi * c[t-1] + rnorm(1, 0, 0.1)  # AR(1) process for c
    }
    
    ## sample h
    ystar = log((y - tau - c)^2 + .0001)  # Include AR component
    result = SVRW(ystar, h, omega2h, Vh)
    h = result[[1]]
    
    ## sample phi (AR coefficient)
    if (sum(c[1:(T-1)]^2) > 0) {
      phi_var = 1 / (1 + sum(c[1:(T-1)]^2))
      phi_mean = phi_var * sum(c[1:(T-1)] * c[2:T])
      phi = rnorm(1, phi_mean, sqrt(phi_var))
      phi = max(-0.99, min(0.99, phi))  # Keep stationary
    }
    
    ## sample omega2tau
    newltau = ltau + sum((tau[2:nrow(tau)] - tau[1:(nrow(tau)-1)])^2)/2
    omega2tau = 1/gamrand(newatau, newltau)
    
    ## sample omega2h
    newlh = lh + sum((h[2:nrow(h)] - h[1:(nrow(h)-1)])^2)/2
    omega2h = 1/gamrand(newah, newlh)
    
    ## save draws
    if (loop > burnin) {
      i = loop - burnin
      store_tau[i,] = as.vector(tau)
      store_h[i,] = as.vector(h)
      store_c[i,] = as.vector(c)
      store_phi[i,] = phi
      store_omega2tau[i,] = omega2tau
      store_omega2h[i,] = omega2h
    }
  }
  
  # Posterior means
  tauhat = colMeans(store_tau)
  hhat   = colMeans(store_h)
  chat   = colMeans(store_c)
  phihat = mean(store_phi)
  
  # Calculate fitted values correctly
  fitted_series = tauhat + chat
  
  # Optionally plot
  if (plot_results) {
    # Convert to time series objects
    ts_y = ts(as.vector(y), start = c(1991, 1), frequency = 12)
    ts_fitted = ts(fitted_series, start = c(1991, 1), frequency = 12)
    ts_tau = ts(tauhat, start = c(1991, 1), frequency = 12)
    ts_c = ts(chat, start = c(1991, 1), frequency = 12)
    
    par(mfrow = c(1, 2))
    #Fitted vs Actual
    plot(ts_y, main = "Fitted vs Actual Values", 
         ylab = "Value", xlab = "Year", col = "cornflowerblue", lwd = 1.5)
    lines(ts_fitted, col = "coral2", lwd = 2)
    legend("topleft", c("Actual", "Fitted"), col = c("cornflowerblue", "coral2"), 
           lty = 1, lwd = c(1.5, 2), bg = "white")
    grid(col = "lightgray", lty = "dotted")
    #Estimated Trend
    plot(ts_tau, main = "Estimated Trend (tau)", 
         ylab = "Value", xlab = "Year", col = "cornflowerblue", lwd = 1.5)
    grid(col = "lightgray", lty = "dotted")
    
    par(mfrow = c(1, 1))
  }
  
  # Return results
  return(list(
    tauhat = tauhat,
    hhat = hhat,
    chat = chat,
    phihat = phihat,
    fitted_values = fitted_series,
    actual = as.vector(y),
    omega2tau = mean(store_omega2tau),
    omega2h = mean(store_omega2h)
  ))
}

####################
### With Factors ###
####################
fit_ucsv_gif = function(y, global_factor, nloop = 1000, burnin = 200, plot_results = TRUE) {
  # y: numeric vector (e.g. us_infl, uk_infl, etc.)
  # global_factor: numeric vector of global factor values
  # nloop, burnin: MCMC settings
  # plot_results: if TRUE, will show plots
  
  if (!is.numeric(y)) stop("Input must be a numeric vector")
  if (!is.numeric(global_factor)) stop("Global factor must be a numeric vector")
  if (length(global_factor) != length(y)) stop("Global factor must be the same length as the series")
  
  # Convert to column matrix
  y = matrix(y, ncol = 1)
  global_factor = matrix(global_factor, ncol = 1)
  T = length(y)
  
  
  # --- Priors ---
  Vtau = 9; Vh = 9
  atau = 10; ltau = .25^2 * (atau - 1)
  ah   = 10; lh   = .2^2  * (ah - 1)
  
  # Prior for global factor coefficient
  Vbeta = 10  # Prior variance for global factor coefficient
  
  # Initialize Markov chain
  omega2tau = .25^2
  omega2h   = .2^2
  h = as.numeric(log(var(y) * .8)) * matrix(rep(1, T))
  
  
  # Initialize tau (trend component)
  tau = rep(mean(y), T)
  
  # Initialize global factor coefficient
  beta_global = 0.5  # Coefficient for global factor
  
  # --- Initialize AR Component ---
  c = numeric(T)
  phi = 0.5
  
  # --- State-Space Set Up ---
  Sp = diag(T)
  d  = matrix(rep(0, T), 1)
  sparse = rbind(d, Sp)
  sparse = sparse[-(T+1), ]
  H = diag(T) - sparse
  
  # --- Storage for MCMC draws ---
  store_omega2tau = matrix(0, nloop - burnin, 1)
  store_omega2h   = matrix(0, nloop - burnin, 1)
  store_tau       = matrix(0, (nloop - burnin), T)
  store_h         = matrix(0, (nloop - burnin), T)
  store_c         = matrix(0, (nloop - burnin), T)
  store_phi       = matrix(0, (nloop - burnin), 1)
  store_beta      = matrix(0, (nloop - burnin), 1)  # Store global factor coefficient
  
  # Precompute constants
  newatau = (T - 1) / 2 + atau
  newah   = (T - 1) / 2 + ah
  
  # --- Main MCMC loop ---
  for (loop in 1:nloop) {
    ## Sample beta (global factor coefficient)
    # y = tau + beta * global_factor + c + epsilon
    residual_y = y - tau - c
    beta_var = 1 / (1/Vbeta + sum(global_factor^2 * exp(-h)))
    beta_mean = beta_var * sum(global_factor * residual_y * exp(-h))
    beta_global = rnorm(1, beta_mean, sqrt(beta_var))
    
    ## sample tau
    #Adjust for global factor in observation equation
    y_adjusted = y - beta_global * global_factor - c
    invOmegatau = diag(T) * c(1/Vtau, 1/omega2tau*rep(1, T-1))
    invSigy = diag(T) * c(exp(-h))
    Ktau = t(H) %*% invOmegatau %*% H + invSigy
    Ctau = t(chol(Ktau))
    tauhat_current = solve(Ktau, (invSigy %*% y_adjusted))
    tau = tauhat_current + solve(t(Ctau), matrix(rnorm(T),,1))
    
    #Update AR component using current estimates
    residual_for_ar = y - tau - beta_global * global_factor
    c[1] = 0  # Initialize first observation
    for (t in 2:T) {
      c[t] = phi * c[t-1] + rnorm(1, 0, 0.1)  # AR(1) process for c
    }
    
    #sample h
    ystar = log((y - tau - beta_global * global_factor - c)^2 + .0001)
    result = SVRW(ystar, h, omega2h, Vh)
    h = result[[1]]
    
    #sample phi (AR coefficient)
    if (sum(c[1:(T-1)]^2) > 0) {
      phi_var = 1 / (1 + sum(c[1:(T-1)]^2))
      phi_mean = phi_var * sum(c[1:(T-1)] * c[2:T])
      phi = rnorm(1, phi_mean, sqrt(phi_var))
      phi = max(-0.99, min(0.99, phi))  # Keep stationary
    }
    
    #sample omega2tau
    newltau = ltau + sum((tau[2:nrow(tau)] - tau[1:(nrow(tau)-1)])^2)/2
    omega2tau = 1/gamrand(newatau, newltau)
    
    #sample omega2h
    newlh = lh + sum((h[2:nrow(h)] - h[1:(nrow(h)-1)])^2)/2
    omega2h = 1/gamrand(newah, newlh)
    
    #save draws
    if (loop > burnin) {
      i = loop - burnin
      store_tau[i,] = as.vector(tau)
      store_h[i,] = as.vector(h)
      store_c[i,] = as.vector(c)
      store_phi[i,] = phi
      store_beta[i,] = beta_global
      store_omega2tau[i,] = omega2tau
      store_omega2h[i,] = omega2h
    }
  }
  
  # Posterior means
  tauhat = colMeans(store_tau)
  hhat   = colMeans(store_h)
  chat   = colMeans(store_c)
  phihat = mean(store_phi)
  betahat = mean(store_beta)
  
  # Calculate fitted values correctly
  fitted_series = tauhat + betahat * as.vector(global_factor) + chat
  
  # Calculate residuals
  residuals = as.vector(y) - fitted_series
  
  if (plot_results) {
    # Convert to time series objects
    ts_y = ts(as.vector(y), start = c(1991, 1), frequency = 12)
    ts_fitted = ts(fitted_series, start = c(1991, 1), frequency = 12)
    ts_tau = ts(tauhat, start = c(1991, 1), frequency = 12)
    ts_c = ts(chat, start = c(1991, 1), frequency = 12)
    
    par(mfrow = c(1, 2))
    #Fitted vs Actual
    plot(ts_y, main = "Fitted vs Actual Values", 
         ylab = "Value", xlab = "Year", col = "cornflowerblue", lwd = 1.5)
    lines(ts_fitted, col = "coral2", lwd = 2)
    legend("topleft", c("Actual", "Fitted"), col = c("cornflowerblue", "coral2"), 
           lty = 1, lwd = c(1.5, 2), bg = "white")
    grid(col = "lightgray", lty = "dotted")
    #Estimated Trend
    plot(ts_tau, main = "Estimated Trend (tau)", 
         ylab = "Value", xlab = "Year", col = "cornflowerblue", lwd = 1.5)
    grid(col = "lightgray", lty = "dotted")
    
   
    
    par(mfrow = c(1, 1))
  }
  
  # Return results
  result = list(
    # Main model components
    tauhat = tauhat,
    hhat = hhat,
    chat = chat,
    betahat = betahat,
    phihat = phihat,
    
    # Model fit
    fitted_values = fitted_series,
    residuals = residuals,
    actual = as.vector(y),
    
    # Model parameters
    omega2tau = mean(store_omega2tau),
    omega2h = mean(store_omega2h),
    
    # Data info (like ARIMA function)
    time_series = ts(as.vector(y), start = c(1991, 1), frequency = 12),
    start_date = as.Date("1991-01-01"),
    end_date = as.Date("2023-12-01"),
    frequency = 12,
    data_length = T,
    global_factor = as.vector(global_factor),
    
    # MCMC storage for further analysis
    mcmc_draws = list(
      tau = store_tau,
      h = store_h,
      c = store_c,
      phi = store_phi,
      beta = store_beta,
      omega2tau = store_omega2tau,
      omega2h = store_omega2h
    )
  )
  
  return(result)
}
