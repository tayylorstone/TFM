# --- Set workspace and source other functions --- #
setwd("/Users/taylorstone/Desktop/TFM")
source("SVRW.r")
source("gamrand.r")

fit_ucsv <- function(y, nloop = 1000, burnin = 200, plot_results = TRUE) {
  # y: numeric vector (e.g. us_infl, uk_infl, etc.)
  # nloop, burnin: MCMC settings
  # plot_results: if TRUE, will show tauhat and hhat plots
  
  if (!is.numeric(y)) stop("Input must be a numeric vector")
  
  # Convert to column matrix
  y <- matrix(y, ncol = 1)
  T <- length(y)
  
  # --- Priors ---
  Vtau <- 9; Vh <- 9
  atau <- 10; ltau <- .25^2 * (atau - 1)
  ah   <- 10; lh   <- .2^2  * (ah - 1)
  
  # Initialize Markov chain
  omega2tau <- .25^2
  omega2h   <- .2^2
  h <- as.numeric(log(var(y) * .8)) * matrix(rep(1, T))
  
  # --- State-Space Set Up ---
  Sp <- diag(T)
  d  <- matrix(rep(0, T), 1)
  sparse <- rbind(d, Sp)
  sparse <- sparse[-(T+1), ]
  H <- diag(T) - sparse
  
  # --- Storage for MCMC draws ---
  store_omega2tau <- matrix(0, nloop - burnin, 1)
  store_omega2h   <- matrix(0, nloop - burnin, 1)
  store_tau       <- matrix(0, (nloop - burnin) * T, T)
  store_h         <- matrix(0, (nloop - burnin) * T, T)
  
  # Precompute constants
  newatau <- (T - 1) / 2 + atau
  newah   <- (T - 1) / 2 + ah
  
  # --- Main MCMC loop ---
  for (loop in 1:nloop) {
    ## sample tau
    invOmegatau <- diag(T) * c(1/Vtau, 1/omega2tau*rep(1, T-1))
    invSigy <- diag(T) * c(exp(-h))
    Ktau <- t(H) %*% invOmegatau %*% H + invSigy
    Ctau <- t(chol(Ktau))
    tauhat <- solve(Ktau, (invSigy %*% y))
    tau <- tauhat + solve(t(Ctau), matrix(rnorm(T),,1))
    
    ## sample h
    ystar <- log((y - tau)^2 + .0001)
    result <- SVRW(ystar, h, omega2h, Vh)
    h <- result[[1]]
    
    ## sample omega2tau
    newltau <- ltau + sum((tau[2:nrow(tau)] - tau[1:(nrow(tau)-1)])^2)/2
    omega2tau <- 1/gamrand(newatau, newltau)
    
    ## sample omega2h
    newlh <- lh + sum((h[2:nrow(h)] - h[1:(nrow(h)-1)])^2)/2
    omega2h <- 1/gamrand(newah, newlh)
    
    ## save draws
    if (loop > burnin) {
      i <- loop - burnin
      store_tau[i,] <- t(tau)
      store_h[i,] <- t(h)
      store_omega2tau[i,] <- omega2tau
      store_omega2h[i,] <- omega2h
    }
  }
  
  # Posterior means
  tauhat <- matrix(rowMeans(t(store_tau)))
  hhat   <- matrix(colMeans(store_h))
  
  # Optionally plot
  if (plot_results) {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    plot(tauhat, type = 'l', col = 'blue', main="Estimated Trend (tau)")
    
    plot(hhat, type = 'l', col = 'blue', main="Estimated Volatility (h)")
    
    plot(tauhat, main = "Trend vs Actual Values", 
         ylab = "Value", xlab = "Index", col = "blue", lwd = 1.5) 
    lines(y, col = "red", lwd = 2) 
    legend("topleft", c("Actual", "Trend"), col = c("blue", "red"), 
           lty = 1, lwd = c(1.5, 2), bg = "white") 
  }
  par(mfrow=c(2,2))
  
}

