forecast_results <- function(data, forecasts, alpha){
  # postprocess forecasting draws
  #
  # data: T-vector of observed data values
  # forecasts: T x H x M array of M predictive samples forecasting period t 
  #            from model or horizon h
  # alpha: A-vector of coverage probabilities for intervals
  
  T <- length(data)
  A <- length(alpha)
  H <- dim(forecasts)[2]
  
  metrics <- array(0, c(H, 6, A, T))
  colnames(metrics) <- c("MAE", "MSE", "INT-SIZE", "INT-COV", "INT-SCORE", "CRPS")
  dimnames(metrics)[[3]] <- paste("alpha =", alpha)
  
  for(h in 1:H){
    for(t in 1:T){
      obs = data[t]
      fcast_draws = forecasts[t, h, ]
      for(a in 1:A){
        
        alph = alpha[a]
        
        fcast_med = median(fcast_draws)
        fcast_mean = mean(fcast_draws)
        fcast_hdi = hdi(fcast_draws, 1 - alph)
        hdi_low = fcast_hdi[1]
        hdi_hi = fcast_hdi[2]
        
        metrics[h, 1, a, t] = abs(fcast_med - obs)
        metrics[h, 2, a, t] = (fcast_mean - obs) ^ 2
        metrics[h, 3, a, t] = hdi_hi - hdi_low
        metrics[h, 4, a, t] = (hdi_low <= obs) & (obs <= hdi_hi)
        metrics[h, 5, a, t] = interval_score(obs, hdi_low, hdi_hi, 100 * (1 - alph))
        metrics[h, 6, a, t] = scoringRules::crps_sample(obs, fcast_draws)
        
        # JZ: am I handling interval bounds correctly for discrete variables?
      }
    }
  }
  
  return( round(apply(metrics, c(1, 2, 3), FUN = mean), digits = 3) )
  
}

forecast_results_discrete <- function(data, forecasts, alpha){
  # postprocess forecasting draws
  #
  # data: T-vector of observed data values
  # forecasts: T x H x M array of M predictive samples forecasting period t 
  #            from model or horizon h
  # alpha: A-vector of coverage probabilities for intervals
  
  T <- length(data)
  A <- length(alpha)
  H <- dim(forecasts)[2]
  
  metrics <- array(0, c(H, 6, A, T))
  colnames(metrics) <- c("MAE", "MSE", "INT-SIZE", "INT-COV", "INT-SCORE", "CRPS")
  dimnames(metrics)[[3]] <- paste("alpha =", alpha)
  
  for(h in 1:H){
    for(t in 1:T){
      obs = data[t]
      fcast_draws = forecasts[t, h, ]
      for(a in 1:A){
        
        alph = alpha[a]
        
        fcast_med = median(fcast_draws)
        fcast_mean = mean(fcast_draws)
        #fcast_hdi = hdi(fcast_draws, 1 - alph)
        hdi_low = quantile(fcast_draws, alph / 2)#fcast_hdi[1]
        hdi_hi = quantile(fcast_draws, 1 - alph/2)#fcast_hdi[2]
        
        metrics[h, 1, a, t] = abs(fcast_med - obs)
        metrics[h, 2, a, t] = (fcast_mean - obs) ^ 2
        metrics[h, 3, a, t] = hdi_hi - hdi_low
        metrics[h, 4, a, t] = (hdi_low <= obs) & (obs <= hdi_hi)
        metrics[h, 5, a, t] = interval_score(obs, hdi_low, hdi_hi, 100 * (1 - alph))
        metrics[h, 6, a, t] = scoringRules::crps_sample(obs, fcast_draws)
        
        # JZ: am I handling interval bounds correctly for discrete variables?
      }
    }
  }
  
  return( round(apply(metrics, c(1, 2, 3), FUN = mean), digits = 3) )
  
}

forecast_results <- function(data, forecasts, alpha){
  # postprocess forecasting draws
  #
  # data: T-vector of observed data values
  # forecasts: T x H x M array of M predictive samples forecasting period t 
  #            from model or horizon h
  # alpha: A-vector of coverage probabilities for intervals
  
  T <- length(data)
  A <- length(alpha)
  H <- dim(forecasts)[2]
  
  metrics <- array(0, c(H, 6, A, T))
  colnames(metrics) <- c("MAE", "MSE", "INT-SIZE", "INT-COV", "INT-SCORE", "CRPS")
  dimnames(metrics)[[3]] <- paste("alpha =", alpha)
  
  for(h in 1:H){
    for(t in 1:T){
      obs = data[t]
      fcast_draws = forecasts[t, h, ]
      for(a in 1:A){
        
        alph = alpha[a]
        
        fcast_med = median(fcast_draws)
        fcast_mean = mean(fcast_draws)
        fcast_hdi = hdi(fcast_draws, 1 - alph)
        hdi_low = fcast_hdi[1]
        hdi_hi = fcast_hdi[2]
        
        metrics[h, 1, a, t] = abs(fcast_med - obs)
        metrics[h, 2, a, t] = (fcast_mean - obs) ^ 2
        metrics[h, 3, a, t] = hdi_hi - hdi_low
        metrics[h, 4, a, t] = (hdi_low <= obs) & (obs <= hdi_hi)
        metrics[h, 5, a, t] = interval_score(obs, hdi_low, hdi_hi, 100 * (1 - alph))
        metrics[h, 6, a, t] = scoringRules::crps_sample(obs, fcast_draws)
        
        # JZ: am I handling interval bounds correctly for discrete variables?
      }
    }
  }
  
  return( round(apply(metrics, c(1, 2, 3), FUN = mean), digits = 3) )
  
}

# Add a method for processing discrete stuff

# 
get_variable_fcasts <- function(i, fcasts){
  # fcasts: H x n x ndraw x T
  # new_fcasts: T, H, ndraw
  H = dim(fcasts)[1]
  ndraw = dim(fcasts)[3]
  T = dim(fcasts)[4]
  new_fcasts = array(0, c(T, H, ndraw))
  for(t in 1:T){
    for(h in 1:H){
      # do an if statement to make sure t - h is not less than 1
      if(t - h >= 1){
        new_fcasts[t, h, ] = fcasts[h, i, , t - h]
      }
    }
  }
  return(new_fcasts)
}