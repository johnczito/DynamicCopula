# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

set.seed(8675309)

# ==============================================================================
# Pick DGP and simulate fake data
# ==============================================================================



Y = as.matrix(read.csv("datasets/ccm_2016_JBES_4var_1964Q1_2013Q4.csv", header = FALSE))

T <- nrow(Y)
n <- ncol(Y)

# ==============================================================================
# Plot
# ==============================================================================

par(mfcol = c(2, n - 1), mar = c(2, 2, 2, 2))
labels = c("GDP growth", "Unemployment", "Inflation", "Federal funds rate")
dates = seq(1964, 2013.75, by = 0.25)

for( i in 1:3){
  plot(dates, Y[, i], type = "l", main = labels[i])
  hist(Y[, i], breaks = "Scott", main = "")
}

# ==============================================================================
# DLM hyperparameters
# ==============================================================================

m0 = 0
C0 = 1
shape.y = 1
rate.y = 1
shape.theta = 1
rate.theta = 1

# ==============================================================================
# Sampling settings
# ==============================================================================

H = 10
t0 = 50

dlm_ndraw  = 1000
bvar_ndraw = 1000
dfc_ndraw  = 1000
nburn = 0
nthin = 1

# ==============================================================================
# Preallocate storage
# ==============================================================================

fcast_dlm  = array(0, c(H, n, dlm_ndraw, T))
fcast_bvar = array(0, c(H, n, bvar_ndraw, T))
fcast_dfc  = array(0, c(H, n, dfc_ndraw, T))

# ==============================================================================
# Generate forecasts
# ==============================================================================

for(t in t0:T){
  
  # ---------------------------
  # DLM
  # ---------------------------
  
  for(i in 1:n){
    
    outGibbs2 = dlmGibbsDIG(Y[1:t, i],
                            dlmModPoly(order = 1, m0 = 0, C0 = C0),
                            shape.y = shape.y,
                            rate.y = rate.y,
                            shape.theta = shape.theta,
                            rate.theta = rate.theta,
                            n.sample = dlm_ndraw,
                            thin = nthin,
                            ind = 1,
                            save.states = TRUE,
                            progressBar = FALSE)
    
    
    for(l in 1:dlm_ndraw){
      new_s = numeric(H)
      s0 = outGibbs2$theta[t + 1, 1, l]
      for(h in 1:H){
        new_s[h] = rnorm(1, mean = s0, sd = sqrt(outGibbs2$dW[l, 1]))
        s0 = new_s[h]
      }
      fcast_dlm[, i, l, t] = rnorm(H, mean = new_s, sd = sqrt(outGibbs2$dV[l]))
    }
    
  }
  
  # ---------------------------
  # BVAR
  # ---------------------------
  
  my_bvar_model <- bvar(Y[1:t, ], lags = 1, verbose = FALSE,
                        n_draw = bvar_ndraw * nthin + nburn, n_burn = nburn, n_thin = nthin)
  fcast_bvar[, , , t] <- aperm(predict(my_bvar_model, horizon = H)$fcast, c(2, 3, 1))
  
  # ---------------------------
  # DGFC
  # ---------------------------
  
  draws = DGFC.mcmc(Y[1:t, ], ndraw = dfc_ndraw, burn = nburn, thin = nthin)
  fcast_dfc[, , , t] = DGFC.forecast(H, draws)
  
  message(paste("Stage ", t, " of ", T, " done!", sep = ""))
  
}



# ==============================================================================
# Get results
# ==============================================================================

par(mfcol = c(5, n))

dfc_col = "red"
var_col = "blue"
dlm_col = "green"
a = alpha[3]
alpha_level = paste("alpha =", a)
box_lwd = 1.5
common_margin = 0
left_mar_1 = 4
left_mar_2 = 2
right_mar = 1
H = 5

for(i in 1:n){
  dgp = c("gdp", "unemployment", "inflation", "interest rate")[i]
  alpha = c(0.01, 0.05, 0.1, 0.25, 0.5)
  
  dfc_stuff = get_variable_fcasts(i, fcast_dfc)
  var_stuff = get_variable_fcasts(i, fcast_bvar)
  dlm_stuff = get_variable_fcasts(i, fcast_dlm)
  
  dfc_results = forecast_results(Y[(t0 + 1):T, i], dfc_stuff[(t0 + 1):T, , ], alpha)
  var_results = forecast_results(Y[(t0 + 1):T, i], var_stuff[(t0 + 1):T, , ], alpha)
  dlm_results = forecast_results(Y[(t0 + 1):T, i], dlm_stuff[(t0 + 1):T, , ], alpha)
  
  # ==============================================================================
  # Plot things
  # ==============================================================================
  
  if(i == 1){
    side_mar = left_mar_1
  }else{
    side_mar = left_mar_2
  }
  
  # -------------------------
  # Point forecasts
  # -------------------------
  
  par(mar = c(common_margin, side_mar, 2, right_mar))
  
  dfc_mse = dfc_results[1:H, "MSE", alpha_level]
  var_mse = var_results[1:H, "MSE", alpha_level]
  dlm_mse = dlm_results[1:H, "MSE", alpha_level]
  
  L = min(dfc_mse, var_mse, dlm_mse)
  U = max(dfc_mse, var_mse, dlm_mse)
  
  plot(1:H, dfc_mse, type = "l", col = dfc_col, ylim = c(L, U), xaxt = "n",
       xlab = "forecast horizon", ylab = "MSFE", main = dgp)
  lines(1:H, var_mse, col = var_col)
  lines(1:H, dlm_mse, col = dlm_col)
  
  points(1:H, dfc_mse, col = dfc_col, pch = 19)
  points(1:H, var_mse, col = var_col, pch = 19)
  points(1:H, dlm_mse, col = dlm_col, pch = 19)
  
  if(i == 2){
    legend("bottomright", c("BVAR", "DLM", "DGFC"), lty = 1, bty = "n",
           col = c(var_col, dlm_col, dfc_col))
  }
  
  box(lwd = box_lwd)
  
  # -------------------------
  # Interval coverage
  # -------------------------
  
  par(mar = c(common_margin, side_mar, common_margin, right_mar))
  
  dfc_int = dfc_results[1:H, "INT-COV", alpha_level]
  var_int = var_results[1:H, "INT-COV", alpha_level]
  dlm_int = dlm_results[1:H, "INT-COV", alpha_level]
  
  L = min(dfc_int, var_int, dlm_int)
  U = max(dfc_int, var_int, dlm_int)
  
  plot(1:H, dfc_int, type = "l", col = dfc_col, ylim = c(L, U), xaxt = "n",
       xlab = "forecast horizon", ylab = "Interval coverage")
  lines(1:H, var_int, col = var_col)
  lines(1:H, dlm_int, col = dlm_col)
  
  
  points(1:H, dfc_int, col = dfc_col, pch = 19)
  points(1:H, var_int, col = var_col, pch = 19)
  points(1:H, dlm_int, col = dlm_col, pch = 19)
  
  abline(h = 1 - a, lty = 2, col = "grey")
  
  box(lwd = box_lwd)
  
  # -------------------------
  # Interval size
  # -------------------------
  
  par(mar = c(common_margin, side_mar, common_margin, right_mar))
  
  dfc_int = dfc_results[1:H, "INT-SIZE", alpha_level]
  var_int = var_results[1:H, "INT-SIZE", alpha_level]
  dlm_int = dlm_results[1:H, "INT-SIZE", alpha_level]
  
  L = min(dfc_int, var_int, dlm_int)
  U = max(dfc_int, var_int, dlm_int)
  
  plot(1:H, dfc_int, type = "l", col = dfc_col, ylim = c(L, U), xaxt = "n",
       xlab = "forecast horizon", ylab = "Average interval size")
  lines(1:H, var_int, col = var_col)
  lines(1:H, dlm_int, col = dlm_col)
  
  
  points(1:H, dfc_int, col = dfc_col, pch = 19)
  points(1:H, var_int, col = var_col, pch = 19)
  points(1:H, dlm_int, col = dlm_col, pch = 19)
  
  box(lwd = box_lwd)
  
  #paste("Score of", alpha_level, "HPD interval")
  
  # -------------------------
  # Interval size
  # -------------------------
  
  par(mar = c(common_margin, side_mar, common_margin, right_mar))
  
  dfc_sc = dfc_results[1:H, "INT-SCORE", alpha_level]
  var_sc = var_results[1:H, "INT-SCORE", alpha_level]
  dlm_sc = dlm_results[1:H, "INT-SCORE", alpha_level]
  
  L = min(dfc_sc, var_sc, dlm_sc)
  U = max(dfc_sc, var_sc, dlm_sc)
  
  plot(1:H, dfc_sc, type = "l", col = dfc_col, ylim = c(L, U), xaxt = "n",
       xlab = "forecast horizon", ylab = "Average interval score")
  lines(1:H, var_sc, col = var_col)
  lines(1:H, dlm_sc, col = dlm_col)
  
  
  points(1:H, dfc_sc, col = dfc_col, pch = 19)
  points(1:H, var_sc, col = var_col, pch = 19)
  points(1:H, dlm_sc, col = dlm_col, pch = 19)
  
  box(lwd = box_lwd)
  
  # -------------------------
  # Density forecasts
  # -------------------------
  
  par(mar = c(4, side_mar, common_margin, right_mar))
  
  dfc_crps = dfc_results[1:H, "CRPS", alpha_level]
  var_crps = var_results[1:H, "CRPS", alpha_level]
  dlm_crps = dlm_results[1:H, "CRPS", alpha_level]
  
  L = min(dfc_crps, var_crps, dlm_crps)
  U = max(dfc_crps, var_crps, dlm_crps)
  
  plot(1:H, dfc_crps, type = "l", col = dfc_col, ylim = c(L, U),
       xlab = "forecast horizon", ylab = "CRPS")
  lines(1:H, var_crps, col = var_col)
  lines(1:H, dlm_crps, col = dlm_col)
  
  points(1:H, dfc_crps, col = dfc_col, pch = 19)
  points(1:H, var_crps, col = var_col, pch = 19)
  points(1:H, dlm_crps, col = dlm_col, pch = 19)
  
  box(lwd = box_lwd)
  
}

