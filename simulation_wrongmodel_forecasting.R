# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

set.seed(8675309)

# ==============================================================================
# exercise settings
# ==============================================================================

# data settings

T <- 300
n <- 2
H = 10
t0 = 10
nrep = 10

# Which DGP? VARMA, VARCH, VARMA Copula
dgp <- "varma_copula"

# sampling settings

dlm_ndraw  = 1000
dlm_burn = 0
dlm_thin = 1

bvar_ndraw = 1000
bvar_burn = 0
bvar_thin = 1

dfc_ndraw  = 1000
dfc_burn = 0
dfc_thin = 1

# ==============================================================================
# Pick DGP and simulate fake data
# ==============================================================================

if(dgp == "varma"){
  
  p = 3
  q = 6
  m = rnorm(n)
  AR = simulate_nonexplosive_var_params(n, p, numeric(n*n*p), 0.1*diag(n*n*p))
  ARarray = array(c(AR), dim = c(n, n, p))
  MAarray = array(rnorm(n*n*q), dim = c(n, n, q))
  S = riwish(n + 1, diag(n))
  
  Y = simulate_stationary_varma(T, m, ARarray, MAarray, S)
  
  abs(eigen(companion(AR))$values)
  
}else if(dgp == "varch"){
  
  v = n + 2 + rexp(1)
  A = riwish(n + 1, diag(n))
  Y = simulate_varch(T, v, A)
  
}else if(dgp == "varma_copula"){
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # get stationary distribution
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Finv1 <- function(x){qpois(x, lambda = 5)}
  Finv2 <- function(x){qst(x, alpha = 2, nu = 3)}
  
  Finv <- c(Finv1, Finv2)
  
  F1 <- function(x){ppois(x, lambda = 5)}
  F2 <- function(x){pst(x, alpha = 2, nu = 3)}
  
  F <- c(F1, F2)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # model settings
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  n = length(Finv)
  p = 1
  q = 1
  
  band_col = rgb(1, 0.6, 0, 0.1)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # generate random parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  AR = simulate_nonexplosive_var_params(n, p, numeric(n*n*p), 0.1*diag(n*n*p))
  ARarray = array(c(AR), dim = c(n, n, p))
  MAarray = array(rnorm(n*n*q), dim = c(n, n, q))
  S = riwish(n + 1, diag(n))
  
  # how persistent is the time series process?
  
  abs(eigen(companion(AR))$values)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulate fake data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Y = simulate_varma_copula(T, Finv, ARarray, MAarray, S)
  
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
  
  mystart <- Sys.time()
  
  for(i in 1:n){
    
    outGibbs2 = dlmGibbsDIG(Y[1:t, i],
                            dlmModPoly(order = 1, m0 = 0, C0 = C0),
                            shape.y = shape.y,
                            rate.y = rate.y,
                            shape.theta = shape.theta,
                            rate.theta = rate.theta,
                            n.sample = dlm_ndraw,
                            thin = dlm_thin,
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
  
  myend <- Sys.time()
  
  print(myend - mystart)
  
  # ---------------------------
  # BVAR
  # ---------------------------
  
  mystart <- Sys.time()
  
  my_bvar_model <- bvar(Y[1:t, ], lags = 1, verbose = FALSE,
                        n_draw = bvar_ndraw * bvar_thin + bvar_burn, n_burn = bvar_burn, n_thin = bvar_thin)
  fcast_bvar[, , , t] <- aperm(predict(my_bvar_model, horizon = H)$fcast, c(2, 3, 1))
  
  myend <- Sys.time()
  
  print(myend - mystart)
  
  # ---------------------------
  # DGFC
  # ---------------------------
  
  mystart <- Sys.time()
  
  draws = DGFC.mcmc(Y[1:t, ], ndraw = dfc_ndraw, burn = dfc_burn, thin = dfc_thin)
  fcast_dfc[, , , t] = DGFC.forecast(H, draws)
  
  myend <- Sys.time()
  
  print(myend - mystart)
  
  message(paste("Stage ", t, " of ", T, " done!", sep = ""))
  
}



# ==============================================================================
# Get results
# ==============================================================================

i = 2
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

par(mfrow = c(4, 1))

dfc_col = "red"
var_col = "blue"
dlm_col = "green"
a = alpha[2]
alpha_level = paste("alpha =", a)
box_lwd = 1.5
common_margin = 0
side_mar = 2

# -------------------------
# Point forecasts
# -------------------------

par(mar = c(common_margin, side_mar, 2, 2))

dfc_mse = dfc_results[, "MSE", alpha_level]
var_mse = var_results[, "MSE", alpha_level]
dlm_mse = dlm_results[, "MSE", alpha_level]

L = min(dfc_mse, var_mse, dlm_mse)
U = max(dfc_mse, var_mse, dlm_mse)

plot(1:H, dfc_mse, type = "l", col = dfc_col, ylim = c(L, U), xaxt = "n",
     xlab = "forecast horizon", ylab = "MSFE", main = paste("Forecasting", dgp, "data"))
lines(1:H, var_mse, col = var_col)
lines(1:H, dlm_mse, col = dlm_col)

points(1:H, dfc_mse, col = dfc_col, pch = 19)
points(1:H, var_mse, col = var_col, pch = 19)
points(1:H, dlm_mse, col = dlm_col, pch = 19)

legend("bottomright", c("BVAR", "DLM", "DGFC"), lty = 1, bty = "n",
       col = c(var_col, dlm_col, dfc_col))

box(lwd = box_lwd)

# -------------------------
# Interval coverage
# -------------------------

par(mar = c(common_margin, side_mar, common_margin, 2))

dfc_int = dfc_results[, "INT-COV", alpha_level]
var_int = var_results[, "INT-COV", alpha_level]
dlm_int = dlm_results[, "INT-COV", alpha_level]

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

par(mar = c(common_margin, side_mar, common_margin, 2))

dfc_int = dfc_results[, "INT-SIZE", alpha_level]
var_int = var_results[, "INT-SIZE", alpha_level]
dlm_int = dlm_results[, "INT-SIZE", alpha_level]

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
# Density forecasts
# -------------------------

par(mar = c(4, side_mar, common_margin, 2))

dfc_crps = dfc_results[, "CRPS", alpha_level]
var_crps = var_results[, "CRPS", alpha_level]
dlm_crps = dlm_results[, "CRPS", alpha_level]

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



# ==============================================================================
# Plot for slides
# ==============================================================================

png(paste("_images/", dgp, "_forecast_slides.png", sep = ""), 
    width = 8, height = 2, units = "in", res = 1000)

par(mfcol = c(1, 4))

dfc_col = "red"
var_col = "blue"
dlm_col = "orange"
a = alpha[2]
alpha_level = paste("alpha =", a)
box_lwd = 1.5
common_margin = 0
side_mar = 2

# -------------------------
# Point forecasts
# -------------------------

par(mar = c(2, 2, 2, 2))

dfc_mse = dfc_results[, "MSE", alpha_level]
var_mse = var_results[, "MSE", alpha_level]
dlm_mse = dlm_results[, "MSE", alpha_level]

L = min(dfc_mse, var_mse, dlm_mse)
U = max(dfc_mse, var_mse, dlm_mse)

plot(1:H, dfc_mse, type = "l", col = dfc_col, ylim = c(L, U),
     xlab = "forecast horizon", ylab = "MSFE", main = "MSE", cex.main = 2)
lines(1:H, var_mse, col = var_col)
lines(1:H, dlm_mse, col = dlm_col)

points(1:H, dfc_mse, col = dfc_col, pch = 19)
points(1:H, var_mse, col = var_col, pch = 19)
points(1:H, dlm_mse, col = dlm_col, pch = 19)

#legend("bottomright", c("BVAR", "DLM", "DGFC"), lty = 1, bty = "n",
#       col = c(var_col, dlm_col, dfc_col))

box(lwd = box_lwd)

# -------------------------
# Interval coverage
# -------------------------

par(mar = c(2, 2, 2, 2))

dfc_int = dfc_results[, "INT-COV", alpha_level]
var_int = var_results[, "INT-COV", alpha_level]
dlm_int = dlm_results[, "INT-COV", alpha_level]

L = min(dfc_int, var_int, dlm_int)
U = max(dfc_int, var_int, dlm_int)

plot(1:H, dfc_int, type = "l", col = dfc_col, ylim = c(L, U),
     xlab = "forecast horizon", main = "Interval coverage", cex.main = 2)
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

par(mar = c(2, 2, 2, 2))

dfc_int = dfc_results[, "INT-SIZE", alpha_level]
var_int = var_results[, "INT-SIZE", alpha_level]
dlm_int = dlm_results[, "INT-SIZE", alpha_level]

L = min(dfc_int, var_int, dlm_int)
U = max(dfc_int, var_int, dlm_int)

plot(1:H, dfc_int, type = "l", col = dfc_col, ylim = c(L, U), 
     xlab = "forecast horizon", main = "Interval size", cex.main = 2)
lines(1:H, var_int, col = var_col)
lines(1:H, dlm_int, col = dlm_col)


points(1:H, dfc_int, col = dfc_col, pch = 19)
points(1:H, var_int, col = var_col, pch = 19)
points(1:H, dlm_int, col = dlm_col, pch = 19)

box(lwd = box_lwd)

#paste("Score of", alpha_level, "HPD interval")

# -------------------------
# Density forecasts
# -------------------------

par(mar = c(2, 2, 2, 2))

dfc_crps = dfc_results[, "CRPS", alpha_level]
var_crps = var_results[, "CRPS", alpha_level]
dlm_crps = dlm_results[, "CRPS", alpha_level]

L = min(dfc_crps, var_crps, dlm_crps)
U = max(dfc_crps, var_crps, dlm_crps)

plot(1:H, dfc_crps, type = "l", col = dfc_col, ylim = c(L, U),
     xlab = "forecast horizon", main = "CRPS", cex.main = 2)
lines(1:H, var_crps, col = var_col)
lines(1:H, dlm_crps, col = dlm_col)

points(1:H, dfc_crps, col = dfc_col, pch = 19)
points(1:H, var_crps, col = var_col, pch = 19)
points(1:H, dlm_crps, col = dlm_col, pch = 19)

box(lwd = box_lwd)

dev.off()