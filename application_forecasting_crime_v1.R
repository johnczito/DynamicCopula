# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

set.seed(8675309)

# ==============================================================================
# load data
# https://www.bocsar.nsw.gov.au/Pages/bocsar_datasets/Offence.aspx
# ==============================================================================

crime_data = read.csv("datasets/full_nsw_crime_dataset.csv", header = TRUE)
crime_data$date <- as.Date(crime_data$date, format = "%m/%d/%Y")

myvars = c("murder", "attempted", "accessory_murder", "manslaughter", "assault_nondomestic",
           "assault_police", "abduction", "blackmail", "Escape_custody", "resist_officer")

var_names = myvars

Y <- as.matrix(crime_data[, myvars])

T <- nrow(Y)
n <- ncol(Y)

# ==============================================================================
# load data
# https://www.bocsar.nsw.gov.au/Pages/bocsar_datasets/Offence.aspx
# ==============================================================================

#maxY <- max(Y)
#var_names <- colnames(crime_data[, 2:ncol(crime_data)])

#empirical_probs <- matrix(0, maxY + 1, n)

#for(i in 1:n){
#  empirical_probs[, i] = table(factor(Y[, i], levels = 0:maxY)) / T
#}

# ==============================================================================
# plot data
# ==============================================================================

#par(mfcol = c(2, n), mar = c(2, 2, 2, 1))

#for(i in 1:n){
#  plot(crime_data$Date, Y[, i], type = "l", 
#       main = paste(var_names[i], "cases in NSW"),
#       ylim = c(0, max(Y)),
#       yaxt = "n")
#  if(i == 1){axis(2)}
#  barplot(empirical_probs[, i], names.arg = 0:maxY, 
#          ylim = c(0, max(empirical_probs)), yaxt = "n")
#  if(i == 1){axis(2)}
#  legend("topright", 
#         c(paste("mean:", round(mean(Y[, i]), 3)), 
#           paste("var:", round(var(Y[, i]), 3))),
#         bty = "n", cex = 1.25)
#}

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

H = 5
t0 = 12 * 5

dlm_ndraw  = 500

pois_ndraw = 2000
nbinom_ndraw = 2000

bvar_ndraw = 1000
bvar_nburn = 1000
bvar_nthin = 5

dfc_ndraw  = 1000
dfc_nburn = 0
dfc_nthin = 1

# ==============================================================================
# Preallocate storage
# ==============================================================================

fcast_dlm  = array(0, c(H, n, dlm_ndraw, T))
fcast_bvar = array(0, c(H, n, bvar_ndraw, T))
fcast_dfc  = array(0, c(H, n, dfc_ndraw, T))
fcast_pois = array(0, c(H, n, pois_ndraw, T))
fcast_nbinom  = array(0, c(H, n, nbinom_ndraw, T))

# ==============================================================================
# which model to run
# ==============================================================================

run_DLM = FALSE
run_BVAR = TRUE
run_DGFC = TRUE
run_POISSON = TRUE
run_NBINOM = TRUE

# ==============================================================================
# Generate forecasts
# ==============================================================================

for(t in t0:T){
  
  # ---------------------------
  # DLM
  # ---------------------------
  
  if(run_DLM){
    
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
    
  }
  
  # ---------------------------
  # Poisson DGLM
  # ---------------------------
  
  mystart <- Sys.time()
  
  fcast_pois[, , , t] <- get_poisson_mvforecasts(Y[1:t, ], H, pois_ndraw)
  
  myend <- Sys.time()
    
  print(myend - mystart)
  
  # ---------------------------
  # Negative binomial DGLM
  # ---------------------------
  
  if(run_NBINOM){
    
    mystart <- Sys.time()
    
    fcast_nbinom[, , , t] <- get_nbinom_mvforecasts(Y[1:t, ], H, nbinom_ndraw)
    
    myend <- Sys.time()
    
    print(myend - mystart)
    
  }
  
  # ---------------------------
  # BVAR
  # ---------------------------
  
  mystart <- Sys.time()
  
  my_bvar_model <- bvar(Y[1:t, ], lags = 1, verbose = FALSE,
                        n_draw = bvar_ndraw * bvar_nthin + bvar_nburn, n_burn = bvar_nburn, n_thin = bvar_nthin)
  fcast_bvar[, , , t] <- aperm(predict(my_bvar_model, horizon = H)$fcast, c(2, 3, 1))
  
  myend <- Sys.time()
  
  print(myend - mystart)
  
  # ---------------------------
  # DGFC
  # ---------------------------
  
  mystart <- Sys.time()
  
  draws = DGFC.mcmc(Y[1:t, ], ndraw = dfc_ndraw, burn = dfc_nburn, thin = dfc_nthin)
  fcast_dfc[, , , t] = DGFC.forecast(H, draws)
  
  myend <- Sys.time()
  
  print(myend - mystart)
  
  message(paste("Stage ", t, " of ", T, " done!", sep = ""))
  
}



# ==============================================================================
# Get results
# ==============================================================================

n = 3

par(mfcol = c(4, n))

dfc_col = "red"
var_col = "blue"
dlm_col = "orange"
nbinom_col = "green"
alpha = c(0.01, 0.05, 0.1, 0.25, 0.5)
a = alpha[3]
alpha_level = paste("alpha =", a)
box_lwd = 1.5
common_margin = 0
left_mar_1 = 4
left_mar_2 = 2
right_mar = 1
Hend = 6

for(i in c(4, 5, 6)){
  dgp = var_names[i]
  
  dfc_stuff = get_variable_fcasts(i, fcast_dfc)
  var_stuff = get_variable_fcasts(i, fcast_bvar)
  dlm_stuff = get_variable_fcasts(i, fcast_pois)
  nbinom_stuff = get_variable_fcasts(i, fcast_nbinom)
  
  dfc_results = forecast_results_discrete(Y[(t0 + 1):T, i], dfc_stuff[(t0 + 1):T, , ], alpha)
  var_results = forecast_results_discrete(Y[(t0 + 1):T, i], var_stuff[(t0 + 1):T, , ], alpha)
  dlm_results = forecast_results_discrete(Y[(t0 + 1):T, i], dlm_stuff[(t0 + 1):T, , ], alpha)
  nbinom_results = forecast_results_discrete(Y[(t0 + 1):T, i], nbinom_stuff[(t0 + 1):T, , ], alpha)
  
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
  
  dfc_mse = dfc_results[1:Hend, "MAE", alpha_level]
  var_mse = var_results[1:Hend, "MAE", alpha_level]
  dlm_mse = dlm_results[1:Hend, "MAE", alpha_level]
  nbinom_mse = nbinom_results[1:Hend, "MAE", alpha_level]
  
  L = min(dfc_mse, var_mse, dlm_mse, nbinom_mse)
  U = max(dfc_mse, var_mse, dlm_mse, nbinom_mse)
  
  plot(1:Hend, dfc_mse, type = "l", col = dfc_col, ylim = c(0.95*L, 1.05*U), xaxt = "n",
       xlab = "forecast horizon", ylab = "MAFE", main = dgp)
  lines(1:Hend, var_mse, col = var_col)
  lines(1:Hend, dlm_mse, col = dlm_col)
  #lines(1:Hend, nbinom_mse, col = nbinom_col)
  
  points(1:Hend, dfc_mse, col = dfc_col, pch = 19)
  points(1:Hend, var_mse, col = var_col, pch = 19)
  points(1:Hend, dlm_mse, col = dlm_col, pch = 19)
  #points(1:Hend, nbinom_mse, col = nbinom_col, pch = 19)
  
  if(i == 2){
    #legend("bottomright", c("BVAR", "DLM", "DGFC"), lty = 1, bty = "n",
    #       col = c(var_col, dlm_col, dfc_col))
  }
  
  box(lwd = box_lwd)
  
  # -------------------------
  # Interval coverage
  # -------------------------
  
  par(mar = c(common_margin, side_mar, common_margin, right_mar))
  
  dfc_int = dfc_results[1:Hend, "INT-COV", alpha_level]
  var_int = var_results[1:Hend, "INT-COV", alpha_level]
  dlm_int = dlm_results[1:Hend, "INT-COV", alpha_level]
  nbinom_int = nbinom_results[1:Hend, "INT-COV", alpha_level]
  
  L = min(dfc_int, var_int, dlm_int, nbinom_int)
  U = max(dfc_int, var_int, dlm_int, nbinom_int)
  
  plot(1:Hend, dfc_int, type = "l", col = dfc_col, ylim = c(0.95*L, 1.05*U), xaxt = "n",
       xlab = "forecast horizon", ylab = "Interval coverage")
  lines(1:Hend, var_int, col = var_col)
  lines(1:Hend, dlm_int, col = dlm_col)
  #lines(1:Hend, nbinom_int, col = nbinom_col)
  
  
  points(1:Hend, dfc_int, col = dfc_col, pch = 19)
  points(1:Hend, var_int, col = var_col, pch = 19)
  points(1:Hend, dlm_int, col = dlm_col, pch = 19)
  #points(1:Hend, nbinom_int, col = nbinom_col, pch = 19)
  
  abline(h = 1 - a, lty = 2, col = "grey")
  
  box(lwd = box_lwd)
  
  # -------------------------
  # Interval size
  # -------------------------
  
  par(mar = c(common_margin, side_mar, common_margin, right_mar))
  
  dfc_int = dfc_results[1:Hend, "INT-SIZE", alpha_level]
  var_int = var_results[1:Hend, "INT-SIZE", alpha_level]
  dlm_int = dlm_results[1:Hend, "INT-SIZE", alpha_level]
  nbinom_int = nbinom_results[1:Hend, "INT-SIZE", alpha_level]
  
  L = min(dfc_int, var_int, dlm_int, nbinom_int)
  U = max(dfc_int, var_int, dlm_int, nbinom_int)
  
  plot(1:Hend, dfc_int, type = "l", col = dfc_col, ylim = c(0.95*L, 1.05*U), xaxt = "n",
       xlab = "forecast horizon", ylab = "Average interval size")
  lines(1:Hend, var_int, col = var_col)
  lines(1:Hend, dlm_int, col = dlm_col)
  #lines(1:Hend, nbinom_int, col = nbinom_col)
  
  
  points(1:Hend, dfc_int, col = dfc_col, pch = 19)
  points(1:Hend, var_int, col = var_col, pch = 19)
  points(1:Hend, dlm_int, col = dlm_col, pch = 19)
  #points(1:Hend, nbinom_int, col = nbinom_col, pch = 19)
  
  box(lwd = box_lwd)
  
  #paste("Score of", alpha_level, "HPD interval")
  
  # -------------------------
  # Interval size
  # -------------------------
  
  #par(mar = c(common_margin, side_mar, common_margin, right_mar))
  
  #dfc_sc = dfc_results[1:H, "INT-SCORE", alpha_level]
  #var_sc = var_results[1:H, "INT-SCORE", alpha_level]
  #dlm_sc = dlm_results[1:H, "INT-SCORE", alpha_level]
  
  #L = min(dfc_sc, var_sc, dlm_sc)
  #U = max(dfc_sc, var_sc, dlm_sc)
  
  #plot(1:H, dfc_sc, type = "l", col = dfc_col, ylim = c(L, U), xaxt = "n",
  #     xlab = "forecast horizon", ylab = "Average interval score")
  #lines(1:H, var_sc, col = var_col)
  #lines(1:H, dlm_sc, col = dlm_col)
  
  
  #points(1:H, dfc_sc, col = dfc_col, pch = 19)
  #points(1:H, var_sc, col = var_col, pch = 19)
  #points(1:H, dlm_sc, col = dlm_col, pch = 19)
  
  #box(lwd = box_lwd)
  
  # -------------------------
  # Density forecasts
  # -------------------------
  
  par(mar = c(4, side_mar, common_margin, right_mar))
  
  dfc_crps = dfc_results[1:Hend, "CRPS", alpha_level]
  var_crps = var_results[1:Hend, "CRPS", alpha_level]
  dlm_crps = dlm_results[1:Hend, "CRPS", alpha_level]
  nbinom_crps = nbinom_results[1:Hend, "CRPS", alpha_level]
  
  L = min(dfc_crps, var_crps, dlm_crps, nbinom_crps)
  U = max(dfc_crps, var_crps, dlm_crps, nbinom_crps)
  
  plot(1:Hend, dfc_crps, type = "l", col = dfc_col, ylim = c(0.95*L, 1.05*U),
       xlab = "forecast horizon", ylab = "RPS")
  lines(1:Hend, var_crps, col = var_col)
  lines(1:Hend, dlm_crps, col = dlm_col)
  #lines(1:Hend, nbinom_crps, col = nbinom_col)
  
  points(1:Hend, dfc_crps, col = dfc_col, pch = 19)
  points(1:Hend, var_crps, col = var_col, pch = 19)
  points(1:Hend, dlm_crps, col = dlm_col, pch = 19)
  #points(1:Hend, nbinom_crps, col = nbinom_col, pch = 19)
  
  box(lwd = box_lwd)
  
}



##################

# for slide





i = 1

dfc_stuff = get_variable_fcasts(i, fcast_dfc)
var_stuff = get_variable_fcasts(i, fcast_bvar)
dlm_stuff = get_variable_fcasts(i, fcast_pois)
nbinom_stuff = get_variable_fcasts(i, fcast_nbinom)

dfc_results = forecast_results_discrete(Y[(t0 + 1):T, i], dfc_stuff[(t0 + 1):T, , ], alpha)
var_results = forecast_results_discrete(Y[(t0 + 1):T, i], var_stuff[(t0 + 1):T, , ], alpha)
dlm_results = forecast_results_discrete(Y[(t0 + 1):T, i], dlm_stuff[(t0 + 1):T, , ], alpha)
nbinom_results = forecast_results_discrete(Y[(t0 + 1):T, i], nbinom_stuff[(t0 + 1):T, , ], alpha)


png(paste("_images/", var_names[i], "_forecast_slides.png", sep = ""), 
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

dfc_mse = dfc_results[, "MAE", alpha_level]
var_mse = var_results[, "MAE", alpha_level]
dlm_mse = dlm_results[, "MAE", alpha_level]

L = min(dfc_mse, var_mse, dlm_mse)
U = max(dfc_mse, var_mse, dlm_mse)

plot(1:H, dfc_mse, type = "l", col = dfc_col, ylim = c(L, U),
     xlab = "forecast horizon", ylab = "MSFE", main = "MAE", cex.main = 2)
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
