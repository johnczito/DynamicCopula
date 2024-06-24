# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

set.seed(8675309)

# ==============================================================================
# load full data
# https://www.bocsar.nsw.gov.au/Pages/bocsar_datasets/Offence.aspx
# ==============================================================================

crime_data = read.csv("datasets/full_nsw_crime_dataset.csv", header = TRUE)

# ==============================================================================
# subset of data for forecasting
# ==============================================================================

dates <- as.Date(crime_data$date, format = "%m/%d/%Y")

var_names <- c("murder", "attempted", "accessory_murder", "manslaughter", 
               "assault_nondomestic", "assault_police", "abduction", 
               "blackmail", "Escape_custody", "resist_officer")

my_labs <- c("Murder", "Attempted", "Accessory", "Manslaughter", "Non-domestic assault",
             "Officer assault", "Abduction", "Blackmail", "Escape custody", "Resist officer")

Y <- as.matrix(crime_data[, var_names])

T <- nrow(Y)
n <- ncol(Y)

# JZ: cross-section of variables that are vaguely normal, sorta continuous,
#     low counts, zero-inflated, asymmetric, heavy-tail
  
# ==============================================================================
# plot time series and histograms for each chosen variable
# ==============================================================================

png(paste("_images/crime_data_", n, "_var.png", sep = ""), 
    width = 6.5, height = 8, units = "in", res = 650)

par(mfrow = c(5, 4))

for(i in c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10)){
  y = Y[, i]
  minY = min(y)
  maxY = max(y)
  
  par(mar = c(2, 2, 2, 0.5))
  plot(dates, y, type = "l", main = my_labs[i])
  mtext(var_names[i], outer = TRUE, line = 1)
  
  par(mar = c(2, 0.5, 2, 1))
  probs = table(factor(y, levels = minY:maxY)) / T
  barplot(probs, names.arg = minY:maxY, col = "blue", border = NA, yaxt = "n")
}

dev.off()
  
# JZ: need to have the variable name bestride both plots

# ==============================================================================
# sampling settings
# ==============================================================================

pois_ndraw = 1000
nbinom_ndraw = 1000

bvar_ndraw = 1000
bvar_nburn = 1000
bvar_nthin = 1

dfc_ndraw  = 1000
dfc_nburn = 0
dfc_nthin = 1

# ==============================================================================
# forecast settings
# ==============================================================================

H = 5
Tstart = 12 * 5
Tstop = T

# ==============================================================================
# model settings
# ==============================================================================

bvar_lag = 1
factor_dim = 3

col_bvar = "orange"
col_dfc = "blue"
col_pois = "red"
col_nbinom = "darkgreen"

# ==============================================================================
# preallocate storage
# ==============================================================================

fcast_bvar = array(0, c(H, n, bvar_ndraw, T))
fcast_dfc  = array(0, c(H, n, dfc_ndraw, T))
fcast_pois = array(0, c(H, n, pois_ndraw, T))
fcast_nbinom  = array(0, c(H, n, nbinom_ndraw, T))

# ==============================================================================
# which model to run
# ==============================================================================

run_BVAR = TRUE
run_DGFC = TRUE
run_POISSON = TRUE
run_NBINOM = TRUE

# ==============================================================================
# run it hot!
# ==============================================================================

for(t in Tstart:Tstop){
  # ---------------------------
  # Poisson DGLM
  # ---------------------------
  
  if(run_POISSON){
    mystart <- Sys.time()
    
    fcast_pois[, , , t] <- get_poisson_mvforecasts(Y[1:t, ], H, pois_ndraw)
    
    myend <- Sys.time()
    
    print(myend - mystart)
  }
  
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
  
  if(run_BVAR){
    mystart <- Sys.time()
    
    my_bvar_model <- bvar(Y[1:t, ], lags = bvar_lag, verbose = FALSE,
                          n_draw = bvar_ndraw * bvar_nthin + bvar_nburn, n_burn = bvar_nburn, n_thin = bvar_nthin)
    fcast_bvar[, , , t] <- aperm(predict(my_bvar_model, horizon = H)$fcast, c(2, 3, 1))
    
    myend <- Sys.time()
    
    print(myend - mystart)
  }
  
  # ---------------------------
  # BVAR
  # ---------------------------
  
  if(run_DGFC){
    mystart <- Sys.time()
    
    draws = DGFC.mcmc(Y[1:t, ], k.star = factor_dim, 
                      ndraw = dfc_ndraw, burn = dfc_nburn, thin = dfc_nthin)
    fcast_dfc[, , , t] = DGFC.forecast(H, draws)
    
    myend <- Sys.time()
    
    print(myend - mystart)
  }
  
  message(paste("Stage ", t, " of ", Tstop, " done!", sep = ""))
}

# ==============================================================================
# h-step-ahead forecast table
# ==============================================================================

alpha = 0.05
period = (Tstart + H):Tstop
h = 1
metrics = c("MAE", "INT-COV", "INT-SIZE", "CRPS")
nmod = 4
nmet = length(metrics)
mastertable = matrix(0, 0, nmet)
colnames(mastertable) = metrics

for(i in 1:n){
  # T x H x ndraw
  dfc_stuff = get_variable_fcasts(i, fcast_dfc)
  var_stuff = get_variable_fcasts(i, fcast_bvar)
  pois_stuff = get_variable_fcasts(i, fcast_pois)
  nbinom_stuff = get_variable_fcasts(i, fcast_nbinom)
  
  # H x 6 x A x T
  dfc_metrics = get_forecast_metrics(Y[period, i], dfc_stuff[period, , ], alpha)
  var_metrics = get_forecast_metrics(Y[period, i], var_stuff[period, , ], alpha)
  pois_metrics = get_forecast_metrics(Y[period, i], pois_stuff[period, , ], alpha)
  nbinom_metrics = get_forecast_metrics(Y[period, i], nbinom_stuff[period, , ], alpha)
  
  # H x 6 x A
  dfc_results = apply(dfc_metrics, c(1, 2, 3), FUN = mean)
  var_results = apply(var_metrics, c(1, 2, 3), FUN = mean)
  pois_results = apply(pois_metrics, c(1, 2, 3), FUN = mean)
  nbinom_results = apply(nbinom_metrics, c(1, 2, 3), FUN = mean)
  
  summ <- round(rbind(dfc_results[h, metrics, 1], var_results[h, metrics, 1], 
                      pois_results[h, metrics, 1], nbinom_results[h, metrics, 1]), digits = 2)
  rownames(summ) <- c("DGFC", "BVAR", "Pois-DGLM", "NB-DGLM")
  
  mastertable <- rbind(mastertable, summ)
}

#mastertable <- cbind(mastertable[1:(nmod * n / 2), ], mastertable[((nmod * n / 2) + 1):(nmod*n), ])

#write.csv(mastertable, "crime_fcasts.csv", row.names = FALSE)

my_metric = "MAE"

plot(1:length(period), dfc_metrics[h, my_metric, 1, ], type = "l")
lines(1:length(period), var_metrics[h, my_metric, 1, ])
lines(1:length(period), pois_metrics[h, my_metric, 1, ])


# ==============================================================================
# inspect forecasts
# ==============================================================================

i = 1
alpha = c(0.01, 0.05, 0.1, 0.25, 0.5)
period = (Tstart + H):Tstop
h = 1
metric = "INT-SIZE"

# T x H x ndraw
dfc_stuff = get_variable_fcasts(i, fcast_dfc)
var_stuff = get_variable_fcasts(i, fcast_bvar)
pois_stuff = get_variable_fcasts(i, fcast_pois)
nbinom_stuff = get_variable_fcasts(i, fcast_nbinom)

# H x 6 x A x T
dfc_metrics = get_forecast_metrics(Y[period, i], dfc_stuff[period, , ], alpha)
var_metrics = get_forecast_metrics(Y[period, i], var_stuff[period, , ], alpha)
pois_metrics = get_forecast_metrics(Y[period, i], pois_stuff[period, , ], alpha)
nbinom_metrics = get_forecast_metrics(Y[period, i], nbinom_stuff[period, , ], alpha)

# H x 6 x A
dfc_results = apply(dfc_metrics, c(1, 2, 3), FUN = mean)
var_results = apply(var_metrics, c(1, 2, 3), FUN = mean)
pois_results = apply(pois_metrics, c(1, 2, 3), FUN = mean)
nbinom_results = apply(nbinom_metrics, c(1, 2, 3), FUN = mean)

par(mar = c(2, 2, 2, 2))

dfc_mse = dfc_results[1:H, metric, 1]
var_mse = var_results[1:H, metric, 1]
pois_mse = pois_results[1:H, metric, 1]
nbinom_mse = nbinom_results[1:H, metric, 1]

L = min(dfc_mse, var_mse, pois_mse, nbinom_mse)
U = max(dfc_mse, var_mse, pois_mse, nbinom_mse)

plot(1:H, dfc_mse, type = "l", col = col_dfc, ylim = c(0.95*L, 1.05*U), xaxt = "n",
     xlab = "forecast horizon", ylab = metric)
lines(1:H, var_mse, col = col_bvar)
lines(1:H, pois_mse, col = col_pois)
lines(1:H, nbinom_mse, col = col_nbinom)











h = 1
t = Tstart
alpha = 0.1
yo = Y[t + h, i]

bvar_draws = round(fcast_bvar[h, i, , t])
dfc_draws = fcast_dfc[h, i, , t]
pois_draws = fcast_pois[h, i, , t]
nbinom_draws = fcast_nbinom[h, i, , t]

dfc_med = median(dfc_draws)
dfc_int = c(quantile(dfc_draws, alpha / 2), quantile(dfc_draws, 1 - alpha/2))
dfc_crps = scoringRules::crps_sample(yo, dfc_draws)

pois_med = median(pois_draws)
pois_int = c(quantile(pois_draws, alpha / 2), quantile(pois_draws, 1 - alpha/2))
pois_crps = scoringRules::crps_sample(yo, pois_draws)

nbinom_med = median(nbinom_draws)
nbinom_int = c(quantile(nbinom_draws, alpha / 2), quantile(nbinom_draws, 1 - alpha/2))
nbinom_crps = scoringRules::crps_sample(yo, nbinom_draws)


L = min(c(dfc_draws, pois_draws, nbinom_draws))
U = max(c(dfc_draws, pois_draws, nbinom_draws))

par(mfrow = c(1, 1))
barplot(table(factor(dfc_draws, levels = L:U)) / dfc_ndraw, names.arg = L:U, border = NA, col = rgb(0, 0, 1, 0.5))
barplot(table(factor(pois_draws, levels = L:U)) / pois_ndraw, names.arg = L:U, border = NA, add = TRUE, col = rgb(0, 1, 0, 0.5))
barplot(table(factor(nbinom_draws, levels = L:U)) / nbinom_ndraw, names.arg = L:U, border = NA, add = TRUE, col = rgb(1, 0, 0, 0.5))
#barplot(table(factor(bvar_draws, levels = L:U)) / bvar_ndraw, names.arg = L:U, border = NA, add = TRUE, col = rgb(0, 0, 0, 0.5))

