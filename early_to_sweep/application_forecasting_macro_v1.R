# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

set.seed(8675309)

continuous_forecast_metrics <- function(obs, alph, fcast_draws){
  fcast_mean <- mean(fcast_draws)
  fcast_hdi <- hdi(fcast_draws, 1 - alph)
  hdi_low = fcast_hdi[1]
  hdi_hi = fcast_hdi[2]
  err <- (fcast_mean - obs)^2
  size <- hdi_hi - hdi_low
  cover <- (hdi_low <= obs) & (obs <= hdi_hi)
  crpscore <- scoringRules::crps_sample(obs, fcast_draws)
  return(c(err, cover, size, crpscore))
}

# ==============================================================================
# get real-time data set up
# ==============================================================================

var_names <- c("gdpgrowth", "pcegrowth", "bfigrowth", "resinvgrowth", "ipgrowth", "capu", 
               "employ", "hours", "ur", "gdppigrowth", "pcepigrowth", "ffr", "spread", "realSPgrowth")

unrevised <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_ur_ffr_spread_1964Q1_2014Q1.csv", header = FALSE)
bfigrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_bfigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
capu <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_capu_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
employ <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_employ_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
gdpgrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_gdpgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
gdppigrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_gdppigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
hours <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_hours_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
ipgrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_ipgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
pcegrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_pcegrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
pcepigrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_pcepigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
realSPgrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_realSPgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
resinvgrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_resinvgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)

dates = seq(1964.00, 2014.00, by = 0.25)
vintages <- seq(1985.00, 2014.25, by = 0.25)

n <- length(var_names)
Tmax <- nrow(unrevised)
V <- length(vintages)

real_time_macro_data <- array(0, c(Tmax, n, V))

for(v in 1:V){
  vintage <- vintages[v]
  end_period <- vintage - 0.25
  T <- which(dates == end_period)
  
  real_time_macro_data[1:T, 01, v] <- gdpgrowth[1:T, v]
  real_time_macro_data[1:T, 02, v] <- pcegrowth[1:T, v]
  real_time_macro_data[1:T, 03, v] <- bfigrowth[1:T, v]
  real_time_macro_data[1:T, 04, v] <- resinvgrowth[1:T, v]
  real_time_macro_data[1:T, 05, v] <- ipgrowth[1:T, v]
  real_time_macro_data[1:T, 06, v] <- capu[1:T, v]
  real_time_macro_data[1:T, 07, v] <- employ[1:T, v]
  real_time_macro_data[1:T, 08, v] <- hours[1:T, v]
  real_time_macro_data[1:T, 09, v] <- unrevised[1:T, 1]
  real_time_macro_data[1:T, 10, v] <- gdppigrowth[1:T, v]
  real_time_macro_data[1:T, 11, v] <- pcepigrowth[1:T, v]
  real_time_macro_data[1:T, 12, v] <- unrevised[1:T, 2]
  real_time_macro_data[1:T, 13, v] <- unrevised[1:T, 3]
  real_time_macro_data[1:T, 14, v] <- realSPgrowth[1:T, v]
}

# ==============================================================================
# sampling settings
# ==============================================================================

bvar_ndraw = 1000
bvar_nburn = 1000
bvar_nthin = 1

dfc_ndraw = 2000
dfc_nburn = 0
dfc_nthin = 1

dlm_ndraw = 2000

# ==============================================================================
# preallocate storage
# ==============================================================================

Fend = V - 2

fcast_bvar = array(0, c(n, bvar_ndraw, Fend))
fcast_dfc  = array(0, c(n, dfc_ndraw, Fend))
fcast_dlm =  array(0, c(n, dlm_ndraw, Fend))
actuals = matrix(0, Fend, n)

# ==============================================================================
# run it hot!
# ==============================================================================

for(v in 1:Fend){
  
  vintage <- vintages[v]
  end_period <- vintage - 0.25
  T <- which(dates == end_period)
  Y <- real_time_macro_data[1:T, , v]
  actuals[v, ] = real_time_macro_data[T + 1, , v + 2]
  
  # ---------------------------
  # DLM
  # ---------------------------
  
  mystart <- Sys.time()
  
  inits <- c(numeric(sum(1:n)), log(diag(cov(Y))), numeric(sum(1:(n-1))))
  
  my_dlm <- SSModel(Y ~ SSMtrend(1, Q = matrix(NA, n, n), type = "distinct"), H = matrix(NA, n, n))
  fitted_model = fitSSM(my_dlm, inits = inits)$model
  filter_output <- KFS(fitted_model)
  
  fcast_mean = filter_output$att[T, ]
  fcast_cov = filter_output$Ptt[, , T] + fitted_model$H[, , 1] + fitted_model$Q[, , 1]
  
  fcast_dlm[, , v] = t(mvrnorm(dlm_ndraw, fcast_mean, fcast_cov))
  
  myend <- Sys.time()
  
  print(myend - mystart)
  
  # ---------------------------
  # BVAR
  # ---------------------------
  
  mystart <- Sys.time()
  
  mybvar = bvar(Y, 
                lags = 4, 
                verbose = FALSE, 
                n_draw = bvar_ndraw * bvar_nthin + bvar_nburn, 
                n_burn = bvar_nburn, 
                n_thin = bvar_nthin)
  
  fcast_bvar[, , v] = t(predict(mybvar, horizon = 1)$fcast[, 1, ])
  
  myend <- Sys.time()
  
  print(myend - mystart)
  
  # ---------------------------
  # DGFC
  # ---------------------------
  
  mystart <- Sys.time()
  
  draws = DGFC.mcmc(Y, 
                    k.star = 3, 
                    ndraw = dfc_ndraw, 
                    burn = dfc_nburn, 
                    thin = dfc_nthin)
  fcast_dfc[, , v] = DGFC.forecast(1, draws, use_spline = TRUE)[1, , ]
  
  myend <- Sys.time()
  
  print(myend - mystart)
  
  message(paste("Stage ", v, " of ", Fend, " done!", sep = ""))
  
}

# ==============================================================================
# post-process
# ==============================================================================

alph = 0.05

dfc_metrics = array(0, c(Fend, 4, n))
var_metrics = array(0, c(Fend, 4, n))
dlm_metrics = array(0, c(Fend, 4, n))

for(v in 1:Fend){
  for(i in 1:n){
    obs = actuals[v, i]
    dfc_draws <- fcast_dfc[i, , v]
    var_draws <- fcast_bvar[i, , v]
    dlm_draws <- fcast_dlm[i, , v]
    dfc_metrics[v, , i] <- continuous_forecast_metrics(obs, alph, dfc_draws)
    var_metrics[v, , i] <- continuous_forecast_metrics(obs, alph, var_draws)
    dlm_metrics[v, , i] <- continuous_forecast_metrics(obs, alph, dlm_draws)
  }
}

dfc_final <- t(apply(dfc_metrics, c(2, 3), FUN = mean))
var_final <- t(apply(var_metrics, c(2, 3), FUN = mean))
dlm_final <- t(apply(dlm_metrics, c(2, 3), FUN = mean))

fcast_comparison <- array(0, c(3, 4, n))

for(i in 1:n){
  fcast_comparison[1, , i] <- dfc_final[i, ]
  fcast_comparison[2, , i] <- var_final[i, ]
  fcast_comparison[3, , i] <- dlm_final[i, ]
}

paper_table = round(rbind(cbind(fcast_comparison[, , 1], fcast_comparison[, , 2]),
      cbind(fcast_comparison[, , 3], fcast_comparison[, , 4]),
      cbind(fcast_comparison[, , 5], fcast_comparison[, , 6]),
      cbind(fcast_comparison[, , 7], fcast_comparison[, , 8]),
      cbind(fcast_comparison[, , 9], fcast_comparison[, , 10]),
      cbind(fcast_comparison[, , 11], fcast_comparison[, , 12]),
      cbind(fcast_comparison[, , 13], fcast_comparison[, , 14])), digits = 2)

write.csv(paper_table, "macro_table.csv", row.names = F)
