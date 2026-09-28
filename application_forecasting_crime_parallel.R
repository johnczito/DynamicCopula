# ==============================================================================
# Crime forecasting exercise, parallelized over training/test splits.
#
# Same computation as the original loop, restructured so each forecast origin
# is an independent task.
#
# 1. As a sequential lapply, this reproduces the original for-loop bit for bit 
# (identical() on all four forecast arrays). That establishes the restructuring is faithful.
#
# 2. As future_lapply, results are reproducible across runs and independent
# of the number of workers.
#
# The second step changes the numbers relative to a sequential run, so each origin 
# now draws from its own RNG stream rather than one shared sequential stream. 
# A full run takes around ~10 minutes on 13 workers.
# ==============================================================================

# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

library(future)
library(future.apply)

REPO <- getwd()

plan(multisession, workers = max(1, availableCores() - 1))
# plan(multisession, workers = 2)


# set.seed(42065)
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

# ==============================================================================
# small subset for in-text
# ==============================================================================

png(paste("_images/crime_subset_intext.png", sep = ""), 
    width = 6.5, height = 8 * (2 / 7), units = "in", res = 650)

par(mfrow = c(2, 4))

for(i in c(4, 7, 8, 9)){
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
# Tstop = Tstart + H + 4   # ten-origin validation run

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
# settings bundle
#
# everything run_origin() needs, passed explicitly. Workers do not share the parent's 
# global environment, so anything the function needs arrives as an argument.
# ==============================================================================

settings <- list(
  H = H,
  n = n,
  Tstop = Tstop,
  repo = REPO,
  pois_ndraw = pois_ndraw,
  nbinom_ndraw = nbinom_ndraw,
  bvar_ndraw = bvar_ndraw,
  bvar_nburn = bvar_nburn,
  bvar_nthin = bvar_nthin,
  bvar_lag = bvar_lag,
  dfc_ndraw = dfc_ndraw,
  dfc_nburn = dfc_nburn,
  dfc_nthin = dfc_nthin,
  factor_dim = factor_dim,
  run_POISSON = run_POISSON,
  run_NBINOM = run_NBINOM,
  run_BVAR = run_BVAR,
  run_DGFC = run_DGFC
)

# ==============================================================================
# one forecast origin
#
# returns a list of four H x n x ndraw arrays
#
# Model order (Poisson, negative binomial, BVAR, DGFC) matches the original
# loop. It no longer affects results under per-task seeding, but keeping it
# means this file can be switched back to a plain lapply and compared against
# the original output directly.
# ==============================================================================

run_origin <- function(t, Y, s) {
  
  # One-time setup per worker. multisession workers persist across the call,
  # so this fires once each, not once per origin.
  if (!exists("DGFC.mcmc", mode = "function")) {
    source(file.path(s$repo, "_packages.R"))
    source(file.path(s$repo, "_helpers", "_helpers.R"))
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
  }
  
  t0 <- Sys.time()
  Yt <- Y[1:t, ]
  out <- list()
  
  # Poisson DGLM
  if (s$run_POISSON) {
    out$pois <- get_poisson_mvforecasts(Yt, s$H, s$pois_ndraw)
  } else {
    out$pois <- array(0, c(s$H, s$n, s$pois_ndraw))
  }
  
  # Negative binomial DGLM
  if (s$run_NBINOM) {
    out$nbinom <- get_nbinom_mvforecasts(Yt, s$H, s$nbinom_ndraw)
  } else {
    out$nbinom <- array(0, c(s$H, s$n, s$nbinom_ndraw))
  }
  
  # BVAR
  if (s$run_BVAR) {
    my_bvar_model <- bvar(Yt, lags = s$bvar_lag, verbose = FALSE,
                          n_draw = s$bvar_ndraw * s$bvar_nthin + s$bvar_nburn,
                          n_burn = s$bvar_nburn,
                          n_thin = s$bvar_nthin)
    out$bvar <- aperm(predict(my_bvar_model, horizon = s$H)$fcast, c(2, 3, 1))
  } else {
    out$bvar <- array(0, c(s$H, s$n, s$bvar_ndraw))
  }
  
  # DGFC
  if (s$run_DGFC) {
    draws <- DGFC.mcmc(Yt, k.star = s$factor_dim,
                       ndraw = s$dfc_ndraw, burn = s$dfc_nburn, thin = s$dfc_nthin)
    out$dfc <- DGFC.forecast(s$H, draws)
  } else {
    out$dfc <- array(0, c(s$H, s$n, s$dfc_ndraw))
  }
  
  message(paste("Stage ", t, " of ", s$Tstop, " done in ",
                round(difftime(Sys.time(), t0, units = "mins"), 2),
                " min", sep = ""))
  
  out
}

# ==============================================================================
# run it hot!
# ==============================================================================

origins <- Tstart:Tstop

run_start <- Sys.time()

# future.seed = TRUE gives each origin its own L'Ecuyer stream, so results are
# reproducible and do not depend on worker count. future.chunk.size = 1
# dispatches one origin at a time: cost grows with t, so the default equal
# blocks would leave most workers idle waiting on whoever got the last block.
results <- future_lapply(
  origins, run_origin, Y = Y, s = settings,
  future.seed       = TRUE,
  future.chunk.size = 1,
  future.packages   = c("MASS", "MCMCpack", "mvtnorm", "sn", "matrixNormal",
                        "KFAS", "BVAR", "dlm", "HDInterval",
                        "scoringRules", "scoringutils")
)

message(paste("Loop total:",
              round(difftime(Sys.time(), run_start, units = "mins"), 2), "min"))

# ==============================================================================
# assemble
#
# results[[i]] corresponds to origins[i], NOT to time index i.
# ==============================================================================

for (i in seq_along(origins)) {
  tt <- origins[i]
  fcast_pois[, , , tt]   <- results[[i]]$pois
  fcast_nbinom[, , , tt] <- results[[i]]$nbinom
  fcast_bvar[, , , tt]   <- results[[i]]$bvar
  fcast_dfc[, , , tt]    <- results[[i]]$dfc
}

rm(results)
gc()

# ==============================================================================
# h-step-ahead forecast table
# ==============================================================================

alpha = 0.05

save.image("par_stage2_after_loop.Rdata")

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

mastertable <- cbind(mastertable[1:(nmod * n / 2), ], mastertable[((nmod * n / 2) + 1):(nmod*n), ])

write.csv(mastertable, "crime_fcasts_par_stage2.csv", row.names = FALSE)

message(paste("Total:", round(difftime(Sys.time(), run_start, units = "mins"), 2), "min"))

# Ethan Carlson