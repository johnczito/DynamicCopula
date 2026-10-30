# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

library(future)
library(future.apply)

#[CHANGED]
REPO <- getwd()
plan(multisession, workers = max(1, availableCores() - 1))

set.seed(8675309)
#set.seed(42065)

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

unrevised <- read.csv("datasets/ccm_2016_JBES_ur_ffr_spread_1964Q1_2014Q1.csv", header = FALSE)
bfigrowth <- read.csv("datasets/ccm_2016_JBES_bfigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
capu <- read.csv("datasets/ccm_2016_JBES_capu_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
employ <- read.csv("datasets/ccm_2016_JBES_employ_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
gdpgrowth <- read.csv("datasets/ccm_2016_JBES_gdpgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
gdppigrowth <- read.csv("datasets/ccm_2016_JBES_gdppigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
hours <- read.csv("datasets/ccm_2016_JBES_hours_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
ipgrowth <- read.csv("datasets/ccm_2016_JBES_ipgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
pcegrowth <- read.csv("datasets/ccm_2016_JBES_pcegrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
pcepigrowth <- read.csv("datasets/ccm_2016_JBES_pcepigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
realSPgrowth <- read.csv("datasets/ccm_2016_JBES_realSPgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
resinvgrowth <- read.csv("datasets/ccm_2016_JBES_resinvgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)

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

#bvar_ndraw = 1000
#bvar_nburn = 1000
bvar_ndraw = 5000
bvar_nburn = 5000
bvar_nthin = 1

#[CHANGED]
bvar_lag = 4
svbvar_ndraw  = 5000
svbvar_nburn  = 5000
svbvar_lag    = 4 # match the existing BVAR
svbvar_stable = TRUE
factor_dim = n # saturated, per Dan

#dfc_ndraw = 2000
#dfc_nburn = 0
dfc_ndraw = 5000
dfc_nburn = 5000
dfc_nthin = 1

dlm_ndraw = 2000

# ==============================================================================
# preallocate storage
# ==============================================================================

Fend = V - 2
#Fend = 10 # For short-run

fcast_bvar = array(0, c(n, bvar_ndraw, Fend))
fcast_dfc  = array(0, c(n, dfc_ndraw, Fend))

fcast_dfc_ecdf <- array(0, c(n, dfc_ndraw, Fend))
fcast_svbvar   <- vector("list", Fend)
svbvar_nstable <- numeric(Fend)
ma_top_q       <- array(0, c(5, n, Fend))

fcast_dlm =  array(0, c(n, dlm_ndraw, Fend))
actuals = matrix(0, Fend, n)

CKPT_DIR <- sprintf("_ckpt_macro_k%d", factor_dim)
dir.create(CKPT_DIR, showWarnings = FALSE)

settings <- list(
  n            = n,
  Fend         = Fend,
  vintages     = vintages,
  var_names    = var_names,
  dates        = dates,
  repo         = REPO,
  ckpt_dir     = CKPT_DIR,
  dlm_ndraw    = dlm_ndraw,
  bvar_ndraw   = bvar_ndraw,
  bvar_nburn   = bvar_nburn,
  bvar_nthin   = bvar_nthin,
  bvar_lag     = bvar_lag,
  svbvar_ndraw  = svbvar_ndraw,
  svbvar_nburn  = svbvar_nburn,
  svbvar_lag    = svbvar_lag,
  svbvar_stable = svbvar_stable,
  dfc_ndraw    = dfc_ndraw,
  dfc_nburn    = dfc_nburn,
  dfc_nthin    = dfc_nthin,
  factor_dim   = factor_dim
)

# ==============================================================================
# run it hot!
# ==============================================================================

#[CHANGED] for loop is now a function
run_vintage <- function(v, RTMD, s) {
  
  ckpt <- file.path(s$repo, s$ckpt_dir, sprintf("vintage_%04d.rds", v))
  if (file.exists(ckpt)) {
    cached <- readRDS(ckpt)
    complete <- is.null(cached$error) &&
      all(c("dlm", "bvar", "svbvar", "dfc", "dfc_ecdf") %in% names(cached))
    if (complete) {
      message("Vintage ", v, " of ", s$Fend, " loaded from checkpoint")
      return(cached)
    }
    message("Vintage ", v, " checkpoint is stale or incomplete; recomputing")
  }
  
  if (!exists("DGFC.mcmc", mode = "function")) {
    source(file.path(s$repo, "_packages.R"))
    source(file.path(s$repo, "_helpers", "_helpers.R"))
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
  }
  
  v0 <- Sys.time()
  
  out <- tryCatch({
  
  end_period <- s$vintages[v] - 0.25
    
  TT <- which(s$dates == end_period)
  Y  <- RTMD[1:TT, , v]
  
  o <- list()
  o$actuals <- RTMD[TT + 1, , v + 2]
  
  # ---------------------------
  # DLM
  # ---------------------------
  
  inits <- c(numeric(sum(1:s$n)), log(diag(cov(Y))), numeric(sum(1:(s$n - 1))))
  
  my_dlm <- SSModel(Y ~ SSMtrend(1, Q = matrix(NA, s$n, s$n), type = "distinct"),
                    H = matrix(NA, s$n, s$n))
  fitted_model  <- fitSSM(my_dlm, inits = inits)$model
  filter_output <- KFS(fitted_model)
  
  fcast_mean <- filter_output$att[TT, ]
  fcast_cov  <- filter_output$Ptt[, , TT] + fitted_model$H[, , 1] + fitted_model$Q[, , 1]
  
  o$dlm <- t(mvrnorm(s$dlm_ndraw, fcast_mean, fcast_cov))
  
  # ---------------------------
  # BVAR
  # ---------------------------
  
  mybvar <- bvar(Y,
                 lags = s$bvar_lag,
                 verbose = FALSE,
                 n_draw = s$bvar_ndraw * s$bvar_nthin + s$bvar_nburn,
                 n_burn = s$bvar_nburn,
                 n_thin = s$bvar_nthin)
  
  o$bvar <- t(predict(mybvar, horizon = 1)$fcast[, 1, ])
  
  # ---------------------------
  # SV BVAR (Gruber & Kastner)
  #
  # Namespace-qualified: the BVAR package also exports bvar().
  # stable = TRUE keeps only stationary draws, matching the restriction the
  # DGFC imposes on its own factor VAR by rejection sampling. The surviving
  # count varies a lot by vintage and is worth recording.
  # ---------------------------
  
  Ysv <- Y
  colnames(Ysv) <- s$var_names
  
  sv_fit <- bayesianVARs::bvar(Ysv,
                               lags    = s$svbvar_lag,
                               draws   = s$svbvar_ndraw,
                               burnin  = s$svbvar_nburn,
                               sv_keep = "last",
                               quiet   = TRUE)
  
  sv_pred <- predict(sv_fit, ahead = 1, each = 1, stable = s$svbvar_stable)
  
  o$svbvar         <- sv_pred$predictions[1, , ]
  o$svbvar_nstable <- dim(sv_pred$predictions)[3]
  
  rm(sv_fit, sv_pred, Ysv)
  
  # ---------------------------
  # DGFC
  # ---------------------------
  
  draws <- DGFC.mcmc(Y,
                     k.star = s$factor_dim,
                     ndraw = s$dfc_ndraw,
                     burn = s$dfc_nburn,
                     thin = s$dfc_nthin)
  
  # o$dfc <- DGFC.forecast(1, draws, use_spline = TRUE)[1, , ]
  # Both arms off the same fit. The RNG state is restored between calls so
  # the latent forward simulations are identical and the only difference is
  # the margin transform. That makes the ablation a paired comparison.
  rng_state <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  
  o$dfc <- DGFC.forecast(1, draws, use_spline = TRUE)[1, , ]
  
  if (!is.null(rng_state)) assign(".Random.seed", rng_state, envir = .GlobalEnv)
  
  o$dfc_ecdf <- DGFC.forecast(1, draws, use_spline = TRUE,
                              use_ecdf = TRUE, Y = Y)[1, , ]
  
  # F at the largest observed value, quantiles over draws. Cheap to store,
  # and it tells you whether the margin-adjusted arm's upper tail came from
  # the sampler or from the padding.
  o$ma_top_q <- sapply(seq_len(s$n), function(j) {
    mm <- draws$ma[[j]]
    quantile(mm[dim(mm)[1], 2, ], c(0, .01, .05, .25, .5))
  })
  
  o
  
  }, error = function(e) {
    list(error = conditionMessage(e), v = v)
  })
  
  message(paste("Vintage ", v, " of ", s$Fend,
                if (is.null(out$error)) " done in " else " FAILED after ",
                round(difftime(Sys.time(), v0, units = "mins"), 2),
                " min", sep = ""))
  
  if (is.null(out$error)) saveRDS(out, ckpt)
  
  out
}

vints <- 1:Fend

run_start <- Sys.time()

#[CHANGED]
#results <- lapply(vints, run_vintage, RTMD = real_time_macro_data, s = settings)
results <- future_lapply(
  vints, run_vintage, RTMD = real_time_macro_data, s = settings,
  future.seed       = TRUE,
  future.chunk.size = 1,
  future.packages   = c("MASS", "MCMCpack", "mvtnorm", "sn", "matrixNormal",
                        "KFAS", "BVAR", "dlm", "HDInterval",
                        "scoringRules", "scoringutils")
)

message(paste("Loop total:",
              round(difftime(Sys.time(), run_start, units = "mins"), 2), "min"))

failed <- which(vapply(results, function(r) !is.null(r$error), logical(1)))
if (length(failed)) {
  message("FAILED vintages: ", paste(vints[failed], collapse = ", "))
  message("First error: ", results[[failed[1]]]$error)
  message("Rerun to retry only these; completed vintages load from checkpoint.")
} else {
  message("All vintages completed.")
}

for (i in seq_along(vints)) {
  vv <- vints[i]
  if (!is.null(results[[i]]$error)) next
  fcast_dlm[, , vv]  <- results[[i]]$dlm
  fcast_bvar[, , vv] <- results[[i]]$bvar
  fcast_dfc[, , vv]  <- results[[i]]$dfc
  fcast_dfc_ecdf[, , vv] <- results[[i]]$dfc_ecdf
  fcast_svbvar[[vv]]   <- results[[i]]$svbvar
  svbvar_nstable[vv]   <- results[[i]]$svbvar_nstable
  ma_top_q[, , vv]       <- results[[i]]$ma_top_q
  actuals[vv, ]      <- results[[i]]$actuals
}

rm(results)
gc()

save.image(sprintf("macro_k%d_after_loop.RData", factor_dim))
# ==============================================================================
# post-process
# ==============================================================================

alph = 0.05

dfc_metrics = array(0, c(Fend, 4, n))
var_metrics = array(0, c(Fend, 4, n))
dlm_metrics = array(0, c(Fend, 4, n))
svbvar_metrics = array(0, c(Fend, 4, n))

for(v in 1:Fend){
  for(i in 1:n){
    obs = actuals[v, i]
    dfc_draws <- fcast_dfc[i, , v]
    var_draws <- fcast_bvar[i, , v]
    dlm_draws <- fcast_dlm[i, , v]
    dfc_metrics[v, , i] <- continuous_forecast_metrics(obs, alph, dfc_draws)
    var_metrics[v, , i] <- continuous_forecast_metrics(obs, alph, var_draws)
    dlm_metrics[v, , i] <- continuous_forecast_metrics(obs, alph, dlm_draws)
    svbvar_metrics[v, , i] <- continuous_forecast_metrics(obs, alph, fcast_svbvar[[v]][i, ])  }
}

dfc_final <- t(apply(dfc_metrics, c(2, 3), FUN = mean))
var_final <- t(apply(var_metrics, c(2, 3), FUN = mean))
dlm_final <- t(apply(dlm_metrics, c(2, 3), FUN = mean))
svbvar_final <- t(apply(svbvar_metrics, c(2, 3), FUN = mean))

fcast_comparison <- array(0, c(4, 4, n))

for(i in 1:n){
  fcast_comparison[1, , i] <- dfc_final[i, ]
  fcast_comparison[2, , i] <- var_final[i, ]
  fcast_comparison[3, , i] <- dlm_final[i, ]
  fcast_comparison[4, , i] <- svbvar_final[i, ]
}

paper_table = round(rbind(cbind(fcast_comparison[, , 1], fcast_comparison[, , 2]),
      cbind(fcast_comparison[, , 3], fcast_comparison[, , 4]),
      cbind(fcast_comparison[, , 5], fcast_comparison[, , 6]),
      cbind(fcast_comparison[, , 7], fcast_comparison[, , 8]),
      cbind(fcast_comparison[, , 9], fcast_comparison[, , 10]),
      cbind(fcast_comparison[, , 11], fcast_comparison[, , 12]),
      cbind(fcast_comparison[, , 13], fcast_comparison[, , 14])), digits = 2)

write.csv(paper_table, sprintf("macro_table_k%d.csv", factor_dim), row.names = F)

# ==============================================================================
# ECDF ablation
# ==============================================================================

ecdf_metrics <- array(0, c(Fend, 4, n))

for(v in 1:Fend){
  for(i in 1:n){
    ecdf_metrics[v, , i] <- continuous_forecast_metrics(actuals[v, i], alph,
                                                        fcast_dfc_ecdf[i, , v])
  }
}

ecdf_final <- t(apply(ecdf_metrics, c(2, 3), FUN = mean))

ablation <- data.frame(
  series   = var_names,
  mse_ma   = round(dfc_final[, 1], 3), mse_ec  = round(ecdf_final[, 1], 3),
  cov_ma   = round(dfc_final[, 2], 3), cov_ec  = round(ecdf_final[, 2], 3),
  size_ma  = round(dfc_final[, 3], 3), size_ec = round(ecdf_final[, 3], 3),
  crps_ma  = round(dfc_final[, 4], 3), crps_ec = round(ecdf_final[, 4], 3),
  size_pct = round(100 * (ecdf_final[, 3] / dfc_final[, 3] - 1), 1),
  crps_pct = round(100 * (ecdf_final[, 4] / dfc_final[, 4] - 1), 1),
  stringsAsFactors = FALSE
)

print(ablation)
write.csv(ablation, sprintf("macro_ablation_k%d.csv", factor_dim), row.names = FALSE)

cat("\n--- F at the largest observed value, averaged over vintages ---\n")
tq <- apply(ma_top_q, c(1, 2), mean)
dimnames(tq) <- list(c("0%", "1%", "5%", "25%", "50%"), var_names)
print(round(t(tq), 4))

cat("\n--- SV BVAR stable draws out of", svbvar_ndraw, "---\n")
print(summary(svbvar_nstable))
cat("vintages below 50% stable:", sum(svbvar_nstable < svbvar_ndraw / 2), "of", Fend, "\n")

