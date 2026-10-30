# ==============================================================================
# Does the spline's upper-tail padding cost CRPS?
#
# On the spline path DGFC.forecast pads a knot out to max(xj) + sd(xj), so a
# few percent of predictive mass lands beyond anything ever observed. CRPS
# sees that mass; the 95% HDI mostly does not. That would explain how the
# DGFC can be sharper and better covered at the band while losing on CRPS.
#
# This fits 20 vintages spread across the sample and forecasts each one twice
# off the SAME fit, with the RNG state restored between calls. The latent
# forward simulations are therefore identical and the only difference is the
# top knot, so it is a paired comparison rather than a noisy one.
#
# Requires the pad_mult argument added to DGFC.forecast. pad_mult = 1 is
# current behaviour; pad_mult = 0 anchors the top knot at max(xj).
#
# Run from the repo root.
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

library(future)
library(future.apply)

REPO <- getwd()
plan(multisession, workers = max(1, availableCores() - 1))

set.seed(8675309)

continuous_forecast_metrics <- function(obs, alph, fcast_draws){
  fcast_mean <- mean(fcast_draws)
  fcast_hdi  <- hdi(fcast_draws, 1 - alph)
  c((fcast_mean - obs)^2,
    (fcast_hdi[1] <= obs) & (obs <= fcast_hdi[2]),
    fcast_hdi[2] - fcast_hdi[1],
    scoringRules::crps_sample(obs, fcast_draws))
}

# ==============================================================================
# data, lifted from the macro script
# ==============================================================================

var_names <- c("gdpgrowth", "pcegrowth", "bfigrowth", "resinvgrowth", "ipgrowth",
               "capu", "employ", "hours", "ur", "gdppigrowth", "pcepigrowth",
               "ffr", "spread", "realSPgrowth")

unrevised    <- read.csv("datasets/ccm_2016_JBES_ur_ffr_spread_1964Q1_2014Q1.csv", header = FALSE)
bfigrowth    <- read.csv("datasets/ccm_2016_JBES_bfigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
capu         <- read.csv("datasets/ccm_2016_JBES_capu_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
employ       <- read.csv("datasets/ccm_2016_JBES_employ_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
gdpgrowth    <- read.csv("datasets/ccm_2016_JBES_gdpgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
gdppigrowth  <- read.csv("datasets/ccm_2016_JBES_gdppigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
hours        <- read.csv("datasets/ccm_2016_JBES_hours_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
ipgrowth     <- read.csv("datasets/ccm_2016_JBES_ipgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
pcegrowth    <- read.csv("datasets/ccm_2016_JBES_pcegrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
pcepigrowth  <- read.csv("datasets/ccm_2016_JBES_pcepigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
realSPgrowth <- read.csv("datasets/ccm_2016_JBES_realSPgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
resinvgrowth <- read.csv("datasets/ccm_2016_JBES_resinvgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)

dates    <- seq(1964.00, 2014.00, by = 0.25)
vintages <- seq(1985.00, 2014.25, by = 0.25)

n    <- length(var_names)
Tmax <- nrow(unrevised)
V    <- length(vintages)
Fend <- V - 2

real_time_macro_data <- array(0, c(Tmax, n, V))

for(v in 1:V){
  T <- which(dates == vintages[v] - 0.25)
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
# settings
# ==============================================================================

NVINT  <- 20                            # vintages to test
PADS   <- c(1, 0)                       # 1 = current, 0 = no extrapolation
vsel   <- unique(round(seq(1, Fend, length.out = NVINT)))

settings <- list(n = n, repo = REPO, dates = dates, vintages = vintages,
                 var_names = var_names, pads = PADS,
                 dfc_ndraw = 5000, dfc_nburn = 5000, dfc_nthin = 1,
                 factor_dim = n)

cat("Testing pad_mult in", paste(PADS, collapse = " and "), "on",
    length(vsel), "vintages\n")
cat("vintages:", paste(vsel, collapse = ", "), "\n\n")

# ==============================================================================
# one vintage, both padding values off the same fit
# ==============================================================================

run_one <- function(v, RTMD, s) {

  if (!exists("DGFC.mcmc", mode = "function")) {
    source(file.path(s$repo, "_packages.R"))
    source(file.path(s$repo, "_helpers", "_helpers.R"))
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
  }

  t0 <- Sys.time()

  out <- tryCatch({

    TT <- which(s$dates == s$vintages[v] - 0.25)
    Y  <- RTMD[1:TT, , v]

    draws <- DGFC.mcmc(Y, k.star = s$factor_dim, ndraw = s$dfc_ndraw,
                       burn = s$dfc_nburn, thin = s$dfc_nthin)

    # Restore the RNG state before each forecast so the latent forward
    # simulations are identical and only the top knot differs.
    rng <- if (exists(".Random.seed", envir = .GlobalEnv))
             get(".Random.seed", envir = .GlobalEnv) else NULL

    fc <- lapply(s$pads, function(pm) {
      if (!is.null(rng)) assign(".Random.seed", rng, envir = .GlobalEnv)
      DGFC.forecast(1, draws, use_spline = TRUE, pad_mult = pm)[1, , ]
    })
    names(fc) <- paste0("pad", s$pads)

    list(v = v, fc = fc, actuals = RTMD[TT + 1, , v + 2])

  }, error = function(e) list(error = conditionMessage(e), v = v))

  message("vintage ", v, if (is.null(out$error)) " done in " else " FAILED after ",
          round(difftime(Sys.time(), t0, units = "mins"), 2), " min")
  out
}

run_start <- Sys.time()

results <- future_lapply(
  vsel, run_one, RTMD = real_time_macro_data, s = settings,
  future.seed = TRUE, future.chunk.size = 1,
  future.packages = c("MASS", "MCMCpack", "mvtnorm", "sn", "matrixNormal",
                      "KFAS", "HDInterval", "scoringRules")
)

message("Total: ", round(difftime(Sys.time(), run_start, units = "mins"), 2), " min")

ok <- vapply(results, function(r) is.null(r$error), logical(1))
if (any(!ok)) {
  message("FAILED: ", paste(vsel[!ok], collapse = ", "))
  message("First error: ", results[[which(!ok)[1]]]$error)
}
results <- results[ok]

# ==============================================================================
# score
# ==============================================================================

alph <- 0.05
nv   <- length(results)

met <- lapply(PADS, function(pm) array(0, c(nv, 4, n)))
names(met) <- paste0("pad", PADS)

for (r in seq_len(nv)) {
  for (i in 1:n) {
    obs <- results[[r]]$actuals[i]
    for (pm in paste0("pad", PADS)) {
      met[[pm]][r, , i] <- continuous_forecast_metrics(obs, alph,
                                                       results[[r]]$fc[[pm]][i, ])
    }
  }
}

fin <- lapply(met, function(m) t(apply(m, c(2, 3), mean)))

cmp <- data.frame(
  series    = var_names,
  crps_pad1 = round(fin$pad1[, 4], 4),
  crps_pad0 = round(fin$pad0[, 4], 4),
  crps_pct  = round(100 * (fin$pad0[, 4] / fin$pad1[, 4] - 1), 2),
  size_pad1 = round(fin$pad1[, 3], 3),
  size_pad0 = round(fin$pad0[, 3], 3),
  cov_pad1  = round(fin$pad1[, 2], 3),
  cov_pad0  = round(fin$pad0[, 2], 3),
  mse_pct   = round(100 * (fin$pad0[, 1] / fin$pad1[, 1] - 1), 2),
  stringsAsFactors = FALSE
)
rownames(cmp) <- NULL

cat("\n--- pad_mult = 0 vs 1, on", nv, "vintages ---\n")
cat("crps_pct negative = removing the padding IMPROVES CRPS\n\n")
print(cmp[order(cmp$crps_pct), ])

cat("\nseries improved by dropping the pad:",
    sum(cmp$crps_pct < 0), "of", n, "\n")
cat("median CRPS change: ", round(median(cmp$crps_pct), 2), "%\n", sep = "")
cat("mean coverage: pad1 ", round(mean(fin$pad1[, 2]), 4),
    "  pad0 ", round(mean(fin$pad0[, 2]), 4), "\n", sep = "")

cat("\nFor reference, the SV BVAR beats the DGFC by a median of about 5%",
    "\non CRPS in the full run, so a median improvement near that size here",
    "\nwould be worth a full rerun. Under 1% and the padding is not the story.\n")

saveRDS(list(cmp = cmp, met = met, fin = fin, vsel = vsel, pads = PADS),
        "pad_mult_test.rds")
