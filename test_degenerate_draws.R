# ==============================================================================
# Do degenerate margin-adjustment draws hurt the DGFC's forecasts?
#
# The margin adjustment at the largest observed value, F-top = max_x Fhat(x),
# should sit near 1. Across draws it has a long left tail: the 5% quantile
# runs 0.78 to 0.96 on macro and the minimum is around 0.52. A draw with
# F-top = 0.52 is saying half the probability mass lies above the largest of
# ~200 observations, which no well-mixed draw should produce.
#
# This asks whether those draws are degrading the predictive distribution.
#
# Design notes:
#
#   - Dropping draws shrinks the sample, and the empirical CRPS estimator is
#     biased upward at smaller m. So every targeted filter is compared against
#     a RANDOM subsample of identical size. The targeted-minus-random
#     difference is the quantity of interest, not the raw change.
#
#   - Draws are scored by their MINIMUM F-top across series, so filtering is
#     global and the joint predictive distribution stays intact.
#
#   - The per-draw max of the latent Z is stored too. F-top is
#     Phi(max_t z_t / sqrt(Zcov[j,j])), so max-z separates a compressed latent
#     series (the data step) from an inflated covariance scale.
#
# Requires nothing beyond the current helpers.
# Run from the repo root.
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

library(future)
library(future.apply)
library(coda)

REPO <- getwd()
plan(multisession, workers = max(1, availableCores() - 1))

set.seed(8675309)

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

NVINT <- 20
vsel  <- unique(round(seq(1, Fend, length.out = NVINT)))
FRACS <- c(0.01, 0.05, 0.10, 0.25)   # fraction of draws to discard
NRAND <- 10                          # random control replicates per level

settings <- list(n = n, repo = REPO, dates = dates, vintages = vintages,
                 dfc_ndraw = 5000, dfc_nburn = 5000, dfc_nthin = 1,
                 factor_dim = n)

cat("Fitting", length(vsel), "vintages, storing per-draw F-top and max-z\n\n")

# ==============================================================================
# one vintage
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

    # forecast draw m uses MCMC draw m, so these line up by index
    fc <- DGFC.forecast(1, draws, use_spline = TRUE)[1, , ]       # n x ndraw

    # F-hat at the largest observed value. The margin adjustment is a CDF,
    # so its maximum over the support is its last row.
    ftop <- t(sapply(seq_len(s$n), function(j) {
      mm <- draws$ma[[j]]
      mm[dim(mm)[1], 2, ]
    }))                                                            # n x ndraw

    maxz <- apply(draws$Z, c(2, 3), max)                           # n x ndraw
    
    # spectral radius of the factor VAR, per draw. Gamma0 blows up as this
    # approaches 1, which would inflate Zcov and squash F-top.
    maxeig <- sapply(seq_len(dim(draws$G)[3]), function(m)
      max(Mod(eigen(draws$G[, , m], only.values = TRUE)$values)))
    
    list(v = v, fc = fc, ftop = ftop, maxz = maxz, maxeig = maxeig,
         actuals = RTMD[TT + 1, , v + 2])

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
nv <- length(results)

# ==============================================================================
# 1. how bad is the left tail, and does it cluster in the chain?
# ==============================================================================

cat("\n", strrep("=", 62), "\n", sep = "")
cat("1. The badness itself\n")
cat(strrep("=", 62), "\n\n")

fmin <- lapply(results, function(r) apply(r$ftop, 2, min))   # per vintage, length ndraw

qs <- t(sapply(fmin, quantile, c(0, .01, .05, .25, .50)))
cat("min-across-series F-top, quantiles over draws, by vintage:\n")
print(round(colMeans(qs), 4))

cat("\ndraws below 0.90:", round(100 * mean(sapply(fmin, function(f) mean(f < 0.90))), 2), "%\n")
cat("draws below 0.80:", round(100 * mean(sapply(fmin, function(f) mean(f < 0.80))), 2), "%\n")

ess <- sapply(fmin, function(f) effectiveSize(as.mcmc(f)))
cat("\nESS of the min-F-top series, out of", settings$dfc_ndraw, "draws:\n")
print(round(summary(ess), 1))
cat("\nLow ESS here means the bad draws arrive in runs, which would point at\n")
cat("the sampler wandering rather than at isolated outliers.\n")

# ==============================================================================
# 2. mechanism: compressed latent series, or inflated scale?
# ==============================================================================

cat("\n", strrep("=", 62), "\n", sep = "")
cat("2. Mechanism\n")
cat(strrep("=", 62), "\n\n")

mech <- do.call(rbind, lapply(results, function(r) {
  f  <- apply(r$ftop, 2, min)
  bad <- f < quantile(f, 0.05)
  # implied sd: F-top = Phi(maxz / sd)  =>  sd = maxz / qnorm(F-top)
  sd_imp <- r$maxz / qnorm(pmin(pmax(r$ftop, 1e-6), 1 - 1e-6))
  c(maxz_bad  = mean(r$maxz[, bad]),   maxz_good  = mean(r$maxz[, !bad]),
    sdimp_bad = mean(sd_imp[, bad]),   sdimp_good = mean(sd_imp[, !bad]))
}))

cat("averaged over vintages, worst 5% of draws vs the rest:\n")
print(round(colMeans(mech), 4))
cat("\nmax-z lower in the bad draws  -> the latent data step is compressing\n")
cat("implied sd higher in the bad draws -> the covariance scale is inflated\n")

# ==============================================================================
# 2b. is the inflated scale a near-unit-root factor VAR?
#
# Gamma0, the stationary covariance of the factor VAR, blows up as the
# spectral radius of G approaches 1. If that is the source, low F-top draws
# should have G sitting near the boundary.
# ==============================================================================

cat("\n", strrep("=", 62), "\n", sep = "")
cat("2b. Is it a near-unit-root factor VAR?\n")
cat(strrep("=", 62), "\n\n")

eig <- do.call(rbind, lapply(results, function(r) {
  f   <- apply(r$ftop, 2, min)
  cut <- quantile(f, 0.05)
  c(cor      = cor(f, r$maxeig),
    eig_bad  = mean(r$maxeig[f <  cut]),
    eig_good = mean(r$maxeig[f >= cut]),
    eig_max  = max(r$maxeig))
}))

cat("averaged over vintages:\n")
print(round(colMeans(eig), 4))

cat("\nper-vintage correlation between min-F-top and spectral radius of G:\n")
print(round(summary(eig[, "cor"]), 3))

cat("\nA strongly negative correlation, with the bad draws' spectral radius\n")
cat("near 1, confirms the inflated scale is a near-unit-root factor VAR\n")
cat("rather than anything wrong with the sampler.\n")

# ==============================================================================
# 3. does filtering help, beyond the sample-size effect?
# ==============================================================================

cat("\n", strrep("=", 62), "\n", sep = "")
cat("3. Targeted filter vs matched random subsample\n")
cat(strrep("=", 62), "\n\n")

crps_of <- function(fc, obs, idx)
  sapply(seq_along(obs), function(i) scoringRules::crps_sample(obs[i], fc[i, idx]))

tab <- NULL

for (fr in FRACS) {

  tgt <- matrix(0, nv, n)
  rnd <- matrix(0, nv, n)
  base <- matrix(0, nv, n)

  for (r in seq_len(nv)) {
    res  <- results[[r]]
    f    <- apply(res$ftop, 2, min)
    m    <- length(f)
    keep <- round((1 - fr) * m)

    base[r, ] <- crps_of(res$fc, res$actuals, seq_len(m))
    tgt[r, ]  <- crps_of(res$fc, res$actuals, order(f, decreasing = TRUE)[seq_len(keep)])

    rr <- replicate(NRAND, crps_of(res$fc, res$actuals, sample.int(m, keep)))
    rnd[r, ] <- rowMeans(rr)
  }

  b <- colMeans(base); t_ <- colMeans(tgt); r_ <- colMeans(rnd)

  tab <- rbind(tab, data.frame(
    dropped     = paste0(100 * fr, "%"),
    targeted_pct = round(100 * mean(t_ / b - 1), 3),
    random_pct   = round(100 * mean(r_ / b - 1), 3),
    difference   = round(100 * mean(t_ / b - r_ / b), 3),
    series_helped = sum(t_ < r_),
    stringsAsFactors = FALSE
  ))
}

cat("CRPS change vs using all draws, averaged across series.\n")
cat("'difference' is targeted minus random. Negative means the F-top filter\n")
cat("beats a same-size random subsample, i.e. the bad draws were hurting.\n\n")
print(tab, row.names = FALSE)

cat("\nRule of thumb: a difference under about 0.5% is noise. The SV BVAR's\n")
cat("edge is a median of roughly 5%, so anything near that size here would\n")
cat("matter and would justify a full rerun.\n")

saveRDS(list(tab = tab, qs = qs, ess = ess, mech = mech, eig = eig,
             vsel = vsel, var_names = var_names),
        "degenerate_draw_test.rds")
