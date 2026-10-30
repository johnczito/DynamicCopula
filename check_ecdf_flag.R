# ==============================================================================
# ECDF ablation: the statistical comparison, done per series
#
# check_ecdf_flag.R establishes that the edit is inert and that the ECDF path
# runs. This script reads the actual effect.
#
# Two things it fixes relative to check 4 of that script:
#
#   - per-series normalization. Crime runs from single-digit counts to
#     thousands, so pooling scales across series makes any ratio meaningless.
#   - the right spline setting. Splining count data changes the interpolant
#     between knots for reasons unrelated to the uncertainty being removed,
#     which contaminates the comparison. Crime is read with the spline OFF,
#     as it actually runs. Macro is read with it ON, likewise.
#
# Both arms run at the same seed, so the latent forward simulations are
# identical and only the margin transform differs. Paired, not noisy.
#
# Run from the repo root.
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

library(coda)

NDRAW <- 5000
NBURN <- 1000
H     <- 5
SEED  <- 1

# ------------------------------------------------------------------------------
# the comparison, for one fit
# ------------------------------------------------------------------------------

compare_arms <- function(Y, draws, H, use_spline, labs = colnames(Y)) {
  
  n <- ncol(Y)
  
  f <- function(...) { set.seed(SEED); DGFC.forecast(H, draws, ...) }
  
  ma <- f(use_spline = use_spline)
  ec <- f(use_spline = use_spline, use_ecdf = TRUE, Y = Y)
  
  wid <- function(X) apply(X, c(1, 2), function(z) diff(quantile(z, c(.025, .975))))
  med <- function(X) apply(X, c(1, 2), median)
  
  w_ma <- wid(ma); w_ec <- wid(ec)
  m_ma <- med(ma); m_ec <- med(ec)
  
  # support granularity: typical spacing between adjacent observed values.
  # A width change smaller than this is below the resolution of the data
  # and must not be read as an effect.
  gran <- sapply(seq_len(n), function(j) {
    xj <- sort(unique(Y[, j]))
    if(length(xj) < 2) NA_real_ else median(diff(xj))
  })
  
  out <- data.frame(
    series        = labs,
    scale         = round(apply(Y, 2, median), 1),
    n_uniq        = sapply(seq_len(n), function(j) length(unique(Y[, j]))),
    width_ma      = round(apply(w_ma, 2, median), 3),
    width_ec      = round(apply(w_ec, 2, median), 3),
    width_pct     = round(100 * (apply(w_ec, 2, median) / apply(w_ma, 2, median) - 1), 1),
    resolvable    = abs(apply(w_ec, 2, median) - apply(w_ma, 2, median)) > gran,
    pct_narrower  = round(100 * colMeans(w_ec < w_ma), 0),
    med_shift_rel = round(apply(abs(m_ec - m_ma), 2, max) / apply(w_ma, 2, median), 3),
    stringsAsFactors = FALSE
  )
  
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------------------
# does the support match between arms?
#
# If the margin adjustment is defined on sort(unique(Y[, j])), then the two
# arms differ ONLY in F and the comparison is clean. If not, the support is
# a second moving part and the numbers below need reading with that in mind.
# ------------------------------------------------------------------------------

support_matches <- function(Y, draws) {
  all(sapply(seq_len(ncol(Y)), function(j)
    isTRUE(all.equal(sort(unique(Y[, j])), draws$ma[[j]][, 1, 1]))))
}

# ==============================================================================
# CRIME, spline off
# ==============================================================================

crime_data <- read.csv("datasets/full_nsw_crime_dataset.csv", header = TRUE)

var_names <- c("murder", "attempted", "accessory_murder", "manslaughter",
               "assault_nondomestic", "assault_police", "abduction",
               "blackmail", "Escape_custody", "resist_officer")

my_labs <- c("Murder", "Attempted", "Accessory", "Manslaughter",
             "Non-domestic assault", "Officer assault", "Abduction",
             "Blackmail", "Escape custody", "Resist officer")

Yc <- as.matrix(crime_data[, var_names])

cat("=== CRIME, spline off ===\n\n")
cat("Fitting, ndraw =", NDRAW, "burn =", NBURN, "k.star = 10 (saturated)\n")

set.seed(8675309)
t0 <- Sys.time()
dc <- DGFC.mcmc(Yc, k.star = ncol(Yc), ndraw = NDRAW, burn = NBURN, thin = 1)
cat("Fit took", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n\n")

cat("support identical between arms:", support_matches(Yc, dc), "\n\n")

tab_c <- compare_arms(Yc, dc, H, use_spline = FALSE, labs = my_labs)
print(tab_c)

cat("\nwidth_pct  : median interval width under ECDF vs margin adjustment.",
    "\n             Negative = narrower, which is the prediction.",
    "\npct_narrower: share of (horizon, series) cells that narrowed.",
    "\nmed_shift_rel: largest median shift, as a fraction of that series'",
    "\n             own median interval width. Small = point predictions held.\n\n")

cat("Summary:",
    sprintf("\n  resolvable series: %d of %d", sum(tab_c$resolvable), nrow(tab_c)),
    sprintf("\n  median width change, resolvable only: %+.1f%%",
            median(tab_c$width_pct[tab_c$resolvable])),
    sprintf("\n  of those, narrowed: %d", sum(tab_c$width_pct[tab_c$resolvable] < 0)),
    sprintf("\n  largest median shift, relative: %.3f", max(tab_c$med_shift_rel)),
    "\n\n")

# ------------------------------------------------------------------------------
# Is the margin adjustment's top value near 1, or is Fxj[length(xj)] = 1
# overwriting a real value and giving the two arms different top knots?
# ------------------------------------------------------------------------------

top <- sapply(seq_len(ncol(Yc)), function(j) {
  m <- apply(dc$ma[[j]][, 2, ], 2, max)
  c(min = min(m), median = median(m), max = max(m))
})
colnames(top) <- my_labs
print(round(t(top), 4))

# the two F functions side by side, for the extreme case and a clean case
for (nm in c("Blackmail", "Escape custody")) {
  j  <- which(my_labs == nm)
  xj <- sort(unique(Yc[, j]))
  Fm <- rowMeans(dc$ma[[j]][, 2, ])
  Fq <- apply(dc$ma[[j]][, 2, ], 1, quantile, c(.025, .975))
  cat("\n---", nm, "---\n")
  print(cbind(x = xj, ecdf = round(ecdf(Yc[, j])(xj), 4),
              ma_mean = round(Fm, 4),
              ma_lo = round(Fq[1, ], 4), ma_hi = round(Fq[2, ], 4)))
}

tq <- sapply(seq_len(ncol(Yc)), function(j)
  quantile(apply(dc$ma[[j]][, 2, ], 2, max), c(0, .001, .01, .05, .25, .5)))
colnames(tq) <- my_labs
round(t(tq), 4)

# ==============================================================================
# MACRO, spline on
#
# Fill in the data block from your macro script. The ECDF must come from the
# same training matrix the MCMC saw at that vintage, not the full series.
# ==============================================================================

# Ym <- <training matrix for one vintage, built the way the macro script builds it>
#
# cat("=== MACRO, spline on ===\n\n")
# set.seed(8675309)
# dm <- DGFC.mcmc(Ym, k.star = ncol(Ym), ndraw = NDRAW, burn = NBURN, thin = 1)
# cat("support identical between arms:", support_matches(Ym, dm), "\n\n")
# tab_m <- compare_arms(Ym, dm, H = 1, use_spline = TRUE)
# print(tab_m)