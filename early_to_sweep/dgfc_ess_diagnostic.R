# ==============================================================================
# DGFC mixing diagnostic
#
# One long fit on the full crime sample, to decide the Stage 4 MCMC settings
# empirically rather than by guessing a multiple of the current ones.
#
# Monitors IDENTIFIED quantities only. Individual entries of Lambda, G and
# Sigma are not identified in a factor model (rotations and sign flips leave
# the likelihood unchanged), so their traces wander even when the chain is
# mixing fine, and their ESS means nothing. What is identified:
#
# - Zcov = Lambda Gam Lambda' + diag(v), the latent covariance
# - C0   = cov2cor(Zcov), the stationary correlation
# - v, idiosyncratic variances
# - the margin adjustment, which is what forecasts actually depend on
#
# Run this from the repo root.
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

library(coda)

set.seed(8675309)

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------

NDRAW <- 5000
NBURN <- 5000
NTHIN <- 1
KSTAR <- 7 # the formula's value for n = 10, per Section 3.4

OUT_PREFIX <- "_images/dgfc_ess"

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

crime_data <- read.csv("datasets/full_nsw_crime_dataset.csv", header = TRUE)

var_names <- c("murder", "attempted", "accessory_murder", "manslaughter",
               "assault_nondomestic", "assault_police", "abduction",
               "blackmail", "Escape_custody", "resist_officer")

my_labs <- c("Murder", "Attempted", "Accessory", "Manslaughter",
             "Non-domestic assault", "Officer assault", "Abduction",
             "Blackmail", "Escape custody", "Resist officer")

Y <- as.matrix(crime_data[, var_names])
n <- ncol(Y)

# ------------------------------------------------------------------------------
# the fit
# ------------------------------------------------------------------------------

cat("Fitting: ndraw =", NDRAW, " burn =", NBURN, " thin =", NTHIN,
    " k.star =", KSTAR, "\n")

t0 <- Sys.time()
draws <- DGFC.mcmc(Y, k.star = KSTAR,
                   ndraw = NDRAW, burn = NBURN, thin = NTHIN)
fit_mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

cat("Done in", round(fit_mins, 1), "min\n\n")

# ------------------------------------------------------------------------------
# rebuild identified quantities, draw by draw
# ------------------------------------------------------------------------------

M <- dim(draws$Lambda)[3]

n_offdiag <- n * (n - 1) / 2
C0_chain  <- matrix(0, M, n_offdiag)
v_chain   <- matrix(0, M, n)

offdiag_labs <- character(n_offdiag)
k <- 1
for (i in 1:(n - 1)) for (j in (i + 1):n) {
  offdiag_labs[k] <- paste0(my_labs[i], " / ", my_labs[j]); k <- k + 1
}

for (m in 1:M) {
  Lam  <- draws$Lambda[, , m]
  Gam  <- stationary_var1_covariance(draws$G[, , m], draws$S[, , m])
  Zcov <- Lam %*% Gam %*% t(Lam) + diag(draws$sig[, m])
  C0   <- cov2cor(Zcov)
  
  C0_chain[m, ] <- C0[lower.tri(C0)]
  v_chain[m, ]  <- draws$sig[, m]
}

colnames(C0_chain) <- offdiag_labs
colnames(v_chain)  <- my_labs

# ------------------------------------------------------------------------------
# margin adjustment at a few quantiles of each series
#
# draws$ma[[j]] is (unique values) x 2 x M: column 1 the y values,
# column 2 the CDF. Track F-hat at the 25th, 50th and 75th percentile
# of each series.
# ------------------------------------------------------------------------------

probs   <- c(0.25, 0.50, 0.75)
ma_cols <- n * length(probs)
ma_chain <- matrix(0, M, ma_cols)
ma_labs  <- character(ma_cols)

k <- 1
for (j in 1:n) {
  yq <- quantile(Y[, j], probs = probs, type = 1)
  for (q in seq_along(probs)) {
    xj  <- draws$ma[[j]][, 1, 1]            # support, constant across draws
    idx <- which.min(abs(xj - yq[q]))
    ma_chain[, k] <- draws$ma[[j]][idx, 2, ]
    ma_labs[k] <- paste0(my_labs[j], " q", probs[q] * 100)
    k <- k + 1
  }
}
colnames(ma_chain) <- ma_labs

# constant columns (discrete series where a quantile sits at a mass point
# with no posterior spread) break effectiveSize; drop them
keep     <- apply(ma_chain, 2, function(z) sd(z) > 0)
ma_chain <- ma_chain[, keep, drop = FALSE]

# ------------------------------------------------------------------------------
# effective sample sizes
# ------------------------------------------------------------------------------

ess_of <- function(X) {
  e <- effectiveSize(as.mcmc(X))
  names(e) <- colnames(X)
  e
}

ess_C0 <- ess_of(C0_chain)
ess_v  <- ess_of(v_chain)
ess_ma <- ess_of(ma_chain)

report <- function(label, e) {
  cat(sprintf("%-28s min %7.1f   median %7.1f   max %7.1f\n",
              label, min(e), median(e), max(e)))
  worst <- sort(e)[1:min(3, length(e))]
  for (nm in names(worst)) cat(sprintf("      %-40s %7.1f\n", nm, worst[nm]))
}

cat("Effective sample size, out of", M, "retained draws\n")
cat(strrep("-", 70), "\n")
report("latent correlations",   ess_C0)
report("idiosyncratic variances", ess_v)
report("margin adjustment",     ess_ma)

all_ess   <- c(ess_C0, ess_v, ess_ma)
worst_ess <- min(all_ess)

cat("\nWorst ESS overall:", round(worst_ess, 1), "of", M, "\n")
cat("Efficiency:", round(worst_ess / M, 4), "effective draws per retained draw\n")

needed <- ceiling(1000 / (worst_ess / M))
cat("\nFor 1000 effective draws on the worst quantity you would need about\n")
cat("  ", needed, "retained draws (", round(needed * fit_mins / M, 1),
    "min per fit at this k.star)\n\n")

# ------------------------------------------------------------------------------
# trace plots for the three worst
# ------------------------------------------------------------------------------

worst3 <- names(sort(all_ess)[1:3])

png(paste0(OUT_PREFIX, "_traces.png"),
    width = 6.5, height = 6, units = "in", res = 300)
par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))

for (nm in worst3) {
  z <- if (nm %in% colnames(C0_chain)) C0_chain[, nm]
  else if (nm %in% colnames(v_chain))  v_chain[, nm]
  else                                  ma_chain[, nm]
  plot(z, type = "l", xlab = "", ylab = "",
       main = paste0(nm, "  (ESS ", round(all_ess[nm], 1), ")"))
}

dev.off()

cat("Trace plots for the three worst written to ", OUT_PREFIX, "_traces.png\n\n",
    sep = "")

saveRDS(draws, "dgfc_full_fit.rds")
saveRDS(list(ess_C0 = ess_C0, ess_v = ess_v, ess_ma = ess_ma,
             M = M, fit_mins = fit_mins,
             settings = c(ndraw = NDRAW, burn = NBURN,
                          thin = NTHIN, k.star = KSTAR)),
        "dgfc_ess_diagnostic.rds")

m <- 100 # any draw
c <- 3 # any positive scale
j <- 7 # Abduction

Lam <- draws$Lambda[, , m]
sig <- draws$sig[, m]
Gam <- stationary_var1_covariance(draws$G[, , m], draws$S[, , m])

# original
Zcov <- Lam %*% Gam %*% t(Lam) + diag(sig)
u    <- pnorm(draws$Z[, j, m], 0, sqrt(Zcov[j, j]))

# rescaled: row j of Lambda by c, v_j by c^2, z_j by c
Lam2 <- Lam; Lam2[j, ] <- c * Lam[j, ]
sig2 <- sig; sig2[j]   <- c^2 * sig[j]
Zcov2 <- Lam2 %*% Gam %*% t(Lam2) + diag(sig2)
u2    <- pnorm(c * draws$Z[, j, m], 0, sqrt(Zcov2[j, j]))

max(abs(u - u2))                           # ~1e-16
max(abs(cov2cor(Zcov) - cov2cor(Zcov2)))   # ~1e-16

# Pasted console results:
# Fitting: ndraw = 5000  burn = 5000  thin = 1  k.star = 7 
# Done in 2.9 min
# 
# Effective sample size, out of 5000 retained draws
# ---------------------------------------------------------------------- 
#   latent correlations          min  2890.0   median  4145.5   max  5000.0
# Accessory / Officer assault               2890.0
# Abduction / Blackmail                     2934.0
# Officer assault / Blackmail               3001.3
# idiosyncratic variances      min     3.5   median    19.4   max   237.5
# Abduction                                    3.5
# Manslaughter                                 4.8
# Murder                                       5.1
# margin adjustment            min     6.9   median   197.8   max  2542.1
# Officer assault q50                          6.9
# Resist officer q50                          16.3
# Blackmail q50                               16.4
# 
# Worst ESS overall: 3.5 of 5000 
# Efficiency: 7e-04 effective draws per retained draw
# 
# For 1000 effective draws on the worst quantity you would need about
# 1431043 retained draws ( 840.5 min per fit at this k.star)
# 
# null device 
# 1 
# Trace plots for the three worst written to _images/dgfc_ess_traces.png
# 
# [1] 1.110223e-16
# [1] 1.665335e-16