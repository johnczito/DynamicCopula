# ==============================================================================
# Does the single-site truncated Gibbs update mix, in the simplest possible
# setting?
#
# Dan's univariate static case: z_i iid N(0,1), y_i = F^{-1}(Phi(z_i)).
# The margin adjustment gives Fhat(y_(k)) = u_(k) = Phi(z_(k)), so the exact
# posterior for F at the data is the k-th order statistic of n uniforms,
# Beta(k, n-k+1). That is also the sampling law of F(y_(k)), so the exact band
# covers at nominal. This is confirmed analytically and is the `exact` column
# below
#
# The DGFC does not sample that posterior directly. It updates each latent
# singly, truncated to an interval set by its neighbouring ranks. This script
# runs exactly that update on the univariate static problem, with no model, no
# factors and no dynamics attached, and asks whether the resulting band matches
# the exact one
#
#   gibbs coverage at nominal, width ratio near 1
#     -> single-site updating mixes fine here, and whatever is happening in the
#        DGFC comes from the factor or dynamic structure rather than from the
#        update itself
#
#   gibbs coverage short of nominal, width ratio below 1, recovering as thin
#   grows
#     -> the mechanism is demonstrated in fifteen lines with no model attached,
#        and the remedy is a blocked or reparameterized update of z
#
# The chain is initialised at a draw from the target, so there is no burn-in
# question. Anything below nominal is autocorrelation in the retained sample
#
# Run from anywhere. Depends only on base R plus future/future.apply
# ==============================================================================

library(future)
library(future.apply)

plan(multisession, workers = max(1, availableCores() - 1))

# ==============================================================================
# settings
# ==============================================================================

NS     <- c(25, 100, 400) # sample sizes
THINS  <- c(2, 20) # retained every `thin`-th sweep
NKEEP  <- 500 # retained draws, fixed, so quantiles compare
NREP   <- 30 # replicate datasets per cell
LEVEL  <- 0.95

cat("cells: ", length(NS) * length(THINS),
    "   replicates each: ", NREP,
    "   retained draws: ", NKEEP, "\n", sep = "")
cat("total sweeps per fit: ", paste(THINS * NKEEP, collapse = ", "), "\n\n",
    sep = "")

# ==============================================================================
# the sampler under test
#
# One sweep updates z_(1) through z_(n) in order, each from a standard normal
# truncated to its neighbours. This mirrors the rank-level loop in DGFC.mcmc,
# which walks rank levels in order and uses the already-updated neighbour below
# and the not-yet-updated neighbour above
# ==============================================================================

gibbs_static <- function(n, nkeep, thin) {
  z <- sort(rnorm(n)) # a draw from the target
  U <- matrix(0, n, nkeep)
  for (m in seq_len(nkeep)) {
    for (s in seq_len(thin)) {
      for (k in seq_len(n)) {
        lo <- if (k == 1L) -Inf else z[k - 1L]
        hi <- if (k == n)   Inf else z[k + 1L]
        z[k] <- qnorm(runif(1, pnorm(lo), pnorm(hi)))
      }
    }
    U[, m] <- pnorm(z)
  }
  U
}

# ==============================================================================
# one replicate of one cell
# ==============================================================================

run_one <- function(task, s) {

  n    <- task$n
  thin <- task$thin
  a    <- (1 - s$level) / 2

  y  <- sort(rnorm(n)) # any continuous F works; standard normal
  Ft <- pnorm(y)  # the truth at the observed points

  U  <- gibbs_static(n, s$nkeep, thin)
  lo <- apply(U, 1, quantile, probs = a,     names = FALSE)
  hi <- apply(U, 1, quantile, probs = 1 - a, names = FALSE)

  k     <- seq_len(n) # exact pointwise posterior, no simulation
  eb_lo <- qbeta(a,     k, n - k + 1)
  eb_hi <- qbeta(1 - a, k, n - k + 1)

  data.frame(n          = n,
             sweeps     = thin * s$nkeep,
             rep        = task$rep,
             cov_gibbs  = mean(Ft >= lo    & Ft <= hi),
             cov_exact  = mean(Ft >= eb_lo & Ft <= eb_hi),
             w_gibbs    = mean(hi - lo),
             w_exact    = mean(eb_hi - eb_lo),
             stringsAsFactors = FALSE)
}

# ==============================================================================
# run
# ==============================================================================

settings <- list(nkeep = NKEEP, level = LEVEL)

tasks <- expand.grid(n = NS, thin = THINS, rep = seq_len(NREP))
tasks <- split(tasks, seq_len(nrow(tasks)))

set.seed(8675309)
t0 <- Sys.time()

res <- future_lapply(tasks, run_one, s = settings,
                     future.seed = TRUE, future.chunk.size = 1)

cat("Total: ", round(difftime(Sys.time(), t0, units = "mins"), 2), " min\n",
    sep = "")

R <- do.call(rbind, res)
saveRDS(R, "gibbs_static_check.rds")

# ==============================================================================
# read-out
# ==============================================================================

agg <- aggregate(cbind(cov_gibbs, cov_exact, w_gibbs, w_exact) ~ n + sweeps,
                 data = R, FUN = mean)
agg <- agg[order(agg$n, agg$sweeps), ]
agg$width_ratio <- round(agg$w_gibbs / agg$w_exact, 3)

out <- data.frame(
  n           = agg$n,
  sweeps      = agg$sweeps,
  `cov gibbs` = round(agg$cov_gibbs, 3),
  `cov exact` = round(agg$cov_exact, 3),
  `width ratio` = agg$width_ratio,
  check.names = FALSE)

cat("\n--- single-site Gibbs against the exact Beta posterior ---\n")
cat("nominal coverage ", LEVEL, "; width ratio is the sampled band over the\n",
    "exact band, so below 1 means the sampler is under-dispersed\n\n", sep = "")
print(out, row.names = FALSE)

cat("\nThe `cov exact` column is a control and should sit at ", LEVEL,
    " everywhere.\n", sep = "")
cat("If it does not, something is wrong with the framework rather than the\n")
cat("sampler, and nothing else in the table means anything.\n\n")

cat("Read `cov gibbs` and `width ratio` together:\n")
cat("  both at nominal / near 1 even at n = 400\n")
cat("    -> single-site updating is not the problem; the DGFC's behaviour\n")
cat("       comes from the factor and dynamic structure.\n")
cat("  ratio below 1, falling with n, recovering as sweeps rise\n")
cat("    -> the mechanism, shown with no model attached, and the remedy is a\n")
cat("       blocked or reparameterized update rather than a longer chain.\n")
