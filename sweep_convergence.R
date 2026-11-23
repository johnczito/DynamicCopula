# ==============================================================================
# Does the margin band converge, or is it still widening?
#
# The 5k-against-20k check on the macro data showed the band growing about 16%
# with four times the sweeps, at both T = 50 and T = 199. That means every band
# coverage figure reported so far is a lower bound on an unknown quantity: we
# cannot say whether the pseudo-posterior is miscalibrated or merely
# under-sampled until the width plateaus
#
# This runs one DGP at T = 200 across three sweep budgets, holding the number of
# RETAINED draws fixed at 500 so the sample quantiles are estimated comparably,
# and using the SAME simulated data at all three budgets within a replicate so
# the comparison is paired. Because the DGP has closed-form margins we can
# measure coverage, not just width
#
#   width plateaus and coverage climbs to 0.95  -> computational throughout
#   width plateaus and coverage settles below   -> the pseudo-posterior is narrow
#
# Run from the repo root
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")
source("_helpers/_sim_dgp.R")
source("_helpers/_sim_metrics.R")

library(future)
library(future.apply)

plan(multisession, workers = max(1, availableCores() - 1))

# ==============================================================================
# settings
# ==============================================================================

#DGP    <- "varch"     # worst band coverage in the full run, 0.773, and continuous
DGPS_SWEEP <- c("varch", "varma", "varma_copula")
TTRAIN <- 200
NREP   <- 6 #10
NDRAW  <- 500         # retained draws, held fixed across budgets
BURNIN <- 200         # simulator burn-in, not MCMC
ALPH   <- 0.05

# (burn, thin) chosen so total sweeps = burn + thin * NDRAW is 10k, 40k, 160k
SWEEPS <- list(
  list(label =  "10k", burn =   7500, thin =  5),
  list(label =  "40k", burn =  30000, thin = 20),
  list(label = "160k", burn = 120000, thin = 80)
)

#RUN_TAG  <- sprintf("sweepconv_%s_T%d_R%d", DGP, TTRAIN, NREP)
RUN_TAG  <- sprintf("sweepconv_all_T%d_R%d", TTRAIN, NREP)
CKPT_DIR <- sprintf("_ckpt_%s", RUN_TAG)
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

settings <- list(repo = getwd(), TT = TTRAIN, ndraw = NDRAW,
                 burnin = BURNIN, alph = ALPH, sweeps = SWEEPS,
                 ckpt_dir = CKPT_DIR)

tasks <- expand.grid(dgp = DGPS_SWEEP, rep = seq_len(NREP),
                     stringsAsFactors = FALSE)
tasks <- split(tasks, seq_len(nrow(tasks)))

cfg_file <- file.path(CKPT_DIR, "_config.rds")
cfg <- settings[setdiff(names(settings), c("repo", "ckpt_dir"))]
if (file.exists(cfg_file)) {
  if (!identical(readRDS(cfg_file), cfg))
    stop("checkpoint dir ", CKPT_DIR, " holds a different configuration")
} else saveRDS(cfg, cfg_file)

cat("run tag: ", RUN_TAG, "\n", sep = "")
cat("tasks: ", length(tasks),
    "   sweep budgets: ",
    paste(vapply(SWEEPS, `[[`, "", "label"), collapse = ", "),
    "   retained draws: ", NDRAW, "\n\n", sep = "")

# ==============================================================================
# one replicate, three budgets, one dataset
# ==============================================================================

run_one <- function(task, s) {
  
  if (!exists("DGFC.mcmc", mode = "function")) {
    source(file.path(s$repo, "_packages.R"))
    source(file.path(s$repo, "_helpers", "_helpers.R"))
    source(file.path(s$repo, "_helpers", "_sim_dgp.R"))
    source(file.path(s$repo, "_helpers", "_sim_metrics.R"))
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
  }

  ckpt <- file.path(s$ckpt_dir, sprintf("%s_rep_%03d.rds", task$dgp, task$rep))
  if (file.exists(ckpt)) return(readRDS(ckpt))

  t0 <- Sys.time()

  out <- tryCatch({

    d <- make_dgp(task$dgp, burnin = s$burnin)
    Y <- d$sim(s$TT)                 # one dataset, reused at every budget

    rows <- list()

    for (sw in s$sweeps) {

      draws <- DGFC.mcmc(Y, ndraw = s$ndraw, burn = sw$burn, thin = sw$thin)

      for (j in seq_len(d$n)) {
        mm <- margin_metrics(draws$ma[[j]], d$F[[j]], 1 - s$alph)
        rows[[length(rows) + 1]] <- data.frame(
          dgp    = task$dgp,
          rep    = task$rep,
          sweeps = sw$label,
          var    = j,
          width  = mm[["ma_width"]],
          cov    = mm[["ma_cov"]],
          l1     = mm[["ma_l1"]],
          stringsAsFactors = FALSE)
      }

      rm(draws); gc(verbose = FALSE)
    }

    do.call(rbind, rows)

  }, error = function(e)
      data.frame(dgp = task$dgp, rep = task$rep,
               error = conditionMessage(e), stringsAsFactors = FALSE))
  
  message(task$dgp, " rep ", task$rep,
          if (is.null(out$error)) " done in " else " FAILED after ",
          round(difftime(Sys.time(), t0, units = "mins"), 2), " min")

  if (is.null(out$error)) saveRDS(out, ckpt)
  out
}

# ==============================================================================
# run
# ==============================================================================

set.seed(8675309)
run_start <- Sys.time()

res <- future_lapply(
  tasks, run_one, s = settings,
  future.seed = TRUE, future.chunk.size = 1,
  future.packages = c("MASS", "MCMCpack", "mvtnorm", "sn", "matrixNormal",
                      "KFAS", "HDInterval", "scoringRules")
)

message("Total: ", round(difftime(Sys.time(), run_start, units = "mins"), 2),
        " min")

ok <- vapply(res, function(z) is.null(z$error), logical(1))
if (any(!ok)) message("FAILED: ", sum(!ok), " of ", length(ok))

R <- do.call(rbind, res[ok])
R$sweeps <- factor(R$sweeps, levels = vapply(SWEEPS, `[[`, "", "label"))
saveRDS(R, sprintf("%s.rds", RUN_TAG))

# ==============================================================================
# read-out
# ==============================================================================

agg <- aggregate(cbind(width, cov, l1) ~ dgp + sweeps, data = R, FUN = mean)
agg <- agg[order(agg$dgp, agg$sweeps), ]

# percentage against the smallest budget, within each DGP
agg$pct_vs_10k <- round(100 * (agg$width /
                                 ave(agg$width, agg$dgp, FUN = function(z) z[1]) - 1), 1)
agg$width <- round(agg$width, 4)
agg$cov   <- round(agg$cov,   3)
agg$l1    <- round(agg$l1,    4)

cat("\n--- T = ", TTRAIN, ", ", NREP, " replicates per DGP ---\n", sep = "")
cat("nominal coverage 0.95\n\n")
print(agg, row.names = FALSE)

cat("\nRead it this way:\n")
cat("  width plateaus and coverage reaches ~0.95  -> computational, not a\n")
cat("    property of the method; A.1 gets rewritten or cut.\n")
cat("  width plateaus and coverage stalls lower   -> the pseudo-posterior is\n")
cat("    genuinely narrow; A.1 stands with corrected numbers.\n")
cat("  width still climbing at 160k               -> inconclusive, and every\n")
cat("    coverage figure so far is a lower bound on an unknown.\n")

cat("\nl1 should be flat across budgets. A posterior mean is unbiased under\n")
cat("autocorrelation, only imprecise, so movement there means something else\n")
cat("is wrong.\n")
