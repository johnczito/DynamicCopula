# ==============================================================================
# Replaces Section 6.2
#
# ECDF against margin adjustment across sample sizes, on the three DGPs from
# Section 6.1. Both arms come off the SAME posterior fit with the RNG state
# restored between forecast calls, so the latent forward simulations are
# identical and only the margin treatment differs. Paired, not noisy
#
# Differences from the old simulation_wrongmodel_forecasting.R:
#   - replicates (its nrep = 10 was declared and never used)
#   - burn = 5000 as in the consistency script, not burn = 0
#   - fixed sample-size grid matching 6.1 instead of averaging over an
#     expanding window
#   - margin accuracy scored against the true F, which the closed-form
#     stationary distributions make possible
#   - no competitor models; this section is about the margins
#
# Requires the use_spline vectorization in DGFC.forecast
#
# Run from the repo root
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")
source("_helpers/_sim_dgp.R")
source("_helpers/_sim_metrics.R")

library(future)
library(future.apply)

REPO <- getwd()
plan(multisession, workers = max(1, availableCores() - 1))

# ==============================================================================
# settings
#
# T_GRID and the MCMC settings come from simulation_wrongmodel_consistency.R,
# NOT from the old forecasting script, which ran burn = 0
# ==============================================================================

# MODE: "full" is the production run, "smoke" a two-minute shakeout, and
# "longchain" a single-cell check on whether the margin-adjustment bias and
# band undercoverage at T = 200 survive a four-times-longer chain. If they
# hold, they are properties of the pseudo-posterior; if they move, the full
# run needs longer chains before anything gets written up

MODE <- "full"

DGPS   <- c("varma", "varch", "varma_copula")
T_GRID <- 25 * 2^(0:3)   # 25 50 100 200, as in 6.1
H      <- 10
NREP   <- 100
BURNIN <- 200
NDRAW  <- 1000
NBURN  <- 5000
NTHIN  <- 5
ALPH   <- 0.05

if (MODE == "smoke") {
  
  T_GRID <- c(25, 50)
  H <- 3
  NREP <- 2
  NDRAW <- 200
  NBURN <- 200
  NTHIN <- 1
  
} else if (MODE == "longchain") {
  
  DGPS <- c("varch", "varma")
  T_GRID <- 200
  H <- 1 # forecasts are not what this check is about
  NREP <- 20
  NDRAW <- 4000 # 40,000 sweeps against the full run's 10,000
  NBURN <- 20000
  NTHIN <- 5
  
}

# k.star is left at its default, ceiling(0.7 * n) = 2 at n = 2, matching 6.1,
# which also omits it

# The full-run tag is left exactly as it was so the existing checkpoint
# directory and rds still match. Other modes carry the draw count, since
# chain length is the only thing that distinguishes this check
RUN_TAG <- if (MODE == "full") {
  sprintf("margincmp_R%d_H%d_T%d", NREP, H, max(T_GRID))
} else {
  sprintf("margincmp_%s_R%d_H%d_T%d_d%d", toupper(MODE),
          NREP, H, max(T_GRID), NDRAW)
}
CKPT_DIR <- sprintf("_ckpt_%s", RUN_TAG)

settings <- list(repo = REPO, T_grid = T_GRID, H = H, burnin = BURNIN,
                 ndraw = NDRAW, nburn = NBURN, nthin = NTHIN, alph = ALPH,
                 ckpt_dir = CKPT_DIR)

# ------------------------------------------------------------------------------
# stamp the config so a settings change can never silently resume stale results
# ------------------------------------------------------------------------------

dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)
cfg_file <- file.path(CKPT_DIR, "_config.rds")
cfg <- settings[setdiff(names(settings), c("repo", "ckpt_dir"))]
if (file.exists(cfg_file)) {
  if (!identical(readRDS(cfg_file), cfg))
    stop("checkpoint dir ", CKPT_DIR, " holds a different configuration; ",
         "delete it or change RUN_TAG")
} else {
  saveRDS(cfg, cfg_file)
}

tasks <- expand.grid(dgp = DGPS, rep = seq_len(NREP), stringsAsFactors = FALSE)
tasks <- split(tasks, seq_len(nrow(tasks)))

cat("run tag: ", RUN_TAG, "\n", sep = "")
cat("tasks: ", length(tasks),
    "  T grid: ", paste(T_GRID, collapse = " "),
    "  H: ", H, "\n\n", sep = "")

# ==============================================================================
# one DGP, one replicate, all sample sizes
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

  ckpt <- file.path(s$ckpt_dir, sprintf("%s_%04d.rds", task$dgp, task$rep))
  if (file.exists(ckpt)) return(readRDS(ckpt))

  t0 <- Sys.time()

  out <- tryCatch({

    d <- make_dgp(task$dgp, burnin = s$burnin)
    Y <- d$sim(max(s$T_grid) + s$H)

    fc_rows <- list()
    ma_rows <- list()

    for (TT in s$T_grid) {

      Ytr <- Y[1:TT, , drop = FALSE]
      Yfu <- Y[(TT + 1):(TT + s$H), , drop = FALSE]

      draws <- DGFC.mcmc(Ytr, ndraw = s$ndraw, burn = s$nburn, thin = s$nthin)

      # Spline on continuous variables only. The VARMA Copula has a Poisson
      # variable and a skew-t variable in one series, so this has to vary by
      # variable; splining a step function is wrong
      sp <- !d$discrete

      rng <- if (exists(".Random.seed", envir = .GlobalEnv))
               get(".Random.seed", envir = .GlobalEnv) else NULL
      fc_ma <- DGFC.forecast(s$H, draws, use_spline = sp)
      if (!is.null(rng)) assign(".Random.seed", rng, envir = .GlobalEnv)
      fc_ec <- DGFC.forecast(s$H, draws, use_spline = sp,
                             use_ecdf = TRUE, Y = Ytr)

      oor <- out_of_range(Ytr, Yfu)

      for (j in seq_len(d$n)) {

        scorer <- if (d$discrete[j]) discrete_forecast_metrics
                  else               continuous_forecast_metrics

        for (h in seq_len(s$H)) {
          for (arm in c("ma", "ecdf")) {
            fcj <- if (arm == "ma") fc_ma[h, j, ] else fc_ec[h, j, ]
            mm  <- scorer(Yfu[h, j], s$alph, fcj)
            fc_rows[[length(fc_rows) + 1]] <- data.frame(
              dgp = task$dgp, rep = task$rep, T = TT, var = j, h = h,
              arm = arm,
              point = mm[1], cov = mm[2], size = mm[3], crps = mm[4],
              oor = oor[j], rho = d$rho, discrete = d$discrete[j],
              stringsAsFactors = FALSE)
          }
        }

        ma_rows[[length(ma_rows) + 1]] <- data.frame(
          dgp = task$dgp, rep = task$rep, T = TT, var = j, arm = "ma",
          as.list(margin_metrics(draws$ma[[j]], d$F[[j]], 1 - s$alph)),
          stringsAsFactors = FALSE)
        ma_rows[[length(ma_rows) + 1]] <- data.frame(
          dgp = task$dgp, rep = task$rep, T = TT, var = j, arm = "ecdf",
          as.list(ecdf_margin_metrics(Ytr[, j], d$F[[j]])),
          stringsAsFactors = FALSE)
      }
    }

    list(fc = do.call(rbind, fc_rows), ma = do.call(rbind, ma_rows),
         dgp = task$dgp, rep = task$rep)

  }, error = function(e) list(error = conditionMessage(e),
                              dgp = task$dgp, rep = task$rep))

  message(task$dgp, " rep ", task$rep,
          if (is.null(out$error)) " done in " else " FAILED after ",
          round(difftime(Sys.time(), t0, units = "mins"), 2), " min",
          if (!is.null(out$error)) paste0(": ", out$error) else "")

  if (is.null(out$error)) saveRDS(out, ckpt) # never checkpoint an error
  out
}

# ==============================================================================
# run
# ==============================================================================

set.seed(8675309)
run_start <- Sys.time()

results <- future_lapply(
  tasks, run_one, s = settings,
  future.seed = TRUE, future.chunk.size = 1,
  future.packages = c("MASS", "MCMCpack", "mvtnorm", "sn", "matrixNormal",
                      "KFAS", "HDInterval", "scoringRules")
)

message("Total: ", round(difftime(Sys.time(), run_start, units = "mins"), 2),
        " min")

ok <- vapply(results, function(r) is.null(r$error), logical(1))
if (any(!ok)) {
  message("FAILED: ", sum(!ok), " of ", length(ok))
  message("First error: ", results[[which(!ok)[1]]]$error)
}

fc <- do.call(rbind, lapply(results[ok], `[[`, "fc"))
ma <- do.call(rbind, lapply(results[ok], `[[`, "ma"))

saveRDS(list(fc = fc, ma = ma, settings = settings, failed = sum(!ok)),
        sprintf("%s.rds", RUN_TAG))

# ==============================================================================
# quick look, so a smoke run tells you something without a separate script
# ==============================================================================

cat("\nforecast rows:", nrow(fc), "  margin rows:", nrow(ma), "\n\n")

# ------------------------------------------------------------------------------
# Paired differences. Both arms come off the same fit, so differencing within
# (dgp, rep, T, var, h) removes the replicate-to-replicate variation that
# dominates a comparison of group means
# ------------------------------------------------------------------------------

key <- c("dgp", "rep", "T", "var", "h")
w <- merge(fc[fc$arm == "ma",   c(key, "crps", "cov", "size", "oor", "rho")],
           fc[fc$arm == "ecdf", c(key, "crps", "cov", "size")],
           by = key, suffixes = c(".ma", ".ec"))
w$d_crps_pct <- 100 * (w$crps.ec / w$crps.ma - 1)

p1 <- w[w$h == 1, ]
m1 <- aggregate(d_crps_pct ~ dgp + T + var, data = p1, FUN = mean)
s1 <- aggregate(d_crps_pct ~ dgp + T + var, data = p1,
                FUN = function(z) sd(z) / sqrt(length(z)))
names(m1)[names(m1) == "d_crps_pct"] <- "mean_pct"
names(s1)[names(s1) == "d_crps_pct"] <- "se_pct"
out1 <- merge(m1, s1)
out1$mean_pct <- round(out1$mean_pct, 1)
out1$se_pct   <- round(out1$se_pct,   1)

cat("--- h = 1 CRPS, paired difference, ecdf against margin adjustment ---\n")
cat("positive means freezing F is worse\n\n")
print(out1[order(out1$dgp, out1$var, out1$T), ], row.names = FALSE)

cat("\n--- CONVERGENCE CHECK: L1 error against the TRUE F must fall with T ---\n")
cat("(CRPS is not a convergence check here. It depends on the realized future\n")
cat(" value, so across-T differences are dominated by which window was drawn.)\n\n")
print(aggregate(ma_l1 ~ dgp + T + arm, data = ma, FUN = mean), row.names = FALSE)

cat("\n--- does the margin's 95% credible band cover the true F? ---\n")
print(aggregate(ma_cov ~ dgp + T, data = ma[ma$arm == "ma", ], FUN = mean),
      row.names = FALSE)

cat("\n--- forecast interval coverage, both arms, h = 1 ---\n")
print(aggregate(cov ~ dgp + T + arm, data = fc[fc$h == 1, ], FUN = mean),
      row.names = FALSE)

cat("\nsaved to ", RUN_TAG, ".rds\n", sep = "")

# ------------------------------------------------------------------------------
# for longchain mode
# ------------------------------------------------------------------------------

if (MODE == "longchain") {
  
  ref_file <- "margincmp_R100_H10_T200.rds"
  
  if (file.exists(ref_file)) {
    
    ref <- readRDS(ref_file)
    
    pick <- function(m) {
      z <- m[m$dgp %in% DGPS & m$T == max(T_GRID), ]
      c(ma_l1_margin = mean(z$ma_l1[z$arm == "ma"]),
        ma_l1_ecdf   = mean(z$ma_l1[z$arm == "ecdf"]),
        ma_cov       = mean(z$ma_cov[z$arm == "ma"]),
        n_rep        = length(unique(z$rep)))
    }
    
    cmp <- rbind(`full run, 10k sweeps`   = pick(ref$ma),
                 `long chain, 40k sweeps` = pick(ma))
    
    cat("\n--- longchain check: ", DGPS, " at T = ", max(T_GRID), " ---\n", sep = "")
    cat("Rows agree  -> the bias and undercoverage are real, report them.\n")
    cat("Long chain improves -> residual mixing, the full run needs longer chains.\n\n")
    print(round(cmp, 5))
    
  } else {
    cat("\n(", ref_file, " not found, cannot compare)\n", sep = "")
  }
}
