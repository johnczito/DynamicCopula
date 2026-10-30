# ==============================================================================
# DGLM diagnostic
#
# The Poisson and negative binomial competitors are not MCMC. from_brian.R
# fits them by maximum likelihood (deterministic) and then draws from the
# predictive distribution by importance sampling, resampling nsim times with
# replacement from nsim weighted signal paths.
#
# If the importance weights are uneven, those draws come from far fewer
# distinct paths than nsim suggests, and different seeds land on different
# subsets. This script measures that.
#
# Run from the repo root. Takes a couple of minutes.
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

set.seed(8675309)

crime_data <- read.csv("datasets/full_nsw_crime_dataset.csv", header = TRUE)

var_names <- c("murder", "attempted", "accessory_murder", "manslaughter",
               "assault_nondomestic", "assault_police", "abduction",
               "blackmail", "Escape_custody", "resist_officer")

Y_full <- as.matrix(crime_data[, var_names])

# from_brian.R calls importMvSim(object, nsim / 4, H) with nsim = 1000,
# so importanceSSM gets 250 base draws and antithetics expand it to 1000.
NSIM_BASE  <- 250
NSIM_TOTAL <- 4 * NSIM_BASE

# ------------------------------------------------------------------------------

diagnose <- function(Y, family) {
  
  mod <- SSModel(Y ~ -1 + SSMtrend(degree = 1, Q = matrix(NA), type = "distinct"),
                 distribution = family)
  
  object <- fitSSM(mod, inits = c(diag(ncol(Y))))$model
  
  imp <- KFAS::importanceSSM(object, "signal",
                             nsim = NSIM_BASE, antithetics = TRUE)
  
  w <- imp$weights / sum(imp$weights)
  
  ess <- 1 / sum(w^2)
  top10 <- sum(sort(w, decreasing = TRUE)[1:10])
  maxw <- max(w)
  
  # the resampling step, as importMvSim does it
  idx <- sample(seq_along(w), size = length(w), replace = TRUE, prob = w)
  distinct <- length(unique(idx))
  
  list(n_draws = length(w),
       ess = ess,
       top10 = top10,
       maxw = maxw,
       distinct = distinct,
       u = unique(as.vector(object$u)))
}

# ------------------------------------------------------------------------------

windows  <- c(100, 200, nrow(Y_full))
families <- c("poisson", "negative binomial")

cat("\nDGLM importance weight diagnostic\n")
cat("nsim =", NSIM_TOTAL, "(", NSIM_BASE, "base draws x 4 antithetics )\n\n")
cat(sprintf("%-20s %5s %8s %8s %8s %10s\n",
            "family", "T", "ESS", "top10 w", "max w", "distinct"))
cat(strrep("-", 64), "\n")

u_report <- list()

for (t in windows) {
  for (fam in families) {
    d <- diagnose(Y_full[1:t, ], fam)
    cat(sprintf("%-20s %5d %8.1f %8.3f %8.4f %10d\n",
                fam, t, d$ess, d$top10, d$maxw, d$distinct))
    u_report[[paste(fam, t)]] <- d$u
  }
}

cat("\n")
cat("Dispersion / offset parameter u (KFAS default is 1):\n")
for (nm in names(u_report)) {
  v <- u_report[[nm]]
  cat(sprintf("  %-26s %s\n", nm,
              if (length(v) == 1) as.character(v)
              else paste0(length(v), " distinct values, range ",
                          signif(min(v), 4), " to ", signif(max(v), 4))))
}

cat("\n")
cat("Reading it:\n")
cat("  ESS near", NSIM_TOTAL, "= weights are close to uniform, importance\n")
cat("    sampling is not the source of the instability.\n")
cat("  ESS in the tens or low hundreds = the predictive draws rest on few\n")
cat("    distinct signal paths, which would explain run-to-run movement.\n")
cat("  'distinct' is how many unique paths the resampler actually picked.\n")
cat("  u fixed at 1 for the negative binomial means dispersion is not being\n")
cat("    estimated, which would widen and destabilise that competitor on its own.\n\n")
