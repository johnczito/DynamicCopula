# ==============================================================================
# Scorers for the margin-comparison simulation (new Section 6.2)
#
# Place in _helpers/ and source after _helpers.R
#
# The two forecast scorers are scalar-at-a-time versions of the paper's
# existing functions in _helpers/forecasting.R, and use identical formulas:
#
#   continuous_forecast_metrics <-> forecast_results
#     MSE       (mean(draws) - obs)^2
#     INT-COV   hdi(draws, 1 - alph)
#     INT-SIZE  hdi upper - hdi lower
#     CRPS      scoringRules::crps_sample
#
#   discrete_forecast_metrics <-> forecast_results_discrete
#     MAE       abs(median(draws) - obs)
#     INT-COV   equal-tailed quantiles, not HDI
#     INT-SIZE  upper - lower
#     CRPS      scoringRules::crps_sample
#
# Scalar versions are used rather than the array functions because this
# simulation produces one forecast per (replicate, T, horizon) rather than the
# T x H x M panel that forecast_results expects, so get_variable_fcasts and its
# t - h alignment do not apply
#
# I chose not to compute INT-SCORE, as it is the only metric needing
# interval_score(), which requires scoringutils pinned to 1.2.2
# ==============================================================================

continuous_forecast_metrics <- function(obs, alph, fcast_draws){
  fcast_mean <- mean(fcast_draws)
  fcast_hdi  <- hdi(fcast_draws, 1 - alph)
  c((fcast_mean - obs)^2,
    (fcast_hdi[1] <= obs) & (obs <= fcast_hdi[2]),
    fcast_hdi[2] - fcast_hdi[1],
    scoringRules::crps_sample(obs, fcast_draws))
}

discrete_forecast_metrics <- function(obs, alph, fcast_draws){
  fcast_med <- median(fcast_draws)
  q <- quantile(fcast_draws, c(alph / 2, 1 - alph / 2), names = FALSE)
  c(abs(fcast_med - obs),
    (q[1] <= obs) & (obs <= q[2]),
    q[2] - q[1],
    scoringRules::crps_sample(obs, fcast_draws))
}

# ==============================================================================
# Margin accuracy against the true F
#
# This is what the referee comment is actually about, and it is only possible
# because every DGP has closed-form stationary margins
#
# draws$ma[[j]] is [length(unique(Y[,j])), 2, ndraw]. The x grid in [, 1, ] is
# sort(unique(Y[, j])) and is identical across draws, which is also the grid the 
# use_ecdf path builds, so the two arms are compared on the same points
# ==============================================================================

margin_metrics <- function(ma_j, Ftrue, level = 0.95) {
  x  <- ma_j[, 1, 1]
  Fh <- ma_j[, 2, ] # npoints x ndraw
  Ft <- Ftrue(x)

  Fbar <- rowMeans(Fh)
  a    <- (1 - level) / 2
  lo   <- apply(Fh, 1, quantile, probs = a,     names = FALSE)
  hi   <- apply(Fh, 1, quantile, probs = 1 - a, names = FALSE)

  c(ma_sup   = max(abs(Fbar - Ft)), # point accuracy of the posterior mean
    ma_l1    = mean(abs(Fbar - Ft)),
    ma_cov   = mean(Ft >= lo & Ft <= hi), # does the band cover the truth
    ma_width = mean(hi - lo))
}

# The frozen arm has no uncertainty at all. Coverage is undefined and width is 
# identically zero. Reporting that is more informative than leaving the row out
ecdf_margin_metrics <- function(y_j, Ftrue) {
  x  <- sort(unique(y_j))
  Fh <- ecdf(y_j)(x)
  Ft <- Ftrue(x)
  c(ma_sup   = max(abs(Fh - Ft)),
    ma_l1    = mean(abs(Fh - Ft)),
    ma_cov   = NA_real_,
    ma_width = 0)
}

# ==============================================================================
# The mechanism
#
# An empirical CDF can place no probability outside the observed range, so this
# is the proximate cause of every ECDF failure seen in the real data.
# ==============================================================================

out_of_range <- function(y_train, y_future) {
  lo <- apply(y_train, 2, min)
  hi <- apply(y_train, 2, max)
  vapply(seq_len(ncol(y_future)),
         function(j) mean(y_future[, j] < lo[j] | y_future[, j] > hi[j]), 0)
}
