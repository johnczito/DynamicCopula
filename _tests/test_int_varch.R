# Does simulated VARCH data match true stationary distribution?

# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

# ==============================================================================
# exercise settings
# ==============================================================================

n = 4
T = 5000

#set.seed(8675309)

# ==============================================================================
# simulate parameter values
# ==============================================================================

v = n + 2 + rexp(1)
A = riwish(n + 1, diag(n))

# ==============================================================================
# calculate stationary distribution
# ==============================================================================

VARCHparams = varch_stationary_params(v, A)
df = VARCHparams$VARCHdf
S = VARCHparams$VARCHscale

# ==============================================================================
# simulate
# ==============================================================================

Y = simulate_varch(T, v, A)

# ==============================================================================
# compare simulation to stationary distribution
# ==============================================================================

# CDF plots

par(mfrow = c(n, 1), mar = c(2, 2, 2, 2))

for(i in 1:n){
  plot(ecdf(Y[, i]), do.points = FALSE, main = paste("Variable", i))
  x_grid = seq(min(Y[, i]), max(Y[, i]), length.out = 500)
  y_vals = pt(x_grid / sqrt(S[i, i]), df = df)
  lines(x_grid, y_vals, col = "red")
}

# Q-Q plots

par(mfrow = c(n, 1), mar = c(2, 2, 2, 2))

a = seq(0.005, 0.995, length.out = 250)

for(i in 1:n){
  samp_q = quantile(Y[, i], probs = a)
  true_q = qt(a, df = df) * sqrt(S[i, i])
  plot(samp_q, true_q, pch = 19, cex = 0.1, main = paste("Variable", i))
  abline(a = 0, b = 1, col = "red", lty = 2)
}

# Compare moments

array(c((df / (df - 2)) * S, cov(Y)), dim = c(n, n, 2))