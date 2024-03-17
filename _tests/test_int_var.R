# Does simulated VAR data match true stationary distribution?

# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

# ==============================================================================
# exercise settings
# ==============================================================================

n = 2
p = 4
T = 25000

#set.seed(8675309)

# ==============================================================================
# simulate parameter values
# ==============================================================================

m = rnorm(n)
AR = simulate_nonexplosive_var_params(n, p, numeric(n*n*p), 0.5 * diag(n*n*p))
ARarray = array(c(AR), dim = c(n, n, p))
S = riwish(n + 1, diag(n))

# ==============================================================================
# calculate stationary distribution
# ==============================================================================

VARparams = var_stationary_params(m, ARarray, S)

# ==============================================================================
# simulate
# ==============================================================================

Y = simulate_stationary_var(T, m, ARarray, S)

#Y = simulate_stationary_var2(T, m, ARarray[, , 1], ARarray[, , 2], S)
#Y = simulate_stationary_var1(T, m, AR, S)

# ==============================================================================
# compare simulation to stationary distribution
# ==============================================================================

# CDF plots

par(mfrow = c(n, 1), mar = c(2, 2, 2, 2))

for(i in 1:n){
  plot(ecdf(Y[, i]), main = paste("Variable", i), do.points = FALSE, col.01line = NULL)
  x_grid = seq(min(Y[, i]), max(Y[, i]), length.out = 500)
  y_vals = pnorm(x_grid, mean = VARparams$VARmean[i], sd = sqrt(VARparams$VARcov[i, i]))
  lines(x_grid, y_vals, col = "red")
}

# Q-Q plots

par(mfrow = c(n, 1), mar = c(2, 2, 2, 2))

a = seq(0.005, 0.995, length.out = 250)

for(i in 1:n){
  samp_q = quantile(Y[, i], probs = a)
  true_q = qnorm(a, mean = VARparams$VARmean[i], sd = sqrt(VARparams$VARcov[i, i]))
  plot(samp_q, true_q, pch = 19, cex = 0.1, main = paste("Variable", i))
  abline(a = 0, b = 1, col = "red", lty = 2)
}

# Compare moments

array(c(VARparams$VARcov, cov(Y)), dim = c(n, n, 2))
cbind(VARparams$VARmean, colMeans(Y))

# (abs) eigenvalues of companion matrix

abs(eigen(companion(AR))$value)