# Does ECDF of simulated VARMA data match true stationary distribution?

# ==============================================================================
# load resources
# ==============================================================================

library(MASS)
library(LaplacesDemon)
library(MTS)

# ==============================================================================
# exercise settings
# ==============================================================================

n = 2
p = 2
q = 2
T = 5000

set.seed(8675309)

# ==============================================================================
# simulate parameter values
# ==============================================================================

m = numeric(n)#rnorm(n)
AR = simulate_nonexplosive_var_params(n, p, numeric(n*n*p), diag(n*n*p))
MA = matrix(rnorm(n*n*q), n, n*q)
S = rinvwishart(n + 1, diag(n))

# ==============================================================================
# calculate stationary distribution
# ==============================================================================

VARMAcov = VARMAcov(Phi = AR, Theta = MA, Sigma = S, lag = 12, trun = 120)$autocov[1:n, 1:n]

# ==============================================================================
# simulate
# ==============================================================================

Y = VARMAsim(T, arlags = p, malags = q, cnst = m, phi = AR, theta = MA, sigma = S)$series

# ==============================================================================
# compare simulation to stationary distribution
# ==============================================================================

par(mfrow = c(n, 1), mar = c(2, 2, 2, 2))

for(i in 1:n){
  plot(ecdf(Y[, i]), do.points = F, main = paste("Variable", i))
  x_grid = seq(min(Y[, i]), max(Y[, i]), length.out = 500)
  y_vals = pnorm(x_grid, mean = 0, sd = sqrt(VARMAcov[i, i]))
  lines(x_grid, y_vals, col = "red")
}

array(c(VARMAcov, cov(Y)), dim = c(n, n, 2))