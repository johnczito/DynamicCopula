# Does simulated DGFC data match true stationary distribution?

# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

# ==============================================================================
# exercise settings
# ==============================================================================

Finv1 <- function(x){
  qnorm(x, mean = 0, sd = 1)
}

F1 <- function(x){
  pnorm(x, mean = 0, sd = 1)
}

Finv2 <- function(x){
  qpois(x, lambda = 10)
}

F2 <- function(x){
  ppois(x, lambda = 10)
}

Finv3 <- function(x){
  qst(x, alpha = 2, nu = 3)
}

F3 <- function(x){
  pst(x, alpha = 2, nu = 3)
}

Finv4 <- function(x){
  qgamma(x, shape = 1, rate = 1)
}

F4 <- function(x){
  pgamma(x, shape = 1, rate = 1)
}

Finv <- c(Finv1, Finv2, Finv3, Finv4)
F <- c(F1, F2, F3, F4)

n = length(Finv)
p = 3
k = ceiling(0.7 * n)
T = 25000

#set.seed(8675309)

# ==============================================================================
# simulate parameter values
# ==============================================================================

m = rinvgamma(n, 1, 1)
M = diag(m)
L = matrix(rnorm(n * k), n, k)
AR = simulate_nonexplosive_var_params(k, p, numeric(k*k*p), 0.5 * diag(k*k*p))
G = array(c(AR), dim = c(k, k, p))
S = riwish(k + 1, diag(k))

# ==============================================================================
# simulate
# ==============================================================================

Y = simulate_dgfc(T, Finv, G, S, L, M)

# ==============================================================================
# compare simulation to stationary distribution
# ==============================================================================

# CDF plots

par(mfrow = c(n, 1), mar = c(2, 2, 2, 2))

for(i in 1:n){
  plot(ecdf(Y[, i]), main = paste("Variable", i), do.points = FALSE, col.01line = NULL)
  x_grid = seq(min(Y[, i]), max(Y[, i]), length.out = 500)
  y_vals = F[[i]](x_grid)
  lines(x_grid, y_vals, col = "red")
}

# Q-Q plots

par(mfrow = c(n, 1), mar = c(2, 2, 2, 2))

a = seq(0.005, 0.995, length.out = 250)

for(i in 1:n){
  samp_q = quantile(Y[, i], probs = a)
  true_q = Finv[[i]](a)
  plot(samp_q, true_q, pch = 19, cex = 0.1, main = paste("Variable", i))
  abline(a = 0, b = 1, col = "red", lty = 2)
}

# (abs) eigenvalues of companion matrix

abs(eigen(companion(AR))$value)
