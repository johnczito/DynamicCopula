# ==============================================================================
# load resources
# ==============================================================================

library(MASS)
library(sn)

source("_helpers/_helpers.R")

# ==============================================================================
# marginal CDFs
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

Finv <- c(Finv1, Finv2, Finv3)
F <- c(F1, F2, F3)

# ==============================================================================
# exercise settings
# ==============================================================================

n <- length(Finv)
k <- 2#floor(0.7 * n)
T <- 100000

# ==============================================================================
# simulate fake data
# ==============================================================================

set.seed(8675309)

m = rinvgamma(n, 1, 1)
M = diag(m)
L = matrix(rnorm(n * k), n, k)
S = diag(k)
Gmat = simulate_nonexplosive_var_params(k, 1, c(diag(k)), diag(k * k))
G = array(c(Gmat), dim = c(k, k, 1))

Y = simulate_dgfc(T, Finv, G, S, L, M)

# ==============================================================================
# plot ECDF against margin adjustment
# ==============================================================================

par(mfrow = c(n, 1), mar = c(2, 2, 2, 2))

for(i in 1:n){
  plot(ecdf(Y[, i]), cex = 0.1, col = "black")
  x_vals <- seq(min(Y[, i]) - 2, max(Y[, i]) + 2, length.out = 500)
  y_vals <- F[[i]](x_vals)
  lines(x_vals, y_vals, col = "red")
}
