# Simulate Theorem 1 in Feldman, Kowal (2023) when the Z are correlated

# ==============================================================================
# load resources
# ==============================================================================

library(MASS)
library(MCMCpack)

source("_helpers/_helpers.R")

# ==============================================================================
# Strictly increasing transformation
# ==============================================================================

h <- function(x){return(exp(x))} 

# ==============================================================================
# Simulation settings
# ==============================================================================

sample.sizes <- 10 * 2 ^ (0:6)
M <- length(sample.sizes)
G <- 500
dependent <- TRUE

set.seed(8675309)

# ==============================================================================
# True CDF of Y = h(Z)
# ==============================================================================

x_grid <- seq(0, 5, length.out = G)
y_grid <- plnorm(x_grid)

# ==============================================================================
# Simulate master Z sample to chew over
# ==============================================================================

N <- max(sample.sizes)

if(dependent == TRUE){
  H <- riwish(N + 1, diag(N))
  R <- diag(1/(sqrt(diag(H))))
  C <- R %*% H %*% R
  #C <- 0.5 * matrix(1, N, N) + (1 - 0.5) * diag(N)
}else{
  C <- diag(N)
}

Z <- mvrnorm(1, numeric(N), C)

# ==============================================================================
# GO!
# ==============================================================================

par(mfrow = c(M, 1), mar = c(2, 2, 2, 2))

for(m in 1:M){
  n <- sample.sizes[m]
  z <- Z[1:n]
  y <- h(z)
  F <- sapply(x_grid, margin.adjustment, z, y, 0, 1)
  plot(x_grid, y_grid, type = "l", col = "red", main = paste("Sample size =", n))
  lines(x_grid, F)
}