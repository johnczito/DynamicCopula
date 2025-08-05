# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

# ==============================================================================
# exercise settings
# ==============================================================================

n = 2
p = 1

Tmax <- 200


set.seed(8675309)

# ==============================================================================
# Ground truth params
# ==============================================================================

Finv1 <- function(x){qpois(x, lambda = 1)}
Finv2 <- function(x){qst(x, alpha = 2, nu = 1)}

Finv <- c(Finv1, Finv2)

F1 <- function(x){ppois(x, lambda = 1)}
F2 <- function(x){pst(x, alpha = 2, nu = 1)}

F <- c(F1, F2)


A <- simulate_nonexplosive_var_params(n, p, numeric(n*n*p), 0.1 * diag(n*n*p))
SU <- riwish(n + 1, diag(n))
GS = normalize_var1_params(A, SU)
G = GS$G
S = GS$S
O0 = stationary_var1_covariance(G, S)

abs(eigen(G)$values)

# ==============================================================================
# Simulate fake data
# ==============================================================================

Z = simulate_stationary_var1(Tmax, numeric(n), G, S)

Y = matrix(0, Tmax, n)
for(i in 1:n){
  Y[, i] = Finv[[i]](pnorm(Z[, i], mean = 0, sd = 1))
}

par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))

plot(1:Tmax, Y[, 1], type = "l", main = "Poisson(1) series", ylab = "", xlab = "t")
plot(1:Tmax, Y[, 2], type = "l", main = "Skew-t(2, 1) series", ylab = "", xlab = "t")



