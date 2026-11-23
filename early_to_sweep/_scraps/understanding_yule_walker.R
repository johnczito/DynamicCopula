# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

set.seed(8675309)

n <- 2
p <- 1
T <- 3

G <- simulate_nonexplosive_var_params(n, p, numeric(n*n*p), diag(n*n*p))
S <- riwish(n + 1, diag(n))

S0 = matrix(solve(diag(n*n) - G %x% G, c(S)), n, n)
S1 = G %*% S0
O = cbind(rbind(S0, t(S1)), rbind(S1, S0))
C = cov2cor(O)
C0 = C[1:n, 1:n]
C1 = C[(n+1):(2*n), 1:n]

G2 = C1 %*% solve(C0)
S2 = matrix((diag(n*n) - G2 %x% G2) %*% c(C0), n, n)
