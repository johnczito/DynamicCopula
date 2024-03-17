# Are the MNIW MvReg calculations correct?

# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

# ==============================================================================
# dimensions
# ==============================================================================

n = 3
p = 1
T = 20
inclconst = TRUE

# ==============================================================================
# prior hyperparameters
# ==============================================================================

v0 = n + 2 * rexp(1)
P0 = riwish(n + 1, diag(n))
B0 = matrix(rnorm(n*n*p + n*inclconst), n*p + inclconst, n)
invO0 = rwish(n*p + inclconst + 1, diag(n*p + inclconst))

# ==============================================================================
# fake data
# ==============================================================================

Y = matrix(rnorm(T*n), T, n)
X = varp_design_matrix(Y, p, inclconst)
Y = Y[(p + 1):T, ]

# ==============================================================================
# posterior hyperparameters
# ==============================================================================

params = mvregposterior(Y, X, v0, P0, B0, invO0)
v = params$v
P = params$P
B = params$B
invO = params$invO

# ==============================================================================
# posterior hyperparameters
# ==============================================================================

Btest = matrix(rnorm(n*n*p + n*inclconst), n*p + inclconst, n)
Stest = riwish(n + 1, diag(n))


out1 = mvregloglikelihood(Y, X, Btest, Stest) + dmniw(Btest, Stest, v0, P0, B0, invO0) - dmniw(Btest, Stest, v, P, B, invO)
out2 = mvreglogmdd(n, v0, v, P0, P, invO0, invO)

c(out1, out2)
