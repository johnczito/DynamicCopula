source("_helpers/_helpers.R")

#A = matrix(runif(1, -1, 1), 1, 1)
#SU = matrix(rexp(1), 1, 1)
#GS = normalize_var1_params(A, SU)

#n <- 2
#p <- 1

G <- simulate_nonexplosive_var_params(n, p, numeric(n*n*p), diag(n*n*p))
S <- riwish(n + 1, diag(n))

newstuff = normalize_var1_params(G, S)