# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

# ==============================================================================
# Pick DGP and simulate fake data
# ==============================================================================

T <- 200
n <- 2

# VARMA, VARCH, DGFC
dgp <- "DGFC"

if(dgp == "VARMA"){
  
  p = 3
  q = 6
  m = rnorm(n)
  AR = simulate_nonexplosive_var_params(n, p, numeric(n*n*p), 0.5*diag(n*n*p))
  ARarray = array(c(AR), dim = c(n, n, p))
  MAarray = array(rnorm(n*n*q), dim = c(n, n, q))
  S = riwish(n + 1, diag(n))
  
  Y = simulate_stationary_varma(T, m, ARarray, MAarray, S)
  
}else if(dgp == "VARCH"){
  
  v = n + 2 + rexp(1)
  A = riwish(n + 1, diag(n))
  Y = simulate_varch(T, v, A)
  
}else if(dgp == "DGFC"){
  
  Finv1 <- function(x){qnorm(x, mean = 0, sd = 1)}
  Finv2 <- function(x){qpois(x, lambda = 5)}
  Finv3 <- function(x){qst(x, alpha = 2, nu = 3)}
  Finv4 <- function(x){qgamma(x, shape = 1, rate = 1)}
  
  Finv <- c(Finv1, Finv2, Finv3, Finv4)
  
  n = length(Finv)
  p = 3
  k = ceiling(0.7 * n)
  m = rinvgamma(n, 1, 1)
  M = diag(m)
  L = matrix(rnorm(n * k), n, k)
  AR = simulate_nonexplosive_var_params(k, p, numeric(k*k*p), 0.5 * diag(k*k*p))
  G = array(c(AR), dim = c(k, k, p))
  S = riwish(k + 1, diag(k))
  
  Y = simulate_dgfc(T, Finv, G, S, L, M)
  
}

T <- nrow(Y)
n <- ncol(Y)

# ==============================================================================
# DLM hyperparameters
# ==============================================================================

m0 = 0
C0 = 1
shape.y = 1
rate.y = 1
shape.theta = 1
rate.theta = 1

# ==============================================================================
# Sampling settings
# ==============================================================================

H = 10
t0 = 100

dlm_ndraw  = 1000
bvar_ndraw = 1000
dfc_ndraw  = 1000

# ==============================================================================
# Preallocate storage
# ==============================================================================

fcast_dlm  = array(0, c(H, n, dlm_ndraw, T))
fcast_bvar = array(0, c(H, n, bvar_ndraw, T))
fcast_dfc  = array(0, c(H, n, dfc_ndraw, T))

# ==============================================================================
# Generate forecasts
# ==============================================================================

for(t in t0:T){
  
  my_bvar_model <- bvar(Y[1:t, ], lags = 1, n_draw = bvar_ndraw, n_burn = 0, verbose = FALSE)
  fcast_bvar[, , , t] <- aperm(predict(my_bvar_model, horizon = H)$fcast, c(2, 3, 1))
  
  draws = DGFC.mcmc(Y[1:t, ], ndraw = dfc_ndraw)
  fcast_dfc[, , , t] = DGFC.forecast(H, draws)
  
  message(paste("Stage ", t, " of ", T, " done!", sep = ""))
  
}



# ==============================================================================
# Get results
# ==============================================================================

alpha = c(0.01, 0.05, 0.1, 0.25, 0.5)
