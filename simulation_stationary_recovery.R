# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")
# Hello World?
# ==============================================================================
# exercise settings
# ==============================================================================

n = 2
sample.sizes <- 25 * 2 ^ (0:5)
Tmax <- max(sample.sizes)
Nsamp <- length(sample.sizes)

#set.seed(8675309)

# ==============================================================================
# Pick DGP, simulate fake data, and store true stationary distributions
# ==============================================================================

# VARMA, VARCH, DGFC
dgp <- "DGFC"

if(dgp == "VARMA"){
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # generate random parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  p = 3
  q = 6
  m = rnorm(n)
  AR = simulate_nonexplosive_var_params(n, p, numeric(n*n*p), 0.1*diag(n*n*p))
  ARarray = array(c(AR), dim = c(n, n, p))
  MAarray = array(rnorm(n*n*q), dim = c(n, n, q))
  S = riwish(n + 1, diag(n))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulate fake data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Y = simulate_stationary_varma(Tmax, m, ARarray, MAarray, S)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # get stationary distribution
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  VARMAparams = varma_stationary_params(m, ARarray, MAarray, S)
  
  F <- vector(mode = 'list', length = n)
  
  for(i in 1:n){
    F[[i]] = function(x){pnorm(x, mean = VARMAparams$VARMAmean[i], sd = sqrt(VARMAparams$VARMAcov[i, i]))}
  }
  
}else if(dgp == "VARCH"){
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # generate random parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  v = n + 2 + rexp(1)
  A = riwish(n + 1, diag(n))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulate fake data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Y = simulate_varch(Tmax, v, A)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # get stationary distribution
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  VARCHparams = varch_stationary_params(v, A)
  
  F <- vector(mode = 'list', length = n)
  
  for(i in 1:n){
    F[[i]] = function(x){pt(x / sqrt(VARCHparams$VARCHscale[i, i]), df = VARCHparams$VARCHdf)}
  }
  
}else if(dgp == "DGFC"){
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # generate random parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Finv1 <- function(x){qnorm(x, mean = 0, sd = 1)}
  Finv2 <- function(x){qgamma(x, shape = 1, rate = 1)}
  Finv3 <- function(x){qst(x, alpha = 2, nu = 3)}
  Finv4 <- function(x){qpois(x, lambda = 5)}
  
  Finv <- c(Finv1, Finv2, Finv3)#c(Finv1, Finv2, Finv3, Finv4)
  
  n = length(Finv)
  p = 1
  k = n#ceiling(0.7 * n)
  m = rinvgamma(n, 1, 1)
  M = diag(n)#diag(m)
  L = rbind(diag(k), matrix(0, n - k, k))#matrix(rnorm(n * k), n, k)
  AR = 0.5 * diag(k)#simulate_nonexplosive_var_params(k, p, numeric(k*k*p), 0.1 * diag(k*k*p))
  G = array(c(AR), dim = c(k, k, p))
  S = diag(k)#riwish(k + 1, diag(k))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulate fake data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Y = simulate_dgfc(Tmax, Finv, G, S, L, M)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # get stationary distribution
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  F1 <- function(x){pnorm(x, mean = 0, sd = 1)}
  F2 <- function(x){pgamma(x, shape = 1, rate = 1)}
  F3 <- function(x){pst(x, alpha = 2, nu = 3)}
  F4 <- function(x){ppois(x, lambda = 5)}
  
  F <- c(F1, F2, F3)#c(F1, F2, F3, F4)
  
}

T <- nrow(Y)
n <- ncol(Y)

# ==============================================================================
# run MCMC on an expanding window of data
# ==============================================================================

MAdraws <- vector(mode = 'list', length = Nsamp)

for(l in 1:Nsamp){
  T <- sample.sizes[l]
  MAdraws[[l]] <- DGFC.mcmc(Y[1:T, ], ndraw = 1000, thin = 10)$ma
  message(paste("Stage ", l, " of ", Nsamp, " done!", sep = ""))
}

# ==============================================================================
# plots true stationary CDF, ECDF, and MA credible bands for each sample size
# ==============================================================================

i = 2
a = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75)
x_grid = seq(min(Y[, i]), max(Y[, i]), length.out = 500)
y_vals = F[[i]](x_grid)

#par(mfrow = c(Nsamp, 1), mar = c(2, 2, 2, 2))
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))

for(l in Nsamp:1){
  T <- sample.sizes[l]
  plot(ecdf(Y[1:T, i]), main = paste("Variable", i, "; T =", sample.sizes[l]), 
       do.points = FALSE, col.01line = NULL, xlim = c(min(Y[, i]), max(Y[, i])))
  plot_ma_band(MAdraws[[l]][[i]][, 1, 1], MAdraws[[l]][[i]][, 2, ], a, rgb(0, 0, 1, 0.1))
  lines(x_grid, y_vals, col = "red")
}

# how persistent is the time series process?

abs(eigen(companion(AR))$values)

# ADD TRACE PLOTS!

plot(1:1000, MAdraws[[Nsamp]][[i]][4, 2, ], type ="l")
