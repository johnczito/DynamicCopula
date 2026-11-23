else if(dgp == "DGFC"){
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # generate random parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Finv1 <- function(x){qnorm(x, mean = 0, sd = 1)}
  Finv2 <- function(x){qgamma(x, shape = 1, rate = 1)}
  Finv3 <- function(x){qst(x, alpha = 2, nu = 3)}
  Finv4 <- function(x){qpois(x, lambda = 5)}
  
  Finv <- c(Finv3, Finv4)#c(Finv1, Finv2, Finv3, Finv4)
  
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
  
  F <- c(F3, F4)#c(F1, F2, F3, F4)
  
  band_col = rgb(1, 0.6, 0, 0.1)
  
}