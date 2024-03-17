isexplosive <- function(A){
# VAR(1): y[t] = A * y[t - 1] + e[t]
# Does A have any eigenvalues outside the unit circle?
#
# A: n x n matrix
  
  return(all(abs(eigen(A)$values) < 1) == FALSE)
}

companion <- function(VARparams){
# Calculate the companion form matrix.
#
# VARparams: n x n*p matrix
  
  n = nrow(VARparams)
  p = ncol(VARparams) / n
  return( rbind(VARparams, cbind(diag(n*(p-1)), matrix(0, n*(p-1), n))) )
}

varp2var1 <- function(VARparams, S){
# Calculate the parameters of the VAR(1) representation of an n-dimensional VAR(p)
# Source: Lutkepohl (2005) pg. 15
#
# AR: n x n*p matrix
#  S: n x n covariance matrix
  n = nrow(VARparams)
  np = ncol(VARparams)
  A = companion(VARparams)
  SU = matrix(0, np, np)
  SU[1:n, 1:n] = S
  return(list(A = A, SU = SU))
}

stationary_var1_covariance <- function(A, SU){
# Calculate marginal covariance of Y[t] in stationary VAR(1): 
# Y[t] = A * Y[t - 1] + E[t], E[t] ~ N(0, SU)
# Source: Lutkepohl (2005) Section 2.1.4 (pg. 29)
# Note: matrix inverse should fail if A has eigenvalues outside the unit circle
#
#  A: n x n matrix
# SU: n x n covariance matrix
  
  n = nrow(A)
  return(matrix(solve(diag(n^2) - A %x% A, c(SU)), n, n))
}

simulate_nonexplosive_var_params <- function(n, p, g0, O0){
# Think of VAR(p) parameters concatenated into a matrix G = [G1 G2 ... Gp].
# Simulate vec(G) ~ N(g0, O0) truncated to have eigenvalues in the unit circle.
# 
#  n: integer
#  p: integer
# g0: n*n*p-length vector
# O0: n*n*p x n*n*p matrix
  
  if(n*n*p != length(g0)){
    stop('n*n*p is not equal to length of g0')
  }
  
  g <- mvrnorm(1, g0, O0)
  G <- matrix(g, n, n*p)
  Gcomp <- companion(G)
  while(isexplosive(Gcomp)){
    g <- mvrnorm(1, g0, O0)
    G <- matrix(g, n, n*p)
    Gcomp <- companion(G)
  }
  
  return(G)
  
}

var_stationary_params <- function(m, ARarray, S){
# Compute the mean and covariance matrix of the stationary 
# distribution of an n-variable VAR(p):
# y[t] = m + AR[1] * y[t - 1] + ... + AR[p] * y[t - p] + e[t], e[t] ~ (0, S)
#
#       m: n-vector
# ARarray: n x n x p array
#       S: n x n matrix
  
  # dimensions
  n = length(m)
  p = dim(ARarray)[3]
  
  # flatten array to n x n*p matrix
  AR = t(apply(ARarray, 1, c))
  
  # compute marginal E(y[t]) = inv(I - A1 - A2 - ... - Ap) * m
  # Source: Lutkepohl (2005) pg. 16
  
  VARmean = solve(diag(n) - rowSums(ARarray, dims = 2), m)
  
  # compute parameters of (centered) VAR(1) representation of VARMA
  # read off the marginal covariance matrix of y[t]
  # Source: Lutkepohl (2005) pg. 27
  
  ASU = varp2var1(AR, S)
  L = stationary_var1_covariance(ASU$A, ASU$SU)
  VARcov = L[1:n, 1:n]
  
  return(list(VARmean = VARmean, VARcov = VARcov))
}

simulate_stationary_var <- function(T, m, ARarray, S, Y0 = NULL){
# Simulate forward a stationary VAR(p):
# y[t] = m + AR[1] * y[t - 1] + ... + AR[p] * y[t - p] + e[t], e[t] ~ N(0, S)
#
#       T: integer
#       m: n-vector
# ARarray: n x n x p array
#       S: n x n matrix
#      Y0: p x n matrix
  
  # dimensions
  n = length(m)
  p = dim(ARarray)[3]
  
  # flatten array to n x n*p matrix
  AR = t(apply(ARarray, 1, c))
  
  if(is.null(Y0)){
    params = var_stationary_params(m, ARarray, S)
    Y0 = matrix(mvrnorm(p, params$VARmean, params$VARcov), p, n)
    # NOTE: really these should be correlated, but for now oh well (2/26/2024)
  }
  
  # construct lag vector
  x = c(t(Y0[p:1, ]))
  
  # pre-allocate storage for output and simulate that thing
  Y = matrix(0, T, n)
  
  for(t in 1:T){
    Y[t, ] = m + AR %*% x + mvrnorm(1, numeric(n), S)
    
    if(p > 1){
      x[(n + 1):(n * p)] = x[1:(n * (p - 1))]
    }
    x[1:n] = Y[t, ]
    
  }
  
  return(Y)
}

varp_design_matrix <- function(Y, p, intercept){
  TT = nrow(Y)
  n = ncol(Y)
  XX = matrix(1, TT - p, n * p)
  
  for(l in 1:p){
    stop  = TT - l
    start = stop - TT + p + 1
    
    XX[, (1:n) + n * (l - 1)] = Y[start:stop, ]
  }
  
  if(intercept == TRUE){
    XX = cbind(rep(1, TT - p), XX)
  }
  
  return(XX)
}

# ==============================================================================
# Extra stuff I don't need
# ==============================================================================

simulate_stationary_var1 <- function(T, m, A, S){
  n = length(m)
  params = var_stationary_params(m, array(A, dim = c(n, n, 1)), S)
  y = mvrnorm(1, params$VARmean, params$VARcov)
  Y = matrix(0, T, n)
  for(t in 1:T){
    y = m + A %*% y + mvrnorm(1, numeric(n), S)
    Y[t, ] = y
  }
  return(Y)
}

simulate_stationary_var2 <- function(T, m, A1, A2, S){
  n = length(m)
  params = var_stationary_params(m, array(c(A1, A2), dim = c(n, n, 2)), S)
  y1 = mvrnorm(1, params$VARmean, params$VARcov)
  y2 = mvrnorm(1, params$VARmean, params$VARcov)
  Y = matrix(0, T, n)
  for(t in 1:T){
    Y[t, ] = m + A1 %*% y1 + A2 %*% y2 + mvrnorm(1, numeric(n), S)
    y2 = y1
    y1 = Y[t, ]
  }
  return(Y)
}