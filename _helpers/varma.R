varma2var1 <- function(AR, MA, S){
# Calculate the parameters of the VAR(1) representation of a K-dimensional
# VARMA(p, q): Y[t] = A * Y[t - 1] + U[t], U[t] ~ N(0, SU) 
# Source: Lutkepohl (2005) Section 11.3.2 (pg. 426)
#
# AR: K x K * p matrix
# MA: K x K * q matrix
#  S: K x K covariance matrix
  
  # extract dimensions
  K = nrow(AR)
  p = ncol(AR) / K
  q = ncol(MA) / K
  
  # compute companion form transition matrix
  A11 = rbind(AR, cbind(diag(K * (p-1)), matrix(0, K * (p-1), K)))
  A12 = rbind(MA, matrix(0, K*(p-1), K*q))
  A21 = matrix(0, K*q, K*p)
  A22 = rbind(matrix(0, K, K*q), cbind(diag(K*(q-1)), matrix(0, K*(q-1), K)))
  A = rbind(cbind(A11, A12), cbind(A21, A22))
  
  # compute error covariance matrix
  SU = matrix(0, K*p + K*q, K*p + K*q)
  SU[1:K, 1:K] = S
  SU[K*p + (1:K), K*p + (1:K)] = S
  SU[1:K, K*p + (1:K)] = S
  SU[K*p + (1:K), 1:K] = S
  
  return(list(A = A, SU = SU))
}

varma_stationary_params <- function(m, ARarray, MAarray, S){
# Compute the mean and covariance matrix of the stationary 
# distribution of an n-variable VARMA(p, q):
# y[t] = m + AR[1] * y[t - 1] + ... + AR[p] * y[t - p]
#          + MA[1] * e[t - 1] + ... + MA[q] * e[t - q] + e[t], e[t] ~ (0, S)
#
#       m: n-vector
# ARarray: n x n x p array
# MAarray: n x n x q array 
#       S: n x n matrix

  # dimensions
  n = length(m)
  p = dim(ARarray)[3]
  q = dim(MAarray)[3]
  
  # flatten arrays to n x n*p and n x n*q matrices
  AR = t(apply(ARarray, 1, c))
  MA = t(apply(MAarray, 1, c))
  
  # pad out AR with extra zero matrices to compute VARMA representation
  if(p <= q){
    k = q - p + 1
    AR <- cbind(AR, matrix(0, n, n*k))
  }
  
  # compute marginal E(y[t]) = inv(I - A1 - A2 - ... - Ap) * m
  # Source: Lutkepohl (2005) pg. 424
  
  VARMAmean = solve(diag(n) - rowSums(ARarray, dims = 2), m)
  
  # compute parameters of (centered) VAR(1) representation of VARMA
  # read off the marginal covariance matrix of y[t]
  # Source: Lutkepohl (2005) pg. 430
  
  ASU = varma2var1(AR, MA, S)
  L = stationary_var1_covariance(ASU$A, ASU$SU)
  VARMAcov = L[1:n, 1:n]
  
  return(list(VARMAmean = VARMAmean, VARMAcov = VARMAcov))
}

simulate_stationary_varma <- function(T, m, ARarray, MAarray, S){
# Simulate forward a stationary VARMA(p, q):
# y[t] = m + AR[1] * y[t - 1] + ... + AR[p] * y[t - p]
#          + MA[1] * e[t - 1] + ... + MA[q] * e[t - q] + e[t], e[t] ~ (0, S)
#
#       T: integer
#       m: n-vector
# ARarray: n x n x p array
# MAarray: n x n x q array 
#       S: n x n matrix
  
  # dimensions
  n = length(m)
  p = dim(ARarray)[3]
  q = dim(MAarray)[3]
  
  # flatten arrays to n x n*p and n x n*q matrices
  AR = t(apply(ARarray, 1, c))
  MA = t(apply(MAarray, 1, c))
  
  # simulate initial conditions iid from the stationary distribution
  # NOTE: really these should be correlated, but for now oh well (2/26/2024)
  params = varma_stationary_params(m, ARarray, MAarray, S)
  Y0 = mvrnorm(p, params$VARMAmean, params$VARMAcov)
  
  # simulate MA lags 
  E0 = mvrnorm(q, numeric(n), S)
  
  # construct lag vectors 
  # Note: because the initial conditions are simulated iid, the order doesn't 
  #       matter at this point, but strictly speaking, the rows of matrices are 
  #       indexed in chronological order, and so the BOTTOM row of the matrices
  #       should be at the beginning of the vector, not the end as it is here
  x = c(t(Y0))
  e = c(t(E0))
  
  # pre-allocate storage for output and simulate that thing
  Y = matrix(0, T, n)
  
  for(t in 1:T){
    err = mvrnorm(1, numeric(n), S)
    Y[t, ] = m + AR %*% x + MA %*% e + err
    
    if(p > 1){
      x[(n + 1):(n * p)] = x[1:(n * (p - 1))]
    }
    
    if(q > 1){
      e[(n + 1):(n * q)] = e[1:(n * (q - 1))]
    }
    
    x[1:n] = Y[t, ]
    e[1:n] = err
  }
  
  return(Y)
}

simulate_varma_copula <- function(T, Finv, ARarray, MAarray, S){
  n = length(Finv)
  m = numeric(n)
  Z <- simulate_stationary_varma(T, m, ARarray, MAarray, S)
  W = varma_stationary_params(m, ARarray, MAarray, S)$VARMAcov
  Y = matrix(0, T, n)
  for(i in 1:n){
    Y[, i] = Finv[[i]](pnorm(Z[, i], mean = 0, sd = sqrt(W[i, i])))
  }
  return(Y)
}

# Notes:
# - initialization of VARMA in simulator ought to incorporate the 
#   autocorrelation between the initial conditions.
# - computing the companion form matrices is copy-pasted.
# - convert array to matrix is copy-pasted
# - `simulate_nonexplosive_var_params` takes a stand on how the vectorized 
#   VAR params should be reshaped to get the matrix of params, but I may want 
#   to change this later.
