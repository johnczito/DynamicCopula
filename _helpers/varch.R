simulate_varch <- function(T, v, A, y0 = NULL){
  n = ncol(A)
  if(v <= n + 2){
    stop('Degrees of freedom too small: must have v > n + 2')
  }
  if(is.null(y0)){
    VARCHparams = varch_stationary_params(v, A)
    y0 = c(rmvt(n = 1, delta = VARCHparams$VARCHmean, 
                       sigma = VARCHparams$VARCHscale, 
                       df = VARCHparams$VARCHdf))
  }
  Y = matrix(0, T, n)
  y = y0
  for (t in 1:T){
    y = c(rmvt(n = 1, delta = numeric(n), 
                      sigma = (A + y %*% t(y)) / v, 
                      df = v))
    Y[t, ] = y
  }
  return(Y)
}

varch_stationary_params <- function(v, A){
  if(v <= ncol(A) + 2){
    stop('Degrees of freedom too small: must have v > n + 2')
  }
  return(list(VARCHdf = v - 1, VARCHmean = numeric(nrow(A)), VARCHscale = A / (v - 1)))
}