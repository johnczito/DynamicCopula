plot_ma_band <- function(ma_support, ma_vals, alphas, col){
  # ma_support: M vector
  # ma_vals: M x ndraw matrix
  
  K = length(alphas)
  M = length(ma_support)
  qU = matrix(0, K, M)
  qL = matrix(0, K, M)
  for(j in 1:M){
    qL[, j] = quantile(ma_vals[j, ], prob = a / 2)
    qU[, j] = quantile(ma_vals[j, ], prob = 1 - (a / 2))
  }
  for(k in 1:K){
    polygon(c(ma_support, rev(ma_support)), c(qL[k, ], rev(qU[k, ])), col = col,
            border = NA)
  }
}