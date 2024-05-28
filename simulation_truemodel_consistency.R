# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

# ==============================================================================
# exercise settings
# ==============================================================================

n = 2
p = 1
sample.sizes <- 25 * 2 ^ (0:5)
Tmax <- max(sample.sizes)
Nsamp <- length(sample.sizes)
ndraw = 1000
burn = 5000
thin = 10
band_col = rgb(0, 0, 1, 0.1)

set.seed(8675309)

# ==============================================================================
# Ground truth params
# ==============================================================================

Finv1 <- function(x){qgamma(x, shape = 1, rate = 1)}
Finv2 <- function(x){qst(x, alpha = 2, nu = 3)}
Finv <- c(Finv1, Finv2)

F1 <- function(x){pgamma(x, shape = 1, rate = 1)}
F2 <- function(x){pst(x, alpha = 2, nu = 3)}
F <- c(F1, F2)

A <- simulate_nonexplosive_var_params(n, p, numeric(n*n*p), diag(n*n*p))
SU <- riwish(n + 1, diag(n))
GS = normalize_var1_params(A, SU)
G = GS$G
S = GS$S
O0 = stationary_var1_covariance(G, S)

# ==============================================================================
# Simulate fake data
# ==============================================================================

Z = simulate_stationary_var1(Tmax, numeric(n), G, S)

Y = matrix(0, Tmax, n)
for(i in 1:n){
  Y[, i] = Finv[[i]](pnorm(Z[, i], mean = 0, sd = 1))
}

# ==============================================================================
# how long does the big'n take?
# ==============================================================================

mystart = Sys.time()

#post_draws = DGC.mcmc(Y, ndraw = ndraw, thin = thin, burn = burn)

endstart = Sys.time()

endstart - mystart

# ==============================================================================
# run MCMC on an expanding window of data
# ==============================================================================

MAdraws <- vector(mode = 'list', length = Nsamp)
G_draws = array(0, c(n, n, ndraw, Nsamp))
Sigma_draws = array(0, c(n, n, ndraw, Nsamp))

for(l in 1:Nsamp){
  T <- sample.sizes[l]
  post_draws <- DGC.mcmc(Y[1:T, ], ndraw = ndraw, thin = thin, burn = burn)
  MAdraws[[l]] = post_draws$ma
  G_draws[, , , l] = post_draws$G
  Sigma_draws[, , , l] = post_draws$S
  message(paste("Stage ", l, " of ", Nsamp, " done!", sep = ""))
}

# ==============================================================================
# plots true stationary CDF and MA credible bands for each sample size
# ==============================================================================

# fix this when visualizing the discrete MA and everything

i = 2
a = seq(0.1, 0.9, length.out = 9)#c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75)
x_grid = seq(min(Y[, i]), max(Y[, i]), length.out = 500)
y_vals = F[[i]](x_grid)

#par(mfrow = c(Nsamp, 1), mar = c(2, 2, 2, 2))
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))

for(l in Nsamp:1){
  T <- sample.sizes[l]
  plot(ecdf(Y[1:T, i]), main = paste("Variable", i, "; T =", sample.sizes[l]), 
       do.points = FALSE, col.01line = NULL, xlim = c(min(Y[, i]), max(Y[, i])),
       col = "white")
  #plot_ma_band(MAdraws[[l]][[i]][, 1, 1], MAdraws[[l]][[i]][, 2, ], a, rgb(0, 0, 1, 0.1))
  plot_ma_band_discrete(MAdraws[[l]][[i]][, 1, 1], MAdraws[[l]][[i]][, 2, ], a, band_col)
  lines(x_grid, y_vals, col = "black")
}

# how persistent is the time series process?

abs(eigen(G)$values)

# ADD TRACE PLOTS!

par(mfrow = c(1, 1))

plot(1:ndraw, MAdraws[[Nsamp]][[i]][12, 2, ], type ="l")


# ==============================================================================
# plot posterior convergence for parameters
# ==============================================================================

par(mfrow = c(2 * n, n), mar = c(2, 2, 1, 1))

for(i in 1:n){
  for(j in 1:n){
    
    true_value = G[i, j]
    
    my_draws = matrix(0, ndraw, Nsamp)
    
    for(k in 1:Nsamp){
      my_draws[, k] = G_draws[i, j, , k]
    }
    
    if(i == 1 & j == 1){
      boxplot(my_draws, 
              outline = FALSE,
              names = sample.sizes, 
              xlab = "",
              ylab = "",
              main = "")
    }else{
      boxplot(my_draws, 
              outline = FALSE,
              xaxt = "n", 
              xlab = "",
              ylab = "",
              main = "")
    }
    
    abline(h = true_value, col = "red", lwd = 2)
    legend("bottomright", paste("G[", i, ", ", j, "]", sep = ""), bty = "n")
    
  }
}

for(i in 1:n){
  for(j in 1:n){
    if(i < j){
      plot.new()
    } else {
      true_value = S[i, j]
      
      my_draws = matrix(0, ndraw, Nsamp)
      
      for(k in 1:Nsamp){
        my_draws[, k] = Sigma_draws[i, j, , k]
      }
      
      boxplot(my_draws, 
              outline = FALSE,
              names = NULL, 
              xlab = "",
              ylab = "",
              xaxt = "n",
              main = "")
      abline(h = true_value, col = "red", lwd = 2)
      legend("topright", paste("Sigma[", i, ", ", j, "]", sep = ""), bty = "n")
    }
    
  }
}

par(mfrow = c(4, 1), mar = c(2, 2, 1, 1))

i = 1
a = seq(0.1, 0.9, length.out = 9)#c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75)
x_grid = seq(min(Y[, i]), max(Y[, i]), length.out = 500)
y_vals = F[[i]](x_grid)

#par(mfrow = c(Nsamp, 1), mar = c(2, 2, 2, 2))
par(mfrow = c(4, 1))

for(l in c(1, 2, 3, 4)){
  #if(l == Nsamp){
  #  par(mar = c(2, 1, 0, 1))
  #}else{
  #  par(mar = c(0, 1, 0, 1))
  #}
  
  T <- sample.sizes[l]
  
  if(l == 4){
    plot(ecdf(Y[1:T, i]), main = "", 
         do.points = FALSE, col.01line = NULL, xlim = c(min(Y[, i]), max(Y[, i])),
         col = "white",yaxt = "n")
  }else{
    plot(ecdf(Y[1:T, i]), main = "", 
         do.points = FALSE, col.01line = NULL, xlim = c(min(Y[, i]), max(Y[, i])),
         col = "white", xaxt = "n", yaxt = "n")
  }
  #plot_ma_band(MAdraws[[l]][[i]][, 1, 1], MAdraws[[l]][[i]][, 2, ], a, rgb(0, 0, 1, 0.1))
  lines(x_grid, y_vals, col = "black", lwd = 2)
  plot_ma_band_discrete(MAdraws[[l]][[i]][, 1, 1], MAdraws[[l]][[i]][, 2, ], a, band_col)
  legend("bottomright", paste("T =", T), bty = "n", cex = 2)
  
}

