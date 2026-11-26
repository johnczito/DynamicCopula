# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

# ==============================================================================
# exercise settings
# ==============================================================================

n = 2
sample.sizes <- 25 * 2 ^ (0:3)
Tmax <- max(sample.sizes)
Nsamp <- length(sample.sizes)
ndraw = 1000
burn = 5000
thin = 5

# ==============================================================================
# Pick DGP, simulate fake data, and store true stationary distributions
# ==============================================================================

# varma, varch, varma_copula
dgp <- "varma_copula"

if(dgp == "varma"){
  
  set.seed(1)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # model settings
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  p = 3
  q = 6
  
  band_col = rgb(1, 0, 0, 0.1)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # generate random parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  m = rnorm(n)
  AR = simulate_nonexplosive_var_params(n, p, numeric(n*n*p), 0.1*diag(n*n*p))
  ARarray = array(c(AR), dim = c(n, n, p))
  MAarray = array(rnorm(n*n*q), dim = c(n, n, q))
  S = riwish(n + 1, diag(n))
  
  # how persistent is the time series process?
  
  abs(eigen(companion(AR))$values)
  
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
  
}else if(dgp == "varch"){
  
  set.seed(1)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # generate random parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  v = n + 2 + rexp(1)
  A = riwish(n + 1, diag(n))
  
  band_col = rgb(0, 0, 1, 0.15)
  
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
  
}else if(dgp == "varma_copula"){
  
  set.seed(1)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # get stationary distribution
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Finv1 <- function(x){qpois(x, lambda = 5)}
  Finv2 <- function(x){qst(x, alpha = 2, nu = 3)}
  
  Finv <- c(Finv1, Finv2)
  
  F1 <- function(x){ppois(x, lambda = 5)}
  F2 <- function(x){pst(x, alpha = 2, nu = 3)}
  
  F <- c(F1, F2)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # model settings
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  n = length(Finv)
  p = 1
  q = 1
  
  band_col = rgb(1, 0.6, 0, 0.1)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # generate random parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  AR = simulate_nonexplosive_var_params(n, p, numeric(n*n*p), 0.1*diag(n*n*p))
  ARarray = array(c(AR), dim = c(n, n, p))
  MAarray = array(rnorm(n*n*q), dim = c(n, n, q))
  S = riwish(n + 1, diag(n))
  
  # how persistent is the time series process?
  
  abs(eigen(companion(AR))$values)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulate fake data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Y = simulate_varma_copula(Tmax, Finv, ARarray, MAarray, S)
  
}

# ==============================================================================
# run MCMC on an expanding window of data
# ==============================================================================

MAdraws <- vector(mode = 'list', length = Nsamp)

for(l in 1:Nsamp){
  T <- sample.sizes[l]
  MAdraws[[l]] <- DGFC.mcmc(Y[1:T, ], ndraw = ndraw, thin = thin, burn = burn)$ma
  message(paste("Stage ", l, " of ", Nsamp, " done!", sep = ""))
}

# ==============================================================================
# inspect random set of trace plots
# ==============================================================================

nrows = 3
ncols = 3

par(mfrow = c(nrows, ncols), mar = c(2, 2, 2, 2))

for(q in 1:(nrows*ncols)){
  i <- sample(1:n, 1)
  t <- sample(1:nrow(MAdraws[[Nsamp]][[i]][, , 1]), 1)
  plot(1:ndraw, MAdraws[[Nsamp]][[i]][t, 2, ], 
       type = "l",
       main = paste(i, t, sep = ", "),
       ylim = c(0, 1))
}

# ==============================================================================
# plot true stationary CDF and MA credible bands for each sample size
# ==============================================================================

i = 1
a = seq(0.1, 0.9, length.out = 9)
x_grid = seq(min(Y[, i]) - 1, max(Y[, i]) + 1, length.out = 500)
y_vals = F[[i]](x_grid)

par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))

for(l in Nsamp:1){
  T <- sample.sizes[l]
  plot(ecdf(Y[1:T, i]), main = paste("Variable", i, "; T =", sample.sizes[l]), 
       do.points = FALSE, col.01line = NULL, xlim = c(min(Y[, i]), max(Y[, i])),
       col = "white")
  plot_ma_band_discrete(MAdraws[[l]][[i]][, 1, 1], MAdraws[[l]][[i]][, 2, ], a, band_col)
  lines(x_grid, y_vals, col = "black")
}

# ==============================================================================
# panel of plots for the paper
# ==============================================================================

i = 2
a = seq(0.1, 0.9, length.out = 9)
x_grid = seq(min(Y[, i]) - 1, max(Y[, i]) + 1, length.out = 500)
y_vals = F[[i]](x_grid)
which_sizes = 1:Nsamp
npanel = length(which_sizes)
add_legend = dgp == "varma"
add_yticks = dgp == "varma"
right_panel = 0.1
left_panel = 1.55

png(paste("_images/", dgp, "_consistency_paper.png", sep = ""), 
          width = 1.5, height = 6, units = "in", res = 1000)

par(mfrow = c(npanel, 1))

for(l in which_sizes){
  
  if(l == which_sizes[npanel]){
    par(mar = c(2, left_panel, 0, right_panel))
  }else if(l == which_sizes[1]){
    par(mar = c(0, left_panel, 0.4, right_panel))
  }else{
    par(mar = c(0, left_panel, 0, right_panel))
  }
  
  T <- sample.sizes[l]
  
  plot(ecdf(Y[1:T, i]), main = "", 
       do.points = FALSE, col.01line = NULL, xlim = c(min(Y[, i]), max(Y[, i])),
       col = "white", xaxt = "n", yaxt = "n", yaxs = "i")
  lines(x_grid, y_vals, col = "black", lwd = 2)
  plot_ma_band_discrete(MAdraws[[l]][[i]][, 1, 1], MAdraws[[l]][[i]][, 2, ], a, band_col)
  
  if(add_legend){
    legend("topleft", paste("T =", T), bty = "n", cex = 1.25)
  }
  
  if(l == which_sizes[npanel]){
    axis(1, cex.axis = 1.2)
  }
  
  if(l == which_sizes[1] & add_yticks){
    axis(2, at = c(0, 1), las = 1, cex.axis = 1.1)
  }
  
}

dev.off()

# ==============================================================================
# panel of plots for the slides
# ==============================================================================

i = 2
a = seq(0.1, 0.9, length.out = 9)
x_grid = seq(min(Y[, i]) - 1, max(Y[, i]) + 1, length.out = 500)
y_vals = F[[i]](x_grid)
which_sizes = 1:Nsamp
npanel = length(which_sizes)
add_legend = dgp == "varma"
add_yticks = TRUE
right_panel = 0.1
left_panel = 1.55

png(paste("_images/", dgp, "_consistency_slides.png", sep = ""), 
    width = 6, height = 1.5, units = "in", res = 1000)

par(mfcol = c(1, npanel))

for(l in which_sizes){
  
  if(l == which_sizes[npanel]){
    par(mar = c(0.5, 0, 0.5, 0.1))
  }else if(l == which_sizes[1]){
    par(mar = c(0.5, 2, 0.5, 0))
  }else{
    par(mar = c(0.5, 0, 0.5, 0))
  }
  
  T <- sample.sizes[l]
  
  plot(ecdf(Y[1:T, i]), main = "", 
       do.points = FALSE, col.01line = NULL, xlim = c(min(Y[, i]), max(Y[, i])),
       col = "white", xaxt = "n", yaxt = "n", yaxs = "i")
  lines(x_grid, y_vals, col = "black", lwd = 2)
  plot_ma_band_discrete(MAdraws[[l]][[i]][, 1, 1], MAdraws[[l]][[i]][, 2, ], a, band_col)
  
  if(add_legend){
    legend("topleft", paste("T =", T), bty = "n", cex = 1.25)
  }
  
  if(l == which_sizes[1] & add_yticks){
    axis(2, at = c(0, 1), las = 1, cex.axis = 1.1)
  }
  
}

dev.off()

