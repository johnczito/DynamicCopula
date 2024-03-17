# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

# ==============================================================================
# exercise settings
# ==============================================================================

n = 3
p = 3
q = 6
sample.sizes <- 10 * 2 ^ (1:5)
Tmax <- max(sample.sizes)
Nsamp <- length(sample.sizes)

#set.seed(8675309)

# ==============================================================================
# simulate parameter values
# ==============================================================================

m = rnorm(n)
AR = simulate_nonexplosive_var_params(n, p, numeric(n*n*p), 0.25*diag(n*n*p))
ARarray = array(c(AR), dim = c(n, n, p))
MAarray = array(rnorm(n*n*q), dim = c(n, n, q))
S = riwish(n + 1, diag(n))

# ==============================================================================
# calculate true stationary distribution
# ==============================================================================

VARMAparams = varma_stationary_params(m, ARarray, MAarray, S)

# ==============================================================================
# simulate fake data
# ==============================================================================

Y = simulate_stationary_varma(Tmax, m, ARarray, MAarray, S)

# ==============================================================================
# simulate fake data
# ==============================================================================

# JZ: This is awful. Need to loop through, store the MA output, and then do 
# all the plotting as post-processing

i = 1
a = c(0.01, 0.05, 0.1, 0.25, 0.5)
x_grid = seq(min(Y[, i]), max(Y[, i]), length.out = 500)
y_vals = pnorm(x_grid, mean = VARMAparams$VARMAmean[i], sd = sqrt(VARMAparams$VARMAcov[i, i]))

par(mfrow = c(Nsamp, 1), mar = c(2, 2, 2, 2))

for(l in 1:Nsamp){
  T <- sample.sizes[l]
  draws = DGFC.mcmc(Y[1:T, ])
  
  plot(ecdf(Y[1:T, i]), main = paste("Variable", i, ), do.points = FALSE, col.01line = NULL,
       xlim = c(min(Y[, i]), max(Y[, i])))
  plot_ma_band(draws$ma[[i]][, 1, 1], draws$ma[[i]][, 2, ], a, rgb(0, 0, 1, 0.1))
  lines(x_grid, y_vals, col = "red")
}

