# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

set.seed(8675309)

# ==============================================================================
# load full data
# https://www.bocsar.nsw.gov.au/Pages/bocsar_datasets/Offence.aspx
# ==============================================================================

crime_data = read.csv("datasets/full_nsw_crime_dataset.csv", header = TRUE)

# ==============================================================================
# subset of data for forecasting
# ==============================================================================

dates <- as.Date(crime_data$date, format = "%m/%d/%Y")

var_names <- c("murder", "attempted", "accessory_murder", "manslaughter", 
               "assault_nondomestic", "assault_police", "abduction", 
               "blackmail", "Escape_custody", "resist_officer")

my_labs <- c("Murder", "Attempted", "Accessory", "Manslaughter", "Non-domestic assault",
             "Officer assault", "Abduction", "Blackmail", "Escape custody", "Resist officer")

Y <- as.matrix(crime_data[, var_names])

T <- nrow(Y)
n <- ncol(Y)

# ==============================================================================
# run MCMC
# ==============================================================================

dfc_ndraw = 5000
dfc_nburn = 5000
dfc_nthin = 5

draws = DGFC.mcmc(Y, 
                  k.star = 3, 
                  ndraw = dfc_ndraw, 
                  burn = dfc_nburn, 
                  thin = dfc_nthin)

# ==============================================================================
# plot MA and ECDF
# ==============================================================================

par(mfrow = c(1, 1))

a = seq(0.1, 0.9, length.out = 9)
band_col = rgb(1, 0, 0, 0.15)

for(i in 1:n){
  plot(ecdf(Y[, i]), do.points = FALSE, col = "white", main = my_labs[i],
       yaxs = "i", xlab = "", ylab = "")
  plot_ma_band_discrete(draws$ma[[i]][, 1, 1], draws$ma[[i]][, 2, ], a, band_col)
  lines(ecdf(Y[, i]), col = "black", do.points = FALSE)
}

# ==============================================================================
# plot MA and ECDF
# ==============================================================================

t = 1
k = ncol(draws$G[, , 1])

Cdraws <- array(0, c(n * t, n, dfc_ndraw))

for(m in 1:dfc_ndraw){
  G = draws$G[, , m]
  S = draws$S[, , m]
  L = draws$Lambda[, , m]
  V = diag(draws$sig[, m])
  O = stationary_var1_covariance(G, S)
  C = cov2cor(L %*% O %*% t(L) + V)
  #E = full_var1_covariance(t, G, S)
  #C = cov2cor( (diag(t) %x% L) %*% E %*% (diag(t) %x% t(L)) + (diag(t) %x% V) )
  Cdraws[, , m] = C[1:(t*n), 1:n]
}

Cmean = apply(Cdraws, c(1,2), mean)

col_palette <- colorRampPalette(c("orange", "white", "blue"))(500)

heatmap(Cmean[1:n, ], Rowv = NA, Colv = NA, scale = "none",
        margins = c(3,4), col = col_palette, labCol = my_labs, labRow = (my_labs))

#heatmap(cor(Y), Rowv = NA, Colv = NA, scale = "none", col = col_palette,
#        margins = c(3,4), labCol = var_names, labRow = (var_names))

#heatmap(Cmean, Rowv = NA, Colv = NA, scale = "none",
#        margins = c(3,4), labCol = var_names)#, labRow = var_names, labCol = var_names)
#heatmap(Cmean[(n+1):(2*n), ], Rowv = NA, Colv = NA, scale = "none",
#        margins = c(3,4), labRow = var_names, labCol = var_names)
#heatmap(Cmean[(2*n+1):(3*n), ], Rowv = NA, Colv = NA, scale = "none",
#        margins = c(3,4), labRow = var_names, labCol = var_names)
#heatmap(Cmean[(3*n+1):(4*n), ], Rowv = NA, Colv = NA, scale = "none",
#        margins = c(3,4), labRow = var_names, labCol = var_names)
#heatmap(Cmean[(8*n+1):(9*n), ], Rowv = NA, Colv = NA, scale = "none",
#        margins = c(3,4), labRow = var_names, labCol = var_names)



