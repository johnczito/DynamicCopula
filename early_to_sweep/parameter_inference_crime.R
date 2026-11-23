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

#par(mfrow = c(n, n), mar = c(2, 2, 0, 0))
#for(i in (n+1):(2*n)){
#  for(j in 1:n){
#    if(1 == 1){#j<i
#      hist(Cdraws[i, j, ], breaks = "scott", xlim = c(-1, 1), freq = FALSE,
#           main = "", border = NA, col = "lightblue")
#      abline(v = 0)
#    }else{
#      plot.new()
#    }
#  }
#}


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





t = 2
k = ncol(draws$G[, , 1])

Cdraws <- array(0, c(n * t, n, dfc_ndraw))

for(m in 1:dfc_ndraw){
  G = draws$G[, , m]
  S = draws$S[, , m]
  L = draws$Lambda[, , m]
  V = diag(draws$sig[, m])
  #O = stationary_var1_covariance(G, S)
  #C = cov2cor(L %*% O %*% t(L) + V)
  E = full_var1_covariance(t, G, S)
  C = cov2cor( (diag(t) %x% L) %*% E %*% (diag(t) %x% t(L)) + (diag(t) %x% V) )
  Cdraws[, , m] = C[1:(t*n), 1:n]
}

hpd_includes_zero = matrix(0, n * t, n)

for(i in 1:(2*n)){
  for(j in 1:n){
    hpdi = hdi(Cdraws[i, j, ], 1 - 0.05 / (n*t*n))
    hpd_includes_zero[i, j] = (hpdi[1] < 0) & (0 < hpdi[2])
  }
}

significant_corrs <- !hpd_includes_zero




lags <- c(expression(Mur[t]), expression(Att[t]), expression(Acc[t]), expression(Man[t]),
          expression(NDA[t]), expression(OA[t]), expression(Abd[t]), expression(B[t]),
          expression(EC[t]), expression(RO[t]))
leads <- c(expression(Mur[t+1]), expression(Att[t+1]), expression(Acc[t+1]), expression(Man[t+1]),
           expression(NDA[t+1]), expression(OA[t+1]), expression(Abd[t+1]), expression(B[t+1]),
           expression(EC[t+1]), expression(RO[t+1]))

Cmean = apply(Cdraws, c(1,2), mean)
Cmean_upper = Cmean

for(i in 1:n){
  for(j in 1:n){
    if(j > i){
      Cmean_upper[i, j] = NA
    }
  }
}

col_palette <- colorRampPalette(c("orange", "white", "blue"))(1000)

png("_images/param_inference_crime_corr.png", 
    width = 4, height = 6, units = "in", res = 650)

par(mfcol = c(2, 1), mar = c(3, 3.75, 0.5, 2))

image.plot(Cmean_upper[1:n, ], col = col_palette[ceiling(500*(1-abs(min(Cmean[1:n, ])))):1000], 
           main = "", axes = FALSE, legend.mar = 3.5)

locations = seq(0, 1, length.out = 10)

for(i in 2:n){
  for(j in 1:(i - 1)){
    if(significant_corrs[i, j] == TRUE){
      points(locations[i], locations[j], pch = 8)
    }
  }
}

axis(1, at = seq(0, 1, length.out = 10), labels = lags, las = 2)
axis(2, at = seq(0, 1, length.out = 10), labels = lags, las = 2)

image.plot(Cmean_upper[(n + 1):(2*n), ], col = col_palette[ceiling(500*(1-abs(min(Cmean[(n + 1):(2*n), ])))):floor(500*(1+max(Cmean[(n + 1):(2*n), ])))], 
           main = "", axes = FALSE, legend.shrink = 0.0001,
           legend.mar = 3.5)

cross_sig = significant_corrs[(n + 1):(2*n), ]

for(i in 1:n){
  for(j in 1:n){
    if(cross_sig[i, j] == TRUE){
      points(locations[i], locations[j], pch = 8)
    }
  }
}

#axis(1, at = seq(0, 1, length.out = 10), labels = lags, las = 2)
axis(2, at = seq(0, 1, length.out = 10), labels = leads, las = 1)

dev.off()

png("_images/param_inference_crime_CDFs.png", 
    width = 2, height = 6, units = "in", res = 650)

par(mfcol = c(4, 1), mar = c(2, 2, 2, 2))

a = seq(0.1, 0.9, length.out = 9)
band_col = rgb(1, 0, 0, 0.15)

for(i in c(4, 7, 8, 9)){
  plot(ecdf(Y[, i]), do.points = FALSE, col = "white", main = my_labs[i],
       yaxs = "i", xlab = "", ylab = "", bty = "n", yaxt = "n")
  plot_ma_band_discrete(draws$ma[[i]][, 1, 1], draws$ma[[i]][, 2, ], a, band_col)
  lines(ecdf(Y[, i]), col = "black", do.points = FALSE)
  abline(h = c(0, 1))
  axis(2, at = c(0, 1), las = 1, cex.axis = 1.1, col = "white")
  if(i == 4){
    legend("bottomright", c("ECDF", "MA"), lty = 1, col = c("black", "red"), bty = "n")
  }
}

dev.off()

