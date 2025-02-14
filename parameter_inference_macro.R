# ==============================================================================
# load resources
# ==============================================================================

source("_packages.R")
source("_helpers/_helpers.R")

set.seed(8675309)

# ==============================================================================
# pull macro data
# ==============================================================================

var_names <- c("gdpgrowth", "pcegrowth", "bfigrowth", "resinvgrowth", "ipgrowth", "capu", 
               "employ", "hours", "ur", "gdppigrowth", "pcepigrowth", "ffr", "spread", "realSPgrowth")

fred_handles <- c("GDPC1", "PCE", "PNFI", "PRFI", "INDPRO", "CUMFNS", "PAYEMS",
                  "AWHI", "UNRATE", "GDPCTPI", "PCEPI", "DFF")
level <- c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE)


start_date <- as.Date("1965-01-01")
end_date <- as.Date("2019-12-31")
dates <- seq(1965.00, 2019.75, by = 0.25)

T <- length(dates)
n = 14

Y <- matrix(0, T, n)

for(i in 1:n){
  
  if(var_names[i] == "spread"){
    bondyield <- fredr(
      series_id = "DGS10",
      observation_start = start_date,
      observation_end = end_date,
      frequency = "q", 
      aggregation_method = "avg"
    )$value
    
    billyield <- fredr(
      series_id = "DTB3",
      observation_start = start_date,
      observation_end = end_date,
      frequency = "q", 
      aggregation_method = "avg"
    )$value
    
    raw_data = bondyield - billyield
    
  }else if(var_names[i] == "realSPgrowth"){
    getSymbols("^GSPC", src = "yahoo", from = start_date - 365 / 4, to = end_date)
    
    # Aggregate the data to a quarterly frequency by averaging
    GSPC_quarterly <- apply.quarterly(GSPC, colMeans, na.rm = TRUE)
    
    #sp <- as.vector(GSPC_quarterly$GSPC.Close[145:345]) # AUTOMATE THIS!
    sp <- as.vector(GSPC_quarterly$GSPC.Close)
    
    pcepi <- fredr(
      series_id = "PCEPI",
      observation_start = start_date - 365 / 4,
      observation_end = end_date,
      frequency = "q", 
      aggregation_method = "avg"
    )$value
    
    raw_data <- sp / pcepi
  }else{
    raw_data <- fredr(series_id = fred_handles[i],
                      observation_start = start_date - !level[i]*365/4,
                      observation_end = end_date,
                      frequency = "q", 
                      aggregation_method = "avg")$value
  }
  
  if(level[i] == TRUE){
    Y[, i] = raw_data
  }else{
    Y[, i] = 400 * diff(log(raw_data))
  }
  #plot(dates, Y[, i], type = 'l', main = var_names[i])
  
}

# ==============================================================================
# run MCMC
# ==============================================================================

dfc_ndraw = 5000
dfc_nburn = 5000
dfc_nthin = 3

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
  plot(ecdf(Y[, i]), do.points = FALSE, col = "white", main = var_names[i])
  plot_ma_band_discrete(draws$ma[[i]][, 1, 1], draws$ma[[i]][, 2, ], a, band_col)
  lines(ecdf(Y[, i]), col = "black", do.points = FALSE)
}

# ==============================================================================
# plot MA and ECDF
# ==============================================================================

library(fields)

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

lags <- c(expression(GDP[t]), 
          expression(PCE[t]), 
          expression(BFI[t]),
          expression(RI[t]),
          expression(IP[t]), 
          expression(CU[t]), 
          expression(E[t]), 
          expression(H[t]),
          expression(UR[t]), 
          expression(Infl1[t]),
          expression(Infl2[t]),
          expression(FFR[t]),
          expression(YS[t]),
          expression(SP[t]))
leads <- c(expression(GDP[t+1]), 
          expression(PCE[t+1]), 
          expression(BFI[t+1]),
          expression(RI[t+1]),
          expression(IP[t+1]), 
          expression(CU[t+1]), 
          expression(E[t+1]), 
          expression(H[t+1]),
          expression(UR[t+1]), 
          expression(Infl1[t+1]),
          expression(Infl2[t+1]),
          expression(FFR[t+1]),
          expression(YS[t+1]),
          expression(SP[t+1]))

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

png("_images/param_inference_macro_corr.png", 
    width = 4, height = 6, units = "in", res = 650)

par(mfcol = c(2, 1), mar = c(3, 3.75, 0.5, 2))

image.plot(Cmean_upper[1:n, ], col = col_palette[ceiling(500*(1-abs(min(Cmean[1:n, ])))):1000], 
           main = "", axes = FALSE, legend.mar = 3.5)

locations = seq(0, 1, length.out = n)

for(i in 2:n){
  for(j in 1:(i - 1)){
    if(significant_corrs[i, j] == TRUE){
      points(locations[i], locations[j], pch = 8)
    }
  }
}

axis(1, at = seq(0, 1, length.out = n), labels = lags, las = 2)
axis(2, at = seq(0, 1, length.out = n), labels = lags, las = 2)

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
axis(2, at = seq(0, 1, length.out = n), labels = leads, las = 1)

dev.off()

png("_images/param_inference_macro_CDFs.png", 
    width = 2, height = 6, units = "in", res = 650)

par(mfcol = c(4, 1), mar = c(2, 2, 2, 2))

a = seq(0.1, 0.9, length.out = 9)
band_col = rgb(1, 0, 0, 0.15)

my_labs <- c("Real GDP", "PCE", "BFI", "Residential invesment",
               "Industrial production", "Capacity utilization", "Payroll employment",
               "Hours", "Unemployment rate", "GDP PI inflation", "PCE inflation",
               "Federal funds rate", "Yield spread", "Real S&P 500")

for(i in c(9,10,12,13)){
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




# =======
# FOR THE DEFENSE
# =======


png("_images/param_inference_macro_corr.png", 
    width = 4, height = 6, units = "in", res = 650)

par(mfcol = c(2, 1), mar = c(3, 3.75, 0.5, 2))

image.plot(Cmean_upper[1:n, ], col = col_palette[ceiling(500*(1-abs(min(Cmean[1:n, ])))):1000], 
           main = "", axes = FALSE, legend.mar = 3.5)

locations = seq(0, 1, length.out = n)

for(i in 2:n){
  for(j in 1:(i - 1)){
    if(significant_corrs[i, j] == TRUE){
      points(locations[i], locations[j], pch = 8)
    }
  }
}

axis(1, at = seq(0, 1, length.out = n), labels = lags, las = 2)
axis(2, at = seq(0, 1, length.out = n), labels = lags, las = 2)

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
axis(2, at = seq(0, 1, length.out = n), labels = leads, las = 1)

dev.off()

png("_images/param_inference_macro_CDFs.png", 
    width = 2, height = 6, units = "in", res = 650)

par(mfcol = c(4, 1), mar = c(2, 2, 2, 2))

a = seq(0.1, 0.9, length.out = 9)
band_col = rgb(1, 0, 0, 0.15)

my_labs <- c("Real GDP", "PCE", "BFI", "Residential invesment",
             "Industrial production", "Capacity utilization", "Payroll employment",
             "Hours", "Unemployment rate", "GDP PI inflation", "PCE inflation",
             "Federal funds rate", "Yield spread", "Real S&P 500")

for(i in c(9,10,12,13)){
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

