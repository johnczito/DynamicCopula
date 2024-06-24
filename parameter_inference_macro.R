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

var_names <- c("gdp", "pce", "bfi", "resinv", "ip", "capu", 
               "employ", "hours", "ur", "gdpinfl", "pceinfl", "ffr", "spread", "S&P")

col_palette <- colorRampPalette(c("orange", "white", "blue"))(500)

heatmap(Cmean[1:n, ], Rowv = NA, Colv = NA, scale = "none",
        margins = c(3,4), col = col_palette, labCol = var_names, labRow = (var_names))

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



